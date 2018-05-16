program galtogadget
implicit none
! main program that handles galactics to gadget conversion

! load variables
include 'variables.inc'
include 'units.inc'
include 'rvarrays.inc'

real*8 :: halofrac_comp
real*8 :: getdr
integer*8 :: ipartthing

!! first section is to store in the variables
print*,"--------------"
write (*,*) "running galactics to gadget"
print*,"--------------"
print*,"Output file        = " ,outFile
print*,"Dark halo file     = " ,darkhalofile
print*,"Stellar Disc file  = " ,discfile
print*,"Stellar Bulge file = " ,bulgefile
print*,"lcool              = " ,lcool
print*,"Metallicity        = " ,zmet0
print*,"Star fraction      = " ,starfrac
print*,"--------------"
print*,"--------------------"
print*, "units"
write(*,*) 'munit/g      =',munit
write(*,*) 'lunit/cm        =',lunit
write(*,*) 'tunit/s         =',tunit
write(*,*) 'vunit/(cm/s)    =',vunit
write(*,*) 'nunit/(cm^{-3}) =',nunit
write(*,*) 'eunit/erg       =',eunit
write(*,*) 'Kunit/K         =',Kunit
write(*,*) 'Cunit/(cm^{-3}/s) = ', cunit
print*,"--------------------"
print *, "running procedures now"
print*,"--------------------"


! dont need this
!zmet = zmet0
!call createcool(zmet)

! initial temperature is 10000K


!! LETS LOAD THE DARK HALO
iobj = 0
darkmid = 0

open(unit=43,file=darkhalofile)
read(43,*) ndark
print *,"ndark=",ndark
do loopi=1,ndark
   iobj = iobj + 1
   read(43,*)    r(4  ,iobj),r(1:3,iobj),v(1:3,iobj)
   darkmid     = darkmid     + r(1:3,iobj)
   v(1:3,iobj) = v(1:3,iobj) * vconffactor     ! 100 km/s to internal   
   r(4  ,iobj) = r(4  ,iobj) * mconffactor ! Weird units tointernal
   r(4  ,iobj) = r(4  ,iobj) * (1.0-halofrac)
   itype(iobj) = itypedark

enddo
close(43)
print *, "dark read"
! find the midpoint
darkmid = darkmid/ndark
rmdark  = r(4,iobj) ! mass of a single dm particle

print *, rmdark, ndark

print*, "Total Dark mass       = ", rmdark*ndark, " Msun e10"
print*, "DM particle mass   = ", rmdark, " Msun e10"



eminin = 1.0e4 / Kunit
ehalo = 1.0e6 / Kunit

! mgas = 0.5 dm mass

Mm =  mgas_p_frac * halofrac * rmdark * ndark * 1e10
mgas = halofrac * rmdark * 1e10
! NOTES
!     Call programme to find the positions of the Hot Halo Gas, based upon
!     rho(r) = rho0 * (1 + (r/rc)**2)**-alpha
!       output: ngas = number of hot halo gas particles
!       output: Mt = total mass of gas halo
!       intput: rc = core radius (kpc)
!       intput: rf = outer radius (kpc)
!       intput: alpha = 1.5*beta for the profile of the halo
!       intput: Mm (M_sun) gas mass within a radius of rm (kpc)
!       intput: mgas (M_sun) mass of a gas particle


if (halofrac .gt. 0.000) then

call makehalo(nhalogas, rc, rm, rf, Mm, Mtotal, mgas, alpha,rho0, ipr)
call calculate_h0_fun

open(unit=43,file=darkhalofile)
read(43,*) ndark
print *,"ndark=",ndark
close(43)

print *, "ndark again", ndark
open(unit=43,file='halogas_positions.txt')
read(43, *) nhalogas, Mtotal
close(43)

halofrac_comp = Mtotal * halofrac  !!! silly galactics mass units
print*, 'The fraction of the halo that is being converted from" &
       dark matter to hot gas is ', halofrac_comp
print*, 'rho0 = ',rho0,'M_sun/kpc^3'

rho0cgs = rho0*Msun/(lunit)**3
rccgs = rc*(lunit)
rcint = rc
rcinti2= 1.0/rcint**2

print*, 'rc = ',rc,'kpc = ',rccgs,'cm'
print*, 'rc = ',rcint,'(internal)'
print*, '1/rc^2 = ',rcinti2,'(internal)'


else
nhalogas = 0
endif

Mgal = 0.0

!! Read in disc potential profile
open(unit=43,file=proffile)
!     First two lines are titles, 3rd line is just 0 really


read(43,*)
read(43,*)
read(43,*)

!!     Load in circular velocit
do i=1,npot
   read(43,*) vcrad(i),realdummy,realdummy,realdummy,vcvc(i)
   ! vcrad is already in kpc
   vcvc (i) = vcvc (i) * vconffactor  ! 100 km/s to internal
end do
close(43)


!! create the hot gas halo

! this was created as a DM halo as per Galactics, just using considerably fewer
! particles
! The particle masses will be set to replace the mass we removed above
! h is calculated the same way here as we do below in the gas disc

if (nhalogas > 0) then
open(unit=431, file='find_rh_halo.txt')
open(unit=43,file='halogas_positions.txt')

read(43, *) nhalogas, Mtotal
print *,"nhalogas=",nhalogas
rmgashalo = rmdark*ndark/( (1.0-halofrac) * nhalogas ) * halofrac
do loopi=1,nhalogas
   iobj = iobj + 1
   read(43,*) idummy, r(1:3,iobj),realdummy
   read(431,*) realdummy, realdummy, realdummy, h(iobj), realdummy, realdummy
   ! convert h from kpc to internal (it is the same) and scale by 1/2 since h is defined as the
   ! half radius of a sphere containing the desired number (30-80) particles
   h(iobj) = h(iobj)*0.5
   r(4  , iobj) = rmgashalo !! this will be eventually overwritten
   dn(iobj) = r(4,iobj)
   dnp   (iobj) = dn(iobj)
   sm    (iobj) = 0.0
   dftime(iobj) = 0.0
   sfr   (iobj) = 0.0
   ! This is a temporary itype number; will be reset to itgas later
   itype( iobj) = itypegastemp

enddo
close(43)
close(431)

print*, "hot halo particle mass = ", rmgashalo, " 1e10 Msun"
print*, "hot halo mass   = ",rmgashalo*nhalogas, " 1e10 Msun"
print*, "total halo mass = ",rmdark*ndark + rmgashalo*nhalogas , " 1e10 Msun"
print*, "DM, halo gas particle mass:", rmdark, rmgashalo  , " 1e10 Msun"

endif
!! bulge

open(unit=43,file=bulgefile)
read(43,*) nbulge
print *,"Bulge stars=",nbulge
do loopi=1,nbulge
    iobj        = iobj + 1
    read(43,*)    rmbulge , realdummy, realdummy, realdummy, realdummy, realdummy, realdummy
enddo
close(43)

!! DISK

open(unit=43,file=discfile)
read(43,*) ndisc
print *,"Disc stars=",ndisc

do loopi=1,ndisc
   iobj        = iobj + 1
   read(43,*)    rmdiskstar , realdummy, realdummy, realdummy, realdummy, realdummy, realdummy
end do
close(43)

rmdiskstar = r(4,iobj)
print *, "disk particle mass = ", rmdiskstar,  " 1e10 Msun"
print*, "Stellar disk mass = ",r(4,iobj)*ndisc, " 1e10 Msun"

gasdiscmass   = rmdiskstar/starfrac*(1-starfrac)*ndisc
bulgemass     = rmbulge*nbulge
realdummy = ((1-starfrac)*rmdiskstar * mconffactor + &
          (nbulge*rmbulge)/ndisc)

if (nhalogas > 0) then
ngas = ndiscgas + nhalogas
else
ngas = ndiscgas
endif
print *, "Gas disk particle mass  = ", r(4,iobj), " 10^10 Msun"
print*, "Gas disk mass = ",r(4,iobj) * ndisc, " 10^10 Msun"

!     Total number of objects in the galaxy
nobj=ndark+nhalogas


!! ok, this is the part where we do some moving around and reallocation.
!! and this is when the design stars deviating away from the origonal slightly
!! since array allocations will speed things up


! 1 = x
! 2 = y
! 3 = z
! 4 = mass
! 5 = vx
! 6 = vy
! 7 = vz
! 8 = energy
print *, ndark, "ndark"
print *, "allocating arrays"
print *, nobj
allocate(tempgas(1:8,1:ngas))
print *, "gas"

if (nhalogas > 0) then
allocate(temphalogas(1:8,1:nhalogas))
endif
print *, "halogas"
allocate(tempothergas(1:8,1:ndiscgas))
print *, "othergas"
allocate(tempdark(1:8,1:ndark))
print *, "dark"
allocate(tempbulge(1:8,1:nbulge))
print *, "bulge"
allocate(tempdisk(1:8,1:ndisc))
print *, "disk"
allocate(tempstar(1:8,1:nstar))
print *, "star"

! commenting out
!allocate(tempsink(1:8,1))
!print *, "sink"

print *, "allocating outputs"

allocate(pos(1:3,1:nobj))
allocate(vel(1:3,1:nobj))
allocate(id(1:nobj))
allocate(masses(1:nobj))
allocate(edum(1:nhalogas))
allocate(dndum(1:nhalogas))
allocate(hdum(1:nhalogas))
!allocate(edum(1:nhalogas))
!allocate(dndum(1:nhalogas))
!allocate(hdum(1:nhalogas))



!! move things into the temp arrays for further computation
!! order is
!! 0 to ndark
!! ndark to halogass
!! bulge stars
!! disk stars
!! disk gas
!! black hole

iobj = 0
do i = 1, ndark
        iobj = iobj + 1
        tempdark(1:4,i) = r(:,iobj)
        tempdark(5:8,i) = v(:,iobj)
end do

print *,  tempdark(4,i)

if (nhalogas > 0) then
do i = 1, nhalogas
        iobj = iobj + 1
        temphalogas(1:4,i) = r(:,iobj)
        temphalogas(5:8,i) = v(:,iobj)
 

end do
endif


!!! overwrite halo gas mass
temphalogas(4,:) = realdummy

!iobj = iobj + 1
!tempsink(1:4,1) = r(:,iobj)
!tempsink(5:8,1) = v(:,iobj)




!! gas temperatures
istrt = 1

! parallelise


if (nhalogas > 0) then
write(*,*) "Calculating halo gas temperatures - brute force method 1"
!! computed via hydrostatic equilibirum ... see
!!http://adsabs.harvard.edu/abs/2012MNRAS.421.2170W
print *, "allocating arrays"
allocate(r_gas_sort(1:nhalogas))
allocate(id_gas_sort(1:nhalogas))

allocate(r_dark_sort(1:ndark))
allocate(id_dark_sort(1:ndark))

allocate(r_sort(1:nobj))
allocate(id_sort(1:nobj))

print *, "initial setup"
do i = 1, nobj
   r_sort(i) = sqrt((r(1,i) ** 2.0) + (r(2,i) ** 2.0) + (r(3,i) ** 2.0))
   id_sort(i) = i
enddo 

do i = 1, nhalogas
   r_gas_sort(i) = sqrt((temphalogas(1,i) ** 2.0) + (temphalogas(2,i) ** 2.0) + (temphalogas(3,i) ** 2.0))
   id_gas_sort(i) = i
enddo

do i = 1, ndark
   r_dark_sort(i) = sqrt((tempdark(1,i) ** 2.0) + (tempdark(2,i) ** 2.0) + (tempdark(3,i) ** 2.0))
   id_dark_sort(i) = i
enddo



r_max_val = maxval(r_sort(:))
dr_bin_width = r_max_val / dr_bins

! create the dr bins
print *, "creating dr bins"
do i = 1, dr_bins
   dr_vals(i) = i * dr_bin_width
enddo

! sort both r arrays
dm_vals(:) = 0.0
dd_vals(:) = 0.0
dtemp_vals(:) = 0.0
dr_count_vals(:) = 0

dr_dark_count_vals(:) = 0
dd_dark_vals(:) = 0.0

print *, nhalogas
print *, "sorting array 1"
call SSORTID(r_gas_sort, id_gas_sort,nhalogas,2)
print *, "sorting array 2"
call SSORTID(r_sort,id_sort,nobj,2)

call SSORTID(r_dark_sort,id_dark_sort,ndark,2)

!! fill mass dr array
i_bin_val = 1
dr_bin_val = dr_bin_width

print *, "computing mass in bins"
do i = 1, nobj

    if (r_sort(i) .gt. dr_bin_val) then
       r_temp = r_sort(i) - dr_bin_val
       ! check again
       if (r_sort(i) .gt. r_temp) then
          dr_bin_val = dr_bin_val + dr_bin_width
          i_bin_val = i_bin_val + 1
          ! copy across
          dm_vals(i_bin_val) = dm_vals(i_bin_val - 1)
       else
          ! need to do this more than once xD
          do while (r_temp .lt. dr_bin_val)
             dr_bin_val = dr_bin_val + dr_bin_width
             i_bin_val = i_bin_val + 1
             ! copy across
             dm_vals(i_bin_val) = dm_vals(i_bin_val - 1)
          end do
       end if
    end if

    dm_vals(i_bin_val) = dm_vals(i_bin_val) + r(4,id_sort(i))
    !print *, dm_vals(i_bin_val), r(4,id_sort(i)), id_sort(i), i
enddo

!! cleaning up
do i = 2, dr_bins
   if (dm_vals(i) .eq. 0.0) then
      if (dm_vals(i-1) .gt. 0.0) then
         dm_vals(i) = dm_vals(i-1)
      end if
   end if
enddo 

m_total_check = 0
do i = 1, nobj
   m_total_check = m_total_check + r(4,i)
end do

print *, "checking total mass ",  m_total_check, dm_vals(dr_bins)

i_bin_val = 1
dr_bin_val = dr_bin_width
print *, "computing gas density in bins"
do i = 1, nhalogas
    if (r_gas_sort(i) .gt. dr_bin_val) then
       r_temp = r_gas_sort(i) - dr_bin_val
       ! check again
       if (r_gas_sort(i) .gt. r_temp) then
          dr_bin_val = dr_bin_val + dr_bin_width
          i_bin_val = i_bin_val + 1
       else
          ! need to do this more than once xD
          do while (r_temp .lt. dr_bin_val)
             dr_bin_val = dr_bin_val + dr_bin_width
             i_bin_val = i_bin_val + 1
          end do
       end if
    end if
    !! essentially for a loop counter
    !! need to check if there is even a particle in this bin
    dr_count_vals(i_bin_val) = dr_count_vals(i_bin_val) + 1.0
    !! mass
    dd_vals(i_bin_val) = dd_vals(i_bin_val) + temphalogas(4,id_gas_sort(i))
enddo

i_bin_val = 1
dr_bin_val = dr_bin_width
print *, "computing dark density in bins"
do i = 1, ndark
    if (r_dark_sort(i) .gt. dr_bin_val) then
       r_temp = r_dark_sort(i) - dr_bin_val
       ! check again
       if (r_dark_sort(i) .gt. r_temp) then
          dr_bin_val = dr_bin_val + dr_bin_width
          i_bin_val = i_bin_val + 1
       else
          ! need to do this more than once xD
          do while (r_temp .lt. dr_bin_val)
             dr_bin_val = dr_bin_val + dr_bin_width
             i_bin_val = i_bin_val + 1
          end do
       end if
    end if
    !! essentially for a loop counter
    !! need to check if there is even a particle in this bin
    dr_dark_count_vals(i_bin_val) = dr_dark_count_vals(i_bin_val) + 1.0
    !! mass
    dd_dark_vals(i_bin_val) = dd_dark_vals(i_bin_val) + tempdark(4,id_dark_sort(i))
enddo


!1 convert dd_vals into gas densities

do i = 1, dr_bins
    if (i == 1) then
       vol_area = (4.0/3.0) * pi_val * (dr_vals(i) ** 3.0) * (lunit ** 3)
    else
       vol_area =  (4.0/3.0) * pi_val * ((((dr_vals(i) ** 3.0) *  lunit ** 3.0)) - ((dr_vals(i-1) ** 3.0) * lunit ** 3.0))
    end if
    dd_vals(i) = dd_vals(i) * munit / vol_area
    dd_dark_vals(i) = dd_dark_vals(i) * munit / vol_area     
end do


!! finally do the integral
print *, "integrating"
do i = 1, dr_bins
   !! now remember, the approximation to an integral is a sumation of width dr
   if (dr_count_vals(i) .ne. 0) then
       do j = i, dr_bins
           dtemp_vals(i) = dtemp_vals(i) + ((((dd_vals(j) * dm_vals(j) / (dr_vals(j) ** 2.0)) * dr_bin_width)) * (munit / lunit))
       end do
   
       !! integral prefactor last
       dtemp_vals(i) = dtemp_vals(i) * (tempconst / dd_vals(i) ) / kunit
   else
       dtemp_vals(i) = 0.0
   endif

end do
    
!! now assign temperatures to particles
print *, "assigning temp to particles"
i_bin_val = 1
dr_bin_val = dr_bin_width

i_bin_val = 1
do i = 1, dr_bins
   if (dr_count_vals(i) .ne. 0) then
      do j = 1, dr_count_vals(i)
          temphalogas(8,id_gas_sort(i_bin_val)) = dtemp_vals(i)
          i_bin_val = i_bin_val + 1
      end do
   endif
end do


write(*,*) "Temperature calculation complete"


write(*,*) "Calculating halo gas velocities - brute force method 2"

!! efficient method
do i = 1, dr_bins
   vcirc_arr(i) = sqrt(G*dm_vals(i) * munit/( dr_vals(i) * lunit) ) / vunit
  ! print *, vcirc_arr(i), dm_vals(i), dr_vals(i)
end do

!! compute r_alpha and v_alpha

r_alpha = dr_vals(1500)
v_alpha = vcirc_arr(1500)

!! see that paper in robs email
do i = 1, nhalogas
    rad = sqrt((temphalogas(1,i)**2) + (temphalogas(2,i)**2) + (temphalogas(3,i)**2) )
    phi_vtheta = acos(temphalogas(3,i) / rad)
    v_theta = v_alpha * ((rad / r_alpha) ** alpha_vtheta) * (sin(phi_vtheta) ** chi_vtheta)
    theta_vtheta = atan(temphalogas(2,i) / temphalogas(1,i))
    !! abs(temphalogas(1,i)) / temphalogas(1,i) is either 1 or -1 depending on
    !! the origional value of abs(temphalogas(1,i))
    temphalogas(5,i) = abs(v_theta *sin(theta_vtheta)) * (abs(temphalogas(2,i))/temphalogas(2,i)) 
    temphalogas(6,i) = -abs(v_theta * cos(theta_vtheta)) * (abs(temphalogas(1,i)) / temphalogas(1,i))
    temphalogas(7,i) = 0.0

! overwriting spin for now
    temphalogas(5,i) = 0.0
    temphalogas(6,i) = 0.0
end do    
endif


!! lets compute an analytical vcirc curve

write(*,*) "computing analytical results for vcirc and temp"

do i = 1, 200
        if (i == 1) then
                r_circ_analytical(i) = 0.0
                v_circ_analytical(i) = 0.0
                t_circ_analytical(i) = 0.0
        else
                r_circ_analytical(i) = i - 1 ! in kpc
                rad2 = r_circ_analytical(i) ** 2
                Menc = 0.0
                do j = istrt,iobj 
                        dis2 = (r(1,j))**2 + (r(2,j))**2 + (r(3,j))**2
                        if (dis2.lt.rad2) Menc = Menc + r(4,j)
                end do
                v_circ_analytical(i) = sqrt(G*Menc*munit/( sqrt(rad2) * lunit) ) / vunit
                t_circ_analytical(i) = mum * (v_circ_analytical(i) ** 2 * (vunit **2))  / (2 * kb) 
        end if
end do


write(*,*) "outputting the analytical results"


open (unit=456, file='analytical_results.txt')
do i = 1, 200
      write(456, '(3Es)') r_circ_analytical(i), v_circ_analytical(i), t_circ_analytical(i)
end do
close(456)

if (nhalogas > 0) then
open (unit=123, file='00halotemp.txt')
do i = 1, nhalogas
      rad2   = (temphalogas(1,i))**2 + (temphalogas(2,i))**2 + (temphalogas(3,i))**2
      write(123, '(6Es)') temphalogas(1:4, i), sqrt(rad2), temphalogas(8,i)*(Kunit**2)
      !if (temphalogas(8,i).ge.1e8) write(*, '(a, 5Es)') 'hot', temphalogas(1:4, i), temphalogas(8,i)*(Kunit**2)
end do
close(123)
endif

!! loop over gas
! gas density profile and  dm density profile already computed
! dm vx phi profile r = 180 and x > 0
! dm vx phi profile r = 180 and x > 0

!! loop over gas

!print *, "print angular vx distribution:"

!i_bin_val = 0
!do i = 3000, dr_bins
!       if (dr_count_vals(i) .ne. 0) then
!                i_bin_val = i_bin_val + 1
!       endif
!end do

!print *, i_bin_val, "i_bin_val"
!j = 0

!! count how many vx
!do i = nhalogas - 3000, nhalogas
!        if ((temphalogas(1,id_gas_sort(i)) > 0.0) .and. (ABS(temphalogas(3,id_gas_sort(i))) < 10.0)) then
!                j = j + 1
!        endif
!enddo
!print *, j, "j"
!allocate(vx_theta_vals(j))
!allocate(vy_theta_vals(j))
!allocate(vz_theta_vals(j))

!allocate(theta_vals(j))

!j = 1

! theta
!do i = nhalogas - 3000, nhalogas
!        if (temphalogas(1,id_gas_sort(i)) > 0.0 .AND. ABS(temphalogas(3,id_gas_sort(i))) < 10.0) then
!                vx_theta_vals(j) = (temphalogas(5,id_gas_sort(i)))
!                vy_theta_vals(j) = (temphalogas(6,id_gas_sort(i)))
!                vz_theta_vals(j) = (temphalogas(7,id_gas_sort(i)))
!                theta_vals(j) = atan(temphalogas(2,id_gas_sort(i)) / temphalogas(1,id_gas_sort(i)))
!                j = j + 1
!        endif
!
!enddo

!phi

!j = 0

!! count how many vx
!do i = nhalogas - 3000, nhalogas
!        if (temphalogas(1,id_gas_sort(i)) > 0.0 .AND. ABS(temphalogas(2,id_gas_sort(i))) < 10.0  ) then
!                j = j + 1
!        endif
!enddo

!allocate(vx_phi_vals(j))
!allocate(vy_phi_vals(j))
!allocate(vz_phi_vals(j))

!allocate(phi_vals(j))

!j = 1

! theta
!do i = nhalogas - 3000, nhalogas
!        if (temphalogas(1,id_gas_sort(i)) > 0.0 .AND. ABS(temphalogas(2,id_gas_sort(i))) < 10.0) then
!                vx_phi_vals(j) = (temphalogas(5,id_gas_sort(i)))
!                vy_phi_vals(j) = (temphalogas(6,id_gas_sort(i)))
!                vz_phi_vals(j) = (temphalogas(7,id_gas_sort(i)))
!
!                rad = sqrt((temphalogas(1,id_gas_sort(i))**2) + (temphalogas(2,id_gas_sort(i))**2) + (temphalogas(3,id_gas_sort(i))**2) )
!                phi_vals(j) = acos(temphalogas(3,id_gas_sort(i)) / rad)
!                j = j + 1
!        endif
!
!enddo


!! now do this for dark
!print *, "vx for dark"

!i_bin_val = 0
!do i = 3000, dr_bins
!       if (dr_dark_count_vals(i) .ne. 0) then
!                i_bin_val = i_bin_val + 1
!       endif
!end do

!j = 0

!! count how many vx
!do i = ndark - 9000, ndark
!        if (tempdark(1,id_dark_sort(i)) > 0.0 .AND. ABS(tempdark(3,id_dark_sort(i))) < 10.0  ) then
!                j = j + 1
!        endif
!enddo

!allocate(vx_theta_dark_vals(j))
!allocate(vy_theta_dark_vals(j))
!allocate(vz_theta_dark_vals(j))
!
!allocate(theta_dark_vals(j))

!j = 1

! theta
!do i = ndark - 9000, ndark
!        if (tempdark(1,id_dark_sort(i)) > 0.0 .AND. ABS(tempdark(3,id_dark_sort(i))) < 10.0) then
!                vx_theta_dark_vals(j) = (tempdark(5,id_dark_sort(i)))
!                vy_theta_dark_vals(j) = (tempdark(6,id_dark_sort(i)))
!                vz_theta_dark_vals(j) = (tempdark(7,id_dark_sort(i)))
!
!                theta_dark_vals(j) = atan(tempdark(2,id_dark_sort(i)) / tempdark(1,id_dark_sort(i)))
!                j = j + 1
!        endif
!
!enddo

!!phi

!j = 0

!!! count how many vx
!do i = ndark - 9000, ndark
!        if (tempdark(1,id_dark_sort(i)) > 0.0 .AND. ABS(tempdark(2,id_dark_sort(i))) < 10.0  ) then
!                j = j + 1
!        endif
!enddo

!allocate(vx_phi_dark_vals(j))
!allocate(vy_phi_dark_vals(j))
!allocate(vz_phi_dark_vals(j))


!allocate(phi_dark_vals(j))

!j = 1

!! theta
!do i = ndark - 9000, ndark
!        if (tempdark(1,id_dark_sort(i)) > 0.0 .AND. ABS(tempdark(2,id_dark_sort(i))) < 10.0) then
!                vx_phi_dark_vals(j) = (tempdark(5,id_dark_sort(i)))
!                vy_phi_dark_vals(j) = (tempdark(6,id_dark_sort(i)))
!                vz_phi_dark_vals(j) = (tempdark(7,id_dark_sort(i)))
!                rad = sqrt((tempdark(1,id_dark_sort(i))**2) + (tempdark(2,id_dark_sort(i))**2) + (tempdark(3,id_dark_sort(i))**2) )
!                phi_dark_vals(j) = acos(tempdark(3,id_dark_sort(i)) / rad)
!                j = j + 1
!        endif
!
!enddo



!! write outputs


!write(*,*) "outputting densities"


!open (unit=456, file='density_results.txt')
!do i = 1, dr_bins
!      write(456, '(3Es)') dr_vals(i), dd_dark_vals(i), dd_vals(i)
!end do
!close(456)


!write(*,*) "outputting vx theta gas"

!open (unit=456, file='vx_theta_gas.txt')
!do i = 1, size(vx_theta_vals)
!      print *, theta_vals(i), vx_theta_vals(i)
!      write(456, '(4Es)') theta_vals(i), vx_theta_vals(i),  vy_theta_vals(i), vz_theta_vals(i)
!end do
!close(456)


!write(*,*) "outputting vx phi gas"
!open (unit=456, file='vx_phi_gas.txt')
!do i = 1, size(vx_phi_vals)
!      write(456, '(4Es)') phi_vals(i), vx_phi_vals(i), vy_phi_vals(i), vz_phi_vals(i)
!end do

!close(456)


!write(*,*) "outputting vx theta dark"
!open (unit=456, file='vx_theta_dark.txt')
!do i = 1, size(vx_theta_dark_vals)
!      write(456, '(4Es)') theta_dark_vals(i), vx_theta_dark_vals(i), vy_theta_dark_vals(i), vz_theta_dark_vals(i)
!end do

!close(456)


!write(*,*) "outputting vx phi dark"
!open (unit=456, file='vx_phi_dark.txt')
!do i = 1, size(vx_phi_dark_vals)
!      write(456, '(4Es)') phi_dark_vals(i), vx_phi_dark_vals(i), vy_phi_dark_vals(i), vz_phi_dark_vals(i)
!end do

!close(456)

!deallocate(phi_dark_vals)
!deallocate(vx_phi_dark_vals)
!deallocate(theta_dark_vals)
!deallocate(vx_theta_dark_vals)
!deallocate(id_dark_sort)
!deallocate(id_gas_sort)

!! disk gas velocities

npart_out(1) = nhalogas
npart_out(2) = ndark
npart_out(3) = 0
npart_out(4) = 0
!npart_out(4) = nbulge KEEP THIS COMMENTED
npart_out(5) = 0
npart_out(6) = 0
!npart_out(6) = 0

nall_out = npart_out
ntot_out = sum(npart_out)

print *, "halo gas", nhalogas
print *, "dark", ndark

massarr_out(1) = 0.0
massarr_out(2) = 0.0
massarr_out(3) = 0.0
!massarr_out(1) = rmgashalo
!massarr_out(2) = rmdark
!massarr_out(3) = rmdiskstar
!massarr_out(4) = rmbulge ! KEEP THIS COMMENTED
massarr_out(4) = 0.0
massarr_out(5) = 0.0
massarr_out(6) = 0.0



unused_out(:) = 0

do i = 1, ntot_out
        id(i) = i
end do



! gas cal
if (nhalogas > 0) then
tempgas(:,1:nhalogas) = temphalogas(:,:)
else
tempgas(:,1:ngas) = tempothergas(:,:)
endif


pos(1:3,1:nhalogas) = temphalogas(1:3,:)
vel(1:3,1:nhalogas) = temphalogas(5:7,:)
masses(1:nhalogas) = temphalogas(4,:)
edum(1:nhalogas) = temphalogas(8,:)
dndum(1:nhalogas) = 0.0
hdum(1:nhalogas) = 0.0

print *, "testing gas", pos(1:3,ngas - 2)

print *, "size of final arrays ", SIZE(masses), SIZE(edum), " total part ", ntot_out, " totgas ", ngas

! dark
i = ngas + 1
j = ngas + ndark

pos(1:3,i:j) = tempdark(1:3,:)
vel(1:3,i:j) = tempdark(5:7,:)
masses(i:j) = tempdark(4,:)



!!! temp



!open(666, file="testdisk.txt")
!do i=(ngas + ndark + ndisc + 1),j
!   write(666,'(3Es)'), vel(1,i), vel(2,i), vel(3,i)
!end do
!close(666)

!open(667, file="testdisk.txt")
!do i=(ngas + ndark + 1),(ngas + ndark + ndisc)
!   write(667,'(6Es)'), pos(1,i), pos(2,i), pos(3,i), vel(1,i), vel(2,i), vel(3,i)
!end do
!close(666)



! black hole

!pos(1:3,ntot) = tempsink(1:3,1)
!vel(1:3,ntot) = tempsink(5:7,1)
!masses(ntot) = tempsink(4,1)
!print *, masses

! diagnostics

print *, " position - velocity - mass - temperature?()"
print *, "dark"
print *, tempdark(4,1), tempdark(4,ndark)
print *, "dark pos"
print *, minval(tempdark(1:3,:)), maxval(tempdark(1:3,:))
print *, "dark vel" 
print *, minval(tempdark(5:7,:)), maxval(tempdark(5:7,:))
print *, "dark mass"
print *, minval(tempdark(4,:)), maxval(tempdark(4,:))
print *, "total mass:", sum(tempdark(4,:))
print *, "no of particles:", ndark


if (nhalogas > 0) then
print *, " position - velocity - mass - temperature?()"
print *, "halogas"
print *, temphalogas(4,1), temphalogas(4,nhalogas)
print *, minval(temphalogas(1:3,:)), maxval(temphalogas(1:3,:)) 
print *, minval(temphalogas(5:7,:)), maxval(temphalogas(5:7,:))
print *, minval(temphalogas(4,:)), maxval(temphalogas(4,:))
print *, minval(temphalogas(8,:)), maxval(temphalogas(8,:))
print *, "total mass: ", sum(temphalogas(4,:))
print *, "no of particles:", nhalogas
endif



! put data into the arrays

print *, "-----------------"
print *, "outputs"
print *, npart_out
print *, "------------"
print *, massarr_out
print *, "------------"
print *, aexp, redshift, flagsfr, flagfeedb
print *, "----------"
print *, nall_out
print *, unused_out
print *, "----------------"


print *, "writing output"



!print *, masses

open(90, file="ics_nogal.dat", status="replace", form="unformatted")
write (90) npart_out, massarr_out, aexp, redshift, flagsfr, flagfeedb, nall_out, unused_out

write(90) pos
write(90) vel
write(90) id
write(90) masses
write(90) edum
write(90) dndum

close(90)
print *, "----------"
print *, "done"

stop
end program
