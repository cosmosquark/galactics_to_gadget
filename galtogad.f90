program galtogadget
use IFPORT

implicit none
! main program that handles galactics to gadget conversion

! load variables
include 'variables.inc'
include 'units.inc'
include 'rvarrays.inc'


real*4, allocatable :: rf_gas_arr_one (:,:)
integer*8 nhalogastemp
real*4 :: r_bin_max_val
real*8 :: halofrac_comp
real*8 :: getdr
integer*8 :: ipartthing
real*4 :: realcare
real*4 :: diskgasfrac
integer*4 :: ival
integer*4 :: nbin_rtrunc
integer*4 :: irando

real*4 :: x_dm_dummy, y_dm_dummy, z_dm_dummy
real*4 :: vx_dm_dummy, vy_dm_dummy, vz_dm_dummy
integer*4 :: dark_cut_count, idark_count

integer*4, allocatable :: gas_trunc_bin_id(:)
integer*4, allocatable :: gas_trunc_keep_id(:)
integer*4, allocatable :: dr_bin_count_cum(:)
integer*4, allocatable :: dr_bin_keep_cum(:)
integer*4, allocatable :: dr_bin_count(:)
real*4, allocatable :: dr_bin_val_arr(:)
integer*4, allocatable :: dr_bin_keep(:)
real*4, allocatable :: dr_bin_max_val(:)
real*4, allocatable :: dr_trunc_comp(:)

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
print *, "------------------------------------"
print *, "gas scale height"
print *, disk_scale_height
print *, "------------------------------------"
print *, "-------------------------------------"
print *, "running procedures now"
print*,"--------------------"

!!! new way



! disk

open(unit=43,file=discfile)
read(43,*) ndisc
print *,"Disc stars=",ndisc
do loopi=1,1
   read(43,*) realcare,realdummy,realdummy,realdummy,realdummy,realdummy,realdummy
enddo
close(43)

rmdiskstar = realcare * mconffactor

diskgasfrac = 0.1

rmdiskstar = (1  - diskgasfrac) * realcare * mconffactor
ndiscgas = ndisc
rmdiskgas = realcare * diskgasfrac * mconffactor

print *, "DISK STAR PARTICLES: ", ndisc
print *, "DISK STAR MASS: ", rmdiskstar
print *, "Total stellar disk Mass: ", rmdiskstar * ndisc

print *, "DISK GAS PARTICLES: ", ndiscgas
print *, "DISK GAS MASS: ", rmdiskgas
print *, "Total gas disk Mass: ", rmdiskgas * ndiscgas



! 1) load nparticles and masses

! halogas

eminin = 1.0e4 / Kunit
ehalo = 1.0e6 / Kunit

! mgas = 0.5 dm mass

! realcare is the origional galactics mass
Mm =  0.02 * Mdisk 

!! need to scale for this mass lost.
mgas = rmdiskgas * 1d10

if (halofrac .gt. 0.000) then

rmgashalo = rmdiskgas

!! call makehalo the first time
call makehalo(nhalogastemp, rc, rm, rf, Mm, Mtotal, mgas, alpha,rho0, ipr)
!!call makehalo(nhalogas, rc, rm, rf, Mm, Mtotal, mgas, alpha,rho0, ipr)
!call calculate_h0_fun


! outer radius (kpc)
!real*8, parameter :: rf = 275.0

! truncation width (kpc)
!real*8, parameter :: rf_trunc = 25.0

!real*, parameter :: rf_gas = rf + rf_trunc


print *, "generating gas truncation"
print *, "reading file"
!! load the gas particles
open(unit=43,file='halogas_positions.txt')

read(43, *) nhalogastemp, Mtotal
print *,"nhalogas=",nhalogastemp
allocate(temphalogasinit(1:8,1:nhalogastemp))
do i=1,nhalogastemp
   read(43,*) idummy, temphalogasinit(1:3,i),realdummy

enddo
close(43)
temphalogasinit(4,:) = rmgashalo
temphalogasinit(5,:) = 0.0
temphalogasinit(6,:) = 0.0
temphalogasinit(7,:) = 0.0

!! need to give these ID's
print *, "sorting particles"
allocate(r_gas_sort(nhalogastemp))
allocate(id_gas_sort(nhalogastemp))
do i = 1, nhalogastemp
   r_gas_sort(i) = sqrt((temphalogasinit(1,i) ** 2.0) + (temphalogasinit(2,i) ** 2.0) + (temphalogasinit(3,i) ** 2.0))
   id_gas_sort(i) = i
enddo

print *, "calling sort"
call SSORTID(r_gas_sort, id_gas_sort,nhalogastemp,2)
print *, "ending sort"

!! sort it, will save a lot of hastle
allocate(rf_gas_arr_one(1:8,1:nhalogastemp))

print *, "sorting array"
do i = 1, nhalogastemp
    rf_gas_arr_one(1:8,i) = temphalogasinit(:,id_gas_sort(i))
enddo

!! initial array is now sorted
temphalogasinit(:,:) = rf_gas_arr_one(1:8,:)


!! now r_gas_sort maps to temphaloinit

deallocate(rf_gas_arr_one)

!! lets use temphalogasinit(8,:) as a proxy for removed star particle (give it a
!value of -999 for remove
! value of -888 for keep.
!! temphalogasinit(7,:) for the radial bin

!! compute max_r

print *, "figuring out number of bins"
r_temp = rf_trunc_start 
realdummy = maxval(r_gas_sort(:))
ival = 1
do while (r_temp < realdummy)
   ival = ival +  1
   r_temp = r_temp + rf_gas_width
enddo

nbin_rtrunc = ival

allocate(dr_bin_val_arr(nbin_rtrunc))
allocate(dr_trunc_comp(nbin_rtrunc))
allocate(dr_bin_count(nbin_rtrunc))
allocate(dr_bin_count_cum(nbin_rtrunc))
allocate(dr_bin_max_val(nbin_rtrunc))
allocate(dr_bin_keep(nbin_rtrunc))
allocate(dr_bin_keep_cum(nbin_rtrunc))


print *, "populating bins"
!! the r bins
do i = 1, nbin_rtrunc
   dr_bin_val_arr(i) = rf_trunc_start + (i *  rf_gas_width) - (rf_gas_width / 2.0)
   realdummy = (dr_bin_val_arr(i) - rf) / (sqrt(2.0) * rf_trunc)
   dr_trunc_comp(i) = 0.5 * erfc(realdummy)
   dr_bin_max_val(i) = rf_trunc_start + (i *  rf_gas_width)
   !! creating a count array
   dr_bin_count(i) = 0 
enddo


!! assign values and check

allocate(gas_trunc_bin_id(nhalogastemp))
allocate(gas_trunc_keep_id(nhalogastemp))

!! ival here to record how many gas particles we are missing in the loop
!intiially
!! j to map to the binspace
print *, "counting particles in bins"
j = 1
ival = 0
do i = 1, nhalogastemp
   !! allocate what gas particles to keep.. set to 1 by default
   gas_trunc_keep_id(i) = 1
   if ( r_gas_sort(i) < rf_trunc_start) then
       !! in space of the origional format
       gas_trunc_bin_id(i) = 0
       !! number of values before we start truncating
       ival = ival + 1
   else
       !! allocate the bin
       if (r_gas_sort(i) < dr_bin_max_val(j)) then
          gas_trunc_bin_id(i) = j
          dr_bin_count(j) = dr_bin_count(j) + 1
       else
          j = j + 1
          gas_trunc_bin_id(i) = j
          dr_bin_count(j) = dr_bin_count(j) + 1
       endif
   endif
enddo

!! now we know how many values are in each bin.

!! need to make the cumulatitive array, left aggrogated

nhalogas = ival
print *, "figuring out how many to keep"

do i = 1, nbin_rtrunc
   realdummy = dr_trunc_comp(i) * float(dr_bin_count(i))
   if (dr_trunc_comp(i) > 0.5) then
       dr_bin_keep(i) = ceiling(realdummy)
   else
       dr_bin_keep(i) = floor(realdummy)
   endif
   nhalogas = nhalogas + dr_bin_keep(i)
enddo 


do i = 1, nbin_rtrunc
  if (i == 1) then
    dr_bin_count_cum(i) = ival
    dr_bin_keep_cum(i) = ival
  else
    dr_bin_count_cum(i) = dr_bin_count(i-1) + dr_bin_count_cum(i-1)
    dr_bin_keep_cum(i) = dr_bin_keep(i-1) + dr_bin_keep_cum(i-1)
  endif
enddo


!! now we have sorted the whole thing out and we know how many particles we are
!keeping

!! now time to filter everything

!! start with the easy bit

allocate(temphalogas(1:8,1:nhalogas))

do i = 1, ival
  temphalogas(:,i) = temphalogasinit(:,i)
enddo


!! now bin 1
print *, "truncating"
do i = 1, nbin_rtrunc
!  need the plus 1 to prevent overflowing etc.. and to be able to push things
!  down
  allocate(rf_gas_arr_one(1:8,1:dr_bin_count(i)+1))
!
  rf_gas_arr_one(1:8,:) = temphalogasinit(:,dr_bin_count_cum(i)+1:dr_bin_count_cum(i)+dr_bin_count(i))
!!  !! to make sure we don't break the arrays
  rf_gas_arr_one(1:8,dr_bin_count(i)+1) = 0.0 
!
!  !! pop a value from the array
!
  idummy = dr_bin_count(i)

  if (idummy > dr_bin_keep(i)) then
     do while (idummy > dr_bin_keep(i))
        irando = int(rand()*idummy)+1
        !! shift everything down one from irand.. effectively overwriting irand
        rf_gas_arr_one(1:8,irando:dr_bin_count(i)) = rf_gas_arr_one(1:8,irando+1:dr_bin_count(i)+1)
        idummy = idummy - 1
     enddo
   endif


  !! move the temp array into the main one
  temphalogas(:,dr_bin_keep_cum(i)+1:dr_bin_keep_cum(i)+dr_bin_keep(i)) = rf_gas_arr_one(:,:dr_bin_keep(i)) 
  deallocate(rf_gas_arr_one)
enddo
!enddo


print *, "truncation complete"

open(unit=106, file="halogas_sorted_properties.txt")
do i = 1, nhalogas
!! dr, temperature, gas dens, dark dens, enc mass
 write(106,'(4Es)') temphalogas(1,i), temphalogas(2,i), temphalogas(3,i), sqrt((temphalogas(1,i))**2.0 + (temphalogas(2,i))**2.0 + (temphalogas(3,i))**2.0)
enddo
close(106)

open(unit=106, file="halogastemp_sorted_properties.txt")
do i = 1, nhalogastemp
!! dr, temperature, gas dens, dark dens, enc mass
 write(106,'(4Es)') temphalogasinit(1,i), temphalogasinit(2,i), temphalogasinit(3,i), sqrt((temphalogasinit(1,i))**2.0 + (temphalogasinit(2,i))**2.0 + (temphalogasinit(3,i))**2.0)
enddo
close(106)



!! move things to the temp bin
deallocate(id_gas_sort)
deallocate(r_gas_sort)
deallocate(gas_trunc_bin_id)
deallocate(gas_trunc_keep_id)
deallocate(dr_bin_max_val)
deallocate(dr_bin_count)
deallocate(dr_bin_keep)
deallocate(dr_bin_count_cum)
deallocate(dr_bin_keep_cum)
deallocate(temphalogasinit)
deallocate(dr_bin_val_arr)

ival = 0

!!!! end here


!! logic being
!! 1) run make halo
!! 2) truncate particles
!! 3) add the mass loss due to truncation to the rest of the particles

!! apply the truncation




!! realcare is the raw mass



!open(unit=43,file='halogas_positions.txt')
!read(43, *) nhalogas, Mtotal
!close(43)
!rmgashalo = Mtotal / nhalogas / 1d10

print *, "HOT HALO GAS PARTICLES: ", nhalogas
print *, "GAS HALO MASS: ", rmgashalo
print *, "Total Halo Mass: ", rmgashalo * nhalogas


!!!! Need to do more with this.

!! to get dark mass, need to conver % of DM mass into halo gas mass

! dark
open(unit=43,file=darkhalofile)
read(43,*) ndark
do loopi=1,2
   read(43,*) realcare,realdummy,realdummy,realdummy,realdummy,realdummy,realdummy
enddo
close(43)

!! count the number of dm particles that are within 275 kpc
dark_cut_count = 0
open(unit=43,file=darkhalofile)
read(43,*) ndark
do loopi=1,ndark
   read(43,*) realcare,x_dm_dummy,y_dm_dummy,z_dm_dummy,realdummy,realdummy,realdummy
   if (sqrt((x_dm_dummy**2.0) + (y_dm_dummy**2.0) + (z_dm_dummy**2.0)) .le. rf) then
      dark_cut_count = dark_cut_count + 1
   endif
enddo
close(43)





! fraction between the two total  masses. need this to know how much fraction of
! dark matter particle mass is lost to gas mass
print *, "initial dm mass", realcare
realdummy = rmgashalo * nhalogas / (realcare * mconffactor * dark_cut_count)
! realdummy here is the fraction of total gas halo mass, to total dark matter halo mass
print *, "realdummy", realdummy
! we rescale the mass of dark matter by this mass scale
rmdark = (realcare * mconffactor) * (1.0 - realdummy)


!! keep the variable here for now (in dark_cut_count as the true count of dark
!matter particles)
print *, "DARK MATTER PARTICLES: ", dark_cut_count
print *, "DARK MATTER MASS: ", rmdark
print *, "Total DM Mass: ", rmdark * dark_cut_count


!! uncomment this
!deallocate(r_gas_sort(:))
!deallocate(id_gas_sort(:))

else
nhalogas = 0
endif


! bulge

open(unit=43,file=bulgefile)
read(43,*) nbulge
print *,"Bulge stars=",nbulge
do loopi=1,2
   read(43,*) realcare,realdummy,realdummy,realdummy,realdummy,realdummy,realdummy
enddo
close(43)

rmbulge = realcare * mconffactor
print *, "BULGE PARTICLES: ", nbulge
print *, "BULGE MASS: ", rmbulge
print *, "Total bulge Mass: ", rmbulge * nbulge

open(unit=43,file=discfile)
read(43,*) ndisc
close(43)


ngas = ndiscgas + nhalogas
nobj = ngas + dark_cut_count + nbulge + ndisc

! 2) allocate arrays

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
print *, dark_cut_count, "ndark"
print *, "allocating arrays"
print *, nobj
print *, "ngas=", ngas
print *, "ndisc=", ndisc
print *, "ndiscgas=", ndiscgas
allocate(tempgas(1:8,1:ngas))
print *, "gas"

!if (nhalogas > 0) then
!allocate(temphalogas(1:8,1:nhalogas))
!endif
print *, "halogas"
allocate(tempothergas(1:8,1:ndiscgas))
print *, "othergas"
allocate(tempdark(1:7,1:dark_cut_count))
print *, "dark"
allocate(tempbulge(1:7,1:nbulge))
print *, "bulge"
allocate(tempdisk(1:7,1:ndisc))
print *, "disk"
!allocate(tempstar(1:7,1:nstar))
!print *, "star"

! commenting out
!allocate(tempsink(1:8,1))
!print *, "sink"

print *, "allocating outputs"

allocate(pos(1:3,1:nobj))
allocate(vel(1:3,1:nobj))
allocate(id(1:nobj))
allocate(masses(1:nobj))
allocate(edum(1:ngas))
allocate(dndum(1:ngas))
allocate(hdum(1:ngas))
!allocate(edum(1:nhalogas))
!allocate(dndum(1:nhalogas))
!allocate(hdum(1:nhalogas))

! 2.5) poplate arrays

!!! not needed here
!open(unit=43,file='halogas_positions.txt')

!read(43, *) nhalogas, Mtotal
!print *,"nhalogas=",nhalogas
!do i=1,nhalogas
!   read(43,*) idummy, temphalogas(1:3,i),realdummy

!enddo
!close(43)
temphalogas(4,:) = rmgashalo
temphalogas(5,:) = 0.0
temphalogas(6,:) = 0.0
temphalogas(7,:) = 0.0



! temperature
!! new temp method
!b_htemp = 1587182.80415180
!c_htemp = 8900.00000000000
!a_htemp = 28.275

!do i = 1, nhalogas
!        r_htemp = sqrt((temphalogas(1,i) ** 2.0) + (temphalogas(2,i) ** 2.0) + (temphalogas(3,i) ** 2.0))
!        tanh_x_htemp = r_htemp / a_htemp
!        sech_x_htemp = (r_htemp / a_htemp) - 1.0

       ! temphalogas(8,i) = (b_htemp * tanh(tanh_x_htemp)) + (b_htemp * (1.0 / cosh(sech_x_htemp))) - (c_htemp* r_htemp) - 200000.0
!        temphalogas(8,i) = max((b_htemp * tanh(tanh_x_htemp)) + (b_htemp * (1.0 /cosh(sech_x_htemp))) - (c_htemp* r_htemp),334779.0)
!        temphalogas(8,i) = temphalogas(8,i) / kunit

!enddo
!endif

! dark

idark_count = 0
open(unit=43,file=darkhalofile)
read(43,*) ndark !! dummy
do i=1,ndark
    read(43,*) realdummy, x_dm_dummy, y_dm_dummy, z_dm_dummy, vx_dm_dummy, vy_dm_dummy, vz_dm_dummy
      if (sqrt((x_dm_dummy**2.0) + (y_dm_dummy**2.0) + (z_dm_dummy**2.0)) .le. rf) then
         idark_count = idark_count + 1
         tempdark(1,idark_count) = x_dm_dummy
         tempdark(2,idark_count) = y_dm_dummy
         tempdark(3,idark_count) = z_dm_dummy
         tempdark(5,idark_count) = vx_dm_dummy * vconffactor
         tempdark(6,idark_count) = vy_dm_dummy * vconffactor
         tempdark(7,idark_count) = vz_dm_dummy * vconffactor
         tempdark(4,idark_count) = rmdark
      endif
enddo
close(43)

ndark = idark_count

print *, "DARK CHECK"
print *, rmdark
print *, rmdark * ndark
print *, sum(tempdark(4,:))
print *, "are these the same?"

!print *, "dr count"
!print *, dr_count_vals
!print *, dd_vals


! bulge stars
open(unit=43,file=bulgefile)
read(43,*) nbulge
print *, "nbulge=",nbulge
do i=1,nbulge
   read(43,*) realdummy, tempbulge(1:3,i),tempbulge(5:7,i)
   tempbulge(5:7,i) = tempbulge(5:7,i) * vconffactor     ! 100 km/s to internal   
   tempbulge(4,i) = rmbulge

enddo
close(43)


! disk stars
open(unit=43,file=discfile)
read(43,*) ndisc
print *, "ndisk=",ndisc
do i=1,ndisc
   read(43,*) realdummy, tempdisk(1:3,i),tempdisk(5:7,i)
   darkmid   = darkmid     + r(1:3,i)
   tempdisk(5:7,i) = tempdisk(5:7,i) * vconffactor     ! 100 km/s to internal   
   tempdisk(4,i) = rmdiskstar
enddo
close(43)

! disk gas

print *, "modifying height by a scale height"
print *, disk_scale_height

do i=1,ndisc
   tempothergas(1,i) = - tempdisk(1,i)
   tempothergas(2,i) = - tempdisk(2,i)
   tempothergas(3,i) = tempdisk(3,i) / (1.0 / disk_scale_height )
   tempothergas(4,i) = rmdiskgas
   tempothergas(5,i) = - tempdisk(5,i)
   tempothergas(6,i) = - tempdisk(6,i)
   tempothergas(7,i) = 0.0
   tempothergas(8,i) = eminin

!   rad = sqrt((tempothergas(1,i))**2+(tempothergas(2,i))**2)
!   jr  = 1
!   do while ( jr<npot .and. rad>vcrad(jr) )
!       jr = jr + 1
!   end do
!   if ( jr>=npot ) then
!      print *,"Disc particle outside of allowed range from circular velocity table"
!      stop
!   end if

!   if ( jr==1 ) then
!      vcirc = rad/vcrad(1) * vcvc(1)
!   else
                        !! gotta convert the rad factor from kpc to km
!                        vcirc
!                        =((rad-vcrad(jr-1)))*(vcvc(jr)-vcvc(jr-1))+vcvc(jr-1)
                        !! would this not make more sense?
!       vcirc=((rad/vcrad(jr-1)))*(vcvc(jr)-vcvc(jr-1))+vcvc(jr-1)
!   end if

                ! Apply Plummer softening
!   vcirc  = vcirc * 1./(1.+(rsoft/rad)**2)**1.5

!   tempothergas(5,i) = -(tempothergas(2,i)) * vcirc/rad
!   tempothergas(6,i) =  (tempothergas(1,i)) * vcirc/rad
!   tempothergas(7,i) = 0

enddo





! 4) write arrays

!! finally output

npart_out(1) = ngas
npart_out(2) = ndark
npart_out(3) = ndisc
npart_out(4) = nbulge
!npart_out(4) = nbulge KEEP THIS COMMENTED
npart_out(5) = 0
npart_out(6) = 0
!npart_out(6) = 0

nall_out = npart_out
ntot_out = sum(npart_out)

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
tempgas(:,nhalogas+1:ngas) = tempothergas(:,:)
else
tempgas(:,1:ngas) = tempothergas(:,:)
endif


pos(1:3,1:ngas) = tempgas(1:3,:)
vel(1:3,1:ngas) = tempgas(5:7,:)
masses(1:ngas) = tempgas(4,:)
edum(1:ngas) = tempgas(8,:)
dndum(1:ngas) = 0.0
hdum(1:ngas) = 0.0


print *, "size of final arrays ", SIZE(masses), SIZE(edum), " total part ", ntot_out, " totgas ", ngas

! dark
i = ngas + 1
j = ngas + ndark

pos(1:3,i:j) = tempdark(1:3,:)
vel(1:3,i:j) = tempdark(5:7,:)
masses(i:j) = tempdark(4,:)

print *, vel(1,j), vel(2,j), vel(3,j), "TEST THREEE"


i = ngas+ndark + 1
j = ngas + ndark + ndisc
! disk
pos(1:3,i:j) = tempdisk(1:3,:)
vel(1:3,i:j) = tempdisk(5:7,:)
masses(i:j) = tempdisk(4,:)


! bulge
i = ngas + ndark + ndisc + 1
j = ngas + ndark + ndisc + nbulge

!print *, "sanity check , last value is = ", j, " and there are ", ntot, " particles"
pos(1:3,i:j) = tempbulge(1:3,:)
vel(1:3,i:j) = tempbulge(5:7,:)
!print *, "bulge vel test"
!print *, vel(1:3,i+2)
masses(i:j) = tempbulge(4,:)


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
print *, "total mass: ", sum(temphalogas(4,:))
print *, "no of particles:", nhalogas
endif

print *, " position - velocity - mass - temperature?()"
print *, "bulge"
print *, tempbulge(4,1), tempbulge(4,nbulge)
print *, minval(tempbulge(1:3,:)), maxval(tempbulge(1:3,:)) 
print *, minval(tempbulge(5:7,:)), maxval(tempbulge(5:7,:))
print *, minval(tempbulge(4,:)), maxval(tempbulge(4,:))
print *, "total mass: ", sum(tempbulge(4,:))
print *, "no of particles:", nbulge

print *, " position - velocity - mass - temperature?()"
print *, "disk"
print *, tempdisk(4,1), tempdisk(4,ndisc)
print *, minval(tempdisk(1:3,:)), maxval(tempdisk(1:3,:)) 
print *, minval(tempdisk(5:7,:)), maxval(tempdisk(5:7,:))
print *, minval(tempdisk(4,:)), maxval(tempdisk(4,:))
print *, "total mass: ", sum(tempdisk(4,:))
print *, "no of particles:", ndisc
print *, "diskgas"

print *, tempothergas(4,1), tempothergas(4,ndiscgas), minval(tempothergas(4,:))
print *, minval(tempothergas(1:3,:)), maxval(tempothergas(1:3,:)) 
print *, minval(tempothergas(5:7,:)), maxval(tempothergas(5:7,:))
print *, minval(tempothergas(4,:)), maxval(tempothergas(4,:))
print *, minval(tempothergas(8,:)), maxval(tempothergas(8,:))
print *, "total mass: ", sum(tempothergas(4,:))
print *, "no of particles:", ndiscgas




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

open(90, file="ics_out.dat", status="replace", form="unformatted")
write (90) npart_out, massarr_out, aexp, redshift, flagsfr, flagfeedb, nall_out, unused_out

write(90) pos
write(90) vel
write(90) id
write(90) masses




!print *, minval(temphalogas(8,:)), maxval(temphalogas(8,:))
!deallocate(tempothergas)
!deallocate(tempdark)
!deallocate(tempbulge)
!deallocate(tempdisk)
!deallocate(vel)
!


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
   r_sort(i) = sqrt((pos(1,i) ** 2.0) + (pos(2,i) ** 2.0) + (pos(3,i) ** 2.0))
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


print *, "r_sort"
print *, r_gas_sort(1:100)

r_max_val = maxval(r_sort(:))
print *, r_max_val
dr_bin_width = r_max_val / dr_bins

r_min_val = minval(r_sort(:))

! create the dr bins
print *, "creating dr bins"
do i = 1, dr_bins
   dr_vals(i) = (i * dr_bin_width)
enddo

! sort both r arrays
dm_vals(:) = 0.0
dd_vals(:) = 0.0
dtemp_vals(:) = 0.0
dr_count_vals(:) = 0

!dr_dark_count_vals(:) = 0
!dd_dark_vals(:) = 0.0

print *, "dr vals"
print *, dr_vals(1:100)


print *, nhalogas
print *, "sorting array 1, no need"
!call SSORTID(r_gas_sort, id_gas_sort,nhalogas,2)
print *, "sorting array 2"
call SSORTID(r_sort,id_sort,nobj,2)

call SSORTID(r_dark_sort,id_dark_sort,ndark,2)

print *, "post sort"
print *, r_gas_sort(1:100)

!! fill mass dr array
i_bin_val = 1
dr_bin_val = dr_bin_width


!!!!! NB, this should be looping in reverse
!!!!! since we want the mass  exterior to that radius
print *, "computing TOTAAL mass in bins"
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

    dm_vals(i_bin_val) = dm_vals(i_bin_val) + masses(id_sort(i))
   ! print *, dm_vals(i_bin_val), masses(id_sort(i)), id_sort(i), i
enddo

print *, "dm_vals"
print *, r_sort(1:100)
print *, dm_vals(1:100)

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
   m_total_check = m_total_check + masses(i)
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

print *, "dr count"
print *, dr_count_vals(1:100)
print *, dr_vals(1:100)
print *, "masses"
print *, dm_vals(1:20), dm_vals(2000:2020), maxval(dm_vals)

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
!          ! need to do this more than once xD
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

print *, "dr dark count"
print *, dr_dark_count_vals(1:100)
print *, dd_dark_vals(1:100)
print *, "dr vals"
print *, dr_vals(1:100)
print *, "dd vals"
print *, dd_vals(1:100)
print *, pi_val, lunit, munit

!print *, "dr count dark"
!print *, dr_dark_count_vals
!print *, dd_dark_vals

!1 convert dd_vals into gas densities

do i = 1, dr_bins
    if (i == 1) then
       vol_area = (4.0/3.0) * pi_val * (dr_vals(i) ** 3.0) * (lunit ** 3)
    else
       vol_area =  (4.0/3.0) * pi_val * ((((dr_vals(i) ** 3.0) *  lunit ** 3.0)) - ((dr_vals(i-1) ** 3.0) * lunit ** 3.0))
    end if
!    print *, vol_area
    dd_vals(i) = dd_vals(i) * munit / vol_area
    dd_dark_vals(i) = dd_dark_vals(i) * munit / vol_area
end do

!print *, "dd final"
!print *, dd_vals


print *, dd_vals(1:100)

!! finally do the integral
print *, "integrating"
do i = 1, dr_bins
   !! now remember, the approximation to an integral is a sumation of width dr
   if (dr_count_vals(i) .ne. 0) then
       do j = i, dr_bins
           dtemp_vals(i) = dtemp_vals(i) + ((((dd_vals(j) * dm_vals(j) / (dr_vals(j) ** 2.0)) * dr_bin_width)) * (munit / lunit))
       end do

       !! integral prefactor last
       dtemp_vals(i) = max(eminin,dtemp_vals(i) * (tempconst / dd_vals(i) ) / kunit)
   else
       dtemp_vals(i) = 0.0
   endif

end do

open(unit=106, file="halogas_radial_properties.txt")
do ival = 1, dr_bins
!! dr, temperature, gas dens, dark dens, enc mass
 write(106,'(5Es, 2I)') dr_vals(ival), dtemp_vals(ival) * kunit,  dd_vals(ival), dd_dark_vals(ival), dm_vals(ival), dr_count_vals(ival), dr_dark_count_vals(ival)
enddo
close(106)

!! now assign temperatures to particles
print *, "assigning temp to particles"
i_bin_val = 1
dr_bin_val = dr_bin_width

i_bin_val = 1
do i = 1, dr_bins
   if (dr_count_vals(i) .ne. 0) then
      do j = 1, dr_count_vals(i)
          !! can't really have the temperature colder than the disk
          temphalogas(8,id_gas_sort(i_bin_val)) = dtemp_vals(i)
          i_bin_val = i_bin_val + 1
      end do
   endif
end do

print *, "halo temp"
print *, temphalogas(8,1:100)


print *, "halo temp ranges"
print *, minval(temphalogas(8,:)), maxval(temphalogas(8,:))

edum(1:nhalogas) = temphalogas(8,:)


write(90) edum
write(90) dndum

close(90)
print *, "----------"
print *, "done"

stop
end program
