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
real*4 :: realcare
real*4 :: diskgasfrac
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

mgas = rmdiskgas * 1d10

if (halofrac .gt. 0.000) then

call makehalo(nhalogas, rc, rm, rf, Mm, Mtotal, mgas, alpha,rho0, ipr)
!call calculate_h0_fun



!! realcare is the raw mass
!rmdark = realcare * (1.0-halofrac) * mconffactor
!rmgashalo = realcare * halofrac * mconffactor



open(unit=43,file='halogas_positions.txt')
read(43, *) nhalogas, Mtotal
close(43)

! get back into correct units
rmgashalo = Mtotal / nhalogas / 1d10

print *, "HOT HALO GAS PARTICLES: ", nhalogas
print *, "GAS HALO MASS: ", rmgashalo
print *, "Total Halo Mass: ", rmgashalo * nhalogas

!! to get dark mass, need to conver % of DM mass into halo gas mass

! dark
open(unit=43,file=darkhalofile)
read(43,*) ndark
do loopi=1,2
   read(43,*) realcare,realdummy,realdummy,realdummy,realdummy,realdummy,realdummy
enddo
close(43)

! fraction between the two particle masses
print *, "initial dm mass", realcare
realdummy = rmgashalo / (realcare * mconffactor)
print *, "realdummy", realdummy
rmdark = (realcare * mconffactor) * (1.0 - realdummy)


print *, "DARK MATTER PARTICLES: ", ndark
print *, "DARK MATTER MASS: ", rmdark
print *, "Total DM Mass: ", rmdark * ndark



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
nobj = ngas + ndark + nbulge + ndisc

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
print *, ndark, "ndark"
print *, "allocating arrays"
print *, nobj
print *, "ngas=", ngas
print *, "ndisc=", ndisc
print *, "ndiscgas=", ndiscgas
allocate(tempgas(1:8,1:ngas))
print *, "gas"

if (nhalogas > 0) then
allocate(temphalogas(1:8,1:nhalogas))
endif
print *, "halogas"
allocate(tempothergas(1:8,1:ndiscgas))
print *, "othergas"
allocate(tempdark(1:7,1:ndark))
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

if (nhalogas > 0) then

open(unit=43,file='halogas_positions.txt')

read(43, *) nhalogas, Mtotal
print *,"nhalogas=",nhalogas
do i=1,nhalogas
   read(43,*) idummy, temphalogas(1:3,i),realdummy

enddo
close(43)
temphalogas(4,:) = rmgashalo
temphalogas(5,:) = 0.0
temphalogas(6,:) = 0.0
temphalogas(7,:) = 0.0


! temperature
!! new temp method
b_htemp = 1587182.80415180
c_htemp = 8900.00000000000
a_htemp = 28.275

do i = 1, nhalogas
        r_htemp = sqrt((temphalogas(1,i) ** 2.0) + (temphalogas(2,i) ** 2.0) + (temphalogas(3,i) ** 2.0))
        tanh_x_htemp = r_htemp / a_htemp
        sech_x_htemp = (r_htemp / a_htemp) - 1.0

       ! temphalogas(8,i) = (b_htemp * tanh(tanh_x_htemp)) + (b_htemp * (1.0 / cosh(sech_x_htemp))) - (c_htemp* r_htemp) - 200000.0
        temphalogas(8,i) = max((b_htemp * tanh(tanh_x_htemp)) + (b_htemp * (1.0 /cosh(sech_x_htemp))) - (c_htemp* r_htemp),334779.0)
        temphalogas(8,i) = temphalogas(8,i) / kunit

enddo
endif

! dark

open(unit=43,file=darkhalofile)
read(43,*) ndark
print *,"ndark=",ndark
do i=1,ndark
   read(43,*) realdummy, tempdark(1:3,i),tempdark(5:7,i)
   darkmid   = darkmid     + tempdark(1:3,i)
   tempdark(5:7,i) = tempdark(5:7,i) * vconffactor     ! 100 km/s to internal   
   tempdark(4,i) = rmdark

enddo
close(43)


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
print *, minval(temphalogas(8,:)), maxval(temphalogas(8,:))
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
write(90) edum
write(90) dndum

close(90)
print *, "----------"
print *, "done"

stop
end program
