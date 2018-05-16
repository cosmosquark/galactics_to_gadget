!     This routine runs through a list of particles and, for every
!     particle, finds the radius which encompases nh_min & nh_max
!     neighbours
!     There are many text output files that document
!     necessary/interesting results
!     This is a revised version of findhagain2.f, which is a parallel
!     (openmp) version of 'findhagain1.f'
!     Note: runtimes depend on values of nlim & dlnlim
!     Note: this will be looped through once or twice, depending on 
!           existence of gas halo
!     Compile as > ifort calculate_h0.f ssort.f -openmp
! 
      subroutine calculate_h0_fun
      implicit none
! 
      integer :: i            , j            , k            , Nmax, &
            ngas         , nh_min       , nh_max, &
            nlim0        , nlimn        , nlimx        , dnlim, &
            imin         , imax         , imax1        , imax2, &
            itmp         , jm1          , jp1, &
            count0       , count1       , ctrate0      , ctrate1, &
            ynr3         , hyn          , comploop     , idummy
      parameter(Nmax=1000000)
      real*8 ::  rdisk        , ndummy       , h0           , R0, &
            hmin0        , hmax0        , rmin         , rmax, &
            htmp         , hval         , rval         , timesec, &
            r1_11        , r2_11        , r3_11, &
            r1_16        , r2_16        , r3_16        , tol
      real*8 ::  r1  (Nmax)   , r1b  (Nmax)  , r1c  (Nmax), &
            r2  (Nmax)   , r3   (Nmax)  , r    (Nmax)  , rb   (Nmax), &
            hmin(Nmax)   , hmax (Nmax)  , htmpn(Nmax)  , htmpx(Nmax), &
            dr2 (Nmax)   , rnum (Nmax)  , rnum2(Nmax)
      character :: halo_yn*1
      character :: cdummy*128 , discfile*128 , gashalofile*128
! ----------------------------------------------------------------------
!     To time the routine
      call system_clock(count0, ctrate0)
      write(*,1100)
      write(*,1000) 
! 
! ----------------------------------------------------------------------
!     Get the filenames, and determine how many times we need to do this
!
      hyn = 2
!
!     Get relevant filenames
      discfile = "disk"
      gashalofile = "halogas_positions.txt"
! ----------------------------------------------------------------------
!     Initial parameters
!     The former two are the min and max h's desired
!     The latter two are the initial distance and extended distances to
!     check for neighbouring particles
      nh_min =  30
      nh_max =  80
      nlim0  = 5000
      dnlim  = nlim0
!     The radius of the disk as input into GalactICs
      rdisk  = 17.29d0
! 
      do comploop = 1,hyn
!       Open files
        if (comploop.eq.1) then 
          print *, "h0 for disk"
          open(unit=1, file=discfile) !! opening
          open(unit=6, file="find_rh.txt") !! writing to
        else
          print *, "h0 for gas"
          open(unit=1, file=gashalofile)  !! opening
          open(unit=6, file="find_rh_halo.txt") !! writing to
        end if
        !! writing to the rest
        open(unit=3, file="find_hvals.txt")        
        open(unit=4, file="find_log.txt")
        open(unit=5, file="find_hvals2.txt")
        open(unit=8, file="find_original_r.txt")
!
!       Read in the number of particles
        read(1,*) ngas, ndummy
        if (comploop.eq.1) then 
          write(*,3016) ngas
        else
          write(*,3017) ngas
        end if
!
!       Read stuff into internal arrays
!       Note that arrays get rearranged in ssort, thus we need the
!       unsorted r1b to sort r3 after r1 gets rearranged after sorting r2
        do i=1,ngas
          if (hyn.eq.1) then
            read (1,   *) ndummy, r1(i),r2(i),r3(i), ndummy,ndummy,ndummy
            write(8,3003)         r1(i),r2(i),r3(i)
          else
            read(1,    *) idummy, r1(i),r2(i),r3(i), ndummy
            write(8,3003)         r1(i),r2(i),r3(i)
          end if
          r1b  (i) = r1(i)
          r1c  (i) = r1(i)
          rnum (i) = real(i)
          rnum2(i) = rnum(i)
        end do
        close(1)
        close(8)
 
        write(*,1001)
!       Sorts arrays via ssort; done in duplicate to also sort r2 & r3
        call ssort(r1 , r2  , ngas, 2)
        call ssort(r1b, r3  , ngas, 2)
        call ssort(r1c, rnum, ngas, 2)
        write(*,1002)
        print *, "begin paralellisation"
!
!!!$OMP PARALLELDO PRIVATE(i, j, k, jm1, jp1,nlimn, nlimx, dr2) & 
!!!$OMP & SHARED(ngas,nh_min,nh_max,r,rb,r1,r2,r3,nlim0,dnlim,hmin,hmax)
        do i=1,ngas
!         Finding the particle position for later (array in duplicate)
          r (i) = sqrt( r1(i)**2 + r2(i)**2 + r3(i)**2 )
          rb(i) = r(i)
!         reinitialising position from particle
          nlimn = nlim0
          nlimx = nlim0
100       continue
          k = 0
!         This top loop works, but faster with the second one + the 
!         if-statements after the 'call ssort'
!         do j=1,ngas 
          do j=max(1,i-nlimn), min(i+nlimx,ngas)
!           Find the distance to the nlim nearest particles in 
!           the r1-direction
            k      = k + 1
            dr2(k) = ( r1(i)-r1(j) )**2 & 
                   + ( r2(i)-r2(j) )**2 &
                   + ( r3(i)-r3(j) )**2
!           Since we don't want to include the primary particle in its
!           list of neighbours, reset the 'separation' to be really big
            if (i.eq.j) dr2(k) = 1.0d30
          end do
          call ssort(dr2, dr2, k, 1)
!         Confirming that the nearest particles have been mapped
!         (these two if-statements are the same, just in opposite order)
          jm1 = max(1        ,i-nlimn-1)
          jp1 = min(i+nlimx+1,ngas     )
          if ( ( (r1(i)-r1(jm1))**2 .lt. dr2(nh_max) ) .and.  &
                (jm1.ne.1) )                           then
            nlimn = nlimn + dnlim
            if (  ( (r1(i)-r1(jp1))**2 .lt. dr2(nh_max) ) .and. &
                  (jp1.ne.ngas) ) nlimx = nlimx + dnlim

            goto 100
          end if
          if (  ( (r1(i)-r1(jp1))**2 .lt. dr2(nh_max) ) .and. &
                (jp1.ne.ngas) )                         then
            nlimx = nlimx + dnlim
            if ( ( (r1(i)-r1(jm1))**2 .lt. dr2(nh_max) ) .and. &
                  (jm1.ne.1) )  nlimn = nlimn + dnlim
            goto 100
          end if
          hmin(i) = sqrt(dr2(nh_min))
          hmax(i) = sqrt(dr2(nh_max))
        end do
!!!$OMP END PARALLEL DO
!
!       Write to file outside of main loop so parallelisation of main
!       loop is possible
        do i = 1,ngas
          write(3,3015) i, r1(i),  r2(i) , r3(i), hmin(i), hmax(i)
        end do
        close(3)
! 
        write(*,1003)
!       Sorting h's according to position and writing to file
!       Sorts arrays via ssort
        do i=1,ngas
          htmpn(i) = hmin(i)
          htmpx(i) = hmax(i)
        end do
        call ssort(r , htmpn, ngas, 2)
        call ssort(rb, htmpx, ngas, 2)
!       Print results to file
        do i=1,ngas
          write(5,3003) r(i), htmpn(i), htmpx(i)
        end do
        close(5)
! 
        write(*,1004)
!       Finding the desired h's
!       Visual inspection shows that there can be anomalous values, so for
!       best results, lets choose the third smallest h for nh_max.  
!       Likewise, the galaxy's radius (rdisk) is set in GalactICs, thus we
!       will chose nh_min at this radius; note also that visual inspection 
!       shows that 98.4% of all particles are within this radius.
!       This is to find the third smallest h
        hmax0 = 1.0d30
        do i=1,ngas
          if (hmax(i).lt.hmax0) then
            hmax0 = hmax(i)
            imax2 = i
          end if
        end do
        hmax0 = 1.0d30
        do i=1,ngas
          if ( (hmax(i).lt.hmax0).and.(i.ne.imax2) ) then
            hmax0 = hmax(i)
            imax1 = i
          end if
        end do
        hmax0 = 1.0d30
        do i=1,ngas
          if( (hmax(i).lt.hmax0).and.(i.ne.imax2).and.(i.ne.imax1) )then
            hmax0 = hmax(i)
            imax  = i
          end if
        end do
!       Calculate the h furthese from the centre
        if (comploop.eq.1) then 
!         This is to find the h at r closest to rdisk; doing this in reverse 
!         order since rdisk is close to the end of the list
          do i=ngas, 1, -1
            if (r(i).le.rdisk) goto 200
          end do
200       continue
          imin  = i + 1
        else
!         For the halo, just taking the point at the, say, 98% percentile
          imin = int(0.98 * ngas) 
        end if
        hmin0 = htmpn(imin) 
!       To find the actual location of the particle; currently only 
!       necessary for a print statement, but.....
!       NOTE: This name must match that of unit 3
        open(unit=7, file="find_hvals.txt")
        do i = 1,ngas
          read(7,3015)  itmp, ndummy,ndummy,ndummy,htmp,ndummy
          if (htmp.eq.hmin0) then
            imin = itmp
            goto 300
          end if
        end do
300     continue
        close(7)
! 
        write(*,1005)
!       The final values and their locations
        rmin = sqrt( r1(imin)**2 + r2(imin)**2 + r3(imin)**2 )
        rmax = sqrt( r1(imax)**2 + r2(imax)**2 + r3(imax)**2 )
        write(*,2001) nh_min, r1(imin), r2(imin), r3(imin), rmin, hmin0
        write(*,2001) nh_max, r1(imax), r2(imax), r3(imax), rmax, hmax0
        write(4,2001) nh_min, r1(imin), r2(imin), r3(imin), rmin, hmin0
        write(4,2001) nh_max, r1(imax), r2(imax), r3(imax), rmax, hmax0
!
!       Determining the coefficients
        R0 = (rmin - rmax) / log(hmin0/hmax0)
        h0 = hmin0 * exp(-rmin/R0)
        write(*,2002) h0, R0
        write(4,2002) h0, R0
        close(4)
! 
        write(*,1007)
!       writing h in the initial order.  Here, we choose h based upon the
!       equation h = h0 * exp(r/R0) if hmin<h<hmax; otherwise we chose
!       hmin or hmax, depending on where h is.
        call ssort(rnum, rnum2, ngas, 2)
        do i=1,ngas
          j    = int(rnum2(i))
          rval = sqrt( r1(j)**2 + r2(j)**2 + r3(j)**2 )
          hval = h0 * exp(rval/R0)
          if (hval.gt.hmax(j)) hval = hmax(j)
          if (hval.lt.hmin(j)) hval = hmin(j)
          write(6,3006) r1(j), r2(j), r3(j), hval, hmin(j), hmax(j)
        end do
        close(6)
! 
        write(*,1008)
!       A check to make sure that orders are correct 
!       ie: that the order of r in the initial file (unit 1) is the same 
!           as the final file (unit 6)
        open(unit=11, file="find_original_r.txt")
        if (comploop.eq.1) then 
          open(unit=16, file="find_rh.txt")
        else
          open(unit=16, file="find_rh_halo.txt")
        end if
        tol = 1.0d-7
        do i=1,ngas
          read(11,3003) r1_11,r2_11,r3_11
          read(16,3006) r1_16,r2_16,r3_16
          if (r1_11.ne.r1_16) write(*,2004) 'r1', i, 1, r1_11, 1, r1_16
          if (r2_11.ne.r2_16) write(*,2004) 'r2', i, 2, r2_11, 2, r2_16
          if (r3_11.ne.r3_16) write(*,2004) 'r3', i, 3, r3_11, 3, r3_16
        end do
        close(11)
        close(16)
! 
!       Delete all unnecessary files
        call system("rm find_hvals.txt")        
        call system("rm find_log.txt")
        call system("rm find_hvals2.txt")
        call system("rm find_original_r.txt")
!       We are now done the loop; repeat for the halo, or exit
        write(*,1100)
      end do
!
!
!
!     To time the routine
      call system_clock(count1, ctrate1)
      ctrate0 = 0.5*(ctrate0+ctrate1)
      timesec = (dble(count1-count0)/dble(ctrate0))
      write(*,2003) timesec
!     Closing print statment   
      write(*,1006)
      write(*,1100)
! 
1100  format('CALCULATEh0: ')
1000  format('CALCULATEh0: Initialising programme and opening files')
1010  format('CALCULATEh0: The input file for the gas disk is:')
1012  format('CALCULATEh0: The input file for the gas halo is:')
1011  format('CALCULATEh0: ',a)
1009  format('CALCULATEh0: Resetting r3 = 0 before calculation')
1001  format('CALCULATEh0: Sorting by ascending r1 (i.e. x)')
1002  format('CALCULATEh0: Finding h for all particles')
1003  format('CALCULATEh0: Sorting h by particle distance from origin')
1004  format('CALCULATEh0: Find extreme values of h')
1005  format('CALCULATEh0: Finding points for an exponential fit')
1007  format('CALCULATEh0: Writing h per particle in original order')
1008  format('CALCULATEh0: Verifying initial and final particle orders',&
            ' agree')
1006  format('CALCULATEh0: Program complete')
2001  format('CALCULATEh0: For',I5,' neighbours, we use the particleat', &
            F,F,F,' (with distance from centre ',F,') whose h is ',F)
2002  format('CALCULATEh0: The coeffients are h0 = ',F,' and R0 = ',F)
2004  format('CALCULATEh0: WARNING: Initial and final ',a,' disagreeat', &
             I,' with r',I1.1,'_i = ',F,' and r',I1.1,'_f = ',F)
2003  format('CALCULATEh0: Time since run began:',F7.2,' seconds.')
3003  format(  3F)
3004  format(  4F)
3006  format(  6F)
3007  format(  7F)
3011  format(I, F)
3015  format(I,5F)
3016  format('CALCULATEh0: Analysing gas disk with ',I,' particles')
3017  format('CALCULATEh0: Analysing gas halo with ',I,' particles')
 
      return
      end



