!      Makes a hot gas halo diustribution with the density profile
!      rho(r) = rho_0 * (1+(r/r_c)^2)^(-1.5beta)
! 
      subroutine makehalo(ngas,rc,rm,rf,Mm,totmass,mgas,alpha,rho0,ipe)
! 
      implicit none
      integer :: i, j, k, n, ngas, dn, cell, Nmax, ds, nremain,nremain0,ctr,ipe
      !! high res run
!      parameter(n    = 20000)
      !! low res run
      parameter(n    = 2000)
      parameter(Nmax = 15000000)
      real*8 :: pi, ri,rm,rf,alpha,Mm,rc, rho0, menc, mgas, totmass, tol, &
            dM, dr, x0, x1, Mass, dMass,Mtop,Mbot, rinner, zbqluab, &
            r00, phi, theta, z, dMpart, rho02
      integer :: np(n)
      real*8 :: r(n), rho(n),mshell(n), rx(n),r0(Nmax), pos(3,Nmax),rpos(n)
      parameter(pi = 3.14159265)
!     Maybe set n such that rf/n = sft0 in hydra simulations
! 
! ----------------------------------------------------------------------
! INITIAL PARAMETERS
! Define Parameters (r in kpc; Mm in solar masses; alpha is unitless)
      ri    =   0.0
!      rc    =   1.75
!      rm    =  40.0
!      rf    = 275.0
       alpha =   1.0
!      Mm    =   0.02 * 4.1031467394787498d10
!      mgas  =   2.576d5
      write(*,*) &
      'We are normaling the density such that there is ', &
      Mm,'M_sun within the inner ',rm,'kpc'
      print*, ' '
!     initialise the random number generator
      call ZBQLINI(1)
!
! ----------------------------------------------------------------------
! DETERMINE THE NORMALISATION
      menc = 0.0
      rho0 = 1.0
!      print *, "Mm ", Mm
      call totalmass(n,ri,rc,rm,alpha,rho0,menc)
!      print *, menc, rho0
!      print *, menc, rho0
      !!!! not sure about this, it turns the density into a mass ratio.
      rho0 = Mm / menc
      write(*,'(a,F)') 'The normalisation constant is ', rho0
      print*, ' '
!
! ----------------------------------------------------------------------
! DETERMINE THE TOTAL MASS 
      call totalmass(n,ri,rc,rf,alpha,rho0,totmass)
      write(*,'(a,Es)') 'The total hot halo gas mass (for 0<r<r_f) is ', &
      totmass
      print*, ' '
      ngas = int(totmass/mgas + 0.5)
      write(*, *) 'Given a particle mass of ',mgas, &
                             'M_sun, we require ',ngas,' particles'
      print*, ' '
! ----------------------------------------------------------------------
! CREATE SHELLS OF EQUAL MASS
      dM = totmass/float(n)
      call findr(dM, rf, ri, rx, n, rho0,rc,alpha)
!
! ----------------------------------------------------------------------     
! ASSIGN an equal number of particles to each bin
      dn = int(float(ngas)/float(n))
      write(*,*) 'Total number of evenly distributed particles is ',dn
      do i = 1,n
        np(i) = dn
      end do
! Add an additional particle to each cell, starting at the outside; do
! this until we run out of particles
      do i = n, n - (ngas - dn*n-1), -1
        np(i) = np(i) + 1
      end do
      ctr = 0
      do i = 1,n
        ctr = ctr + np(i)
      end do
      print*, 'to verify: total number of particles:', ctr
!
! ----------------------------------------------------------------------
! Distribute the particles within each shell according to the density 
! profile

      rho02 = rho0
      open(unit=101, file="halogas_positions.txt")
      write(101, '(I, Es)') ngas, totmass
      k = 0
      do i = 1,n
!        print *, "death star", i
        dMpart = dM/float(np(i))
        if (i.eq.1) then 
          rinner = 0.0
        else
          rinner = rx(i-1)
        end if
!        print *, "weee", i, i, i
        call findr(dMpart, rx(i), rinner, rpos, np(i), rho02,rc,alpha)
!        print *, "zen"
        do j = 1,np(i)
!          print *, "ken", j
          phi      = zbqluab(dble(0.0),dble(2.0*pi))
!          print *, "ben"
          z        = zbqluab(dble(-rpos(j)),dble(rpos(j)))
!          print *, "ten"
          theta    = asin(z/rpos(j))
          pos(1,j) = rpos(j) * cos(theta)*cos(phi)
          pos(2,j) = rpos(j) * cos(theta)*sin(phi)
          pos(3,j) = rpos(j) * sin(theta)
!          print *, "zen", j
        end do
!       write to file
        do j = 1,np(i)
          k = k + 1
!          print *, "bleh", k
          write(101,'(I,4Es)') k, pos(1,j),pos(2,j),pos(3,j), rpos(j)
!          print *, "derpo", k
        end do
!      print *, "the end", i
      end do
      close(101)
!
! ----------------------------------------------------------------------
! END
      write(*,*) "Exiting makehalo"
      return
      end
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This subroutine will calculate the total mass with ri < r < rf
      subroutine totalmass(n,ri,rc,rf,alpha,rho0,totmass)
      implicit none
      integer :: i, n
      real*8 :: pi, ri, rc,rf,rho0,dr, rhoin, rhoout,alpha,totmass,rho,Mass
      real*8 :: r(n),rhoshell(n),menc(n),mshell(n)
      parameter(pi = 3.1415926535)
!
! Calculate the r-array
      dr = (rf-ri)/float(n)
      do i = 1,n
        r(i) = dr*float(i)
      end do
      if (r(n).ne.rf) then
        print*, 'r(n).ne.rf, thus there is a mistake; abort'
        stop 
      end if
!
! Calculate the density per shell (trapezdoid rule)

      rhoin  = rho(rho0,  0 , rc, alpha)
      rhoout = rho(rho0,r(1), rc, alpha)
      rhoshell(1) = 0.5 * (rhoin + rhoout)
      do i = 2,n
        rhoin  = rho(rho0,r(i-1), rc, alpha)
        rhoout = rho(rho0,r(i  ), rc, alpha)
        rhoshell(i) = 0.5 * (rhoin + rhoout)
      end do

!
! Calculate the mass per shell
      mshell(1) = Mass(r(1),0,rho0,rc, alpha)
      do i = 2,n
        mshell(i) = Mass(r(i),r(i-1),rho0,rc, alpha)
      end do

!
! Calculate the enclosed mass
      menc(1) = mshell(1)
      do i = 2,n
        menc(i) = menc(i-1) + mshell(i)
      end do
      totmass = menc(n)
!
! Write statementes
!      if (rho0.eq.1.0) then
!        open(unit=1, file="00mh_0rm.txt")
!      else
!        open(unit=1, file="00mh_0rf.txt")
!      end if
!      do i=1,n
!        write(1,'(I,4Es)') i, r(i),rhoshell(i),mshell(i),menc(i)
!      end do
!      close(1)
!     
      return
      end subroutine totalmass
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This subroutine will calculate M and dM/dr
      real(kind=8) function Mass(r_i, r_im1, rho0,rc, alpha)
      implicit none
! 
      real*8 :: pi, r_i, r_im1, rho0,rc, alpha, &
          rhoin,rhoout, drho, drhoout, rho
      parameter(pi = 3.1415926535)
! 
!     The two densities
!      print *, "dens1"
      rhoin  = rho(rho0,r_im1, rc, alpha)
!      print *, "dens2"
      rhoout = rho(rho0,r_i  , rc, alpha)
!
!     The mass
!      print *, "the mass"
      mass = 4.0/3.0*pi * 0.5*(rhoin+rhoout) * (r_i**3-r_im1**3)
! 
      end function Mass
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This function will calculate dM/dr
      real(kind=8) function dMass(r_i, r_im1,rho0,rc, alpha)
      implicit none
! 
      real*8 :: pi, r_i, r_im1, rho0,rc, alpha
      real*8 :: rhoin,rhoout, drho, drhoout, rho
      parameter(pi = 3.1415926535)
! 
!     The two densities
      rhoin  = rho(rho0,r_im1, rc, alpha)
      rhoout = rho(rho0,r_i  , rc, alpha)
!     The derivative of the density
      drhoout = drho(rho0,r_i, rc, alpha)
!
!     The derivative of the mass
      dMass = 2.0/3.0*pi * (drhoout * (r_i**3-r_im1**3) & 
           +  3.0*r_i**2 * (rhoin+rhoout) )
! 
      end function dMass
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This function will calculate rho 
      real(kind=8) function rho(rho0,r, rc, alpha)
      implicit none
      real*8 :: rho0,r, rc, alpha, rho02, alpha2
! 
      rho02 = rho0
      alpha2 = alpha
!      print *, "test rho"
!      print *, r
!      print *, rc
!      print *, rho02
!      print *, alpha2                -alpha
      rho = rho02*(1.0 + (r/rc)**2)**(-1.0)
!     
!      print *, "test rho end"
      end function rho
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This function will calculate drho/dr
      real(kind=8) function drho(rho0,r, rc, alpha)
      implicit none
      real*8 :: rho0,r, rc, alpha, rho02, alpha2
! 
      rho02 = rho0
      alpha2 = alpha
      drho = -2.0*alpha2*rho02*r/rc**2 * (1.0 + (r/rc)**2)**(-alpha2-1.0)
! 
      end function drho
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This function will find the borders for shells of equal mass
      subroutine findr(dM, rf, ri, rx, n, rho0,rc,alpha)
      implicit none
      integer :: j, i, n
      real*8 :: dM, rf, ri, rx(n), rho0, rc, alpha, rho02, alpha2
      real*8 :: dr, x0, x1, Mtop, Mbot, Mass, dMass, rinner, tol
!
!     Input tolerance
      tol   =   1.0d-10
!      print *, "starting findr"
      rho02 = rho0
!      print *, "derp"
      alpha2 = 1.0
!      print *, "doop"
!     initial guess for the first loop
      dr = (rf - ri)/float(n)
      x0 = dr + ri
!     The positions will be calculated via Newton's Method
      do i = 1,n
!        print *, i, "newtons"
        if (i.eq.1) then 
!          print *, "i = 1"
          rinner = ri
        else
!          print *, "i = 1 else"
          rinner = rx(i-1)
        end if
!        print *, "weee"
!        print *, rx(i), "rx", rinner, "r inner"
        do j = 1,100000
!          print *, j, "j"
!          print *, x0, "x0"
!          print *, rinner, "rinner"
!          print *, rho02, "rho0"
!          print *, rc, "rc"
!          print *, alpha, "alpha"
!          print *, j, "start j", x0, rinner, rho0, rc, alpha
          Mtop =  Mass(x0, rinner, rho02,rc,alpha2)
!          print *, "herp"
          Mbot = dMass(x0, rinner, rho02,rc,alpha2)
!          print *, "derp"
          x1 = x0 - (Mtop-dM)/Mbot
!          print *, "doop"
          if (max(x1,x0)/min(x1,x0).le.1.0+tol) goto 100
!          print *, "deep"
          x0 = x1
!          print *, j, "end2"
        end do
100     continue
!       write result to array
!        print *, "beep"
        rx(i) = x1
!       guess x0 for next loop
!        print *, "boop"
        x0 = x1 + dr
      end do
!      print *, "baap"
      open(unit=1, file="00mh_0rx.txt")
      do i=1,n
!        print *, "omg lol wtf", i
        write(1,'(I,4Es)') i, rx(i)
      end do
!      print *, "doo"
      close(1)
! 
      end subroutine findr
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! 
!******************************************************************
!*******	FILE: randgen.f				***********
!*******	AUTHORS: Richard Chandler		***********
!*******		 (richard@stats.ucl.ac.uk)	***********
!*******		 Paul Northrop 			***********
!*******		 (northrop@stats.ox.ac.uk)	***********
!*******	LAST MODIFIED: 26/8/03			***********
!*******	See file randgen.txt for details	***********
!******************************************************************

      BLOCK DATA ZBQLBD01
!
!       Initializes seed array etc. for random number generator.
!       The values below have themselves been generated using the
!       NAG generator.
!
      COMMON /ZBQL0001/ ZBQLIX,B,C
      DOUBLE PRECISION  ZBQLIX(43),B,C
      INTEGER I
      DATA (ZBQLIX(I),I=1,43) /8.001441D7,5.5321801D8, &
      1.69570999D8,2.88589940D8,2.91581871D8,1.03842493D8, &
      7.9952507D7,3.81202335D8,3.11575334D8,4.02878631D8, &
      2.49757109D8,1.15192595D8,2.10629619D8,3.99952890D8, &
      4.12280521D8,1.33873288D8,7.1345525D7,2.23467704D8, &
      2.82934796D8,9.9756750D7,1.68564303D8,2.86817366D8, &
      1.14310713D8,3.47045253D8,9.3762426D7 ,1.09670477D8, &
      3.20029657D8,3.26369301D8,9.441177D6,3.53244738D8, &
      2.44771580D8,1.59804337D8,2.07319904D8,3.37342907D8, &
      3.75423178D8,7.0893571D7 ,4.26059785D8,3.95854390D8, &
      2.0081010D7,5.9250059D7,1.62176640D8,3.20429173D8, &
      2.63576576D8/
      DATA B / 4.294967291D9 /
      DATA C / 0.0D0 /
      END
!*****************************************************************
!*****************************************************************
!*****************************************************************
      SUBROUTINE ZBQLINI(SEED)
!*****************************************************************
!       To initialize the random number generator - either
!       repeatably or nonrepeatably. Need double precision
!       variables because integer storage can't handle the
!       numbers involved
!*****************************************************************
!	ARGUMENTS
!	=========
!	SEED	(integer, input). User-input number which generates
!		elements of the array ZBQLIX, which is subsequently used 
!		in the random number generation algorithm. If SEED=0,
!		the array is seeded using the system clock if the 
!		FORTRAN implementation allows it.
!*****************************************************************
!	PARAMETERS
!	==========
!	LFLNO	(integer). Number of lowest file handle to try when
!		opening a temporary file to copy the system clock into.
!		Default is 80 to keep out of the way of any existing
!		open files (although the program keeps searching till
!		it finds an available handle). If this causes problems,
!               (which will only happen if handles 80 through 99 are 
!               already in use), decrease the default value.
!*****************************************************************
      INTEGER LFLNO
      PARAMETER (LFLNO=80)
!*****************************************************************
!	VARIABLES
!	=========
!	SEED	See above
!	ZBQLIX	Seed array for the random number generator. Defined
!		in ZBQLBD01
!	B,C	Used in congruential initialisation of ZBQLIX
!	SS,MM,}	System clock secs, mins, hours and days
!	HH,DD }
!	FILNO	File handle used for temporary file
!	INIT	Indicates whether generator has already been initialised
!
      INTEGER SEED,SS,MM,HH,DD,FILNO,I
      INTEGER INIT
      DOUBLE PRECISION ZBQLIX(43),B,C
      DOUBLE PRECISION TMPVAR1,DSS,DMM,DHH,DDD

      COMMON /ZBQL0001/ ZBQLIX,B,C
      SAVE INIT

!
!	Ensure we don't call this more than once in a program
!
      IF (INIT.GE.1) THEN
       IF(INIT.EQ.1) THEN
        WRITE(*,1)
        INIT = 2
       ENDIF
       RETURN
      ELSE
       INIT = 1
      ENDIF
!
!       If SEED = 0, cat the contents of the clock into a file
!       and transform to obtain ZQBLIX(1), then use a congr.
!       algorithm to set remaining elements. Otherwise take
!       specified value of SEED.
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>	NB FOR SYSTEMS WHICH DO NOT SUPPORT THE  >>>>>>>
!>>>>>>>	(NON-STANDARD) 'CALL SYSTEM' COMMAND,    >>>>>>>
!>>>>>>>	THIS WILL NOT WORK, AND THE FIRST CLAUSE >>>>>>>
!>>>>>>>	OF THE FOLLOWING IF BLOCK SHOULD BE	 >>>>>>>
!>>>>>>>	COMMENTED OUT.				 >>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF (SEED.EQ.0) THEN
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>	COMMENT OUT FROM HERE IF YOU DON'T HAVE  >>>>>>>
!>>>>>>>	'CALL SYSTEM' CAPABILITY ...		 >>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       CALL SYSTEM(' date +%S%M%H%j > zbql1234.tmp')
!
!       Try all file numbers for LFLNO to 999 
!
       FILNO = LFLNO
 10    OPEN(FILNO,FILE='zbql1234.tmp',ERR=11)
       GOTO 12
 11    FILNO = FILNO + 1
       IF (FILNO.GT.999) THEN
        WRITE(*,2)
        RETURN
       ENDIF
       GOTO 10
 12    READ(FILNO,'(3(I2),I3)') SS,MM,HH,DD
       CLOSE(FILNO)
       CALL SYSTEM('rm zbql1234.tmp')
       DSS = DINT((DBLE(SS)/6.0D1) * B)
       DMM = DINT((DBLE(MM)/6.0D1) * B)
       DHH = DINT((DBLE(HH)/2.4D1) * B)
       DDD = DINT((DBLE(DD)/3.65D2) * B)
       TMPVAR1 = DMOD(DSS+DMM+DHH+DDD,B)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<	... TO HERE (END OF COMMENTING OUT FOR 	  <<<<<<<
!<<<<<<<<	USERS WITHOUT 'CALL SYSTEM' CAPABILITY	  <<<<<<<
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ELSE
       TMPVAR1 = DMOD(DBLE(SEED),B)
      ENDIF
      ZBQLIX(1) = TMPVAR1
      DO 100 I = 2,43
       TMPVAR1 = ZBQLIX(I-1)*3.0269D4
       TMPVAR1 = DMOD(TMPVAR1,B)       
       ZBQLIX(I) = TMPVAR1
 100  CONTINUE

 1    FORMAT(//5X,'****WARNING**** You have called routine ZBQLINI ', &
      'more than',/5X,'once. I''m ignoring any subsequent calls.',//)
 2    FORMAT(//5X,'**** ERROR **** In routine ZBQLINI, I couldn''t', &
      ' find an',/5X, &
      'available file number. To rectify the problem, decrease the ', &
      'value of',/5X, &
      'the parameter LFLNO at the start of this routine (in file ', &
      'randgen.f)',/5X, &
      'and recompile. Any number less than 100 should work.')
      end subroutine ZBQLINI
!*****************************************************************
      FUNCTION ZBQLU01(DUMMY)
!
!       Returns a uniform random number between 0 & 1, using
!       a Marsaglia-Zaman type subtract-with-borrow generator.
!       Uses double precision, rather than integer, arithmetic 
!       throughout because MZ's integer constants overflow
!       32-bit integer storage (which goes from -2^31 to 2^31).
!       Ideally, we would explicitly truncate all integer 
!       quantities at each stage to ensure that the double
!       precision representations do not accumulate approximation
!       error; however, on some machines the use of DNINT to
!       accomplish this is *seriously* slow (run-time increased
!       by a factor of about 3). This double precision version 
!       has been tested against an integer implementation that
!       uses long integers (non-standard and, again, slow) -
!       the output was identical up to the 16th decimal place
!       after 10^10 calls, so we're probably OK ...
!
      DOUBLE PRECISION ZBQLU01,DUMMY,B,C,ZBQLIX(43),X,B2,BINV
      INTEGER CURPOS,ID22,ID43

      COMMON /ZBQL0001/ ZBQLIX,B,C
      SAVE /ZBQL0001/
      SAVE CURPOS,ID22,ID43
      DATA CURPOS,ID22,ID43 /1,22,43/

      B2 = B
      BINV = 1.0D0/B
 5    X = ZBQLIX(ID22) - ZBQLIX(ID43) - C
      IF (X.LT.0.0D0) THEN
       X = X + B
       C = 1.0D0
      ELSE
       C = 0.0D0
      ENDIF
      ZBQLIX(ID43) = X
!
!     Update array pointers. Do explicit check for bounds of each to
!     avoid expense of modular arithmetic. If one of them is 0 the others
!     won't be
!
      CURPOS = CURPOS - 1
      ID22 = ID22 - 1
      ID43 = ID43 - 1
      IF (CURPOS.EQ.0) THEN
       CURPOS=43
      ELSEIF (ID22.EQ.0) THEN
       ID22 = 43
      ELSEIF (ID43.EQ.0) THEN
       ID43 = 43
      ENDIF
!
!     The integer arithmetic there can yield X=0, which can cause 
!     problems in subsequent routines (e.g. ZBQLEXP). The problem
!     is simply that X is discrete whereas U is supposed to 
!     be continuous - hence if X is 0, go back and generate another
!     X and return X/B^2 (etc.), which will be uniform on (0,1/B). 
!
      IF (X.LT.BINV) THEN
       B2 = B2*B
       GOTO 5
      ENDIF

      ZBQLU01 = X/B2

      END FUNCTION ZBQLU01
!*****************************************************************
      FUNCTION ZBQLUAB(A,B)
!
!       Returns a random number uniformly distributed on (A,B)
!
      DOUBLE PRECISION A,B,ZBQLU01,ZBQLUAB
      
!
!       Even if A > B, this will work as B-A will then be -ve
!
      IF (A.NE.B) THEN
       ZBQLUAB = A + ( (B-A)*ZBQLU01(0.0D0) )
      ELSE
       ZBQLUAB = A
       WRITE(*,1)
      ENDIF
 1    FORMAT(/5X,'****WARNING**** (function ZBQLUAB) Upper and ', &
      'lower limits on uniform',/5X,'distribution are identical',/)
      END FUNCTION ZBQLUAB
