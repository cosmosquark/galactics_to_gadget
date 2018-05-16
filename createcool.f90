      subroutine createcool(zmetc)
      implicit none

!     interpolates from the Sutherland & Dopita tables to any metallicity

      integer*4 :: ntemp,nz
      parameter(ntemp=151,nz=5)
      real*8 :: temp(ntemp),zmet(nz)
      real*8 :: zmetc
      real*8 :: scale
      real*8 :: told,tnew,tc
      integer*4 :: i,ii,ic
      real*8 :: lfind(ntemp),lamc ! starts with l, but is real
      integer*8, parameter ::  nctab = 151
      real*8 :: ctab0(nctab)

!     Cooling rates at various redshifts from Sutherland and Dopita
      include "coolingtable.inc"              
      include "units.inc"
 
!     conversion to Hydra units (we use (n_e+n_i)**2, not n_e*n_H)
!     n_en_H/n_t^2=0.52*0.44=0.229
      scale=0.229 
                 
      do i=1,ntemp
		  temp(i)=10.**(1.000+0.05*(i-1)) ! From 10K to 10**8.5 K
		  do ii=1,nz  
			  lambda(ii,i)=scale*10.**lambda(ii,i)
		  enddo       
      enddo       
      zmet(1)=0.  
      zmet(2)=0.01
      zmet(3)=0.1 
      zmet(4)=1./3.
      zmet(5)=1.
                 
      ic=0
      do i=1,nz   
        if (zmet(i).le.zmetc) ic=i
      enddo       
      if (ic.eq.0) stop 'metallicity out of range'
      if (zmetc.eq.zmet(ic)) then
		  do  i=1,ntemp
			lfind(i)=lambda(ic,i)
		  enddo       
	  else        
		  do i=1,ntemp
	  lfind(i)=lambda(ic,i)+(lambda(ic+1,i)-lambda(ic,i))* &
                   (zmetc-zmet(ic))/(zmet(ic+1)-zmet(ic))
		  enddo       
      endif       
                  
      told = 10.
      ctab0(1)=0. 
      ic=1        
      write(*,*) "finishing off cooling table"
      do i=2,nctab
         tnew=told*10**0.1
         tc=(tnew+told)/2.
 1       if ((tc.gt.temp(ic+1)).and.((ic+1).lt.ntemp)) then
			 ic=ic+1  
			 goto 1   
         endif    
         lamc=lfind(ic)+(lfind(ic+1)-lfind(ic))* &
             (tc-temp(ic))/(temp(ic+1)-temp(ic))
         ctab0(i)=ctab0(i-1)+1.5*kb*(tnew-told)/lamc/cunit
         write(*,"(I5,E12.5,E12.5)") i,lamc/mum,ctab0(i)*cunit*mum
         told=tnew
      enddo      

      end         

