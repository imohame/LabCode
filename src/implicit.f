      subroutine implicit(y,t0,h,N)

*  ____________________________________________________________________
* |                                                                    |
* |...  ITERATION SCHEME FOR BACKWARD EULER APPROXIMATION OF           |
* |...  AN INITIAL VECTOR VALUED PROBLEM ITERATION SCHEME IS BASED     |
* |...  ON ON A QUASI-NEWTON METHOD IN AN 1984 IMSL SUBROUTINE.        |
* |...  REFERENCE IS THE MINIPACK PACKAGE BASED ON POWELL'S HYBRID     |
* |...  HYBRID ALGORITHM WHICH USES A FINITE DIFFERENCE APPROXIMATION  |
* |...  TO THE JACOBIAN AND TAKES PERCAUTIONS TO AVOID LARGE STEP SIZES|
* |...  OR INCREASING RESIDUALS.                                       |   
* |                                                                    |
* |           x(*) : INITIAL STARTING VECTOR FOR QUASI NEWTON METHOD   |
* |                                                                    |
* |...  THE CLOSER THE GUESS,THE FASTER IT CONVERGES.                  |
* |                                                                    |
* |                                     Modified by .................  |
* |                                                 Waeil M. Ashmawi   |
* |                                                 NC State University|
* |                                                 April 07, 1998     |
* |____________________________________________________________________|
*
 
      parameter (nss = 24)
	  integer n2, N 
      dimension x(N), y(N), wk(N*(3*N+15)/2), par(N+1)
      external fcn

*
* ======================== S T A R T   C O D E ========================
* 
      write(*,*) 'in implicit'
      n2=N+1
      do i = 1, N
         x(i)   = 1.0*y(i)
         par(i) = y(i)
c		 write(*,*) 'x,y',x(i),y(i)
      end do
      par(n2) = h
      nsig    = 3
      itmax   = 200
      do i=1,n2
	    write(*,*) 'par in implicit',par(i)
	  enddo
		
      call zspow(fcn,nsig,N,itmax,par,x,fnorm,wk,ier)

      if (ier .eq. 131) then
	     write(*,*) 'implicit.f ier=131'
         stop
      endif

      do i = 1, N
         y(i) = x(i)
      end do
c      WRITE(27,22)
22    FORMAT('IMPLICIT SCHEME')

*
* ========================== E N D   C O D E ==========================
*

      return
      end
