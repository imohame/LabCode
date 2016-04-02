      Subroutine cgssor3(rhs,sol,n,m,nx)
!
!   Solves via PCG with SSOR preconditioner.
!   Uses sparse implementation on space grid.
!
      implicit none
      real*8 rhs(n,m)
      real*8 sol(n,m)
      real,dimension (0:nx+1,0:nx+1,1:m):: u,p,q,r,rhat
      real,dimension (0:nx+1) :: x,y
      real,dimension (1:m+1) :: oldrho, rho,alpha
      real :: error,dx2,w,time
      integer ::i,j,k,kk,n,nx,m,mcg
      w = 1.7
      u = 0.0
      kk = m
      r = 0.0
      rhat = 0.0
      q = 0.0
      p = 0.0
      do k = 1,kk
      do j = 1,nx
        do i = 1,nx
          r(i,j,k) = rhs(i+(j-1)*(nx),k)
        enddo
      enddo
      enddo
      error = 1. 
      mcg = 0
      rho  = 0.0
  
      do while ((error>.0001).and.(mcg<200))
        mcg = mcg+1
        oldrho = rho
!	Execute SSOR preconditioner
        do k=1,kk
	do j= 1,nx
	  do i = 1,nx
	    rhat(i,j,k) = w*(r(i,j,k)+rhat(i-1,j,k)+rhat(i,j-1,k))*.25
	  enddo
	enddo
	enddo
	rhat(1:nx,1:nx,1:kk) =  ((2.-w)/w)*4.*rhat(1:nx,1:nx,1:kk)
	do k = 1,kk
	do j= nx,1,-1
	  do i = nx,1,-1
	  rhat(i,j,k) = w*(rhat(i,j,k)+rhat(i+1,j,k)+rhat(i,j+1,k))*.25
	  enddo
	enddo
	enddo
!	Find conjugate direction 
	do k = 1,kk      
        rho(k) = sum(r(1:nx,1:nx,k)*rhat(1:nx,1:nx,k))
        If (rho(k).ne.0.) then
        if (mcg.eq.1) then
          p(1:nx,1:nx,k) = rhat(1:nx,1:nx,k)
        else
          p(1:nx,1:nx,k) = rhat(1:nx,1:nx,k) + 
     1     			(rho(k)/oldrho(k))*p(1:nx,1:nx,k)     
        endif
!       Execute matrix product q = Ap
        q(1:nx,1:nx,k)=4.0*p(1:nx,1:nx,k)-p(0:nx-1,1:nx,k)-
     1   		p(2:nx+1,1:nx,k)
     1    		- p(1:nx,0:nx-1,k) - p(1:nx,2:nx+1,k)
!	Find steepest descent
        alpha(k) = rho(k)/sum(p(1:nx,1:nx,k)*q(1:nx,1:nx,k))
        u(1:nx,1:nx,k) = u(1:nx,1:nx,k) + alpha(k)*p(1:nx,1:nx,k)
        r(1:nx,1:nx,k) = r(1:nx,1:nx,k) - alpha(k)*q(1:nx,1:nx,k)
        endif
        enddo
        error = maxval(abs(r(1:nx,1:nx,1:kk)))
      enddo
      do k = 1,kk
      do j = 1,nx
        do i = 1,nx
          sol(i+(j-1)*(nx),k) = u(i,j,k)
        enddo
      enddo
      enddo
!      print*, mcg ,error
      end subroutine
