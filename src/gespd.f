      Subroutine gespd(a,rhs,sol,n,m)
!
!  Solves Ax = d with A a nxn SPD and d a nxm.
!
      implicit none
      real*8, dimension(n,n), intent(inout):: a
      real*8, dimension(n,m), intent(inout):: rhs
      real*8, dimension(n,n+m):: aa
      real*8, dimension(n,m) :: y
      real*8, dimension(n,m),intent(out)::sol
      integer ::k,i,j,l
      integer,intent(in)::n,m
      aa(1:n,1:n)= a
      aa(1:n,(n+1):n+m) = rhs
!
!  Factor A via column version and
!  write over the matrix.
!
      do k=1,n-1
        aa(k+1:n,k) = aa(k+1:n,k)/aa(k,k)
        do j=k+1,n
          do i=k+1,n
            aa(i,j) = aa(i,j) - aa(i,k)*aa(k,j)
          enddo
        enddo
      enddo
!
!  Solve Ly = d via column version and
!  multiple right sides.
!
      do j=1,n-1
        do l =1,m
          y(j,l)=aa(j,n+l)
        enddo
        do i = j+1,n
          do l=1,m
            aa(i,n+l) = aa(i,n+l) - aa(i,j)*y(j,l)
          enddo
        enddo
      enddo
!
!  Solve Ux = y via column version and
!  multiple right sides.
!
      do j=n,2,-1
        do l = 1,m
          sol(j,l) = aa(j,n+l)/aa(j,j)
        enddo
        do i = 1,j-1
          do l=1,m
            aa(i,n+l)=aa(i,n+l)-aa(i,j)*sol(j,l)
          enddo
        enddo
      enddo
      do l=1,m
        sol(1,l) = aa(1,n+l)/a(1,1)
      enddo
      end 