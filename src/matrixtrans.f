      subroutine matrixtrans(A,AT,n,m)
      integer n,m
      real*8 A(n,m),AT(m,n)
      do i=1,n
      do j=1,m
      AT(j,i)=A(i,j)
      enddo
      enddo
      return
      end