      subroutine LU(A,d,x,n,nbwth)
      integer n,nbwth
      real*8 A(n,n),d(n,1),x(n,1)
      do i=1,n
c      write(*,*) d(i,1)
      enddo
c      write(*,*) nbwth

      do k=1,n-1
         do i=k+1,min(k+nbwth,n)
            A(i,k)=A(i,k)/A(k,k)
         enddo
         do j=k+1,min(k+nbwth,n)
            do i=k+1,min(k+nbwth,n)
               A(i,j)=A(i,j)-A(i,k)*A(k,j)
            enddo
         enddo
      enddo
      do i=1,n
c       write(*,*) A(i,111:121)
       enddo
      
      do j=1,n
         do i=j+1,min(j+nbwth,n)
            d(i,1)=d(i,1)-A(i,j)*d(j,1)
         enddo
      enddo
      
      nn=n
      do j=1,n
c      write(*,*) nn
      d(nn,1)=d(nn,1)/A(nn,nn)
         do i=max(1,nn-nbwth),nn-1
            d(i,1)=d(i,1)-A(i,nn)*d(nn,1)
         enddo
      nn=nn-1
      enddo
      x=d
      do i=1,n
c      write(*,*) x(i,1)
      enddo
      return
      end
      
      