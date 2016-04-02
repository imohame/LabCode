      subroutine matrixmult(A,B,C,rowsleft,columnsleft,columnsright)
      integer rowsleft, columnsleft,columnsright
      real*8 A(rowsleft,columnsleft), B(columnsleft,columnsright)
      real*8 C(rowsleft,columnsright)
      
     
      do i=1,rowsleft
      do j=1,columnsright
      C(i,j)=0.
      enddo
      enddo
      do i = 1,rowsleft
      do j = 1,columnsright
          do k = 1, columnsleft
             
                C(i,j)  = C(i,j) + A(i,k)*B(k,j)
             enddo
           enddo
         enddo
         return
         end