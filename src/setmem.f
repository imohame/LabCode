      subroutine setmem (i,n)
c     implicit double precision (a-h,o-z)                                    dp
c
c.... expand memory
c
      common /bk00/ ioff(96)
      nn=ioff(i)+n
      call expndm(nn)
      return
      end
