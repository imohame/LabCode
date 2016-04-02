      subroutine setpnt (j,i,len)
c     implicit double precision (a-h,o-z)                                    dp
      common /bk00/ ioff(96)
      ioff(j)=ioff(i)+len
      return
      end
