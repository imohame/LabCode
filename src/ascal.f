      subroutine ascal(a,n,scal)
c     implicit double precision (a-h,o-z)                                    dp
      dimension a(*)
      do 10 i=1,n
   10 a(i)=a(i)*scal
      return
      end
