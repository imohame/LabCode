      subroutine movvar(a,b,n)
c     implicit double precision (a-h,o-z)                                    dp
      dimension a(*),b(*)
      m=n+1
      do 10 i=1,n
   10 b(m-i)=a(m-i)
      return
      end
