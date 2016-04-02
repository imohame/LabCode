      subroutine inver (a)
c     implicit double precision (a-h,o-z)                                    dp
c
      dimension a(4,*)
c
      do 60 n=1,4
c
      d=a(n,n)
      do 10 j=1,4
   10 a(n,j)=-a(n,j)/d
c
      do 50 i=1,4
      if (n-i) 20,50,20
   20 do 40 j=1,4
      if (n-j) 30,40,30
   30 a(i,j)=a(i,j)+a(i,n)*a(n,j)
   40 continue
   50 a(i,n)=a(i,n)/d
c
      a(n,n)=1.0/d
c
   60 continue
c
      return
      end
