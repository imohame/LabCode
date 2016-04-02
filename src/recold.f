      function recold(x,y,n)
c     implicit double precision (a-h,o-z)                                    dp
c
c     reduce the diagonal term in a column
c
c     note- an assembly language version of this routine
c     is available for the cray-1 and for the cdc-7600.
c
      dimension x(*),y(*)
      common/fissl5/dt(2)

      s=0.
      if(n.eq.0) go to 11
      do 10 i=1,n
      s=s+x(i)*y(i)*x(i)
      x(i)=x(i)*y(i)
   10 continue
   11 y(n+1)=x(n+1)-s
      dt(1)=x(n+1)*y(n+1)
      dt(2)=dim(abs(y(n+1)),.5e-7*abs(x(n+1)+s))
      if(y(n+1).ne.0.) y(n+1)=1./y(n+1)
      x(n+1)=y(n+1)
      recold=x(n+1)
      return
      end
