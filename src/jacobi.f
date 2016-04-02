      subroutine jacobi(a,d,v,b,z,n,ir)
c     implicit double precision (a-h,o-z)                                    dp
      dimension a(*),b(*),d(*),v(n,1),z(*)
      n1=n-1
      call pzero (v,n*n)
      do 10 i=1,n
   10 v(i,i)=1.0
      k=0
      do 20 i=1,n
      k=k+i
   20 b(i)=a(k)
      do 30 i=1,n
      d(i)=b(i)
   30 z(i)=0.0
      ir=0
      do 170 m=1,50
      sm=0.0
      do 40 j=2,n
      k=j*(j-1)/2
      j1=j-1
      do 40 i=1,j1
   40 sm=sm+abs(a(i+k))
      if (sm.le.0.0) return
      thresh=0.0
      if (m.lt.4) thresh=0.20*sm/(n*n)
      do 150 i=1,n1
      ip1=i+1
      do 140 j=ip1,n
      ii=i*(i-1)/2
      jj=j*(j-1)/2
      ij=jj+i
      g=100.*abs(a(ij))
      if ((i.gt.4).and.(abs(d(i))+g.eq.abs(d(i))).and.(abs(d(j))+g.eq
     1 .abs(d(j)))) a(ij)=0.0
      if (abs(a(ij)).le.thresh) go to 140
      h=d(j)-d(i)
      if (abs(h)+g.gt.abs(h)) go to 50
      t=a(ij)/h
      go to 60
   50 theta=0.5*h/a(ij)
      t=1./(abs(theta)+sqrt(1.+theta*theta))
      if (theta.lt.0.0) t=-t
   60 c=1./sqrt(1.+t*t)
      s=t*c
      tau=s/(1.+c)
      h=t*a(ij)
      z(i)=z(i)-h
      z(j)=z(j)+h
      d(i)=d(i)-h
      d(j)=d(j)+h
      a(ij)=0.0
      i1=i-1
      if (i1.le.0) go to 80
cdir$ ivdep
      do 70 k=1,i1
      g=a(k+ii)
      h=a(k+jj)
      a(k+ii)=g-s*(h+g*tau)
   70 a(k+jj)=h+s*(g-h*tau)
   80 i1=i+1
      j1=j-1
      ii=ii+i+i
      if (i1.gt.j1) go to 100
      do 90 k=i1,j1
      g=a(ii)
      h=a(k+jj)
      a(ii)=g-s*(h+g*tau)
      a(k+jj)=h+s*(g-h*tau)
   90 ii=ii+k
  100 j1=j+1
      if (j1.gt.n) go to 120
      ii=ii+j
      jj=jj+j+j
      do 110 k=j1,n
      g=a(ii)
      h=a(jj)
      a(ii)=g-s*(h+g*tau)
      a(jj)=h+s*(g-h*tau)
      ii=ii+k
  110 jj=jj+k
  120 do 130 k=1,n
      g=v(k,i)
      h=v(k,j)
      v(k,i)=g-s*(h+g*tau)
  130 v(k,j)=h+s*(g-h*tau)
      ir=ir+1
  140 continue
  150 continue
      do 160 i=1,n
      d(i)=b(i)+z(i)
      b(i)=d(i)
  160 z(i)=0.0
  170 continue
      return
      end
