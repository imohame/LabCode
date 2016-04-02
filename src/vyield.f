      subroutine vyield(eps,sg,ep,ln,ak,qh)
c     implicit double precision (a-h,o-z)                                    dp
c
c     vectorized subroutine for interpolating tabulated curve
c
c      common/range/mft,mlt,lft,llt,nftm1
c      common/pnts/npoint
c      common/double/iprec,ncpw,unit
c      common/vect16/xmult(128,9),epx(128),
c     1 epsj1(128),epsj(128),sgj1(128),sgj(128),
c     2 factor(128),qhj1(128),qhj(128)
c      dimension eps(*),sg(*),ep(ln,*),ak(*),qh(1)
c
c      do 10 i=3,8
c      n=i-2
c      if (eps(i).eq.0.) go to 20
c   10 n=n+1
c   20 npoint=n+1
c      do 30 i=mft,mlt
c      epx(i)=ep(1,i)
c   30 continue
c
c      do 40 i=mft,mlt
c      epsj1(i)=0.
c      epsj(i)=0.
c      sgj1(i)=0.
c      sgj(i)=0.
c      xmult(i,1)=0.
c   40 xmult(i,n+1)=1.
c      do 60 j=2,n
c      do 50 i=mft,mlt
c   50 xmult(i,j)=.5+sign(.5*unit,eps(j)-epx(i))
c   60 continue
c      do 80 j=1,n
c      nn=n-j+2
c      mm=n-j+1
c      do 70 i=mft,mlt
c   70 xmult(i,nn)=xmult(i,nn)-xmult(i,mm)
c   80 continue
c      do 100 j=1,n
c      do 90 i=mft,mlt
c      epsj1(i)=epsj1(i)+eps(j)*xmult(i,j+1)
c      epsj(i)=epsj(i)+eps(j+1)*xmult(i,j+1)
c      sgj1(i)=sgj1(i)+sg(j)*xmult(i,j+1)
c      sgj(i)=sgj(i)+sg(j+1)*xmult(i,j+1)
c   90 continue
c  100 continue
c      do 110 i=mft,mlt
c      qhj(i)=epsj(i)-epsj1(i)
c      qhj1(i)=1./qhj(i)
c  110 factor(i)=(epx(i)-epsj1(i))*qhj1(i)
c      do 120 i=mft,mlt
c      qhj(i)=sgj(i)-sgj1(i)
c      ak(i)=sgj1(i)+qhj(i)*factor(i)
c  120 qh(i)=qhj(i)*qhj1(i)
      return
      end
