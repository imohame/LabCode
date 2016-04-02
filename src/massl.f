c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine massl (xm,xx,rho,thick)
c     implicit double precision (a-h,o-z)                                    dp
c
c     compute lumped mass matrix
c
      common/bk02/ioofc,iphase,imass,lpar(9)
      common/bk44/h(8),p(2,8),xj(2,2)
      common/bk48/stress(4),strain(4),d(4,4),lst,n,nstate
c
      dimension xm(*),xx(2,8)
      equivalence (lpar(5),ityp2d)
c
      do 10 i=1,4
   10 xm(i)=0.0
      do 60 lst=1,4
      call basis2 (h,p,xj,det,xx)
      if (ityp2d.eq.0) go to 20
      radius=1.0
      go to 40
   20 radius=0.0
      do 30 k=1,4
   30 radius=radius+h(k)*xx(1,k)
   40 fac=radius*det*rho
      do 50 i=1,4
   50 xm(i)=xm(i)+fac*h(i)
   60 continue
c
      return
      end
