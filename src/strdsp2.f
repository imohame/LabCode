      subroutine strdsp2 (x11,x12,x13,x14)
c     implicit double precision (a-h,o-z)                                    dp
c
c     evaluate strain-displacement matrix
c
       use mod_parameters
      common/bk02/ioofc,iphase,imass,lpar(9)
      common/intgrt/nintg
      common/bk44/
     1 h1,h2,h3,h4,p11,p21,p12,p22,p13,p23,p14,p24
      common/range/mft,mlt,lft,llt
      dimension xj11(nelemg),xj12(nelemg),xj21(nelemg),
     > xj22(nelemg),vlinv(nelemg),
     1 xji11(nelemg),xji12(nelemg),xji21(nelemg),xji22(nelemg),
     > radius(nelemg),
     2 radinv(nelemg),x1111(nelemg),x1113(nelemg),x2222(nelemg),
     > x2221(nelemg),
     3 x1221(nelemg),x2111(nelemg),x1222(nelemg),x2113(nelemg)
      common/vect4b/
     1py1(nelemg),py2(nelemg),py3(nelemg),py4(nelemg),
     2pz1(nelemg),pz2(nelemg),pz3(nelemg),pz4(nelemg),
     3ph1(nelemg),ph2(nelemg),ph3(nelemg),ph4(nelemg)
      common/vect13/
     1 x1112(nelemg),x1314(nelemg),x1114(nelemg),x1213(nelemg),
     2 x2122(nelemg),x2324(nelemg),x2124(nelemg),x2223(nelemg)
!!      common/kkk/km,kkjj
c
      dimension x11(*),x12(*),x13(*),x14(*),vlinc(nelemg)
      equivalence (lpar(5),ityp2d)
c


c         write(7777,*) '-- strdsp.f'

      call bssf (h1,p11)
c
c     calculates derivative of coordinates wrt local 
      do 10 i=lft,llt
      xj11(i)=p11*x1112(i)+p13*x1314(i)
      xj12(i)=p11*x2122(i)+p13*x2324(i)
      xj21(i)=p21*x1114(i)+p22*x1213(i)
      xj22(i)=p21*x2124(i)+p22*x2223(i)
      vlinc(i)=xj11(i)*xj22(i)-xj21(i)*xj12(i)
      vlinv(i)=1./vlinc(i)
      xji11(i)=xj22(i)*vlinv(i)
      xji12(i)=xj12(i)*vlinv(i)
      xji21(i)=xj21(i)*vlinv(i)
      xji22(i)=xj11(i)*vlinv(i)
      ph1(i)=0.0
      ph2(i)=0.0
      ph3(i)=0.0
      ph4(i)=0.0
   10 continue
c------------------------------------------------------
c------------------------------------------------------
c      if (km.eq.1) then
c      xji12(1)=-1.0*xji12(1)
c      xji21(1)=-1.0*xji21(1)
c      write(25,15) xji11(1),xji12(1),xji21(1),xji22(1)
c   15 format(//5x,'Jacobian matrix :'//5x,f5.2,2x,f5.2//5x,
c     .       f5.2,2x,f5.2//)
c      xji12(1)=-1.0*xji12(1)
c      xji21(1)=-1.0*xji21(1)
c      endif
c-------------------------------------------------------
c-------------------------------------------------------  

c     if only 1 integration point, go to 50 
      if (nintg.eq.5) go to 50
      do 30 i=lft,llt
      x1111(i)=xji11(i)*p11
      x1113(i)=xji11(i)*p13
      x2222(i)=xji22(i)*p22
      x2221(i)=xji22(i)*p21
      x1221(i)=xji12(i)*p21
      x2111(i)=xji21(i)*p11
      x1222(i)=xji12(i)*p22
      x2113(i)=xji21(i)*p13
      py1(i)  =x1111(i)-x1221(i)
      pz1(i)  =x2221(i)-x2111(i)
      py2(i)  =-x1111(i)-x1222(i)
      pz2(i)  =x2222(i)+x2111(i)
      py3(i)  =x1113(i)+x1222(i)
      pz3(i)  =-x2222(i)-x2113(i)
      py4(i)  =-x1113(i)+x1221(i)
      pz4(i)  =-x2221(i)+x2113(i)
   30 continue
c
c     for 2x2 integration and pl. strain or pl. stress, return
      if (ityp2d.gt.0) return
c
      do 40 i=lft,llt
      radius(i)=h1*x11(i)+h2*x12(i)+h3*x13(i)+h4*x14(i)
      radinv(i)=1./radius(i)
      ph1(i)=h1*radinv(i)
      ph2(i)=h2*radinv(i)
      ph3(i)=h3*radinv(i)
      ph4(i)=h4*radinv(i)
      vlinc(i)=radius(i)*vlinc(i)
   40 continue
      return
	  
c     start of derivative of shape fcns with respect to coordinates
c     for single point integration case
   50 do 60 i=lft,llt
      x1111(i)=xji11(i)*p11
      x2222(i)=xji22(i)*p22
      x2221(i)=xji22(i)*p21
      x1221(i)=xji12(i)*p21
      x2111(i)=xji21(i)*p11
      x1222(i)=xji12(i)*p22
      py1(i)  =x1111(i)-x1221(i)
      pz1(i)  =x2221(i)-x2111(i)
      py2(i)  =-x1111(i)-x1222(i)
      pz2(i)  =x2222(i)+x2111(i)
      py3(i)  =-py1(i)
      pz3(i)  =-pz1(i)
      py4(i)  =-py2(i)
      pz4(i)  =-pz2(i)
   60 continue
c
c     for single point integration and pl. strain or pl. stress, return
      if (ityp2d.gt.0) return
c
      do 70 i=lft,llt
      radius(i)=.25*(x11(i)+x12(i)+x13(i)+x14(i))
      radinv(i)=.25/radius(i)
      ph1(i)=radinv(i)
      ph2(i)=radinv(i)
      ph3(i)=radinv(i)
      ph4(i)=radinv(i)
      vlinc(i)=radius(i)*vlinc(i)
   70 continue
      end
