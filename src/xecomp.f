      subroutine xecomp(icomp,nft,nlt)
c     implicit double precision (a-h,o-z)                                    dp
       use mod_parameters
      common/intgrt/ipt
      common/bks16/lstep,lmstep,xtime,xmtime,xdelt
      common/vect2/
     1 ed11(nelemg),ed21(nelemg),ed12(nelemg),ed22(nelemg),
     2 ed13(nelemg),ed23(nelemg),ed14(nelemg),ed24(nelemg),
     3 fd11(nelemg),fd21(nelemg),fd12(nelemg),fd22(nelemg),
     4 fd13(nelemg),fd23(nelemg),fd14(nelemg),fd24(nelemg)
      common/vect6/
     1 cg11(nelemg),cg21(nelemg),cg12(nelemg),cg22(nelemg),
     2 cg13(nelemg),cg23(nelemg),cg14(nelemg),cg24(nelemg),
     3 ch11(nelemg),ch21(nelemg),ch12(nelemg),ch22(nelemg),
     4 ch13(nelemg),ch23(nelemg),ch14(nelemg),ch24(nelemg),
     5 ci11(nelemg),ci21(nelemg),ci12(nelemg),ci22(nelemg),
     6 ci13(nelemg),ci23(nelemg),ci14(nelemg),ci24(nelemg)
      common/vect20/
     1 sdgt1(nelemg,4),sdgt2(nelemg,4),sdgt3(nelemg,4),
     2 sdgt4(nelemg,4),sdgt5(nelemg,4),dfdet(nelemg,4)
      common/vect91/xargat(nelemg),exx(nelemg,4),eps(nelemg,5),
     > sa(nelemg),sb(nelemg),sc(nelemg),ang1(nelemg),ang2(nelemg)
      if(icomp-100.lt.0)then
        istrn=0
      elseif(icomp-200.lt.0)then
        istrn=1
      elseif(icomp-300.lt.0)then
        istrn=2
      elseif(icomp-400.lt.0)then
        istrn=3
      elseif(icomp-500.lt.0)then
        istrn=4
      endif
      lcomp=icomp-istrn*100
      if(lcomp.lt.19.or.(lcomp.gt.27.and.lcomp.lt.36))then
      call xstran(istrn,nft,nlt)
      do 20 i=nft,nlt
      sa(i)=0.5*(exx(i,1)+exx(i,2))
      sb(i)=0.5*(exx(i,1)-exx(i,2))
      sc(i)=sqrt(sb(i)**2+exx(i,4)**2)
   20 continue
      endif
      if(lcomp.ge.12.and.lcomp.le.15)then
      do 90 i=nft,nlt
      z1=-cg21(i)+cg22(i)+cg23(i)-cg24(i)
      r1=-cg11(i)+cg12(i)+cg13(i)-cg14(i)
      z2=-cg21(i)-cg22(i)+cg23(i)+cg24(i)
      r2=-cg11(i)-cg12(i)+cg13(i)+cg14(i)
      ang1(i)=2.*atan2(z1,r1+1.e-6)
      ang2(i)=2.*atan2(z2,r2+1.e-6)
   90 continue
      endif

      if(lcomp.le.4)then
      do 10 i=nft,nlt
      xargat(i)=exx(i,lcomp)
   10 continue
      elseif(lcomp.eq.5)then
      do 30 i=nft,nlt
      xargat(i)=sa(i)+sc(i)
   30 continue
      elseif(lcomp.eq.6)then
      do 40 i=nft,nlt
      xargat(i)=sa(i)-sc(i)
   40 continue
      elseif(lcomp.eq.7)then
      fac=1.
      if(icomp-100.gt.0)fac=2./3.
      do 50 i=nft,nlt
      s5=sa(i)+sc(i)
      s6=sa(i)-sc(i)
      xargat(i)=fac*sqrt(.5*((s5-s6)**2+(s5-exx(i,3))**2+
     1 (s6-exx(i,3))**2))
   50 continue
      elseif(lcomp.eq.8)then
      fac=1.
      if(icomp-100.lt.0)fac=-1./3.
      do 60 i=nft,nlt
      xargat(i)=fac*(exx(i,1)+exx(i,2)+exx(i,3))
   60 continue
      elseif(lcomp.eq.9)then
      do 70 i=nft,nlt
      xargat(i)=2.*sc(i)
   70 continue
      elseif(lcomp.eq.10)then
      do 80 i=nft,nlt
      xargat(i)=exx(i,1)-exx(i,3)
   80 continue
      elseif(lcomp.eq.11)then
      do 110 i=nft,nlt
      xargat(i)=sc(i)
  110 continue
      elseif(lcomp.eq.12)then
      do 120 i=nft,nlt
      xargat(i)=sa(i)-sb(i)*cos(ang1(i))-exx(i,4)*sin(ang1(i))
  120 continue
      elseif(lcomp.eq.13)then
      do 130 i=nft,nlt
      xargat(i)=sa(i)-sb(i)*cos(ang2(i))-exx(i,4)*sin(ang2(i))
  130 continue
      elseif(lcomp.eq.14)then
      do 140 i=nft,nlt
      xargat(i)=-sb(i)*sin(ang1(i))+exx(i,4)*cos(ang1(i))
  140 continue
      elseif(lcomp.eq.15)then
      do 150 i=nft,nlt
      xargat(i)=-sb(i)*sin(ang2(i))+exx(i,4)*cos(ang2(i))
  150 continue
      elseif(lcomp.eq.16)then
      do 160 i=nft,nlt
      xargat(i)=sa(i)+sc(i)-(exx(i,1)+exx(i,2)+exx(i,3))/3.
  160 continue
      elseif(lcomp.eq.17)then
      do 170 i=nft,nlt
      xargat(i)=sa(i)-sc(i)-(exx(i,1)+exx(i,2)+exx(i,3))/3.
  170 continue
      elseif(lcomp.eq.18)then
      do 180 i=nft,nlt
      xargat(i)=exx(i,3)-(exx(i,1)+exx(i,2)+exx(i,3))/3.
  180 continue
      elseif(lcomp.eq.19)then
      call xghist(19,nft,nlt)
      elseif(lcomp.eq.20)then
      call xghist(20,nft,nlt)
      elseif(lcomp.eq.21)then
      do 280 i=nft,nlt
      xargat(i)=log(dfdet(i,ipt))
  280 continue
      elseif(lcomp.eq.22)then
      do 290 i=nft,nlt
      xargat(i)=.25*(ed11(i)+ed12(i)+ed13(i)+ed14(i))
  290 continue
      elseif(lcomp.eq.23)then
      do 300 i=nft,nlt
      xargat(i)=.25*(ed21(i)+ed22(i)+ed23(i)+ed24(i))
  300 continue
      elseif(lcomp.eq.24)then
      do 310 i=nft,nlt
      fac1=.25*(ed11(i)+ed12(i)+ed13(i)+ed14(i))
      fac2=.25*(ed21(i)+ed22(i)+ed23(i)+ed24(i))
      xargat(i)=sqrt(fac1**2+fac2**2)
  310 continue
      elseif(lcomp.eq.25)then
      do 320 i=nft,nlt
      xargat(i)=.25*(fd11(i)+fd12(i)+fd13(i)+fd14(i))/xdelt
  320 continue
      elseif(lcomp.eq.26)then
      do 330 i=nft,nlt
      xargat(i)=.25*(fd21(i)+fd22(i)+fd23(i)+fd24(i))/xdelt
  330 continue
      elseif(lcomp.eq.27)then
      do 340 i=nft,nlt
      fac1=.25*(fd11(i)+fd12(i)+fd13(i)+fd14(i))/xdelt
      fac2=.25*(fd21(i)+fd22(i)+fd23(i)+fd24(i))/xdelt
      xargat(i)=sqrt(fac1**2+fac2**2)
  340 continue
      elseif(lcomp.eq.28)then
      do 190 i=nft,nlt
      ang3=2.*atan2(cg22(i)-cg21(i),cg12(i)-cg11(i)+1.e-6)
      xargat(i)=sa(i)-sb(i)*cos(ang3)-exx(i,4)*sin(ang3)
  190 continue
      elseif(lcomp.eq.29)then
      do 200 i=nft,nlt
      ang3=2.*atan2(cg23(i)-cg22(i),cg13(i)-cg12(i)+1.e-6)
      xargat(i)=sa(i)-sb(i)*cos(ang3)-exx(i,4)*sin(ang3)
  200 continue
      elseif(lcomp.eq.30)then
      do 210 i=nft,nlt
      ang3=2.*atan2(cg24(i)-cg23(i),cg14(i)-cg13(i)+1.e-6)
      xargat(i)=sa(i)-sb(i)*cos(ang3)-exx(i,4)*sin(ang3)
  210 continue
      elseif(lcomp.eq.31)then
      do 220 i=nft,nlt
      ang3=2.*atan2(cg21(i)-cg24(i),cg11(i)-cg14(i)+1.e-6)
      xargat(i)=sa(i)-sb(i)*cos(ang3)-exx(i,4)*sin(ang3)
  220 continue
      elseif(lcomp.eq.32)then
      do 230 i=nft,nlt
      ang3=2.*atan2(cg22(i)-cg21(i),cg12(i)-cg11(i)+1.e-6)
      xargat(i)=-sb(i)*sin(ang3)+exx(i,4)*cos(ang3)
  230 continue
      elseif(lcomp.eq.33)then
      do 240 i=nft,nlt
      ang3=2.*atan2(cg23(i)-cg22(i),cg13(i)-cg12(i)+1.e-6)
      xargat(i)=-sb(i)*sin(ang3)+exx(i,4)*cos(ang3)
  240 continue
      elseif(lcomp.eq.34)then
      do 250 i=nft,nlt
      ang3=2.*atan2(cg24(i)-cg23(i),cg14(i)-cg13(i)+1.e-6)
      xargat(i)=-sb(i)*sin(ang3)+exx(i,4)*cos(ang3)
  250 continue
      elseif(lcomp.eq.35)then
      do 260 i=nft,nlt
      ang3=2.*atan2(cg21(i)-cg24(i),cg11(i)-cg14(i)+1.e-6)
      xargat(i)=-sb(i)*sin(ang3)+exx(i,4)*cos(ang3)
  260 continue
      elseif(lcomp.eq.36)then
      do 350 i=nft,nlt
      xargat(i)=dfdet(i,ipt)
  350 continue
      elseif(lcomp.eq.37)then
      do 360 i=nft,nlt
      xargat(i)=1./dfdet(i,ipt)-1.
  360 continue
      elseif(lcomp.eq.38)then
      do 370 i=nft,nlt
      xargat(i)=0.
  370 continue
      elseif(lcomp.eq.39)then
      do 380 i=nft,nlt
      xargat(i)=-(exx(i,1)+exx(i,2)+exx(i,3))/3.
  380 continue
      elseif(lcomp.eq.40)then
      call xghist(40,nft,nlt)
      do 390 i=nft,nlt
      xargat(i)=xargat(i)/dfdet(i,ipt)
  390 continue
      elseif(lcomp.gt.40)then
      call xghist(lcomp,nft,nlt)
      endif
      return
      end
