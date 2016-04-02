      subroutine xstran(iopt,nft,nlt)
c     implicit double precision (a-h,o-z)                                    dp
       use mod_parameters
      common/bk01/h4(4,5),p14(4,5),p24(4,5)
      common/bk18/nummat,ityp2d,ako(31)
      common/intgrt/ipt
      common/vect3/
     1 dgi1(nelemg,4),dgi2(nelemg,4),dgi3(nelemg,4),dgi4(nelemg,4),
     > dgi5(nelemg,4),
     2 dgt1(nelemg,4),dgt2(nelemg,4),dgt3(nelemg,4),dgt4(nelemg,4),
     3 f11v(nelemg),f22v(nelemg),f12v(nelemg),f21v(nelemg),
     > dsd5(nelemg),
     4 sig11s(nelemg),sig22s(nelemg),sig33s(nelemg),sig12s(nelemg),
     5 ddp1(nelemg,4),ddp2(nelemg,4),ddp3(nelemg,4),ddp4(nelemg,4),
     > ddp5(nelemg,4)
      common/vect20/
     1 sdgt1(nelemg,4),sdgt2(nelemg,4),sdgt3(nelemg,4),
     2 sdgt4(nelemg,4),sdgt5(nelemg,4),dfdet(nelemg,4)
      common/vect91/xargat(nelemg),exx(nelemg,4),eps(nelemg,5),
     > sa(nelemg),sb(nelemg),sc(nelemg),ang1(nelemg),ang2(nelemg)
      if(iopt.eq.0)then
      do 10 i=nft,nlt
      exx(i,1)=sig11s(i)
      exx(i,2)=sig22s(i)
      exx(i,3)=sig33s(i)
      exx(i,4)=sig12s(i)
   10 continue
      elseif(iopt.eq.1)then
      do 20 i=nft,nlt
      exx(i,1)=sdgt1(i,ipt)-1.
      exx(i,2)=sdgt2(i,ipt)-1.
      exx(i,3)=sdgt5(i,ipt)-1.
      exx(i,4)=.50*(sdgt3(i,ipt)+sdgt4(i,ipt))
   20 continue
      elseif(iopt.eq.2)then
      do 30 i=nft,nlt
      exx(i,1)=.5*(sdgt1(i,ipt)**2+sdgt4(i,ipt)**2-1.)
      exx(i,2)=.5*(sdgt2(i,ipt)**2+sdgt3(i,ipt)**2-1.)
      exx(i,3)=.5*(sdgt5(i,ipt)**2-1.)
      exx(i,4)=.5*(sdgt1(i,ipt)*sdgt3(i,ipt)+
     1              sdgt2(i,ipt)*sdgt4(i,ipt))
   30 continue
      elseif(iopt.eq.3)then
      do 40 i=nft,nlt
      df=1./(sdgt1(i,ipt)*sdgt2(i,ipt)-sdgt3(i,ipt)*sdgt4(i,ipt))
      eps(i,1)=df*sdgt2(i,ipt)
      eps(i,2)=df*sdgt1(i,ipt)
      eps(i,3)=-df*sdgt3(i,ipt)
      eps(i,4)=-df*sdgt4(i,ipt)
      eps(i,5)=1./sdgt5(i,ipt)
   40 continue
      do 50 i=nft,nlt
      exx(i,1)=.5*(1.-eps(i,1)*eps(i,1)-eps(i,4)*eps(i,4))
      exx(i,2)=.5*(1.-eps(i,3)*eps(i,3)-eps(i,2)*eps(i,2))
      exx(i,3)=.5*(1.-eps(i,5)*eps(i,5))
      exx(i,4)=.5*(-eps(i,1)*eps(i,3)-eps(i,4)*eps(i,2))
   50 continue
      elseif(iopt.eq.4)then
      do 60 i=nft,nlt
      exx(i,1)=ddp1(i,ipt)
      exx(i,2)=ddp2(i,ipt)
      exx(i,3)=ddp5(i,ipt)
      exx(i,4)=.50*(ddp3(i,ipt)+ddp4(i,ipt))
   60 continue
      endif
      return
      end
