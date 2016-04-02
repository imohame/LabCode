      subroutine raydmp(raystf)
c     implicit double precision (a-h,o-z)                                    dp
c
       use mod_parameters
      common/bk02/ioofc,iphase,imass,lpar(9)
      common/bk18/nummat,ityp2d,ako(31)
      common/bk26/dt,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
      common/intgp/
     1 c11,c21,c31,c41,c12,c22,c32,c42,c13,c23,c33,c43,
     > c14,c24,c34,c44
      common/range/mft,mlt,lft,llt,nftm1
      common/vect3/
     1 dgi1(nelemg,4),dgi2(nelemg,4),dgi3(nelemg,4),dgi4(nelemg,4),
     > dgi5(nelemg,4),
     2 dgt1(nelemg,4),dgt2(nelemg,4),dgt3(nelemg,4),dgt4(nelemg,4),
     3 f11v(nelemg),f22v(nelemg),f12v(nelemg),f21v(nelemg),
     > dsd5(nelemg),
     4 sig11s(nelemg),sig22s(nelemg),sig33s(nelemg),sig12s(nelemg),
     5 ddp1(nelemg,4),ddp2(nelemg,4),ddp3(nelemg,4),ddp4(nelemg,4),
     > ddp5(nelemg,4)
      common/vect8/
     1 dsave(4,4,nelemg)
      common/vect15/
     1 d1(nelemg),d2(nelemg),d3(nelemg),d4(nelemg)
      dimension raystf(*)
c
      if(raystf(1).eq.0.)return
      if(ityp2d.le.1)then
      fac1=raystf(1)/dt
      fac2=1./((1.+raystf(3))*(1.-2.*raystf(3)))
      c1=raystf(2)/2.*fac2*fac1
      c2=raystf(2)*raystf(3)*fac2*fac1
      c3=raystf(2)/(2.*(1.+raystf(3)))*fac1
      do 10 i=mft,mlt
      sig11s(i)=sig11s(i)+c1*d1(i)+c2*(d2(i)+d3(i))
      sig22s(i)=sig22s(i)+c1*d2(i)+c2*(d1(i)+d3(i))
      sig33s(i)=sig33s(i)+c1*d3(i)+c2*(d1(i)+d2(i))
      sig12s(i)=sig12s(i)+c3*d4(i)
      dsave(1,1,i)=dsave(1,1,i)+c1
      dsave(1,2,i)=dsave(1,2,i)+c2
      dsave(1,3,i)=dsave(1,3,i)+c2
      dsave(2,1,i)=dsave(2,1,i)+c2
      dsave(2,2,i)=dsave(2,2,i)+c1
      dsave(2,3,i)=dsave(2,3,i)+c2
      dsave(3,1,i)=dsave(3,1,i)+c2
      dsave(3,2,i)=dsave(3,2,i)+c2
      dsave(3,3,i)=dsave(3,3,i)+c1
      dsave(4,4,i)=dsave(4,4,i)+c3
   10 continue
      else
      fac1=raystf(1)/dt
      c1=raystf(2)/(1.-raystf(3)**2)*fac1
      c2=c1*raystf(3)*fac1
      c3=raystf(2)/(2.*(1.+raystf(3)))*fac1
      do 20 i=mft,mlt
      sig11s(i)=sig11s(i)+c1*d1(i)+c2*d2(i)
      sig22s(i)=sig22s(i)+c2*d1(i)+c1*d2(i)
      sig12s(i)=sig12s(i)+c3*d4(i)
      dsave(1,1,i)=dsave(1,1,i)+c1
      dsave(1,2,i)=dsave(1,2,i)+c2
      dsave(2,1,i)=dsave(2,1,i)+c2
      dsave(2,2,i)=dsave(2,2,i)+c1
      dsave(4,4,i)=dsave(4,4,i)+c3
   20 continue
      endif
c
      return
      end
