      subroutine s2mn (prop,propc,sig,ener,thick,ln)
c   Orthotropic Elastic
c     implicit double precision (a-h,o-z)                                    dp
c
       use mod_parameters
      common/bk02/ioofc,iphase,imass,lpar(9)
      common/bk18/nummat,ityp2d,ako(31)
      common/bk44/
     1 h1,h2,h3,h4,p11,p21,p12,p22,p13,p23,p14,p24
      common/intgp/
     1 d11,d21,d31,d41,d12,d22,d32,d42,d13,d23,d33,d43,d14,d24,d34,d44
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
      common/vect5/
     1 gp11(nelemg),gp21(nelemg),gp31(nelemg),gp41(nelemg),
     2 gp12(nelemg),gp22(nelemg),gp32(nelemg),gp42(nelemg),
     3 gp13(nelemg),gp23(nelemg),gp33(nelemg),gp43(nelemg)
      common/vect6/
     1 cg11(nelemg),cg21(nelemg),cg12(nelemg),cg22(nelemg),
     2 cg13(nelemg),cg23(nelemg),cg14(nelemg),cg24(nelemg),
     3 ch11(nelemg),ch21(nelemg),ch12(nelemg),ch22(nelemg),
     4 ch13(nelemg),ch23(nelemg),ch14(nelemg),ch24(nelemg),
     5 ci11(nelemg),ci21(nelemg),ci12(nelemg),ci22(nelemg),
     6 ci13(nelemg),ci23(nelemg),ci14(nelemg),ci24(nelemg)
      common/vect8/
     1 dsave(4,4,nelemg)
      common/vect14/c11(nelemg),c12(nelemg),c13(nelemg),
     > c14(nelemg),c22(nelemg),
     1 c23(nelemg),c24(nelemg),c33(nelemg),c34(nelemg),
     > c44(nelemg),sp1(nelemg),sp2(nelemg)
      common/vect15/
     1 d1(nelemg),d2(nelemg),d3(nelemg),d4(nelemg)
      common/vect16/
     1 u11(nelemg),u12(nelemg),u14(nelemg),u21(nelemg),
     > u22(nelemg),u24(nelemg),
     2 u41(nelemg),u42(nelemg),u44(nelemg),q1(nelemg),
     > q2(nelemg),q3(nelemg),q4(nelemg),
     3 t11(nelemg),t12(nelemg),t14(nelemg),t41(nelemg),
     > t44(nelemg)
      common/vect17/beta(nelemg)
c
      dimension prop(*),propc(4,*),c(4,4),sig(ln,*),ener(ln,*),
     1 thick(ln,*)
      equivalence (d11,c)
c
      call rotat1 (sig,ener,ln)
c
      do 10 i=1,4
      do 10 j=1,4
   10 c(i,j)=propc(i,j)
      iopt=prop(8)+1
      go to (20,40,60), iopt
   20 do 30 i=mft,mlt
   30 q3(i)=atan2(ci22(i)-ci21(i),ci12(i)-ci11(i))+beta(i)
      go to 80
   40 call bssf (h1,p11)
      do 50 i=mft,mlt
      q1(i)=h1*ci11(i)+h2*ci12(i)+h3*ci13(i)+h4*ci14(i)-prop(9)
      q2(i)=h1*ci21(i)+h2*ci22(i)+h3*ci23(i)+h4*ci24(i)-prop(10)
   50 q3(i)=atan2(q2(i),q1(i))
      go to 80
   60 do 70 i=mft,mlt
   70 q3(i)=prop(11)
   80 do 90 i=mft,mlt
      q1(i)=sin(q3(i))
   90 q2(i)=cos(q3(i))
      d1112=d11-d12
      d1222=d12-d22
      d1323=d13-d23
c
c     transformation and constitutive matrix
c
      do 100 i=mft,mlt
      t11(i)=q2(i)**2
      t12(i)=q1(i)**2
      t14(i)=q2(i)*q1(i)
      t41(i)=-t14(i)*2.0
      t44(i)=t11(i)-t12(i)
      u11(i)=t11(i)*d11+t12(i)*d12
      u12(i)=t11(i)*d12+t12(i)*d22
      c13(i)=t11(i)*d13+t12(i)*d23
      u14(i)=t41(i)*d44
      u21(i)=t12(i)*d11+t11(i)*d12
      u22(i)=t12(i)*d12+t11(i)*d22
      c23(i)=t12(i)*d13+t11(i)*d23
      u24(i)=-t41(i)*d44
      c33(i)=d33
      u41(i)=t14(i)*d1112
      u42(i)=t14(i)*d1222
      c34(i)=t14(i)*d1323
  100 u44(i)=t44(i)*d44
c
      do 110 i=mft,mlt
      c11(i)=u11(i)*t11(i)+u12(i)*t12(i)+u14(i)*t41(i)
      c12(i)=u11(i)*t12(i)+u12(i)*t11(i)-u14(i)*t41(i)
      c14(i)=u11(i)*t14(i)-u12(i)*t14(i)+u14(i)*t44(i)
      c22(i)=u21(i)*t12(i)+u22(i)*t11(i)-u24(i)*t41(i)
      c24(i)=u21(i)*t14(i)-u22(i)*t14(i)+u24(i)*t44(i)
  110 c44(i)=u41(i)*t14(i)-u42(i)*t14(i)+u44(i)*t44(i)
c
c     rotated cauchy stress
c
      do 130 i=mft,mlt
      sig(1,i)=sig(1,i)+c11(i)*d1(i)+c12(i)*d2(i)+c13(i)*d3(i)+
     > c14(i)*d4(i)
      sig(2,i)=sig(2,i)+c12(i)*d1(i)+c22(i)*d2(i)+c23(i)*d3(i)+
     > c24(i)*d4(i)
      sig(3,i)=sig(3,i)+c13(i)*d1(i)+c23(i)*d2(i)+c33(i)*d3(i)+
     > c34(i)*d4(i)
  130 sig(4,i)=sig(4,i)+c14(i)*d1(i)+c24(i)*d2(i)+c34(i)*d3(i)+
     > c44(i)*d4(i)
c
      call rotat2 (sig,ener,ln)
c
      do 140 i=mft,mlt
      sig11s(i)=sig(1,i)
      sig22s(i)=sig(2,i)
      sig33s(i)=sig(3,i)
  140 sig12s(i)=sig(4,i)
c
      if (iphase.eq.3) return
c
      do 170 i=mft,mlt
      dsave(1,1,i)=c11(i)
      dsave(1,2,i)=c12(i)
      dsave(1,3,i)=c13(i)
      dsave(1,4,i)=c14(i)
      dsave(2,1,i)=c12(i)
      dsave(2,2,i)=c22(i)
      dsave(2,3,i)=c23(i)
      dsave(2,4,i)=c24(i)
      dsave(3,1,i)=c13(i)
      dsave(3,2,i)=c23(i)
      dsave(3,3,i)=c33(i)
      dsave(3,4,i)=c34(i)
      dsave(4,1,i)=c14(i)
      dsave(4,2,i)=c24(i)
      dsave(4,3,i)=c34(i)
      dsave(4,4,i)=c44(i)
  170 continue
c
c     transform constitutive model to n+1/2 configuration
c
      do 180 i=mft,mlt
  180 call tranfc(i)
c
      return
      end
