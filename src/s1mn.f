      subroutine s1mn (prop,sig,ener,thick,ln)
c     implicit double precision (a-h,o-z)                                    dp
c
       use mod_parameters
      common/bk02/ioofc,iphase,imass,lpar(9)
      common/bk08/kprint,nstep,ite,ilimit,newstf
      common/bk18/nummat,ityp2d,ako(31)
      common/bk32/nsref,nequit,time,timep,lprint,nprint
      common/intgp/
     1 c11,c21,c31,c41,c12,c22,c32,c42,c13,c23,c33,c43,c14,
     > c24,c34,c44
      common/range/mft,mlt,lft,llt,nftm1
      common/hokao/lst
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
c
      dimension prop(4,*),sig(ln,*),ener(ln,*),thick(ln,*),c(4,4)
      equivalence (c11,c)
c
      do 10 i=1,4
      do 10 j=1,4
   10 c(i,j)=prop(i,j)
c   
c    The stress is updated with the Green-Naghdi rate. 
c     this cab changed for other objective rates using a similar 
c     analysis as the one for crystal plasticity
      call rotat1 (sig,ener,ln)
c
      do 30 i=mft,mlt
      sig(1,i)=sig(1,i)+c11*d1(i)+c12*d2(i)+c13*d3(i)
      sig(2,i)=sig(2,i)+c21*d1(i)+c22*d2(i)+c23*d3(i)
      sig(3,i)=sig(3,i)+c31*d1(i)+c32*d2(i)+c33*d3(i)
   30 sig(4,i)=sig(4,i)+c44*d4(i)
c
      call rotat2 (sig,ener,ln)
c
      do 40 i=mft,mlt
      sig11s(i)=sig(1,i)
      sig22s(i)=sig(2,i)
      sig33s(i)=sig(3,i)
   40 sig12s(i)=sig(4,i)
c
      if (iphase.eq.3) return
c
      do 80 i=1,4
      do 80 j=1,4
      do 70 l=mft,mlt
   70 dsave(i,j,l)=c(i,j)
   80 continue
c
      return
      end
