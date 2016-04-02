      subroutine s0mn (sig,ener,thick,ln)
c     implicit double precision (a-h,o-z)                                    dp
c
       use mod_parameters
      common/bk02/ioofc,iphase,imass,lpar(9)
      common/bk18/nummat,ityp2d,ako(31)
      common/bk49/bulkmx,ncon(30)
      common/intgp/
     1 c11,c21,c31,c41,c12,c22,c32,c42,c13,c23,c33,c43,c14,
     > c24,c34,c44
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
      real k
c
      dimension sig(ln,*),ener(ln,*),thick(ln,*),c(4,4)
      equivalence (c11,c)
c
      e=.001*bulkmx
      g=e/2.6
      k=e/1.2
      c11=k+4.*g/3.
      c12=k-2.*g/3.
      c13=c12
      c14=0.0
      c21=c12
      c22=c11
      c23=c21
      c24=0.0
      c31=c13
      c32=c23
      c33=c11
      c34=0.0
      c41=0.0
      c42=0.0
      c43=0.0
      c44=g
      do 30 i=mft,mlt
      sig(1,i)=0.0
      sig(2,i)=0.0
      sig(3,i)=0.0
   30 sig(4,i)=0.0
c
      do 40 i=mft,mlt
      sig11s(i)=0.0
      sig22s(i)=0.0
      sig33s(i)=0.0
   40 sig12s(i)=0.0
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
