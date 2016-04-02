      subroutine nonlin(s,ityp2d)
c     implicit double precision (a-h,o-z)                                    dp
       use mod_parameters
      common/range/mft,mlt,lft,llt,nftm1
      common/vect4/
     1 py1(nelemg),py2(nelemg),py3(nelemg),py4(nelemg),
     2 pz1(nelemg),pz2(nelemg),pz3(nelemg),pz4(nelemg),
     3 ph1(nelemg),ph2(nelemg),ph3(nelemg),ph4(nelemg)
      common/vect5/
     1 fac1(nelemg),fac2(nelemg),fac3(nelemg),fac4(nelemg),
     > cb1(nelemg),cb2(nelemg),
     2 cb31(nelemg),cb32(nelemg),cb33(nelemg)
      common/vect7/
     1 diavg(nelemg),dfavg(nelemg),volnd(nelemg),volcd(nelemg),
     > vlinc(nelemg),volgp(nelemg),
     2 dimod(nelemg),sigs11(nelemg),sigs22(nelemg),sigs33(nelemg)
     > ,sigs12(nelemg)
      common/vect19/
     1 sdy1(nelemg),sdy2(nelemg),sdy3(nelemg),sdy4(nelemg),
     2 sdz1(nelemg),sdz2(nelemg),sdz3(nelemg),sdz4(nelemg)
      common/vect21/sig11(nelemg),sig22(nelemg),sig33(nelemg),
     > sig12(nelemg)
      dimension s(44,1)
      
!$OMP PARALLEL DO       
      do 10 i=lft,llt
      sig11(i)=volgp(i)*sig11(i)
      sig22(i)=volgp(i)*sig22(i)
      sig33(i)=volgp(i)*sig33(i)
      sig12(i)=volgp(i)*sig12(i)
      sdy1(i) =py1(i)*sig11(i)+pz1(i)*sig12(i)
      sdz1(i) =pz1(i)*sig22(i)+py1(i)*sig12(i)
      sdy2(i) =py2(i)*sig11(i)+pz2(i)*sig12(i)
      sdz2(i) =pz2(i)*sig22(i)+py2(i)*sig12(i)
      sdy3(i) =py3(i)*sig11(i)+pz3(i)*sig12(i)
      sdz3(i) =pz3(i)*sig22(i)+py3(i)*sig12(i)
      sdy4(i) =py4(i)*sig11(i)+pz4(i)*sig12(i)
      sdz4(i) =pz4(i)*sig22(i)+py4(i)*sig12(i)
!   10 continue
!      do 20 i=lft,llt
      fac1(i)=py1(i)*sdy1(i)+pz1(i)*sdz1(i)
      fac2(i)=py2(i)*sdy1(i)+pz2(i)*sdz1(i)
      fac3(i)=py3(i)*sdy1(i)+pz3(i)*sdz1(i)
      fac4(i)=py4(i)*sdy1(i)+pz4(i)*sdz1(i)
      s(1,i)=s(1,i)+fac1(i)
      s(3,i)=s(3,i)+fac1(i)
      s(4,i)=s(4,i)+fac2(i)
      s(8,i)=s(8,i)+fac2(i)
      s(11,i)=s(11,i)+fac3(i)
      s(17,i)=s(17,i)+fac3(i)
      s(22,i)=s(22,i)+fac4(i)
      s(30,i)=s(30,i)+fac4(i)
!   20 continue
!      do 30 i=lft,llt
      fac2(i)=py2(i)*sdy2(i)+pz2(i)*sdz2(i)
      fac3(i)=py3(i)*sdy2(i)+pz3(i)*sdz2(i)
      fac4(i)=py4(i)*sdy2(i)+pz4(i)*sdz2(i)
      s(6,i)=s(6,i)+fac2(i)
      s(10,i)=s(10,i)+fac2(i)
      s(13,i)=s(13,i)+fac3(i)
      s(19,i)=s(19,i)+fac3(i)
      s(24,i)=s(24,i)+fac4(i)
      s(32,i)=s(32,i)+fac4(i)
!   30 continue
!      do 40 i=lft,llt
      fac3(i)=py3(i)*sdy3(i)+pz3(i)*sdz3(i)
      fac4(i)=py4(i)*sdy3(i)+pz4(i)*sdz3(i)
      s(15,i)=s(15,i)+fac3(i)
      s(21,i)=s(21,i)+fac3(i)
      s(26,i)=s(26,i)+fac4(i)
      s(34,i)=s(34,i)+fac4(i)
!   40 continue
!      do 50 i=lft,llt
      fac4(i)=py4(i)*sdy4(i)+pz4(i)*sdz4(i)
      s(28,i)=s(28,i)+fac4(i)
      s(36,i)=s(36,i)+fac4(i)
!   50 continue
   10 continue
c
      if (ityp2d.ne.0) return
c
!$OMP PARALLEL DO       
      do 60 i=lft,llt
      cb31(i)=sig33(i)*ph1(i)
      cb32(i)=sig33(i)*ph2(i)
      cb33(i)=sig33(i)*ph3(i)
      s(1,i)=s(1,i)-cb31(i)*ph1(i)
      s(4,i)=s(4,i)-cb31(i)*ph2(i)
      s(6,i)=s(6,i)-cb32(i)*ph2(i)
      s(11,i)=s(11,i)-cb31(i)*ph3(i)
      s(13,i)=s(13,i)-cb32(i)*ph3(i)
      s(15,i)=s(15,i)-cb33(i)*ph3(i)
      s(22,i)=s(22,i)-cb31(i)*ph4(i)
      s(24,i)=s(24,i)-cb32(i)*ph4(i)
      s(26,i)=s(26,i)-cb33(i)*ph4(i)
      s(28,i)=s(28,i)-sig33(i)*ph4(i)*ph4(i)
   60 continue
      return
      end
