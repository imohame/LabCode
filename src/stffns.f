      subroutine stffns(s,d)
!c     implicit double precision (a-h,o-z)                                    
       use mod_parameters
      common/range/mft,mlt,lft,llt,nftm1
      common/intgrt/nintg
      common/vect4/
     1 py1(nelemg),py2(nelemg),py3(nelemg),py4(nelemg),
     2 pz1(nelemg),pz2(nelemg),pz3(nelemg),pz4(nelemg),
     3 ph1(nelemg),ph2(nelemg),ph3(nelemg),ph4(nelemg)
      common/vect5/
     1 cb11(nelemg),cb21(nelemg),cb31(nelemg),cb41(nelemg),
     2 cb12(nelemg),cb22(nelemg),cb32(nelemg),cb42(nelemg),
     3 cb13(nelemg),cb23(nelemg),cb33(nelemg),cb43(nelemg),
     4 cb14(nelemg),cb24(nelemg),cb34(nelemg),cb44(nelemg),
     5 cb15(nelemg),cb25(nelemg),cb35(nelemg),cb45(nelemg),
     6 cb16(nelemg),cb26(nelemg),cb36(nelemg),cb46(nelemg),
     7 cb17(nelemg),cb27(nelemg),cb37(nelemg),cb47(nelemg),
     8 cb18(nelemg),cb28(nelemg),cb38(nelemg),cb48(nelemg)
      common/vect10/volume(nelemg),
     1 pavgy1(nelemg),pavgy2(nelemg),pavgy3(nelemg),pavgy4(nelemg),
     2 pavgz1(nelemg),pavgz2(nelemg),pavgz3(nelemg),pavgz4(nelemg),
     3 pavgh1(nelemg),pavgh2(nelemg),pavgh3(nelemg),pavgh4(nelemg)
      common/vect11/row1(nelemg),row2(nelemg),row3(nelemg),
     > row4(nelemg),
     1 bmy1(nelemg),bmy2(nelemg),bmy3(nelemg),bmy4(nelemg),
     2 bmz1(nelemg),bmz2(nelemg),bmz3(nelemg),bmz4(nelemg),
     3 bmh1(nelemg),bmh2(nelemg),bmh3(nelemg),bmh4(nelemg)
      dimension s(44,1),d(4,4,1)


!        write(7777,*) '-- stffns.f'

      third=-1./3.
      do 10 i=lft,llt
      row1(i)=d(1,1,i)+d(1,2,i)+d(1,3,i)
      row2(i)=d(2,1,i)+d(2,2,i)+d(2,3,i)
      row3(i)=d(3,1,i)+d(3,2,i)+d(3,3,i)
      row4(i)=d(4,1,i)+d(4,2,i)+d(4,3,i)
!   10 continue
!      do 20 i=lft,llt
      cb11(i)=row1(i)*bmy1(i)+d(1,3,i)*ph1(i)+d(1,1,i)*py1(i)+d(1,4,i)
     1 *pz1(i)
      cb21(i)=row2(i)*bmy1(i)+d(2,3,i)*ph1(i)+d(2,1,i)*py1(i)+d(2,4,i)
     1 *pz1(i)
      cb31(i)=row3(i)*bmy1(i)+d(3,3,i)*ph1(i)+d(3,1,i)*py1(i)+d(3,4,i)
     1 *pz1(i)
      cb41(i)=row4(i)*bmy1(i)+d(4,3,i)*ph1(i)+d(4,1,i)*py1(i)+d(4,4,i)
     1 *pz1(i)
      cb12(i)=row1(i)*bmz1(i)+d(1,2,i)*pz1(i)+d(1,4,i)*py1(i)
      cb22(i)=row2(i)*bmz1(i)+d(2,2,i)*pz1(i)+d(2,4,i)*py1(i)
      cb32(i)=row3(i)*bmz1(i)+d(3,2,i)*pz1(i)+d(3,4,i)*py1(i)
      cb42(i)=row4(i)*bmz1(i)+d(4,2,i)*pz1(i)+d(4,4,i)*py1(i)
!   20 continue
!      do 30 i=lft,llt
      cb13(i)=row1(i)*bmy2(i)+d(1,3,i)*ph2(i)+d(1,1,i)*py2(i)+d(1,4,i)
     1 *pz2(i)
      cb23(i)=row2(i)*bmy2(i)+d(2,3,i)*ph2(i)+d(2,1,i)*py2(i)+d(2,4,i)
     1 *pz2(i)
      cb33(i)=row3(i)*bmy2(i)+d(3,3,i)*ph2(i)+d(3,1,i)*py2(i)+d(3,4,i)
     1 *pz2(i)
      cb43(i)=row4(i)*bmy2(i)+d(4,3,i)*ph2(i)+d(4,1,i)*py2(i)+d(4,4,i)
     1 *pz2(i)
      cb14(i)=row1(i)*bmz2(i)+d(1,2,i)*pz2(i)+d(1,4,i)*py2(i)
      cb24(i)=row2(i)*bmz2(i)+d(2,2,i)*pz2(i)+d(2,4,i)*py2(i)
      cb34(i)=row3(i)*bmz2(i)+d(3,2,i)*pz2(i)+d(3,4,i)*py2(i)
      cb44(i)=row4(i)*bmz2(i)+d(4,2,i)*pz2(i)+d(4,4,i)*py2(i)
!   30 continue
!      do 40 i=lft,llt
      cb15(i)=row1(i)*bmy3(i)+d(1,3,i)*ph3(i)+d(1,1,i)*py3(i)+d(1,4,i)
     1 *pz3(i)
      cb25(i)=row2(i)*bmy3(i)+d(2,3,i)*ph3(i)+d(2,1,i)*py3(i)+d(2,4,i)
     1 *pz3(i)
      cb35(i)=row3(i)*bmy3(i)+d(3,3,i)*ph3(i)+d(3,1,i)*py3(i)+d(3,4,i)
     1 *pz3(i)
      cb45(i)=row4(i)*bmy3(i)+d(4,3,i)*ph3(i)+d(4,1,i)*py3(i)+d(4,4,i)
     1 *pz3(i)
      cb16(i)=row1(i)*bmz3(i)+d(1,2,i)*pz3(i)+d(1,4,i)*py3(i)
      cb26(i)=row2(i)*bmz3(i)+d(2,2,i)*pz3(i)+d(2,4,i)*py3(i)
      cb36(i)=row3(i)*bmz3(i)+d(3,2,i)*pz3(i)+d(3,4,i)*py3(i)
      cb46(i)=row4(i)*bmz3(i)+d(4,2,i)*pz3(i)+d(4,4,i)*py3(i)
!   40 continue
!      do 50 i=lft,llt
      cb17(i)=row1(i)*bmy4(i)+d(1,3,i)*ph4(i)+d(1,1,i)*py4(i)+d(1,4,i)
     1 *pz4(i)
      cb27(i)=row2(i)*bmy4(i)+d(2,3,i)*ph4(i)+d(2,1,i)*py4(i)+d(2,4,i)
     1 *pz4(i)
      cb37(i)=row3(i)*bmy4(i)+d(3,3,i)*ph4(i)+d(3,1,i)*py4(i)+d(3,4,i)
     1 *pz4(i)
      cb47(i)=row4(i)*bmy4(i)+d(4,3,i)*ph4(i)+d(4,1,i)*py4(i)+d(4,4,i)
     1 *pz4(i)
      cb18(i)=row1(i)*bmz4(i)+d(1,2,i)*pz4(i)+d(1,4,i)*py4(i)
      cb28(i)=row2(i)*bmz4(i)+d(2,2,i)*pz4(i)+d(2,4,i)*py4(i)
      cb38(i)=row3(i)*bmz4(i)+d(3,2,i)*pz4(i)+d(3,4,i)*py4(i)
      cb48(i)=row4(i)*bmz4(i)+d(4,2,i)*pz4(i)+d(4,4,i)*py4(i)
!   50 continue
!      do 60 i=lft,llt
      s(1,i)=s(1,i)+py1(i)*cb11(i)+pz1(i)*cb41(i)+ph1(i)*cb31(i)
      s(2,i)=s(2,i)+pz1(i)*cb21(i)+py1(i)*cb41(i)
      s(3,i)=s(3,i)+pz1(i)*cb22(i)+py1(i)*cb42(i)
      s(4,i)=s(4,i)+py2(i)*cb11(i)+pz2(i)*cb41(i)+ph2(i)*cb31(i)
      s(5,i)=s(5,i)+py2(i)*cb12(i)+pz2(i)*cb42(i)+ph2(i)*cb32(i)
      s(6,i)=s(6,i)+py2(i)*cb13(i)+pz2(i)*cb43(i)+ph2(i)*cb33(i)
      s(7,i)=s(7,i)+pz2(i)*cb21(i)+py2(i)*cb41(i)
      s(8,i)=s(8,i)+pz2(i)*cb22(i)+py2(i)*cb42(i)
      s(9,i)=s(9,i)+pz2(i)*cb23(i)+py2(i)*cb43(i)
      s(10,i)=s(10,i)+pz2(i)*cb24(i)+py2(i)*cb44(i)
!   60 continue
!      do 70 i=lft,llt
      s(11,i)=s(11,i)+py3(i)*cb11(i)+pz3(i)*cb41(i)+ph3(i)*cb31(i)
      s(12,i)=s(12,i)+py3(i)*cb12(i)+pz3(i)*cb42(i)+ph3(i)*cb32(i)
      s(13,i)=s(13,i)+py3(i)*cb13(i)+pz3(i)*cb43(i)+ph3(i)*cb33(i)
      s(14,i)=s(14,i)+py3(i)*cb14(i)+pz3(i)*cb44(i)+ph3(i)*cb34(i)
      s(15,i)=s(15,i)+py3(i)*cb15(i)+pz3(i)*cb45(i)+ph3(i)*cb35(i)
      s(16,i)=s(16,i)+pz3(i)*cb21(i)+py3(i)*cb41(i)
      s(17,i)=s(17,i)+pz3(i)*cb22(i)+py3(i)*cb42(i)
      s(18,i)=s(18,i)+pz3(i)*cb23(i)+py3(i)*cb43(i)
      s(19,i)=s(19,i)+pz3(i)*cb24(i)+py3(i)*cb44(i)
      s(20,i)=s(20,i)+pz3(i)*cb25(i)+py3(i)*cb45(i)
!   70 continue
!      do 80 i=lft,llt
      s(21,i)=s(21,i)+pz3(i)*cb26(i)+py3(i)*cb46(i)
      s(22,i)=s(22,i)+py4(i)*cb11(i)+pz4(i)*cb41(i)+ph4(i)*cb31(i)
      s(23,i)=s(23,i)+py4(i)*cb12(i)+pz4(i)*cb42(i)+ph4(i)*cb32(i)
      s(24,i)=s(24,i)+py4(i)*cb13(i)+pz4(i)*cb43(i)+ph4(i)*cb33(i)
      s(25,i)=s(25,i)+py4(i)*cb14(i)+pz4(i)*cb44(i)+ph4(i)*cb34(i)
      s(26,i)=s(26,i)+py4(i)*cb15(i)+pz4(i)*cb45(i)+ph4(i)*cb35(i)
      s(27,i)=s(27,i)+py4(i)*cb16(i)+pz4(i)*cb46(i)+ph4(i)*cb36(i)
      s(28,i)=s(28,i)+py4(i)*cb17(i)+pz4(i)*cb47(i)+ph4(i)*cb37(i)
      s(29,i)=s(29,i)+pz4(i)*cb21(i)+py4(i)*cb41(i)
      s(30,i)=s(30,i)+pz4(i)*cb22(i)+py4(i)*cb42(i)
!   80 continue
!      do 90 i=lft,llt
      s(31,i)=s(31,i)+pz4(i)*cb23(i)+py4(i)*cb43(i)
      s(32,i)=s(32,i)+pz4(i)*cb24(i)+py4(i)*cb44(i)
      s(33,i)=s(33,i)+pz4(i)*cb25(i)+py4(i)*cb45(i)
      s(34,i)=s(34,i)+pz4(i)*cb26(i)+py4(i)*cb46(i)
      s(35,i)=s(35,i)+pz4(i)*cb27(i)+py4(i)*cb47(i)
      s(36,i)=s(36,i)+pz4(i)*cb28(i)+py4(i)*cb48(i)
!   90 continue
!   
!      do 100 i=lft,llt
      cb11(i)=cb11(i)+cb21(i)+cb31(i)
      cb12(i)=cb12(i)+cb22(i)+cb32(i)
      cb13(i)=cb13(i)+cb23(i)+cb33(i)
      cb14(i)=cb14(i)+cb24(i)+cb34(i)
      cb15(i)=cb15(i)+cb25(i)+cb35(i)
      cb16(i)=cb16(i)+cb26(i)+cb36(i)
      cb17(i)=cb17(i)+cb27(i)+cb37(i)
      cb18(i)=cb18(i)+cb28(i)+cb38(i)
!  100 continue
!  
!      do 999 i=lft,llt
!      do 110 i=lft,llt
      s(1,i)=s(1,i)+bmy1(i)*cb11(i)
      s(2,i)=s(2,i)+bmy1(i)*cb12(i)
      s(3,i)=s(3,i)+bmz1(i)*cb12(i)
      s(4,i)=s(4,i)+bmy1(i)*cb13(i)
      s(5,i)=s(5,i)+bmz1(i)*cb13(i)
      s(6,i)=s(6,i)+bmy2(i)*cb13(i)
      s(7,i)=s(7,i)+bmy1(i)*cb14(i)
      s(8,i)=s(8,i)+bmz1(i)*cb14(i)
      s(9,i)=s(9,i)+bmy2(i)*cb14(i)
      s(10,i)=s(10,i)+bmz2(i)*cb14(i)
!  110 continue
!      do 120 i=lft,llt
      s(11,i)=s(11,i)+bmy1(i)*cb15(i)
      s(12,i)=s(12,i)+bmz1(i)*cb15(i)
      s(13,i)=s(13,i)+bmy2(i)*cb15(i)
      s(14,i)=s(14,i)+bmz2(i)*cb15(i)
      s(15,i)=s(15,i)+bmy3(i)*cb15(i)
      s(16,i)=s(16,i)+bmy1(i)*cb16(i)
      s(17,i)=s(17,i)+bmz1(i)*cb16(i)
      s(18,i)=s(18,i)+bmy2(i)*cb16(i)
      s(19,i)=s(19,i)+bmz2(i)*cb16(i)
      s(20,i)=s(20,i)+bmy3(i)*cb16(i)
      s(21,i)=s(21,i)+bmz3(i)*cb16(i)
!  120 continue
!      do 130 i=lft,llt
      s(22,i)=s(22,i)+bmy1(i)*cb17(i)
      s(23,i)=s(23,i)+bmz1(i)*cb17(i)
      s(24,i)=s(24,i)+bmy2(i)*cb17(i)
      s(25,i)=s(25,i)+bmz2(i)*cb17(i)
      s(26,i)=s(26,i)+bmy3(i)*cb17(i)
      s(27,i)=s(27,i)+bmz3(i)*cb17(i)
      s(28,i)=s(28,i)+bmy4(i)*cb17(i)
!  130 continue
!      do 140 i=lft,llt
      s(29,i)=s(29,i)+bmy1(i)*cb18(i)
      s(30,i)=s(30,i)+bmz1(i)*cb18(i)
      s(31,i)=s(31,i)+bmy2(i)*cb18(i)
      s(32,i)=s(32,i)+bmz2(i)*cb18(i)
      s(33,i)=s(33,i)+bmy3(i)*cb18(i)
      s(34,i)=s(34,i)+bmz3(i)*cb18(i)
      s(35,i)=s(35,i)+bmy4(i)*cb18(i)
      s(36,i)=s(36,i)+bmz4(i)*cb18(i)
!  140 continue
   10 continue
      return
      end
