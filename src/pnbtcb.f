      subroutine pnbtcb (s,d)
c     implicit double precision (a-h,o-z)                                    dp
       use mod_parameters
      common/range/mft,mlt,lft,llt,nftm1
	  common /stressflag/ strflag1(nume),ctr_flag(nume)
c	  common /overlapping/ intersec(4, nume), area_coeff(nume), update_flag
c	  real area_coeff
	  integer strflag1
c      integer update_flag
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
      dimension s(44,1),d(4,4,1)
!$OMP PARALLEL DO       
      do 10 i=lft,llt
      cb11(i)=d(1,1,i)*py1(i)+d(1,4,i)*pz1(i)
      cb21(i)=d(2,1,i)*py1(i)+d(2,4,i)*pz1(i)
      cb41(i)=d(4,1,i)*py1(i)+d(4,4,i)*pz1(i)
      cb12(i)=d(1,2,i)*pz1(i)+d(1,4,i)*py1(i)
      cb22(i)=d(2,2,i)*pz1(i)+d(2,4,i)*py1(i)
      cb42(i)=d(4,2,i)*pz1(i)+d(4,4,i)*py1(i)
!   10 continue
!      do 20 i=lft,llt
      cb13(i)=d(1,1,i)*py2(i)+d(1,4,i)*pz2(i)
      cb23(i)=d(2,1,i)*py2(i)+d(2,4,i)*pz2(i)
      cb43(i)=d(4,1,i)*py2(i)+d(4,4,i)*pz2(i)
      cb14(i)=d(1,2,i)*pz2(i)+d(1,4,i)*py2(i)
      cb24(i)=d(2,2,i)*pz2(i)+d(2,4,i)*py2(i)
      cb44(i)=d(4,2,i)*pz2(i)+d(4,4,i)*py2(i)
!   20 continue
!      do 30 i=lft,llt
      cb15(i)=d(1,1,i)*py3(i)+d(1,4,i)*pz3(i)
      cb25(i)=d(2,1,i)*py3(i)+d(2,4,i)*pz3(i)
      cb45(i)=d(4,1,i)*py3(i)+d(4,4,i)*pz3(i)
      cb16(i)=d(1,2,i)*pz3(i)+d(1,4,i)*py3(i)
      cb26(i)=d(2,2,i)*pz3(i)+d(2,4,i)*py3(i)
      cb46(i)=d(4,2,i)*pz3(i)+d(4,4,i)*py3(i)
!   30 continue
!      do 40 i=lft,llt
      cb17(i)=d(1,1,i)*py4(i)+d(1,4,i)*pz4(i)
      cb27(i)=d(2,1,i)*py4(i)+d(2,4,i)*pz4(i)
      cb47(i)=d(4,1,i)*py4(i)+d(4,4,i)*pz4(i)
      cb18(i)=d(1,2,i)*pz4(i)+d(1,4,i)*py4(i)
      cb28(i)=d(2,2,i)*pz4(i)+d(2,4,i)*py4(i)
      cb48(i)=d(4,2,i)*pz4(i)+d(4,4,i)*py4(i)
!   40 continue
!      do 50 i=lft,llt
      s(1,i)=s(1,i)+py1(i)*cb11(i)+pz1(i)*cb41(i)
      s(2,i)=s(2,i)+pz1(i)*cb21(i)+py1(i)*cb41(i)
      s(3,i)=s(3,i)+pz1(i)*cb22(i)+py1(i)*cb42(i)
      s(4,i)=s(4,i)+py2(i)*cb11(i)+pz2(i)*cb41(i)
      s(5,i)=s(5,i)+py2(i)*cb12(i)+pz2(i)*cb42(i)
      s(6,i)=s(6,i)+py2(i)*cb13(i)+pz2(i)*cb43(i)
      s(7,i)=s(7,i)+pz2(i)*cb21(i)+py2(i)*cb41(i)
      s(8,i)=s(8,i)+pz2(i)*cb22(i)+py2(i)*cb42(i)
      s(9,i)=s(9,i)+pz2(i)*cb23(i)+py2(i)*cb43(i)
      s(10,i)=s(10,i)+pz2(i)*cb24(i)+py2(i)*cb44(i)
!   50 continue
!      do 60 i=lft,llt
      s(11,i)=s(11,i)+py3(i)*cb11(i)+pz3(i)*cb41(i)
      s(12,i)=s(12,i)+py3(i)*cb12(i)+pz3(i)*cb42(i)
      s(13,i)=s(13,i)+py3(i)*cb13(i)+pz3(i)*cb43(i)
      s(14,i)=s(14,i)+py3(i)*cb14(i)+pz3(i)*cb44(i)
      s(15,i)=s(15,i)+py3(i)*cb15(i)+pz3(i)*cb45(i)
      s(16,i)=s(16,i)+pz3(i)*cb21(i)+py3(i)*cb41(i)
      s(17,i)=s(17,i)+pz3(i)*cb22(i)+py3(i)*cb42(i)
      s(18,i)=s(18,i)+pz3(i)*cb23(i)+py3(i)*cb43(i)
      s(19,i)=s(19,i)+pz3(i)*cb24(i)+py3(i)*cb44(i)
      s(20,i)=s(20,i)+pz3(i)*cb25(i)+py3(i)*cb45(i)
!   60 continue
!      do 70 i=lft,llt
      s(21,i)=s(21,i)+pz3(i)*cb26(i)+py3(i)*cb46(i)
      s(22,i)=s(22,i)+py4(i)*cb11(i)+pz4(i)*cb41(i)
      s(23,i)=s(23,i)+py4(i)*cb12(i)+pz4(i)*cb42(i)
      s(24,i)=s(24,i)+py4(i)*cb13(i)+pz4(i)*cb43(i)
      s(25,i)=s(25,i)+py4(i)*cb14(i)+pz4(i)*cb44(i)
      s(26,i)=s(26,i)+py4(i)*cb15(i)+pz4(i)*cb45(i)
      s(27,i)=s(27,i)+py4(i)*cb16(i)+pz4(i)*cb46(i)
      s(28,i)=s(28,i)+py4(i)*cb17(i)+pz4(i)*cb47(i)
      s(29,i)=s(29,i)+pz4(i)*cb21(i)+py4(i)*cb41(i)
      s(30,i)=s(30,i)+pz4(i)*cb22(i)+py4(i)*cb42(i)
!   70 continue
!      do 80 i=lft,llt
      s(31,i)=s(31,i)+pz4(i)*cb23(i)+py4(i)*cb43(i)
      s(32,i)=s(32,i)+pz4(i)*cb24(i)+py4(i)*cb44(i)
      s(33,i)=s(33,i)+pz4(i)*cb25(i)+py4(i)*cb45(i)
      s(34,i)=s(34,i)+pz4(i)*cb26(i)+py4(i)*cb46(i)
      s(35,i)=s(35,i)+pz4(i)*cb27(i)+py4(i)*cb47(i)
      s(36,i)=s(36,i)+pz4(i)*cb28(i)+py4(i)*cb48(i)
!   80 continue
   10 continue
!c     update stiffness for cracked element
!c      do i=lft,llt
!c        ink=i+nftm1
!c		if(strflag1(ink)==1) then
!c		  do j=1,36
!c		    s(j,i)=area_coeff(ink)*s(j,i)
!c	      end do
!c		end if
!c	  end do		   
!c      if (update_flag==1) then
!      
!c      do i=lft,llt
!c	 
!c	    ink=i+nftm1
!c	    if(strflag1(ink)==2) then
!c	    open (unit = 11, file='stiff.dat', status = 'old')
!c	    do j=1,36
!c         read (11,*) s(j,i)
!c		end do
!c		close (unit=11)
!c		end if
!c		
!c      end do
!	  	 
!c	  end if
	  
      return
      end
