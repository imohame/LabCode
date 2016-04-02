      subroutine hgstfb3 (hgs,s)
c     implicit double precision (a-h,o-z)                                    dp
       use mod_parameters
      common/range/mft,mlt,lft,llt,nftm1
      common/bk26/dt,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
      common/vect13b/
     1 hgmod1(nelemg),hgmod2(nelemg),hsv1(nelemg),hsv2(nelemg),
     > hsv3(nelemg),hsv4(nelemg),
     2 hg1(nelemg),hg2(nelemg),hsv11(nelemg),hsv12(nelemg),
     > hsv13(nelemg),hsv14(nelemg),
     3 hsv22(nelemg),hsv23(nelemg),hsv24(nelemg),hsv33(nelemg),
     > hsv34(nelemg),hsv44(nelemg)
      common/vect0/
     1 q11(nelemg),q22(nelemg),q12(nelemg),q21(nelemg),
     2 r11(nelemg),r22(nelemg),r12(nelemg),r21(nelemg),
     3 s11(nelemg),s22(nelemg),s12(nelemg),s21(nelemg)
      common/vect7/
     1 diavg(nelemg),dfavg(nelemg),volnd(nelemg),volcd(nelemg),
     > vlinc(nelemg),volgp(nelemg),
     2 dimod(nelemg),sig11(nelemg),sig22(nelemg),sig33(nelemg),
     > sig12(nelemg)
c      common/hourglassk/khg11(128),khg12(128),khg13(128),khg14(128)
c     1  ,khg22(128),khg23(128),khg24(128),khg33(128),khg34(128),
c     2   khg44(128)
      dimension hgs(*),s(44,1)
!$OMP PARALLEL DO       
      do 10 i=lft,llt
      hsv11(i)=hgs(i)*hsv1(i)*hsv1(i)*volgp(i)
      hsv12(i)=hgs(i)*hsv1(i)*hsv2(i)*volgp(i)
      hsv13(i)=hgs(i)*hsv1(i)*hsv3(i)*volgp(i)
      hsv14(i)=hgs(i)*hsv1(i)*hsv4(i)*volgp(i)
      hsv22(i)=hgs(i)*hsv2(i)*hsv2(i)*volgp(i)
      hsv23(i)=hgs(i)*hsv2(i)*hsv3(i)*volgp(i)
      hsv24(i)=hgs(i)*hsv2(i)*hsv4(i)*volgp(i)
      hsv33(i)=hgs(i)*hsv3(i)*hsv3(i)*volgp(i)
      hsv34(i)=hgs(i)*hsv3(i)*hsv4(i)*volgp(i)
      hsv44(i)=hgs(i)*hsv4(i)*hsv4(i)*volgp(i)
!   10 continue
!      do 20 i=lft,llt
      s(1,i) =s(1,i) +(r11(i)**2+r12(i)**2)*hsv11(i)			!1x1x (11)
	  s(2,i) =s(2,i) +(r11(i)*r21(i)+r12(i)*r22(i))*hsv11(i)	!1x1y (12)
      s(3,i) =s(3,i) +(r21(i)**2+r22(i)**2)*hsv11(i)			!1y1y (22)
      s(4,i) =s(4,i) +(r11(i)**2+r12(i)**2)*hsv12(i)			!1x2x (13)
	  s(5,i) =s(5,i) +(r11(i)*r21(i)+r12(i)*r22(i))*hsv12(i)	!1y2x (23)!!!
      s(6,i) =s(6,i) +(r11(i)**2+r12(i)**2)*hsv22(i)			!2x2x (33) 
	  s(7,i) =s(7,i) +(r11(i)*r21(i)+r12(i)*r22(i))*hsv12(i)	!1x2y (14)
      s(8,i) =s(8,i) +(r21(i)**2+r22(i)**2)*hsv12(i)			!1y2y (24)
	  s(9,i) =s(9,i) +(r11(i)*r21(i)+r12(i)*r22(i))*hsv22(i)	!2x2y (34)
      s(10,i)=s(10,i)+(r21(i)**2+r22(i)**2)*hsv22(i)			!2y2y (44)
      s(11,i)=s(11,i)+(r11(i)**2+r12(i)**2)*hsv13(i)			!1x3x (15)
      s(12,i)=s(12,i)+(r11(i)*r21(i)+r12(i)*r22(i))*hsv13(i)	!1y3x (25)
      s(13,i)=s(13,i)+(r11(i)**2+r12(i)**2)*hsv23(i)	        !2x3x (35)
	  s(14,i)=s(14,i)+(r11(i)*r21(i)+r12(i)*r22(i))*hsv23(i)	!2y3x (45)
      s(15,i)=s(15,i)+(r11(i)**2+r12(i)**2)*hsv33(i)			!3x3x (55)
	  s(16,i)=s(16,i)+(r11(i)*r21(i)+r12(i)*r22(i))*hsv13(i)	!1x3y (16)
      s(17,i)=s(17,i)+(r21(i)**2+r22(i)**2)*hsv13(i)			!1y3y (26)
	  s(18,i)=s(18,i)+(r11(i)*r21(i)+r12(i)*r22(i))*hsv23(i)	!2x3y (36)
      s(19,i)=s(19,i)+(r21(i)**2+r22(i)**2)*hsv23(i)			!2y3y (46)
	  s(20,i)=s(20,i)+(r11(i)*r21(i)+r12(i)*r22(i))*hsv33(i)	!3x3y (56)
      s(21,i)=s(21,i)+(r21(i)**2+r22(i)**2)*hsv33(i)			!3y3y (66)
      s(22,i)=s(22,i)+(r11(i)**2+r12(i)**2)*hsv14(i)			!1x4x (17)
	  s(23,i)=s(23,i)+(r11(i)*r21(i)+r12(i)*r22(i))*hsv14(i)	!1y4x (18)
      s(24,i)=s(24,i)+(r11(i)**2+r12(i)**2)*hsv24(i)			!2x4x (37)
	  s(25,i)=s(25,i)+(r11(i)*r21(i)+r12(i)*r22(i))*hsv24(i)	!2y4x (47)
      s(26,i)=s(26,i)+(r11(i)**2+r12(i)**2)*hsv34(i)			!3x4x (57)
	  s(27,i)=s(27,i)+(r11(i)*r21(i)+r12(i)*r22(i))*hsv34(i)	!3y4x (67)
      s(28,i)=s(28,i)+(r11(i)**2+r12(i)**2)*hsv44(i)			!4x4x (77)
	  s(29,i)=s(29,i)+(r11(i)*r21(i)+r12(i)*r22(i))*hsv14(i)	!1x4y (18)
      s(30,i)=s(30,i)+(r21(i)**2+r22(i)**2)*hsv14(i)			!1y4y (28)
	  s(31,i)=s(31,i)+(r11(i)*r21(i)+r12(i)*r22(i))*hsv24(i)	!2x4y (38)
      s(32,i)=s(32,i)+(r21(i)**2+r22(i)**2)*hsv24(i)			!2y4y (48)
	  s(33,i)=s(33,i)+(r11(i)*r21(i)+r12(i)*r22(i))*hsv34(i)	!3x4y (58)
      s(34,i)=s(34,i)+(r21(i)**2+r22(i)**2)*hsv34(i)			!3y4y (68)
	  s(35,i)=s(35,i)+(r11(i)*r21(i)+r12(i)*r22(i))*hsv44(i)	!4x4y (78)
      s(36,i)=s(36,i)+(r21(i)**2+r22(i)**2)*hsv44(i)			!4y4y (88)
   10 continue
      return
      end
