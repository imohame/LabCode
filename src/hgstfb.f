      subroutine hgstfb (hgs,s)
c     implicit double precision (a-h,o-z)                                    dp
       use mod_parameters
      common/range/mft,mlt,lft,llt,nftm1
      common/vect13/
     1 x1112(nelemg),x1314(nelemg),x1114(nelemg),x1213(nelemg),
     2 x2122(nelemg),x2324(nelemg),x2124(nelemg),x2223(nelemg),
     3 hgmod1(nelemg),hgmod2(nelemg),hsv1(nelemg),hsv2(nelemg),
     > hsv3(nelemg),hsv4(nelemg),
     4 hg1(nelemg),hg2(nelemg),hsv11(nelemg),hsv12(nelemg),
     > hsv13(nelemg),hsv14(nelemg),
     5 hsv22(nelemg),hsv23(nelemg),hsv24(nelemg),hsv33(nelemg),
     > hsv34(nelemg),hsv44(nelemg)
      
      dimension hgs(*),s(44,1)
!$OMP PARALLEL DO       
      do 10 i=lft,llt
      hsv11(i)=hgs(i)*hsv1(i)*hsv1(i)
      hsv12(i)=hgs(i)*hsv1(i)*hsv2(i)
      hsv13(i)=hgs(i)*hsv1(i)*hsv3(i)
      hsv14(i)=hgs(i)*hsv1(i)*hsv4(i)
      hsv22(i)=hgs(i)*hsv2(i)*hsv2(i)
      hsv23(i)=hgs(i)*hsv2(i)*hsv3(i)
      hsv24(i)=hgs(i)*hsv2(i)*hsv4(i)
      hsv33(i)=hgs(i)*hsv3(i)*hsv3(i)
      hsv34(i)=hgs(i)*hsv3(i)*hsv4(i)
      hsv44(i)=hgs(i)*hsv4(i)*hsv4(i)
!   10 continue
!      do 20 i=lft,llt
      s(1,i) =s(1,i) +hsv11(i)
      s(3,i) =s(3,i) +hsv11(i)
      s(4,i) =s(4,i) +hsv12(i)
      s(6,i) =s(6,i) +hsv22(i)
      s(8,i) =s(8,i) +hsv12(i)
      s(10,i)=s(10,i)+hsv22(i)
      s(11,i)=s(11,i)+hsv13(i)
      s(13,i)=s(13,i)+hsv23(i)
      s(15,i)=s(15,i)+hsv33(i)
      s(17,i)=s(17,i)+hsv13(i)
      s(19,i)=s(19,i)+hsv23(i)
      s(21,i)=s(21,i)+hsv33(i)
      s(22,i)=s(22,i)+hsv14(i)
      s(24,i)=s(24,i)+hsv24(i)
      s(26,i)=s(26,i)+hsv34(i)
      s(28,i)=s(28,i)+hsv44(i)
      s(30,i)=s(30,i)+hsv14(i)
      s(32,i)=s(32,i)+hsv24(i)
      s(34,i)=s(34,i)+hsv34(i)
      s(36,i)=s(36,i)+hsv44(i)
   10 continue
      return
      end
