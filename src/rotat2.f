      subroutine rotat2(sig,ener,ln)
c     implicit double precision (a-h,o-z)                                    dp
       use mod_parameters
      common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
      common/range/mft,mlt,lft,llt,nftm1
      common/vect0/
     1 q11(nelemg),q22(nelemg),q12(nelemg),q21(nelemg),
     2 r11(nelemg),r22(nelemg),r12(nelemg),r21(nelemg),
     3 s11(nelemg),s22(nelemg),s12(nelemg),s21(nelemg)
      common/vect5/
     1  g11(nelemg),g12(nelemg),g21(nelemg),g22(nelemg),
     > a11(nelemg),a12(nelemg),
     2  a21(nelemg),a22(nelemg)
      common/vect15/
     1 d1(nelemg),d2(nelemg),d3(nelemg),d4(nelemg)
      common/vect18/
     1 sclo(nelemg),volmd(nelemg),volmd4(nelemg,4)
      common/vect21/sig11(nelemg),sig22(nelemg),sig33(nelemg),
     > sig12(nelemg)
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      common/berg/dels11(nelemg),dels22(nelemg),dels33(nelemg),
     > dels12(nelemg)
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      dimension sig(ln,*),ener(ln,*)
c      if (idump.eq.0) go to 20
      do 10 i=mft,mlt
      ener(1,i)=ener(1,i)+volmd(i)*(d1(i)*sig(1,i)+d2(i)*sig(2,i)
     1         +d3(i)*sig(3,i)+d4(i)*sig(4,i))/2.0
   10 continue
   20 do 30 i=mft,mlt
      a11(i)=sig(1,i)*s11(i)+sig(4,i)*s12(i)
      a12(i)=sig(1,i)*s21(i)+sig(4,i)*s22(i)
      a21(i)=sig(4,i)*s11(i)+sig(2,i)*s12(i)
      a22(i)=sig(4,i)*s21(i)+sig(2,i)*s22(i)
      sig(1,i)=s11(i)*a11(i)+s12(i)*a21(i)
      sig(2,i)=s21(i)*a12(i)+s22(i)*a22(i)
      sig(4,i)=s11(i)*a12(i)+s12(i)*a22(i)
      sig11(i)=sig(1,i)
      sig22(i)=sig(2,i)
      sig33(i)=sig(3,i)
      sig12(i)=sig(4,i)
   30 continue
c&&&&&&&&&&&&&&&&&&&&&&&&&&
      do 40 i=mft,mlt
      dels11(i)=sig(1,i)-dels11(i)
      dels22(i)=sig(2,i)-dels22(i)
      dels33(i)=sig(3,i)-dels33(i)
      dels12(i)=sig(4,i)-dels12(i)
   40 continue
c&&&&&&&&&&&&&&&&&&&&&&&&
      return
      end
