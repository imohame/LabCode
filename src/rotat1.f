      subroutine rotat1(sig,ener,ln)
c     implicit double precision (a-h,o-z)                                    dp
c      this a way of updating the stress rate based on the Green-Naghdi rate
c      unrotated and rotated configuratuions are used
c      CAUTION should be exercised not to confuse rates and configurations
c
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
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      common/berg/dels11(nelemg),dels22(nelemg),dels33(nelemg),
     > dels12(nelemg)
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      dimension sig(ln,*),ener(ln,*)
c&&&&&&&&&&&&&&&&&&&&&&&&&&
      do 5 i=mft,mlt
ck      write(27,15) sig(1,i),sig(2,i),sig(3,i),sig(4,i)
ck   15 format(5x,'s1=',f13.2,2x,'s2=',f13.2,2x,'s3=',f13.2,2x,
ck     .'s4=',f13.2/)
      dels11(i)=sig(1,i)
      dels22(i)=sig(2,i)
      dels33(i)=sig(3,i)
      dels12(i)=sig(4,i)
    5 continue
c&&&&&&&&&&&&&&&&&&&&&&
      do 10 i=mft,mlt
      a11(i)=r11(i)*sig(1,i)+r21(i)*sig(4,i)
      a12(i)=r12(i)*sig(1,i)+r22(i)*sig(4,i)
      a21(i)=r11(i)*sig(4,i)+r21(i)*sig(2,i)
      a22(i)=r12(i)*sig(4,i)+r22(i)*sig(2,i)
      sig(1,i)=r11(i)*a11(i)+r21(i)*a21(i)
      sig(2,i)=r12(i)*a12(i)+r22(i)*a22(i)
      sig(4,i)=r11(i)*a12(i)+r21(i)*a22(i)
   10 continue
c      if (idump.eq.0) return
      do 20 i=mft,mlt
      ener(1,i)=ener(1,i)+volmd(i)*(d1(i)*sig(1,i)+d2(i)*sig(2,i)
     1         +d3(i)*sig(3,i)+d4(i)*sig(4,i))/2.0
   20 continue
      return
      end
