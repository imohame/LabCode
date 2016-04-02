      subroutine rott2s(sig1,sig2,sig4,i)
c     implicit double precision (a-h,o-z)                                    dp
c
c     scalar version of rotat2
c
c      common/vect0/
c     1 q11(128),q22(128),q12(128),q21(128),
c     2 r11(128),r22(128),r12(128),r21(128),
c     3 s11(128),s22(128),s12(128),s21(128)
c
c      a11=sig1*s11(i)+sig4*s12(i)
c      a12=sig1*s21(i)+sig4*s22(i)
c      a21=sig4*s11(i)+sig2*s12(i)
c      a22=sig4*s21(i)+sig2*s22(i)
c      sig1=s11(i)*a11+s12(i)*a21
c      sig2=s21(i)*a12+s22(i)*a22
c      sig4=s11(i)*a12+s12(i)*a22
      return
      end
