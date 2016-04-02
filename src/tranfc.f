      subroutine tranfc(i)
c     implicit double precision (a-h,o-z)                                    dp
c      common/range/mft,mlt,lft,llt,nftm1
c      common/vect0/
c     1 q11(128),q22(128),q12(128),q21(128),
c     2 r11(128),r22(128),r12(128),r21(128),
c     3 s11(128),s22(128),s12(128),s21(128)
c      common/vect8/
c     1 d(4,4,128)
c      t11=s11(i)**2
c      t12=s12(i)**2
c      t14=2.*s11(i)*s12(i)
c      t21=s21(i)**2
c      t22=s22(i)**2
c      t24=2.*s21(i)*s22(i)
c      t41=s11(i)*s21(i)
c      t42=s12(i)*s22(i)
c      t44=s11(i)*s22(i)+s12(i)*s21(i)
c      e11=t11*d(1,1,i)+t12*d(2,1,i)+t14*d(4,1,i)
c      e12=t11*d(1,2,i)+t12*d(2,2,i)+t14*d(4,2,i)
c      e13=t11*d(1,3,i)+t12*d(2,3,i)+t14*d(4,3,i)
c      e14=t11*d(1,4,i)+t12*d(2,4,i)+t14*d(4,4,i)
c      e21=t21*d(1,1,i)+t22*d(2,1,i)+t24*d(4,1,i)
c      e22=t21*d(1,2,i)+t22*d(2,2,i)+t24*d(4,2,i)
c      e23=t21*d(1,3,i)+t22*d(2,3,i)+t24*d(4,3,i)
c      e24=t21*d(1,4,i)+t22*d(2,4,i)+t24*d(4,4,i)
c      e41=t41*d(1,1,i)+t42*d(2,1,i)+t44*d(4,1,i)
c      e42=t41*d(1,2,i)+t42*d(2,2,i)+t44*d(4,2,i)
c      e43=t41*d(1,3,i)+t42*d(2,3,i)+t44*d(4,3,i)
c      e44=t41*d(1,4,i)+t42*d(2,4,i)+t44*d(4,4,i)
c      d(1,1,i)=e11*t11+e12*t12+e14*t14
c      d(1,2,i)=e11*t21+e12*t22+e14*t24
c      d(1,3,i)=e13
c      d(1,4,i)=e11*t41+e12*t42+e14*t44
c      d(2,1,i)=d(1,2,i)
c      d(2,2,i)=e21*t21+e22*t22+e24*t24
c      d(2,3,i)=e23
c      d(2,4,i)=e21*t41+e22*t42+e24*t44
c      d(3,1,i)=e13
c      d(3,2,i)=d(2,3,i)
c      d(3,4,i)=e43
c      d(4,1,i)=d(1,4,i)
c      d(4,2,i)=d(2,4,i)
c      d(4,3,i)=e43
c      d(4,4,i)=e41*t41+e42*t42+e44*t44
      return
      end
