      subroutine rotstr(dp1,dp2,dp3,dp4,dp5)
c     implicit double precision (a-h,o-z)                                    dp
       use mod_parameters
      common/range/mft,mlt,lft,llt,nftm1
      common/vect0/
     1 q11(nelemg),q22(nelemg),q12(nelemg),q21(nelemg),
     2 r11(nelemg),r22(nelemg),r12(nelemg),r21(nelemg),
     3 s11(nelemg),s22(nelemg),s12(nelemg),s21(nelemg)
      common/vect5/
     1  g11(nelemg),g12(nelemg),g21(nelemg),g22(nelemg),
     > a11(nelemg),a12(nelemg),
     2  a21(nelemg),a22(nelemg),eig1(nelemg),eig2(nelemg),
     > b(nelemg),c(nelemg),
     3  ff11(nelemg),ff12(nelemg),ff22(nelemg),
     4  u11(nelemg),u12(nelemg),u22(nelemg),ui11(nelemg),
     > ui12(nelemg),ui22(nelemg)
      common/vect15/
     1 d1(nelemg),d2(nelemg),d3(nelemg),d4(nelemg)
      common/vect91/dpp(nelemg)
      dimension dp1(*),dp2(*),dp3(*),dp4(*),dp5(*)
!$OMP PARALLEL DO       
      do 10 i=lft,llt
      dpp(i)=.50*(dp3(i)+dp4(i))
   10 continue
!$OMP PARALLEL DO       
      do 20 i=lft,llt
      a11(i)=q11(i)*dp1(i)+q21(i)*dpp(i)
      a12(i)=q12(i)*dp1(i)+q22(i)*dpp(i)
      a21(i)=q11(i)*dpp(i)+q21(i)*dp2(i)
   20 a22(i)=q12(i)*dpp(i)+q22(i)*dp2(i)
!$OMP PARALLEL DO       
      do 30 i=lft,llt
      d1(i)=q11(i)*a11(i)+q21(i)*a21(i)
      d2(i)=q12(i)*a12(i)+q22(i)*a22(i)
      d3(i)=dp5(i)
      d4(i)=2.*(q11(i)*a12(i)+q21(i)*a22(i))
   30 continue
c
      return
      end
