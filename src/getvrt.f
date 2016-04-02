      subroutine getvrt(q11,q22,q12,q21)
c     implicit double precision (a-h,o-z)                                    dp
       use mod_parameters
      common/range/mft,mlt,lft,llt,nftm1
      common/vect3/
     1 dgi1(nelemg,4),dgi2(nelemg,4),dgi3(nelemg,4),
     > dgi4(nelemg,4),dgi5(nelemg,4),
     2 dgt1(nelemg,4),dgt2(nelemg,4),dgt3(nelemg,4),dgt4(nelemg,4),
     3 f11(nelemg),f22(nelemg),f12(nelemg),f21(nelemg)
      common/vect5/
     1  g11(nelemg),g12(nelemg),g21(nelemg),g22(nelemg),p11(nelemg),
     > p12(nelemg),
     2  p21(nelemg),p22(nelemg),eig1(nelemg),eig2(nelemg),b(nelemg),
     > c(nelemg),
     3  alen(nelemg),f1122(nelemg),f1221(nelemg)
c
      dimension q11(*),q22(*),q12(*),q21(*)
c
!$OMP PARALLEL DO       
      do 10 i=lft,llt
      f1122(i)=f11(i)+f22(i)
      f1221(i)=f12(i)-f21(i)
      alen(i) =1./sqrt(f1122(i)*f1122(i)+f1221(i)*f1221(i))
!   10 continue
!      do 20 i=lft,llt
      q11(i)=f1122(i)*alen(i)
      q12(i)=f1221(i)*alen(i)
      q22(i)= q11(i)
      q21(i)=-q12(i)
   10 continue
c
      return
      end
