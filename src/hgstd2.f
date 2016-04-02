      subroutine hgstd2 (hgs,s)
c     implicit double precision (a-h,o-z)                                    dp
      common/range/mft,mlt,lft,llt,nftm1
      dimension hgs(*),s(44,1)
!$OMP PARALLEL DO       
      do 10 i=lft,llt
      s(1,i) =s(1,i) +hgs(i)
      s(3,i) =s(3,i) +hgs(i)
      s(4,i) =s(4,i) -hgs(i)
      s(6,i) =s(6,i) +hgs(i)
      s(8,i) =s(8,i) -hgs(i)
      s(10,i)=s(10,i)+hgs(i)
      s(11,i)=s(11,i)+hgs(i)
      s(13,i)=s(13,i)-hgs(i)
      s(15,i)=s(15,i)+hgs(i)
      s(17,i)=s(17,i)+hgs(i)
      s(19,i)=s(19,i)-hgs(i)
      s(21,i)=s(21,i)+hgs(i)
      s(22,i)=s(22,i)-hgs(i)
      s(24,i)=s(24,i)+hgs(i)
      s(26,i)=s(26,i)-hgs(i)
      s(28,i)=s(28,i)+hgs(i)
      s(30,i)=s(30,i)-hgs(i)
      s(32,i)=s(32,i)+hgs(i)
      s(34,i)=s(34,i)-hgs(i)
      s(36,i)=s(36,i)+hgs(i)
   10 continue
      return
      end
