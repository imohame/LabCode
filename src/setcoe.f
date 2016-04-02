      subroutine setcoe (hgs,s)
c     implicit double precision (a-h,o-z)                                    dp
      common/range/mft,mlt,lft,llt,nftm1
      common/bk16/maxint,hgc
      dimension idiag(8),hgs(*),s(44,1)
      data idiag/1,3,6,10,15,21,28,36/
c
      do 10 i=lft,llt
      hgs(i)=0.0
   10 continue
      do 30 k=1,8
      j=idiag(k)
      do 20 i=lft,llt
      hgs(i)=max(s(j,i),hgs(i))
   20 continue
   30 continue
      do 40 i=lft,llt
   40 hgs(i)=hgc*hgs(i)
      return
      end
