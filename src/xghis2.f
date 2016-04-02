      subroutine xghis2(temp,nodes,nft,nlt)
c     implicit double precision (a-h,o-z)                                    dp
       use mod_parameters
      common/vect91/xargat(nelemg),exx(nelemg,4),eps(nelemg,5),
     > sa(nelemg),sb(nelemg),sc(nelemg),ang1(nelemg),ang2(nelemg)
      dimension temp(*),nodes(4,*)
      do 10 i=nft,nlt
      xargat(i)=0.
   10 continue
      do 30 k=1,4
      do 20 i=nft,nlt
      xargat(i)=xargat(i)+.25*temp(nodes(k,i))
   20 continue
   30 continue
      return
      end
