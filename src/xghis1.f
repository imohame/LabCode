      subroutine xghis1(a,ln,nft,nlt)
c     implicit double precision (a-h,o-z)                                    dp
       use mod_parameters
      common/vect91/xargat(nelemg),exx(nelemg,4),eps(nelemg,5),
     > sa(nelemg), sb(nelemg),sc(nelemg),ang1(nelemg),ang2(nelemg)
      dimension a(ln,*)
      do 10 i=nft,nlt
      xargat(i)=a(1,i)
   10 continue
      return
      end
