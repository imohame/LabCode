      subroutine intemp
c     implicit double precision (a-h,o-z)                                    dp
      common/bk02/ioofc,iphase,imass,lpar(9)
      equivalence (lpar(2),numel)
      do 10 nel=1,numel
   10 call tmpint (nel)
      return
      end
