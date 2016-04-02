      subroutine s3int (aux,n)
c     implicit double precision (a-h,o-z)                                    dp
      common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
      common/bk16/maxint
      common/bk18/nummat,ityp2d,ako(31)
      common/bk46/nodm(8),reftem,thicknes
c
      dimension aux(n,1)
c
      do 10 j=1,maxint
      do 10 i=1,n
   10 aux(i,j)=0.0
c
      return
      end
