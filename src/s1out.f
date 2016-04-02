      subroutine s1out (sig,thick)
c     implicit double precision (a-h,o-z)                                    dp
c
c.... extract stresses
c
      common/bk18/nummat,ityp2d,ako(31)
      common/bk48/stress(4),strain(4),d(4,4),ipt,nel,nstate
      dimension sig(*)
c
      do 10 i=1,4
   10 stress(i)=sig(i)
      nstate=1
      strain(1)=0.0
c
      return
      end
