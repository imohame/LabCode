      subroutine prinst
c     implicit double precision (a-h,o-z)                                    dp
c
c.... compute quantities for hsp printing
c
      common/bk18/nummat,ityp2d,ako(31)
      common/bk48/stress(4),ft,p1,p2,ag,d(4,4),ipt,nel,nstate
c
      cc=(stress(1)+stress(2))*0.5
      bb=(stress(1)-stress(2))*0.5
      cr=sqrt(bb**2+stress(4)**2)
      p1=cc+cr
      p2=cc-cr
      ag=0.0
c
      strss3=stress(3)
ckk      if (ityp2d.eq.2) strss3=0.
      if (ft.eq.0.0) ft=sqrt(.50*((stress(1)-stress(2))**2+(stress(2)
     1 -strss3)**2+(strss3-stress(1))**2)+3.*stress(4)**2)
      if (abs(bb).lt.1.0e-8) return
      ag=28.648*atan2(stress(4),bb)
c
      return
      end
