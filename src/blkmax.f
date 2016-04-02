      subroutine blkmax(matype,mt,prop,bulkmx,raystf)
c     implicit double precision (a-h,o-z)                                    dp
      common/bk14/lfna(15),lfnt(6)
      common/bk18/nummat,ityp2d,ako(31)
      dimension matype(*),prop(*),raystf(*)
c
      lc=48*(mt-1)+1
      mtype=matype(mt)
ck      go to (10,20,30,40,50,60,20,70,80,90,30,94,90,30,99,96,40,96,
ck     1 30,94,97,98,30,40,30),mtype
      go to (10,20,30,50,50,50,50,50),mtype
      write(lfnt(2),200)mtype
      call bye(2)
   10 ym=prop(lc+16)
      pr=prop(lc+17)
      go to 100
   20 ym=.33333*(prop(lc)+prop(lc+1)+prop(lc+2))
      pr=.33333*(prop(lc+3)+prop(lc+4)+prop(lc+5))
      go to 100
   30 ym=prop(lc)
      pr=prop(lc+1)
      go to 100
   40 ym=prop(lc+8)
      pr=prop(lc+16)
      go to 100
   50 ym=prop(lc+25)
      pr=prop(lc+26)
      go to 100
c   50 ym=3.*prop(lc+1)
c      pr=0.0
c      go to 100
c   60 ym=3.*prop(lc)
c      pr=0.0
c      go to 100
c   70 ym=3.*prop(lc+16)
c      pr=0.0
c      go to 100
c   80 ym=2.926*prop(lc)
c      pr=.463
c      go to 100
c   90 ym=prop(lc)
c      pr=prop(lc+1)
c      go to 100
c   94 ym=prop(lc+8)
c      pv=prop(lc+16)
c      go to 100
c   96 ym=9*prop(33)/(3.*prop(33)*prop(1)/prop(9)+1.)
c      pr=(3.*prop(33)-ym)/(6.*prop(33))
c      go to 100
c   97 ym=prop(lc+4)
c      pr=0.0
c      go to 100
c   98 ym=6.*(prop(lc+6)+prop(lc+1)/3.)
c      pr=0.0
c      goto 100
c   99 ym=9.*prop(lc)*prop(lc+1)/(3.*prop(lc)+prop(lc+1))
c      pr=(3.*prop(lc)-2.*prop(lc+1))/(2.*(3.*prop(lc)+prop(lc+1)))
  100 blkm=ym/(3.*(1.-2.*pr))
      bulkmx=max(blkm,bulkmx)
      raystf(2)=ym
      raystf(3)=pr
      return
  200 format(' **fatal error** bulk modulus computation',/,
     1       ' for material',i5,' not coded (blkmax)')
      end
