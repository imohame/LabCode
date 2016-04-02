      subroutine header
c     implicit double precision (a-h,o-z)                                    dp
c
      real*8 hed                                                        vax750
      character*8 vn,cdate
      character*72 tx1,tx2
      common/bk07/mbfc,nelpg,hed(12)
      common/bk14/lfna(15),lfnt(6)
      common/vrsn/vn,cdate,tx1,tx2
c

ck      write(lfnt(2),10) hed,vn,cdate
      return
c
   10 format('1',12a6,/,
     1 t20,'NIKE2D VERSION ',a8,'  COMPILED ',a8)
      end
