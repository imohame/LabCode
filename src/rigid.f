      subroutine rigid(ncn,nrcc)
c     implicit double precision (a-h,o-z)                                    dp
c
      common/bk14/lfna(15),lfnt(6)
      dimension ncn(3,*)
      character*80 txts,mssg
c
      if (nrcc.eq.0) return
c
      call header
      write(lfnt(2),40)
      do 10 i=1,nrcc
      call gttxsg (txts,lcount)
      read(unit=txts,fmt=20,err=12) (ncn(j,i),j=1,3)
      if (ncn(3,j).eq.0) ncn(3,j)=3
      write(lfnt(2),30) i,(ncn(j,i),j=1,3)
   10 continue
c
      return
c
   12 write (unit=mssg,fmt=50) i
      call termin (txts,mssg,lcount,1)
c
   20 format(3i5)
   30 format('constraint card',i5,' ties nodes',i5,' and',
     1i5,' together'/,'constraint type=',i1,'  y=1, z=2, y&z=3')
   40 format(///'n o d a l   c o n s t r a i n t   c a r d s'/)
   50 format(' error reading rigid constraint card #',i4)
      end
