      subroutine termin (txts,mssg,lcount,iprint)
c     implicit double precision (a-h,o-z)                                    dp
      character*80 txts,mssg
      common/bk14/lfna(15),lfnt(6)
c
      if (iprint.eq.0) then
      write (lfnt(2),30) lcount,lcount,txts
      write (*,30) lcount,lcount,txts
      else
      write (lfnt(2),20) lcount,lcount,txts,mssg
      write (*,20) lcount,lcount,txts,mssg
      endif
      call bye (2)
      return
   20 format(///
     1 '     LINE NUMBER',i7,' CONTAINS IMPROPERLY FORMATTED DATA',//
     2,'************************************************line#',i7
     3,/,a80,/
     4,'************************************************************',/,
     5   a80,/)
   30 format(///
     1 '     LINE NUMBER',i7,' CONTAINS IMPROPERLY FORMATTED DATA',//
     2,'************************************************line#',i7
     3,/,a80,/
     4,'************************************************************',/)
      end
