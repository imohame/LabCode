      subroutine gttxsg(txts,lcount)
      character*80 txts,mssg
      integer count
      common/bk14/lfna(15),lfnt(6)
      common/cn0/iconv
      save count
      data count/0/

   10 read (lfnt(3),20) txts
      count=count+1
      if (txts(1:1).eq.'*'.or.txts(1:1).eq.'$') then
         write (lfna(10),30) txts
         go to 10
      else
ck      if(iconv.eq.1)write(lfnt(4),20)txts
         lcount=count

c         write(7777,*) 'gttxsg.f : --'
c         write(7777,*) '  txts: ',txts
c         write(7777,*) '  lcount: ',lcount
c         write(7777,*) '------------ end '

         return
      endif
   20 format (a80)
   30 format (//,a80,//)
      end
