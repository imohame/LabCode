c&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine strind
      common/bk02/ioofc,iphase,imass,model,lpar(8)
      common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
      common/bk14/lfna(15),lfnt(6)
      common/bk16/maxint
      common/newz1/ibase,iadd5,iadd6,icnt6,locst6
      common/main_block/ a(1)

      mark=igtpnt(14)
      lmark=maxint*lpar(8)*numelt
      call rdabsg(lfna(11),a(mark),lmark,iadd5,ioerr,bb,0)
      call riosta(lfna(11))
      iadd5=iadd5+lmark
      return
      end
