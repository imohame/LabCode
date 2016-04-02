      subroutine blkcpi (irwflg,a,b,idiskl,nw)
c     implicit double precision (a-h,o-z)                                    dp
      common/bk02/ioofc,iphase,imass,lpar(9)
      common/bk14/lfna(15),lfnt(6)
      common/double/iprec,ncpw,unit
      dimension a(*),b(*)

      if (ioofc.eq.0) go to 10
      if (irwflg.eq.0) call blkcpy(a,b(idiskl),nw)
      if (irwflg.eq.1) call blkcpy(a(idiskl),b,nw)
      return

   10 idsklc=idiskl-igtpnt(54)
      if (irwflg.eq.0) call wrabsf(lfna(5),a,nw*iprec,idsklc*iprec)
      if (irwflg.eq.1) call rdabsf(lfna(5),b,nw*iprec,
     1 idsklc*iprec,ioerr)

      call riosta  (lfna(5))

      return
      end
