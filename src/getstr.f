      subroutine getstr(a)
c     implicit double precision (a-h,o-z)                                    dp
c
c.... call appropriate routine to extract stresses
c
      common/bk02/ioofc,iphase,imass,model,numel,lpar(7)
      common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
      common/bk14/lfna(15),lfnt(6)
      common/bk16/maxint
      common/bk48/stress(4),strain(4),d(4,4),lst,nel,nstate
      dimension a(*)
c
      ipt=min(maxint,lst)
      mtp=locdbl(a(igtpnt(8)),nel-1)
      nm=igtpnt(12)+48*(mtp-1)
      ii=igtpnt(14)+maxint*lpar(7)*(nel-1)
      nn=ii+lpar(7)*(ipt-1)
      nt=nn+lpar(7)-1-idump
c
      goto(10,10,10,10,10,10,10,10,10,10,10,10), model
c
      write(lfnt(2),60)model
      call bye(2)
   10 call s1out (a(nn),a(nt))
      return
   60 format(' **fatal error** orion database write for'
     1 ' material',i5,' not coded (getstr)')
      end
