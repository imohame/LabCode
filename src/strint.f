      subroutine strint
c     implicit double precision (a-h,o-z)                                    dp
c
      common/bk02/ioofc,iphase,imass,model,lpar(8)
      common/bk16/maxint
      common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
      common/bk48/stress(4),strain(4),d(4,4),ipt,nel,nstate
      common/main_block/ a(1)
c
      kk=igtpnt(6)+nel-1
      nn=igtpnt(14)+maxint*lpar(8)*(nel-1)
      k08=igtpnt(8)
      mtp=locdbl(a(k08),nel-1)
      nm =igtpnt(12)+48*(mtp-1)
      k63=igtpnt(63)
      k64=igtpnt(64)
      k81=igtpnt(81)
c
      if(model.eq.0) go to 10
      goto(10,10,10,10,10,10,10,10,10,10,10,10), model
c
   10 call s3int (a(nn),lpar(8))
      return
      end
