      subroutine tmpint(nel)
c     implicit double precision (a-h,o-z)                                    dp
c
      common/bk02/ioofc,iphase,imass,model,lpar(8)
      common/bk16/maxint
      common/main_block/ a(1)
c
      k08=igtpnt(8)
      mtp=locdbl(a(k08),nel-1)
      model=locdbl(a(1),mtp)
c
      go to (10,10,10,20,10,10,20,10,10,10,10), model
c
   10 return
c
   20 nm =igtpnt(12)+48*(mtp-1)
      kk =igtpnt(6)+nel-1
      nn =igtpnt(14)+maxint*lpar(8)*(nel-1)
      k82=igtpnt(82)
      call ritemp (a(nn),a(nn),a(k82+1),a(kk),lpar(8))
      return
c
      end
