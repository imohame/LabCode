      subroutine xghist(lcomp,nft,nlt)
c     implicit double precision (a-h,o-z)                                    dp
       use mod_parameters
      common/bk02/ioofc,iphase,imass,lpar(9)
      common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
      common/bk14/lfna(15),lfnt(6)
      common/bk16/maxint,hgc
      common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
      common/bk18/nummat,ityp2d,ako(31)
      common/intgp/d(4,4),ipt,nel,nelsub
      common/vect91/xargat(nelemg),exx(nelemg,4),eps(nelemg,5),
     > sa(nelemg),sb(nelemg),sc(nelemg),ang1(nelemg),ang2(nelemg)
      common/main_block/ a(1)

      ln=maxint*lpar(9)
      ii=igtpnt(14)+ln*(nel-1)
      nn=ii+lpar(9)*(ipt-1)
      ne=nn+lpar(9)-1
      nt=ne-idump
      mm=igtpnt(2)+4*(nel-1)
      k11=igtpnt(11)
      k82=igtpnt(82)
      if(lcomp.eq.19)then
ckk      if(ityp2d.eq.2)then
ckk      call xghis1(a(nt),ln,nft,nlt)
ckk      else
      do 10 i=nft,nlt
      xargat(i)=0.
   10 continue
ckk      endif
      elseif(lcomp.eq.20)then
      if(idump.eq.1.and.itmpop.eq.0)then
      call xghis1(a(ne),ln,nft,nlt)
      elseif(idump.eq.1.and.itmpop.eq.1)then
      call xghis2(a(k82),a(mm),nft,nlt)
      else
      do 20 i=nft,nlt
      xargat(i)=0.
   20 continue
      endif
      elseif(lcomp.eq.40)then
      call xghis1(a(k11),1,nft,nlt)
      elseif(lcomp.gt.40)then
      call xghis1(a(nn+lcomp-37),ln,nft,nlt)
      endif
      return
      end
