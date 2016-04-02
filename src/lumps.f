      subroutine lumps(gm,id,nod,xmass,neq)
c     implicit double precision (a-h,o-z)                                    dp
      common/bk02/ioofc,iphase,imass,lpar(9)
      common/bk03/numdc,imassn,idampn,irller,penstf
      dimension gm(*),id(2,*),nod(*),xmass(2,*)
      do 10 i=1,neq
   10 gm(i)=0.
      if (imassn.eq.0) return
      do 20 i=1,imassn
      node=nod(i)
      ly=id(1,node)
      lz=id(2,node)
      if (ly.ne.0) gm(ly)=gm(ly)+xmass(1,i)
      if (lz.ne.0) gm(lz)=gm(lz)+xmass(2,i)
   20 continue
      return
      end
