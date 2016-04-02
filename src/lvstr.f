      subroutine lvstr (ir,nd,lr,ka)
c     implicit double precision (a-h,o-z)                                    dp
c***********************************************************************
c     this routine builds a level structure rooted at node ir
c***********************************************************************
      dimension nd(2,*), lr(*), ka(*)
      common /sabr0/ maxnd
      common /sabr1/ mind, maxd, ibw, nbw, ipr, npr, nnc, nzn, nnp
      common /gps1/ lw, ld
      do 10 j=1,nnp
      lr(j)=0
   10 continue
      ln=1
      lw=1
      lr(ir)=ln
      go to 30
   20 lw=max(lw,klw)
   30 klw=0
      ld=ln
      ln=ld+1
      do 50 j=1,nnp
      if (lr(j).ne.ld) go to 50
      k=nd(2,j)
      l=nd(1,j)
      do 40 i=1,l
      node=ka(k+i)
      if (lr(node).ne.0) go to 40
      lr(node)=ln
      klw=klw+1
   40 continue
   50 continue
      if (klw.ne.0) go to 20
      return
      end
