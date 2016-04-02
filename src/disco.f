      subroutine disco (ir,nd,ls,ka,kc,kk)
c     implicit double precision (a-h,o-z)                                    dp
c
c***********************************************************************
c     this routine finds the nodes in disjoint connected set kk        *
c***********************************************************************
      dimension nd(2,*), ls(*), ka(*), kc(*)
      common /sabr0/ maxnd
      common /sabr1/ mind, maxd, ibw, nbw, ipr, npr, nnc, nzn, nnp
      common /gps1/ lw, ld
      ls(ir)=kk
      lw=1
      kc(lw)=ir
   10 continue
      klw=0
      do 30 j=1,nnp
      if (ls(j).ne.kk) go to 30
      k=nd(2,j)
      l=nd(1,j)
      do 20 i=1,l
      node=ka(k+i)
      if (ls(node).ne.0) go to 20
      ls(node)=kk
      lw=lw+1
      kc(lw)=node
      klw=klw+1
   20 continue
   30 continue
      if (klw.ne.0) go to 10
      return
      end
