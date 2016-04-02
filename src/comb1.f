      subroutine comb1 (nd,nr,lv,ka)
c     implicit double precision (a-h,o-z)                                    dp
c
c***********************************************************************
c     this routine calculates the bandwidth with                       *
c     the nodal ordering given in the vector nr                        *
c***********************************************************************
      dimension nd(2,*), nr(*), lv(*), ka(*)
      common /sabr0/ maxnd
      common /sabr1/ mind, maxd, ibw, nbw, ipr, npr, nnc, nzn, nnp
      common /sabr2/ kbw, kpr
      do 10 j=1,nnp
   10 lv(nr(j))=j
      kbw=0
      kpr=0
      do 30 j=1,nnp
      k=nd(2,j)
      l=nd(1,j)
      kkh=0
      do 20 i=1,l
      node=ka(k+i)
      kkh=max(kkh,lv(j)-lv(node))
   20 continue
      kpr=kpr+kkh
      kbw=max(kbw,kkh)
   30 continue
      return
      end
