      subroutine numrl (i1,i2,nr,nd,lv,ka)
c     implicit double precision (a-h,o-z)                                    dp
c***********************************************************************
c     this routine numbers the nodes in level l
c***********************************************************************
      dimension nr(*), nd(2,*), lv(*), ka(*)
      common /sabr0/ maxnd
      common /gps2/ inbr, lnbr
      ki=i1
   10 ki=ki+1
      j=nr(ki)
      k=nd(2,j)
      l=nd(1,j)
      do 20 i=1,l
      node=ka(k+i)
      if (lv(node).ne.lnbr) go to 20
      inbr=inbr+1
      nr(inbr)=node
      lv(node)=0
   20 continue
      if (ki.lt.min(inbr,i2)) go to 10
      return
      end
