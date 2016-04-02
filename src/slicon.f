      subroutine slicon(nd,ka,nsn,nmn,iloc,nsv,msr)
c     implicit double precision (a-h,o-z)                                    dp
      common/bk39/ip(9)
      dimension nd(*),ka(*),iloc(*),nsv(*),msr(*)
      do 20 i=1,nsn
      nc=iloc(i)
      ileft=max(nc-1,1)
      irite=min(nc+1,nmn)
      nn=1
      ip(1)=nsv(i)
      do 10 j=ileft,irite
      nn=nn+1
      ip(nn)=msr(j)
   10 continue
      call sink (nn,nd,ka,ip)
   20 continue
      return
      end
