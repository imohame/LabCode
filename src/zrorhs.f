      subroutine zrorhs(id,nodes,r,idir,n)
c     implicit double precision (a-h,o-z)                                    dp
      dimension id(2,*),nodes(*),r(*),idir(*)
      do 10 i=1,n
      ndf=idir(i)
      if (ndf.ne.3000000) then
      r(ndf)=0.0
      else
      nodnum=nodes(i)
      r(id(1,nodnum))=0.0
      r(id(2,nodnum))=0.0
      endif
   10 continue
      return
      end
