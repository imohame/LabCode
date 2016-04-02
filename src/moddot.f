      subroutine moddot (id,nodes,idir,numdc,v1,v2,change)
c     implicit double precision (a-h,o-z)                                    dp
      dimension id(2,*),nodes(*),idir(*),v1(*),v2(*)
      do 10 i=1,numdc
      ndf=idir(i)
      if (ndf.ne.3000000) then
      change=change+v1(ndf)*v2(ndf)
      else
      numnod=nodes(i)
      ndf1=id(1,numnod)
      ndf2=id(2,numnod)
      change=change+v1(ndf1)*v2(ndf1)+v1(ndf2)*v2(ndf2)
      endif
   10 continue
      return
      end
