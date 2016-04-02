      subroutine swtche(id,nodes,a,b,idir,n)
c     implicit double precision (a-h,o-z)                                    dp
      dimension id(2,*),nodes(*),a(*),b(*),idir(*)
      do 10 i=1,n
      ndf=idir(i)
      if (ndf.ne.3000000) then
      b(ndf)=a(ndf)
      a(ndf)=0.
      else
      nodnum=nodes(i)
      ndf1=id(1,nodnum)
      ndf2=id(2,nodnum)
      b(ndf1)=a(ndf1)
      b(ndf2)=a(ndf2)
      a(ndf1)=0.0
      a(ndf2)=0.0
      endif
   10 continue
      return
      end
