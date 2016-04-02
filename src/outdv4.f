      subroutine outdv4(u,udt,udtt,id,numnp,b,ln)                            nk
c     implicit double precision (a-h,o-z)                                    dp
c     subroutine outdv4(u,udt,udtt,id,numnp,b,ln,temp)                       pl
c
c     write displacements and velocities into plot file
c
      common/bk05/ifil,iadd,maxsiz,head(12)
      common/bk14/lfna(15),lfnt(6)
      common/bk20/ntotal
      common/bk34/bb(1)
c
      dimension u(*),udt(*),udtt(*),id(2,*),b(ln,*)                          nk
c     dimension u(*),udt(*),udtt(*),id(2,*),b(ln,*),temp(*)                  pl
c
      do 10 n=1,numnp
      node1=id(1,n)
      node2=id(2,n)
      n20=2*n
      n21=n20-1
      if (node1.ne.0) then
      b(n21,1)=u(node1)
      b(n21,2)=udt(node1)
      b(n21,3)=udtt(node1)
      else
      b(n21,1)=0.
      b(n21,2)=0.
      b(n21,3)=0.
      endif
      if (node2.ne.0) then
      b(n20,1)=u(node2)
      b(n20,2)=udt(node2)                                                    nk
      b(n20,3)=udtt(node2)
      else
      b(n20,1)=0.
      b(n20,2)=0.
      b(n20,3)=0.
      endif
c     b(n20,2)=temp(n)                                                       pl
   10 continue
ck      call wrabsg(lfna(7),b,3*ln,iadd,bb(ntotal+1),0)
ck      call riosta(lfna(7))
ck      iadd=iadd+3*ln
      return
      end
