      subroutine nsetka(neltyp,nd,ka,ix,matp,ncn,nrcc,nsln,nmln,islt,
     1 ilocs,ilocm,nsv,msr,nsl)
c     implicit double precision (a-h,o-z)                                    dp
c
c     this routine  constructs the connectivity and nodal degree
c     arrays needed for the gps bandwidth minimization algorithm.
c
      dimension nd(*),ka(*),ix(4,*),matp(*),ncn(3,*),nsln(*),nmln(*),
     1  islt(*),ilocs(*),ilocm(*),nsv(*),msr(*)
      common/bk39/ip(9)
c
      do 40 ne=1,neltyp
      if (matp(ne).eq.0) go to 40
      ip(1)=ix(1,ne)
      ip(2)=ix(2,ne)
      ip(3)=ix(3,ne)
      ip(4)=ix(4,ne)
      if (nrcc.eq.0) go to 30
      do 20 i=1,4
      node=ip(i)
      do 10 j=1,nrcc
      if (node.ne.ncn(2,j)) go to 10
      if (ncn(3,j).lt.3) go to 10
      ip(i)=ncn(1,j)
      go to 20
   10 continue
   20 continue
   30 call sink (4,nd,ka,ip)
   40 continue
      if (nsl.eq.0) return
      ls=1
      lm=1
      do 60 nm=1,nsl
      if (islt(nm).eq.6) go to 50
      call slicon (nd,ka,nsln(nm),nmln(nm),ilocs(ls),nsv(ls),msr(lm))
      if (islt(nm).lt.3) go to 50
      call slicon (nd,ka,nmln(nm),nsln(nm),ilocm(lm),msr(lm),nsv(ls))
   50 ls=ls+nsln(nm)
   60 lm=lm+nmln(nm)
      return
      end
