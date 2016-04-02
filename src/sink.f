c***********************************************************************
      subroutine sink (nn,nd,ka,ip)
c     implicit double precision (a-h,o-z)                                    dp
c***********************************************************************
c     this routine stores new node numbers in the connectivity
c     array and updates the nodal degree as required
c***********************************************************************
      dimension nd(2,*), ka(*), ip(*)
      common /sabr0/ maxnd
      np=nn
      if (ip(np-1).eq.ip(np)) np=np-1
      do 40 n=1,np
      i=ip(n)
      l=nd(1,i)
      kk=maxnd*(i-1)
      ll=0
      do 30 m=1,np
      if (m.eq.n) go to 30
      if (l.eq.0) go to 20
      do 10 k=1,l
      if (ka(kk+k).eq.ip(m)) go to 30
   10 continue
   20 ll=ll+1
      if (ll+l.gt.maxnd) go to 50
      ka(kk+ll+l)=ip(m)
   30 continue
      if (ll.eq.0) go to 40
      nd(1,i)=nd(1,i)+ll
   40 continue
      return
   50 write (*,60) i,maxnd
      call bye (2)
   60 format (/' error in input data --- connectivity for node ',i4,' ex
     1ceeds ',i3)
      end
