      subroutine ritemp(aux,iaux,told,tref,n)
c     implicit double precision (a-h,o-z)                                    dp
c
      common/bk16/maxint
      dimension aux(*),iaux(*),told(*)
c
      tref=0.0
      n2=n-2
      do 10 k=1,4
      kk=iaux(k)
   10 tref=tref+.25*told(kk)
      do 20 lst=1,maxint
      loc=13+n2*(lst-1)
   20 aux(loc)=tref
c
      return
      end
