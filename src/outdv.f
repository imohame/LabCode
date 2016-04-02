      subroutine outdv(u,id,numnp)
c     implicit double precision (a-h,o-z)                                    dp
c
c     write mode shapes into thor plotfile
c
      common/bk05/ifil,iadd,maxsiz,head(12)
      common/bk14/lfna(15),lfnt(6)
      common/bk20/ntotal
      common/bk34/b(1)
c
      dimension u(*),id(2,*),d(2)
c
      loc=0
      do 10 n=1,numnp
      d1=0.
      d2=0.
      k1=id(1,n)
      k2=id(2,n)
      if (k1.ne.0) d1=u(k1)
      if (k2.ne.0) d2=u(k2)
      loc=loc+1
      b(loc)=d1
      loc=loc+1
      b(loc)=d2
      if (loc+2.lt.ntotal) go to 10
ck      call wrabsg(lfna(7),b,loc,iadd,b(ntotal+1),0)
ck      call riosta(lfna(7))
ck      iadd=iadd+loc
      loc=0
   10 continue
ck      if (loc.eq.0) go to 20
ck      call wrabsg(lfna(7),b,loc,iadd,b(ntotal+1),0)
ck      call riosta(lfna(7))
ck      iadd=iadd+loc
ck   20 iadd=iadd+4*numnp
      return
      end
