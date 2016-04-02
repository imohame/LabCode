      subroutine actred(ia,a1,ja,a2,ie,ic,d)
c     implicit double precision (a-h,o-z)                                    dp
c
c     reduce block a2 using block a1
c
      dimension ia(*),a1(*),ja(*),a2(*),d(*)
      common/fissn0/neq,mxw,no,n1,nfissl(3)
      common/fissn3/ifissl,kfissl(3)
      common/fissl4/je,jc
      common/fissl5/dt(2)
c
      ije=je-ie
      il=min(ije,ic)
      j2=jc
      do 400 j=1,jc
      j1=j2+1
      j2=ja(j)
      jh=j2-j1
c.... skip if column is too short
      if(jh.le.1) go to 301
      ke=ije+j
      is=max(1,ke-jh)
      if(il.lt.is) go to 301
c.... reduce terms in column j (except diagonal term)
      k=j2-min(ke-1,jh)
      call recol(k-j1,is,il,ia,a1,a2(k),d)
c.... reduce diagonal term
  301 if(ije.ne.0) go to 400
      il=j
      k=je+j
      if(recold(a2(j1),d(j-jh),jh).ne.0.) go to 351
      kfissl(1)=kfissl(1)+1
      if(nfissl(1).lt.kfissl(1)) go to 400
      write(no,1001)k
      go to 400
  351 if(dt(1).ge.0.) go to 352
      kfissl(2)=kfissl(2)+1
      if(nfissl(2).lt.kfissl(2)) go to 352
      write(no,1002)k
  352 if(dt(2).ne.0.) go to 400
      kfissl(3)=kfissl(3)+1
      if(nfissl(3).lt.kfissl(3)) go to 400
      write(no,1003)k
  400 continue
      return

 1001 format(5x,'**warning** zero pivot for equation ',i5
     1,' detected during factorization')
 1002 format(5x,'**warning** change in sign of equation ',i5
     1,' detected during factorization')
 1003 format(5x,'**warning** equation ',i5
     1,' has lost at least 7 significant digits during factorization')
      end
