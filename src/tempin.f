      subroutine tempin(tbase,tmult)
c     implicit double precision (a-h,o-z)                                    dp
c
c     read and print nodal point data
c
      common/bk14/lfna(15),lfnt(6)
      common/cn1/numati,numnpi,numeli,nblk1,nslidi,ntslvi,ntmsri,
     1           nnpbi,nepbi,ncnpi
      common/cn3/ibbari,intgi,nmbfi,ithopi,ithcri,ithini,iengri,
     1           ijinti
      dimension tbase(*),tmult(*)
      character*80 txts,mssg
c
      kn0=0
      nold=0
   20 kn=kn0
      call gttxsg (txts,lcount)
      read(unit=txts,fmt=90,err=75) n,tbase(n),tmult(n),kn0
      if (kn.eq.0) kn=1
      if (nold.eq.0) go to 60
      num=(n-nold)/kn
      numn=num-1
      if (numn.lt.1) go to 60
      xnum=num
      dbase=(tbase(n)-tbase(nold))/xnum
      dmult=(tmult(n)-tmult(nold))/xnum
      k=nold
      do 50 j=1,numn
      kk=k
      k=k+kn
      tbase(k)=tbase(kk)+dbase
      tmult(k)=tmult(kk)+dmult
   50 continue
   60 nold=n
      if (n.ne.numnpi) go to 20
c
      write(lfnt(2),80)
      do 70 i=1,numnpi
      write(lfnt(2),81)i,tbase(i),tmult(i)
   70 continue
      return
c
   75 nold=nold+1
      write (unit=mssg,fmt=110) nold
      call termin (txts,mssg,lcount,1)
c
   80 format(///' nodal temperature data ',/,
     1          ' ---------------------- '//,
     2 t5,'node',t20,'base-temp',t35,'mult-temp'/)
   81 format(t5,i5,t20,e12.5,t35,e12.5)
   90 format(i5,2(e10.1),i5)
  110 format(' error reading nodal temperature',/,
     1       ' cards, probably for node #',i6)
      end
