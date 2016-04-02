      subroutine ebfin(freep)
c     implicit double precision (a-h,o-z)                                    dp
c
c     read and print nodal point data
c
      common/bk07/mbfc,nelpg,hed(12)
      common/bk14/lfna(15),lfnt(6)
      common/cn1/numati,numnpi,numeli,nblk1,nslidi,ntslvi,ntmsri,
     1           nnpbi,nepbi,ncnpi
      common/cn3/ibbari,intgi,nmbfi,ithopi,ithcri,ithini,iengri,
     1           ijinti
      dimension freep(5,*),idum(5)
      character*80 txts,mssg
c
      call gttxsg(txts,lcount)
      read(unit=txts,fmt=130,err=140)mbfc,(idum(j),j=1,5)
      idu=idum(1)+idum(2)+idum(3)+idum(4)+idum(5)
      if(idu.ne.0)goto 140
      kgen=0
      kn0=0
      nold=0
   20 kn=kn0
      call gttxsg(txts,lcount)
      read(unit=txts,fmt=90,err=75) n,freep(1,n),freep(2,n),kn0
      kgen=kgen+1
      if (kn.eq.0) kn=1
      if (nold.eq.0) go to 60
      num=(n-nold)/kn
      numn=num-1
      if (numn.lt.1) go to 60
      xnum=num
      dfrp1=(freep(1,n)-freep(1,nold))/xnum
      dfrp2=(freep(2,n)-freep(2,nold))/xnum
      k=nold
      do 50 j=1,numn
      kgen=kgen+1
      kk=k
      k=k+kn
      freep(1,k)=freep(1,kk)+dfrp1
      freep(2,k)=freep(2,kk)+dfrp2
   50 continue
   60 nold=n
      if(kgen.lt.nmbfi)goto 20
c
      write(lfnt(2),80)
      do 70 i=1,numeli
      write(lfnt(2),81)i,freep(1,i),freep(2,i)
   70 continue
      return
c
   75 nold=nold+1
      write (unit=mssg,fmt=110) nold
      call termin (txts,mssg,lcount,1)
  140 write(lfnt(2),120)
      write(*,120)
      call bye(2)
c
   80 format(///' element body force data',/,
     1          ' ---------------------- '//,
     2 t5,'element',t20,'force-x',t35,'force-y'/)
   81 format(t5,i5,t20,e12.5,t35,e12.5)
   90 format(i5,2(e10.0),i5)
  110 format(' error reading element body force',/,
     1       ' cards, probably for element #',i6)
  120 format(' new format for element body force loads'/
     1       ' is different than described in 91 manual'/
     2       ' first card: load curve # in cols 1-5'/
     3       ' card 2 to nbfl+1 same as manual')
  130 format(6i5)
      end
