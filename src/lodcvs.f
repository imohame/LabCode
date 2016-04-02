      subroutine lodcvs(rv,timv,npc,pld,nlcur)
c     implicit double precision (a-h,o-z)                                    dp
c
      common/bk27/ndummy,nptst,nthpy,nthpz,nthps,xmy,xmz,xms,nload,nptm
c
      common/bk14/lfna(15),lfnt(6)
      dimension rv(*),timv(*),npc(*),pld(*)
      character*80 txts,mssg
c
      if (nlcur.eq.0) return
      call header
      write(lfnt(2),110)
      npc(1)=1
      do 40 l=1,nlcur
      call gttxsg (txts,lcount)
      read(unit=txts,fmt=50,err=42) ll,npts
      write(lfnt(2),90) ll,npts
      if (npts.le.nptm) go to 10
      write(lfnt(2),70) nptm
      call bye (2)
   10 if (l.eq.ll) go to 20
      write(lfnt(2),80)
      call bye (2)
   20 continue
      do 24 i=1,npts
      call gttxsg (txts,lcount)
      read (unit=txts,fmt=60,err=42) timv(i),rv(i)
   24 continue
      write(lfnt(2),100) (timv(i),rv(i),i=1,npts)
      npc(l+1)=npc(l)+2*npts
      i=npc(l)
      do 30 kk=1,npts
      pld(i)=timv(kk)
      i=i+1
      pld(i)=rv(kk)
   30 i=i+1
   40 continue
      nptst=npc(nlcur+1)-1
      return
c
   42 write (unit=mssg,fmt=120) l
      call termin (txts,mssg,lcount,1)
c
   50 format(2i5)
   60 format(2e10.0)
   70 format(/' fatal error - number of points (npts) exceeds maximum ('
     1,i5,')')
   80 format(/' fatal error - load curves out of order  ')
   90 format(//' load curve',i6,8x,'npts=',i5//6x,'time',10x,'value')
  100 format(e14.6,4x,e14.6)
  110 format(///' l o a d   c u r v e s ')
  120 format (' error read load curve data for load curve #',i3)
      end
