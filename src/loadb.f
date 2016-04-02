      subroutine loadb
c     implicit double precision (a-h,o-z)                                    dp
c
c     generate and stores body force time-histories
c
      common/bk14/lfna(15),lfnt(6)
      common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
      common/bk27/nlcur,nptst,nthpy,nthpz,nthps,xmy,xmz,xms,nload,nptm
      character*80 txts,mssg
c
      if (nthpy+nthpz+abs(nthps).eq.0) return
c
      call header
      write(lfnt(2),70)
c
      if (nthpy.eq.0) go to 10
c
c     read in acceleration-time history in radial direction
c
      mssg=' error reading base acceleration card in y-direction'
      call gttxsg (txts,lcount)
      read(unit=txts,fmt=30,err=22) nthpy,xmy
      if (xmy.eq.0.0) xmy=1.0
      write(lfnt(2),40) nthpy,xmy
   10 if (nthpz.eq.0) go to 20
c
c     read in accerleration-time history in z-direction
c
      mssg=' error reading base acceleration card in z-direction'
      call gttxsg (txts,lcount)
      read(unit=txts,fmt=30,err=22) nthpz,xmz
      if (xmz.eq.0.0) xmz=1.0
      write(lfnt(2),50) nthpz,xmz
   20 if (nthps.eq.0) return
c
c     read in spin time history about z-axis
c
      ispny=0
      if (nthps.lt.0) ispny=1
      mssg=' error reading angular velocity card'
      call gttxsg (txts,lcount)
      read(unit=txts,fmt=30,err=22) nthps,xms
      if (xms.eq.0.0) xms=1.0
      write(lfnt(2),60) nthps,xms
      if (ispny.eq.1) nthps=-nthps
c
      return
c
   22 call termin (txts,mssg,lcount,1)
c
   30 format(i5,e10.0)
   40 format(//' y - a c c e l e r a t i o n '//5x,
     1       'load curve number=',i3,5x,'scale factor=',e10.3)
   50 format(//' z - a c c e l e r a t i o n '//5x,
     1       'load curve number=',i3,5x,'scale factor=',e10.3)
   60 format(//' a n g u l a r   v e l o c i t y'//5x,
     1       'load curve number=',i3,5x,'scale factor=',e10.3)
   70 format(///' b o d y  f o r c e  l o a d s')
      end
