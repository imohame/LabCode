      subroutine dbcint(id,nodes,idir,lc,amag,numdc)
c     implicit double precision (a-h,o-z)                                    dp
c
c     read and write displacement boundary condition cards
c
      common/bk14/lfna(15),lfnt(6)
      common/bk18/nummat,ityp2d,ako(31)
      dimension id(2,*),nodes(*),idir(*),lc(*),amag(8,1)
      character*80 txts,mssg
c
      if (numdc.eq.0) return
c
      call header
      write(lfnt(2),30)
c
      do 10 n=1,numdc
      call gttxsg (txts,lcount)
      read(unit=txts,fmt=20,err=14) nodes(n),idir(n),lc(n),amag(1,n),
     1 amag(2,n),amag(3,n),amag(4,n)
      if (amag(1,n).eq.0.0) amag(1,n)=1.0
      if (amag(2,n).eq.0.0) amag(2,n)=1.0e+20
      if (idir(n).ne.3) then
      write(lfnt(2),40) nodes(n),idir(n),lc(n),amag(1,n),amag(2,n)
      else
      if (ityp2d.eq.0) then
      write(*,70)
      write(lfnt(2),70)
      call bye(2)
      endif
      write(lfnt(2),60) nodes(n),lc(n),(amag(i,n),i=1,4)
      endif
      idr=idir(n)
      nnd=nodes(n)
      if (idr.ne.3) then
      idir(n)=id(idr,nnd)
      else
      idir(n)=3000000
      endif
   10 continue
c
      return
c
   14 write (unit=mssg,fmt=50) n
      call termin (txts,mssg,lcount,1)
c
   20 format(3i5,4e10.0)
   30 format(///' b o u n d a r y   d i s p l a c e m e n t  '
     1        ,'  c a r d s '//4x,
     2        ' node    direction    load curve   load curve multipl',
     3'ier   removal time',/'     center of rotation')
   40 format(3x,i6,4x,i7,9x,i5,8x,e12.4,9x,e12.4)
   50 format(' error reading displacement boundary condition card #',i5)
   60 format(3x,i6,3x,'x-rotation',7x,i5,8x,e12.4,9x,e12.4,
     1 /5x,'y=',e12.4,3x,'z=',e12.4)
   70 format(
     1 ' FATAL x-rotation is not permitted in axisymmetric problems',/
     2 ' please switch the analysis type to plane stress or strain.'//)
      end
