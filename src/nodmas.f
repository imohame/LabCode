      subroutine nodmas(y,z,nodm,xmass,yz,nodd,xdamp)
c     implicit double precision (a-h,o-z)                                    dp
c
c     read in concentrated nodal masses and dampers
c
      common/bk03/numdc,imassn,idampn,irller,penstf
      common/bk14/lfna(15),lfnt(6)
      dimension y(*),z(*),nodm(*),xmass(2,*),yz(2,*),nodd(*),xdamp(2,*)
      character*80 txts,mssg
c
c     n o d a l   m a s s e s
c
      if (imassn.eq.0) go to 20
c
      call header
      write(lfnt(2),50)
      do 10 in=1,imassn
      call gttxsg (txts,lcount)
      read(unit=txts,fmt=40,err=34) n,(xmass(i,in),i=1,2)
      yz(1,in)=y(n)
      yz(2,in)=z(n)
      nodm(in)=n
   10 write(lfnt(2),60) n,(xmass(i,in),i=1,2)
c
   20 if (idampn.eq.0) return
c
      call header
      write(lfnt(2),70)
      do 30 in=1,idampn
      call gttxsg (txts,lcount)
      read(unit=txts,fmt=40,err=36) n,(xdamp(i,in),i=1,2)
      nodd(in)=n
   30 write(lfnt(2),60) n,(xdamp(i,in),i=1,2)
c
      return
c
   34 write (unit=mssg,fmt=80) in
      call termin (txts,mssg,lcount,1)
   36 write (unit=mssg,fmt=90) in
      call termin (txts,mssg,lcount,1)
c
   40 format(i5,2e10.0)
   50 format(///' n o d a l   m a s s   d a t a   '//
     1          '   node     y-mass     z-mass')
   60 format(i7,6e11.3)
   70 format(///' n o d a l   d a m p e r   d a t a   '//
     1          '   node   y-damper   z-damper')
   80 format (' error reading nodal mass card #',i5)
   90 format (' error reading discrete damper card #',i5)
      end
