      subroutine nodein (id,y,z,roller,reftem,numnp,ityp2d)
!!!!c     implicit double precision (a-h,o-z)                                    dp
!!!!c
!!!!c     read and print nodal point data
!!!!c
      use CN_Objects_manager
      common/bk03/numdc,imassn,idampn,irller,penstf
      common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
      common/bk14/lfna(15),lfnt(6)
      common/bk23/itemp,itherm,irtin
      common/bk29/numfrq,clengt
      dimension y(*),z(*),id(2,*),roller(*),reftem(*)
      character*80 txts,mssg
!!!!!      common/WMLthermal/thermalflag,thermalconstraint(40000),
!!!!!     >    Tinit(40000),Rqold(40000)
      
!      common/WMLthermalperiod/period(40000)
     
      integer period    !!!!thermalflag, thermalconstraint,
      data nold,ymin,ymax,zmin,zmax/0,1.e20,-1.e20,1.e20,-1.e20/
      real xcoord,ycoord,BC_Mech_flag
      real BC_thermal_t,BC_diffusion_c
      integer BC_thermal_flag,BC_diffusion_flag

       write(7777,*) '-- nodein.f'
      nprnt =0
      irller=0
      bcd=0.
      kn0=0
      do 10 i=1,numnp
      id(1,i)=0
   10 id(2,i)=0
!!!!      id(1:2,1:numnp)=0
!ck----for orion input by m.k.
!c      write(29,115) numnp,numelt,numnp,numelt
!ck----end
!c
   20 kn=kn0
      call gttxsg (txts,lcount)
      read(txts,*,err=75) n,y(n),z(n),BC_Mech_flag,BC_thermal_flag
     1      ,BC_thermal_t,BC_diffusion_flag,BC_diffusion_c   !, Tfl(n),period(n)
!!!!      thermalconstraint(n)=BC_thermal_flag
!!!!      Tinit(n)=BC_thermal_t
      
!----------call the CN manager to set BC for both thermal and diffusion ismail2016-02-17      
      Call CNSetBCdata(n,BC_thermal_t,BC_diffusion_c,BC_thermal_flag
     >                  , BC_diffusion_flag)

      !Thermal modification, WML 82210
!!!!      if (thermalflag.eq.0) then
!!!!      read(unit=txts,fmt=90,err=75) n,bcc,y(n),z(n)!,kn0,tref
!!!!!      y(n) = y(n)*1.
!!!!!	  z(n) = z(n)*1.
!!!!	  elseif (thermalflag > 0) then
!!!!      read(unit=txts,fmt=91,err=75) n,bcc,y(n),z(n),thermalconstraint(n)
!!!!     1      ,Tinit(n,1)   !, Tfl(n),period(n)
!!!!      endif
      
      if(irtin.eq.1) reftem(n)=tref
      roller(n)=BC_Mech_flag
!c
!c     check boundary condition code on axis of symmetry
!c
      if (ityp2d.ne.0) go to 30
      if (y(n).ne.0.0) go to 30
      if (BC_Mech_flag.eq.0.0) BC_Mech_flag=1.0
      if (BC_Mech_flag.eq.2.0) BC_Mech_flag=3.0
   30 if (ymin.gt.y(n)) ymin=y(n)
      if (ymax.lt.y(n)) ymax=y(n)
      if (zmin.gt.z(n)) zmin=z(n)
      if (zmax.lt.z(n)) zmax=z(n)
      if (nprnt.gt.0) go to 40
      nprnt=50
      call header
      write(lfnt(2),80)
   40 nprnt=nprnt-1
      write(lfnt(2),100) n,BC_Mech_flag,y(n),z(n),BC_thermal_flag,
     > BC_thermal_t,BC_diffusion_flag,BC_diffusion_c
      if (BC_Mech_flag.ge.0.0) iflag=0
      if (BC_Mech_flag.lt.0.0) iflag=1
      if (BC_Mech_flag.ne.bcd) iflag=1
      if (BC_Mech_flag.lt.0.0) BC_Mech_flag=-BC_Mech_flag
      
      if (BC_Mech_flag.eq.1.0) id(1,n)=1
      if (BC_Mech_flag.eq.2.0) id(2,n)=1
      if (BC_Mech_flag.eq.3.0) id(1,n)=1
      if (BC_Mech_flag.eq.3.0) id(2,n)=1
      
      if (kn.eq.0) kn=1
      if (nold.eq.0) go to 60
      num=(n-nold)/kn
      numn=num-1
      if (numn.lt.1) go to 60
      xnum=num
      dy=(y(n)-y(nold))/xnum
      dz=(z(n)-z(nold))/xnum
      if(irtin.eq.1) dt=(reftem(n)-reftem(nold))/xnum
      k=nold
      do 50 j=1,numn
      kk=k
      k=k+kn
      y(k)=y(kk)+dy
      z(k)=z(kk)+dz
      if(irtin.eq.1) reftem(k)=reftem(kk)+dt
      roller(k)=roller(kk)
      if (iflag.eq.1) go to 50
      id(1,k)=id(1,kk)
      id(2,k)=id(2,kk)
   50 continue
c
   60 nold=n
      bcd=BC_Mech_flag
      if (n.ne.numnp) go to 20
c
!!!!!!    ----------------------------testing  
!!!!        write(*,*)id(1,1:25)
!!!!        write(*,*)id(2,1:25)
      do 70 i=1,numnp
      absrol=abs(roller(i))
      if(absrol.eq.1..or.absrol.eq.2..or.absrol.eq.3.) roller(i)=0.
      if(roller(i).ne.0.) irller=1
   70 continue
c
      clengt=sqrt((zmax-zmin)**2+(ymax-ymin)**2)
      penstf=10000.*clengt*penstf/sqrt(float(numelt))
c
      return
c
   75 nold=nold+1
      write (unit=mssg,fmt=110) nold
      call termin (txts,mssg,lcount,1)
c
   80 format(///' n o d a l   p o i n t   d a t a '/
     > /' node',5x,' b.c.',10x,'y',10x,'z',10x,'BC',
     > 5x,'Tini',5x,'BC',5x,'Cini')
   90 format(i5,f5.1,2e13.5,i5,e13.5)
   91 format(i5,f5.1,2e13.5,i5,2e13.5,i5)
  100 format(i5,5x,f5.0,5x,e12.4,8x,e12.4,6x,i4,
     > 5x,e12.4,5x,i4,5x,e12.4)
  110 format(' error reading nodal cards, probably for node #',i6)
  115 format(i5,2x,i5,2x,i5,2x,i5)
      end
