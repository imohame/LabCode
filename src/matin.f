      subroutine matin (matype,den,thick,temmat,prop,idump,md18fl,
     >                  mtball,tbarr,mtbend,raystf)


c     implicit double precision (a-h,o-z)                                    dp
c
c.... read in material property cards
c
      common/bk02/ioofc,iphase,imass,lpar(9)
      common/bk03/numdc,imassn,idampn,irller,penstf
      common/bk11/cnwmk(2),iequit,iprint,isref
      integer iequit,iprint,isref
      common/bk14/ lfna(15),lfnt(6)
      common/bk16/ maxint
      common/bk18/ nummat,ityp2d,ako(31)
      common/bk23/ itemp,itherm,irtin
      common/bk49/ bulkmx,ncon(30)
      character*6  headng
      dimension    matype(*),den(*),thick(*),temmat(*),prop(48,*),
     >             mtball(*),tbarr(*),raystf(3,*),headng(12)
      character*80 txts,mssg
      
      common /wblock3/  density_ms, density_ims,thermalEnthalpy(1000)
!!!!!!!!      common/WMLthermal/thermalflag,thermalconstraint(40000),
!!!!!!!!     >    Tinit(40000),Rqold(40000)
      common/WMLthermal3/thermalki(1000),thermalhi(1000),
     >      thermalRoi(1000),thermalcpi(1000),thermalxi(1000)
     >     ,thermalDi(1000)
!!!!!	  common/couplinganalysis/ TDflag
!!!!!      integer TDflag
	  common/hydroembrittle/critfrac(1000), sigfrac0(40000), 
     >       sigfrac(40000),decfrac(40000)
	  common/hydroembrittle110/critfrac110(1000), sigfrac0110(40000)

!!!!!!!      integer thermalflag
	  real critfrac, critfrac110
    
      equivalence  (model,lpar(1))

      open(901,file='properties.out',status='unknown')
       ierr=setvbuf3f_local(901,1,100)

      write(7777,*) '-- matin.f -----------------------------'

c
      if (maxint.eq.1) ncon(4) = ncon(4) + 8
      if (maxint.eq.1) ncon(7) = ncon(7) + 8

      call header

      write(lfnt(2),110)
ck      write(lfnt(2),120)
      lpar(9) = 0
      bulkmx  = 0.
      md18fl  = 0
c
c.... loop over materials
c
      do 10 i = 1, nummat
         mssg=' error reading material control card'
         call gttxsg (txts,lcount)
c        write(7777,*) 'In the do ; txts: ',txts

         read(unit=txts,fmt=20,err=19) n,matype(n),den(n)!,thick(n),
!     >                                 temmat(n),raystf(1,n)

         write(901,*) 'txt (n,matype..)=> ',txts
         write(901,*) 'n          : ', n
         write(901,*) 'matype(n)  : ', matype(n)
		 

         mssg=' error reading material header'
         call gttxsg (txts,lcount)
         read(unit=txts,fmt=40,err=19) headng
ckk      if (ityp2d.eq.2.and.thick(n).eq.0.) thick(n)=1.0
c
c.... read in appropriate properties for material model
c--------------------------------
c.....Elastic-Isotropic Materials
c--------------------------------
         if (matype(n).eq.1) then
            do 1 j=1,6
              call gttxsg (txts,lcount)
              read (unit=txts,fmt=70,err=18) prop(j,n)
    1       continue

c----------------------------------
c.....Elastic-Orthotropic Materials
c----------------------------------
         elseif (matype(n).eq.2) then
            call gttxsg (txts,lcount)
            read (unit=txts,fmt=70,err=18) (prop(j,n),j=1,3)
            call gttxsg (txts,lcount)
            read (unit=txts,fmt=70,err=18) (prop(j,n),j=4,6)
            call gttxsg (txts,lcount)
            read (unit=txts,fmt=70,err=18)  prop(7,n)
            call gttxsg (txts,lcount)
            read (unit=txts,fmt=70,err=18)  prop(8,n)
            call gttxsg (txts,lcount)
            read (unit=txts,fmt=70,err=18) (prop(j,n),j=9,10)
            call gttxsg (txts,lcount)
            read (unit=txts,fmt=70,err=18)  prop(11,n)

c----------------------------------------
c.....Elastoplastic Materials (Von Mises)
c----------------------------------------
         elseif (matype(n).eq.3) then
            call gttxsg (txts,lcount)
            read (unit=txts,fmt=70,err=18) prop(1,n)
            call gttxsg (txts,lcount)
            read (unit=txts,fmt=70,err=18) prop(2,n)
            call gttxsg (txts,lcount)
            read (unit=txts,fmt=70,err=18) prop(3,n),prop(22,n)
            call gttxsg (txts,lcount)
            read (unit=txts,fmt=70,err=18) prop(4,n),prop(5,n)
            call gttxsg (txts,lcount)
            read (unit=txts,fmt=70,err=18) (prop(j,n),j=6,13)
            call gttxsg (txts,lcount)
            read (unit=txts,fmt=70,err=18) (prop(j,n),j=14,21)

!!!!!!c-------------------------------------------------
!!!!!!c.....Single Crystal (Double-Slip)/(Multiple-Slip)
!!!!!!c-------------------------------------------------
      elseif (matype(n).eq. 4) then
!!!!!      write(*,*) thermalflag
         write(901,*) '** inside matl#4; i= ',n,' **'
         
!!!         if (thermalflag.eq.0.or.thermalflag.eq.2) then
            call gttxsg (txts,lcount)
            read (unit=txts,fmt=70,err=18) prop(1,n),thermalki(n),
     >          thermalhi(n),thermalEnthalpy(n),thermalcpi(n),
     >          thermalxi(n),dum7,thermalDi(n)       !,etain(n)
!     set the material density to the thermalRoi
            thermalRoi(n)=den(n)
!!!!!!!     --in case of diffusion read the D coeff from the last real*8 in the line
!!!!!!!!!!!!!!!!            write(*,*)'TDflag',TDflag
!!!!!!            if (TDflag==1) then     !diffusion
!!!!!!                thermalki(n)=dum8
!!!!!!            endif
!!!!!!!!!!!!!!!!!            write(*,*) '---- thermalki(n) ',thermalki(n)
            
            write(901,*) '---- txts1 ',txts
            write(901,*) '---- prop1 ',prop(1,n) !!! Young's modulus

            call gttxsg (txts,lcount)
            read (unit=txts,fmt=70,err=18) prop(2,n)    !!!!possion's ratio
            write(901,*) '---- txts2 ',txts
            write(901,*) '---- prop2 ',prop(2,n)

            call gttxsg (txts,lcount)
            read (unit=txts,fmt=70,err=18) prop(3,n), critfrac(n), 
     >                                  critfrac110(n) !!!Fy,Sfrac,Sfrac on 110
            write(901,*) '---- txts3 ',txts
            write(901,*) '---- prop3 ',prop(3,n)

            call gttxsg (txts,lcount)  
            read (unit=txts,fmt=70,err=18) prop(4,n),dumxi  !!! m=sensitivity parameter ,nu -- thermal_factor = (tempr/temp)**rnu in matpoly
            write(901,*) '---- txts4 ',txts
            write(901,*) '---- m ',prop(4,n)
            write(901,*) '---- temp/tempr^ ',dumxi
                
            call gttxsg (txts,lcount)
            read (unit=txts,fmt=70,err=18) prop(5,n)    !!!
            write(7777,*) '---- txts5 ',txts
            write(7777,*) '---- prop5 ',prop(5,n)

            call gttxsg (txts,lcount)
            read (unit=txts,fmt=70,err=18) (prop(j,n),j=6,9)    !!!!--,---,---,mat. output interval
!!!!            write(*,*)dumxi
!!!!!---------I put in 40 b/c from 10 to 25 are used, see printm and setse1 routines
            prop(40,n)=dumxi
!            --- this is done to overwrite the mat. output interval with the iprint  ismail 2016 01 01
            prop(9,n)=iprint
            do 7878 io=6,9
7878            write(901,*) '---- prop',io,' ',prop(io,n)
                
!!            do j=1,9
!!                prop97(97*(n-1)+i)=prop(i,n)
!!            enddo
!!            prop97(48*(n-1)+10)=dumxi
!!!!        elseif (thermalflag.eq.1) then
!!!!            call gttxsg (txts,lcount)
!!!!            read (unit=txts,fmt=70,err=18) prop(1,n),thermalki(n),
!!!!     >          thermalhi(n),etain(n),ecin(n)
!!!!
!!!!            write(901,*) '---- txts1 ',txts
!!!!            write(901,*) '---- prop1 ',prop(1,n)
!!!!
!!!!
!!!!            call gttxsg (txts,lcount)
!!!!            read (unit=txts,fmt=70,err=18) prop(2,n)
!!!!            write(901,*) '---- txts2 ',txts
!!!!            write(901,*) '---- prop2 ',prop(2,n)
!!!!
!!!!            call gttxsg (txts,lcount)
!!!!            read (unit=txts,fmt=70,err=18) prop(3,n), critfrac(n),
!!!!     >                                  critfrac110(n)
!!!!            write(901,*) '---- txts3 ',txts
!!!!            write(901,*) '---- prop3 ',prop(3,n)
!!!!
!!!!            call gttxsg (txts,lcount)  
!!!!            read (unit=txts,fmt=70,err=18) prop(4,n)
!!!!            write(901,*) '---- txts4 ',txts
!!!!            write(901,*) '---- prop4 ',prop(4,n)
!!!!
!!!!            call gttxsg (txts,lcount)
!!!!            read (unit=txts,fmt=70,err=18) prop(5,n)
!!!!            write(901,*) '---- txts5 ',txts
!!!!            write(901,*) '---- prop5 ',prop(5,n)
!!!!
!!!!            call gttxsg (txts,lcount)
!!!!            read (unit=txts,fmt=70,err=18) (prop(j,n),j=6,9)
!!!!            do 7879 io=6,9
!!!!7879            write(901,*) '---- prop',io,' ',prop(io,n)
!!!!        endif


c------------------------------------------------
c.....Polycrystal (G.B.) --- Double-Slip
c------------------------------------------------
      elseif (matype(n).eq.5) then
      call gttxsg (txts,lcount)
      read (unit=txts,fmt=70,err=18) prop(1,n)
      call gttxsg (txts,lcount)
      read (unit=txts,fmt=70,err=18) prop(2,n)
      call gttxsg (txts,lcount)
      read (unit=txts,fmt=70,err=18) prop(3,n)
      call gttxsg (txts,lcount)
      read (unit=txts,fmt=70,err=18) prop(4,n)
      call gttxsg (txts,lcount)
      read (unit=txts,fmt=70,err=18) prop(5,n)
      call gttxsg (txts,lcount)
      read (unit=txts,fmt=70,err=18) (prop(j,n),j=6,9)
c----------------------------------------------
c.....Bicrystal (Sigma9 G.B.) --- Multiple-Slip
c----------------------------------------------
      elseif (matype(n).eq.6) then
      call gttxsg (txts,lcount)
      read (unit=txts,fmt=70,err=18) prop(1,n)
      call gttxsg (txts,lcount)
      read (unit=txts,fmt=70,err=18) prop(2,n)
      call gttxsg (txts,lcount)
      read (unit=txts,fmt=70,err=18) prop(3,n)
      call gttxsg (txts,lcount)
      read (unit=txts,fmt=70,err=18) prop(4,n)
      call gttxsg (txts,lcount)
      read (unit=txts,fmt=70,err=18) prop(5,n)
      call gttxsg (txts,lcount)
      read (unit=txts,fmt=70,err=18) (prop(j,n),j=6,9)

c-----------------------------------------------
c.....Bicrystal ---> Double-Slip
c-----------------------------------------------
      ELSEIF (matype(n).eq.7) THEN
       CALL gttxsg (txts,lcount)
       read (unit=txts,fmt=70,err=18) prop(1,n)
       call gttxsg (txts,lcount)
       read (unit=txts,fmt=70,err=18) prop(2,n)
       call gttxsg (txts,lcount)
       read (unit=txts,fmt=70,err=18) prop(3,n)
       call gttxsg (txts,lcount)
       read (unit=txts,fmt=70,err=18) prop(4,n)
       call gttxsg (txts,lcount)
       read (unit=txts,fmt=70,err=18) prop(5,n)
       call gttxsg (txts,lcount)
      read (unit=txts,fmt=70,err=18) (prop(j,n),j=6,8)

c------------------------------------------------
c.....Polycrystal --- Double-Slip
c------------------------------------------------
         elseif (matype(n).eq.8) then
            call gttxsg (txts,lcount)
            read (unit=txts,fmt=70,err=18) prop(1,n)
            call gttxsg (txts,lcount)
            read (unit=txts,fmt=70,err=18) prop(2,n)
            call gttxsg (txts,lcount)
            read (unit=txts,fmt=70,err=18) prop(3,n)
            call gttxsg (txts,lcount)
            read (unit=txts,fmt=70,err=18) prop(4,n)
            call gttxsg (txts,lcount)
            read (unit=txts,fmt=70,err=18) prop(5,n)
            call gttxsg (txts,lcount)
            read (unit=txts,fmt=70,err=18) (prop(j,n),j=6,9)

            
c.....END...............................................
c
         else
            write(lfnt(2),200)matype(n)
            call bye(2)
         endif
c
c.... update maximum number of history variables etc. required
         lpar(1)=matype(n)
         lpar(9)=max(ncon(model)+idump+max(0,ityp2d-1),lpar(9))

c.... print material properties and define constants (e.g elast. coef.)
      call printm (n,den(n),thick(n),temmat(n),raystf(1,n),
     1 prop(1,n),prop(1,n+nummat),headng)

c.... find equivalent bulk modulus for material n
      call blkmax(matype,n,prop,bulkmx,raystf(1,n))

   10 continue
c


*----------------------------------------------------------------------
*.... Checking  properties
*----------------------------------------------------------------------
      do j = 1, 4
         do i = 1, 48
            write(901, 2021) i, j, prop(i,j)
         end do
      end do

 2021 format(2x,'prop(',i3,',',i3,' ) = ',f16.9)  
*---------------------------------------------------------------------

       
      close(901)
      penstf=bulkmx
      return
c
   18 write (unit=mssg,fmt=130) i,matype(n)

   19 call termin (txts,mssg,lcount,1)
c

* ---------------------------------------------------------------

   20 format(2i7,4e15.5)
   40 format(12a6)
   70 format(8e10.1)
ck   71 format(16i5)
  110 format(///' m a t e r i a l   d e f i n i t i o n s'//
     1        ' material models                         '/
     2        '     eq.1: isotropic                     '/
     3        '     eq.2: orthotropic                   '/
     4        '     eq.3: elastoplastic ( von mises )   '/
     5        '     eq.4: single crystal (double-slip)  '/
     6        '     eq.5: bicrystal (sigma33a g.b.)     '/
     7        '     eq.6: bicrystal (sigma9 g.b.)       '/
     8        '     eq.7: bicrystal (double-slip)       '/
     9        '     eq.8: polycrystal (double-slip)     ')
  130 format(' error reading material constants for material #',i3,
     1' type=',i3)
  200 format(' **fatal error** material property input for',/,
     1       '  material',i5,' not coded (matin)')
      end









