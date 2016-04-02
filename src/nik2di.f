      subroutine nik2di
      use CN_Objects_manager
!c     implicit double precision (a-h,o-z)                                    dp
!c
      real*8 hed                                                        vax750
      common/bk00/
     1k01,k02,k03,k04,k05,k06,k07,k08,k09,k10,k11,k12,
     2k13,k14,k15,k16,k17,k18,k19,k20,k21,k22,k23,k24,
     3k25,k26,k27,k28,k29,k30,k31,k32,k33,k34,k35,k36,
     4k37,k38,k39,k40,k41,k42,k43,k44,k45,k46,k47,k48,
     5k49,k50,k51,k52,k53,k54,k55,k56,k57,k58,k59,k60,
     6k61,k62,k63,k64,k65,k66,k67,k68,k69,k70,k71,k72,
     7k73,k74,k75,k76,k77,k78,k79,k80,k81,k82,k83,k84,
     8k85,k86,k87,k88,k89,k90,k91,k92,k93,k94,k95,k96
      common/bk02/ioofc,iphase,imass,lpar(9)
      common/bk03/numdc,imassn,idampn,irller,penstf
      common/bk05/ifil,iadd,maxsiz,head(12)
      common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
      common/bk07/mbfc,nelpg,hed(12)
      common/bk08/kprint,nstep,ite,ilimit,newstf
      common/bk09/maxref,rhsn,rhsvn,cvtl,iteref,ectl,tolls
      common/bk10/npb,nodep(2,8)
      common/bk11/cnwmk(2),iequit,iprint,isref
      common/bk12/ntlen
      common/bk14/lfna(15),lfnt(6)
      common/bk16/maxint,hgc
      common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
      common/bk18/nummat,ityp2d,ako(31)
      common/bk20/ntotal
      common/bk23/itemp,itherm,irtin
      common/bk26/dt,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
      common/bk27/nlcur,nptst,nthpy,nthpz,nthps,xmy,xmz,xms,nload,nptm
      common/bk29/numfrq,clengt
      common/bk30/numlp,numpc,h22(2,2),pl2(2,2),h33(3,2),pl3(3,2)
      common/bk33/irfreq,krfreq,iress
      common/bk34/bb(1)
      common/slar3/nsl,nsntl,nmntl,nslnmx,sltol,slhrd
      common/fissn0/maxneq,mwspac,ntpe0,ntpe1,nfissl(3)
      common/fissn1/melemt,nnns,ntpe2,n2g,llls
      common/fissn2/nwpblk,numblk,mwsusd,mxnepb,maxch,matpr,mench,ifa(2)
      common/fissn3/ifissl,kfissl(3)
      common/riksw1/sp1,ds0,cs01,cs02,rrsp,linsch,igso,irco,idamp
      common/riksw2/rlnew,alfa0,dsx,iteopt,idctrl,riksf,numspu,mthunl
      common/automt/dtmin,dtmax,mxback,termtm
      logical rezone
      common/rezone/rezone,nrzn,nctr,irzcnt,nrezon
c     common/psa1/ihimem                                                     pl
      common/xcom0/imeth,interq,imess
      common/cn0/iconv,lpb,nrcc,icon,iband,idirw,ftlst
      common/cn3/ibbari,intgi,nmbfi,ithopi,ithcri,ithini,iengri,
     1           ijinti
      common/taux1/itopaz,ithadd
      common/main_block/ a(1)
!!!!      common/WMLthermal/thermalflag
	  common/mbsize/numelt2, numnp2, neq2
	  common/meshnum/ numnpo, numelto
      integer  numnpo, numelto !!!thermalflag,
      dimension istor(7),ipelem(2,8)
      logical nit
      character*8 namef
      character*2 nt
      character*4 auto
      character*80 txts,mssg
      
      ithadd=0
      




c
c     read control information
c
c.... header card
      iconv=0
      mssg=' error reading header card.'
      call gttxsg(txts,lcount)
	  
	  ! khalil 03/04/2007 modified this read statement to directly
	  ! suit matlab output "MAZOUT" file
      read (unit=txts,fmt=160,err=150) nt
	  hed =  8*0.
	  auto = 'K_ELK'
	  
c.... check for input conversion
      call getnam(lfnt(3),namef)
      if(namef.eq.'convert')then
        iconv=1
      endif
c
c.... new crystal2d format
      if(nt.eq.'91')then
         nit=.true.
      endif
	  
      call newnki(nt)
c
      numnpo=numnp
	  numelto=numelt
      nmelt5=5*numelt
	  numelt2=2*numelt
	  numnp2=2*numnp
      cvtl=cvtl*cvtl
c
c.... set time integration coefficients in case of dynamic problem
      if (dt.ne.0.0) call chgint
c
c     fixed element data
c
      call setpnt ( 1, 1,1)                                                  nk
c     ifield=ihimem                                                          pl
c     call setpnt(1,1,ihimem+1)                                              pl
      call setpnt ( 2, 1,numelt2   )
      call setpnt ( 3, 2,4*numelt2 )
      call setpnt ( 4, 3,numnp2    )
      call setpnt ( 5, 4,numnp2    )
      call setpnt ( 6, 5,numelt2   )
      call setpnt ( 7, 6,numelt2+(1-sign(1,maxint-3))*numelt )
      call setpnt (8,7,5*numelt2)                                             nk
c     call setpnt(8,7,6*numelt)                                              pl
      call setpnt ( 9, 8,2*numelt2 )
      call setpnt (10, 9,4*numelt2 )
      call setpnt (11,10,nummat   )
c.... allocate storage for den,thick,temmat,raystf)
      call setpnt (12,11,6*nummat)
c.... allocate storage for tabulated data pointers
      call setpnt (13,12,97*nummat)
      call setmem (13,20)
      mtbend=0

* ----------------------------------------->>>>>>>>>>>>>>>>>>
       write(7777,*) '-- nik2di.f'
       write(911,*) 'k01 : ', k01
       write(911,*) 'k02 : ', k02
       write(911,*) 'k03 : ', k03
       write(911,*) 'k04 : ', k04
       write(911,*) 'k05 : ', k05
       write(911,*) 'k06 : ', k06
       write(911,*) 'k07 : ', k07
       write(911,*) 'k08 : ', k08
       write(911,*) 'k09 : ', k09
       write(911,*) 'k10 : ', k10
       write(911,*) 'k11 : ', k11
       write(911,*) 'k12 : ', k12
       write(911,*) 'k13 : ', k13
       write(911,*) 'k14 : ', k14
       write(911,*) 'k15 : ', k15
       write(911,*) 'k16 : ', k16
       write(911,*) 'k17 : ', k17       
       write(911,*) 'k18 : ', k18
       write(911,*) 'k19 : ', k19
       write(911,*) 'k20 : ', k20
       write(911,*) 'k21 : ', k21
       write(911,*) 'k22 : ', k22
       write(911,*) 'k23 : ', k23
       write(911,*) 'k24 : ', k24
       write(911,*) 'k25 : ', k25
       write(911,*) 'k26 : ', k26
       write(911,*) 'k27 : ', k27
       write(911,*) 'k28 : ', k28
       write(911,*) 'k29 : ', k29
       write(911,*) 'k30 : ', k30
       write(911,*) 'k31 : ', k31
       write(911,*) 'k32 : ', k32
       write(911,*) 'k33 : ', k33
       write(911,*) 'k34 : ', k34
       write(911,*) 'k35 : ', k35
       write(911,*) 'k36 : ', k36
       write(911,*) 'k37 : ', k37
       write(911,*) 'k38 : ', k38
       write(911,*) 'k39 : ', k39
       write(911,*) 'k40 : ', k40
       
       close(911)
* -----------------------------------------<<<<<<<<<<<<<<<<<<

c
c     read in material properties
c
      call matin (a(k10),a(k11),a(k11+nummat),a(k11+2*nummat),
     1 a(k12),idump,itopaz,a(k12+96*nummat),a(k12+97*nummat),mtbend,
     2 a(k11+3*nummat))
!!!!!!!    matin (matype,den   ,thick        ,temmat         
!!!!!!,prop  ,idump,md18fl,mtball          ,tbarr           ,mtbend,raystf)

       write(*,*) '-----after  matin'
!!c
!!c     set element buffer length
!!c
      nwebuf=max(maxint*lpar(9)*numelt2,27*(nsntl+nmntl))
      nwebuf=max(nwebuf,6*numnp2)
      if (numfrq.ne.0) nwebuf=max(nwebuf,2*numnp2*(numfrq+8))
      nwebuf=max(nwebuf,2*(nlcur+1)*nptm+8*nload+1+nlcur)
c
      call setpnt (13,12,97*nummat+mtbend)
      call setpnt (14,13,0)
      call setpnt (15,14,nwebuf)
      call setpnt (16,13,nwebuf )
      call setpnt (17,16,0      )
      call setpnt (18,17,2*numnp2)
      call setpnt (19,18,3*nrcc )
      call setpnt (20,19,numnp2  )
      call setpnt (21,20,irtin*numnp2)
      call setpnt (22,21,0      )
      call setmem (22,2*numnp2+20)
!c
!c     read nodal data
!c
      call nodein (a(k17),a(k03),a(k04),a(k19),a(k20),numnp,ityp2d)
!!!!       nodein (id    ,y      ,z    ,roller,reftem,numnp,ityp2d)
       write(*,*) '-----after  nodein'
!!!!!!    ----------------------------testing  
!!!!       write(*,*)int(a(k17:k17+24))
!!!!       write(*,*)a(k03:k03+24)
!!!!       write(*,*)a(k04:k04+24)
       
!c
!c     read and store element data
!c
      lpar(7)=4

      call elemin (numelt,a(k02),a(k06),a(k07),a(k08),a(k10),a(k17),
     1 a(k22),a(k09),a(k20),a(k11+2*nummat),nit,nt)
!!!!       elemin(numel  ,ix    ,betan ,freep ,matp  ,matype,id    ,ifree,
!!!!    tref,reftem ,temmat,nit            ,nt)
       write(*,*) '-----after  elemin'
!!!!!!    ----------------------------testing  
!!!!        write(*,*)int(a(k08:k08+15))
       
        
     
!!!!      if (thermalflag > 0) then
!!!!!      call to read and apply any thermal load curve ...IM20150817
!!!!!      also it has to be called before thermalBC bc it modifies the Tinit to account for the first point in the curve
      call ReadThermalLoad()
!!!!!c     call program to set thermal BC's if necessary
!!!!!!!!!!!!!      call thermalBC(numelt,numnp)
!!!!!      ---- read the diff coeff tables if any
      call DiffCoeffTableRead()
!----------call the CN manager to set up the necessary data ismail2016-02-17
      call CN_UpdateElemConnet(a(k02))
      call CNSet_T_init()
!!!!      endif
       write(*,*) '-----after  thermalBC'
!c
!c     write heading, control data, and coordinates into plotfile
!c
      istor(1)=numnp*sign(1,-idump)
      istor(2)=0
      istor(3)=-numelt
      istor(4)=nummat-1100000
      istor(5)=ityp2d
      if (istor(5).eq.2) istor(5)=10001
      istor(6)=999999
      istor(7)=ntime+numfrq+120
c
ck      call wrabsg(lfna(7),hed,12,iadd,bb(ntotal+1),1)
ck      call riosta(lfna(7))
ck      iadd=iadd+12
ck      call wrabsg(lfna(7),istor,7,iadd,bb(ntotal+1),1)
ck      call riosta(lfna(7))
ck      iadd=iadd+7
ck      call wrabsg(lfna(7),a(k11),nummat,iadd,bb(ntotal+1),0)
ck      call riosta(lfna(7))
ck      iadd=iadd+nummat
ck      call wrabsg(lfna(7),hed,12,iadd,bb(ntotal+1),1)
ck      call riosta(lfna(7))
ck      iadd=iadd+17
      call coorot (a(k03),a(k04),a(k22),numnp)
ck      call dumpic (a(k10),a(k02),a(k08))
c
      call setpnt (53,22,0)
      if (nsl.eq.0) go to 20
c
      call setpnt (23,22,0)
      call setpnt (24,23,nsl)
      call setpnt (25,24,nsl)
      call setpnt (26,25,nsl)
      call setpnt (27,26,nsl)
c
c     read slideline control cards
c
c      call slcntl (a(k23),a(k24),a(k25),a(k26),nit)
c
      call setpnt (28,27,nsntl   )
      call setpnt (29,28,nmntl   )
      call setpnt (30,29,nsntl   )
      call setpnt (31,30,nmntl   )
      call setpnt (32,31,nsntl   )
      call setpnt (33,32,nmntl   )
      call setpnt (34,33,nsntl   )
      call setpnt (35,34,nmntl   )
      call setpnt (36,35,nsntl   )
      call setpnt (37,36,nmntl   )
      call setpnt (38,37,nsntl   )
      call setpnt (39,38,nmntl   )
      call setpnt (40,39,nsntl   )
      call setpnt (41,40,nmntl   )
      call setpnt (42,41,nsntl   )
      call setpnt (43,42,nmntl   )
      call setpnt (44,43,3*nsntl )
      call setpnt (45,44,3*nmntl )
      call setpnt (46,45,2*nsntl )
      call setpnt (47,46,2*nmntl )
      call setpnt (48,47,nsntl   )
      call setpnt (49,48,nmntl   )
      call setpnt (50,49,4*nslnmx)
      call setpnt (51,50,nsntl   )
      call setpnt (52,51,nmntl   )
      call setpnt (53,52,k47-k43 )
c
c     expand memory
c
      nn=20
      if (iband.eq.2) nn=nn+numnp
      call setmem (53,nn)
c
c     read and initialize slideline data
c
c      call slint0 (a(k23),a(k24),a(k25),a(k26),a(k27),a(k28),a(k29),a
c     1 (k30),a(k31),a(k32),a(k33),a(k34),nsl,a(k03),a(k04),a(k02),a
c     2 (k08),a(k10),a(k12),nummat,numelt,a(k45),a(k46),a(k47),a(k48),
c     3 a(k50),a(k51),a(k43),a(k44),a(k11+nummat))
c      call slint1 (a(k03),a(k04),a(k29),a(k30),a(k35),a(k36),a(k37),a
c     1 (k38),nsntl,nmntl)
c
c     read nodal print information
c
   20 continue
      call header
      write(lfnt(2),300)
      if (npb.ne.0) then
      call gttxsg (txts,lcount)
      read(unit=txts,fmt=170,err=150) nodep
      go to 30
      endif
      npb=1
      nodep(1,1)=1
      nodep(2,1)=numnp
   30 continue
c
      do 40 l=1,npb
   40 write(lfnt(2),260) l,nodep(1,l),nodep(2,l)
c
c     read element print information
c
      call header
      write(lfnt(2),310)
      if (lpb.ne.0) then
      call gttxsg (txts,lcount)
      read(unit=txts,fmt=170,err=150) ipelem
      go to 50
      endif
      lpb=1
      ipelem(1,1)=1
      ipelem(2,1)=numelt
      go to 60
   50 call epflag (a(k10),a(k05),ipelem,lpb,numelt)
   60 do 70 l=1,lpb
   70 write(lfnt(2),270) l,ipelem(1,l),ipelem(2,l)
c
c     read in nodal constraint cards
c
      call rigid (a(k18),nrcc)
c
c     minimize bandwidth if iband=1 and number equations
c
      call eqnum (a(k17),a(k03),a(k04),a(k02),a(k08),a(k53),numnp
     1 ,numelt,neq,iband,a(k18),nrcc,a(k19),idctrl,idirw)
c
      neq2=neq+numnp2
      k23old=k23
      k17old=k17
      k19old=k19
      k20old=k20
c
      call setpnt (54,16,nwebuf  )
      call setpnt (17,54,ioofc*2*ilimit*neq2)
      if (k17.eq.k54.and.numfrq.ne.0) call setpnt (17,54,2*ilimit*neq2)
ck      if (k17.eq.k54.and.mthsol.eq.6) call setpnt (17,54,neq2)
      call setpnt (18,17,neq2)
      call setpnt (19,18,neq2)
      call setpnt (20,19,neq2)
      call setpnt (21,20,neq2)
      call setpnt (22,21,neq2)
      call setpnt (55,22,neq2)
      call setpnt (55,55,neq2)
ck      if (mthsol.ge.8) call setpnt (55,55,3*neq)
      call setpnt (56,55,max(numnp2,neq2))
      call setpnt (57,56,max(2*numnp2,neq2))
      call setpnt (58,57,2*numnp2+2*nlcur)
      call setpnt (59,58,4*nload    )
      call setpnt (60,59,nload)
      call setpnt (61,60,nload)
      call setpnt (62,61,nload)
      call setpnt (63,62,nload)
      call setpnt (64,63,nlcur+1 )
      call setpnt (65,64,2*nlcur*nptm)
      call setpnt (66,65,3*numpc )
      call setpnt (67,66,3*numpc )
      call setpnt (68,67,4*numpc )
      call setpnt (69,68,3*numpc )
      call setpnt (70,69,numpc)
      call setpnt (71,70,numpc)
      call setpnt (72,71,irller*numnp)
      call setpnt (73,72,numdc)
      call setpnt (74,73,numdc)
      call setpnt (75,74,numdc)
      call setpnt (76,75,8*numdc)
      call setpnt (77,76,imassn  )
      call setpnt (78,77,2*imassn)
      call setpnt (79,78,2*imassn)
      call setpnt (80,79,idampn  )
      call setpnt (81,80,2*idampn)
      if(ithopi.lt.0.or.ithopi.eq.1)then
      call setpnt(82,81,numnp+1+itopaz*numnp)
      call setpnt(83,82,numnp+1+itopaz*numnp)
      elseif(ithopi.eq.2)then
      call setpnt(82,81,numnp+1+itopaz*numnp)
      call setpnt(83,82,numnp+1+(itopaz+2)*numnp)
      else
      call setpnt(83,81,0)
      endif
c.... j-integral allocation
      if(ijinti.eq.0)then
      call setpnt(23,83,0)
      else
      call setpnt(84,83,ijinti+1)
      call setpnt(85,84,ijinti)
      call setpnt(86,85,numnp)
      call setpnt(87,86,numelt)
      call setpnt(88,87,30*ijinti)
      call setpnt(89,88,ijinti)
      call setpnt(90,89,4*numelt)
      call setpnt(23,90,0)
      endif
   80 continue
c
c     storage for slidelines
c
      call setpnt (24,23,nsl     )
      call setpnt (25,24,nsl     )
      call setpnt (26,25,nsl     )
      call setpnt (27,26,nsl     )
      call setpnt (28,27,nsntl   )
      call setpnt (29,28,nmntl   )
      call setpnt (30,29,nsntl   )
      call setpnt (31,30,nmntl   )
      call setpnt (32,31,nsntl   )
      call setpnt (33,32,nmntl   )
      call setpnt (34,33,nsntl   )
      call setpnt (35,34,nmntl   )
      call setpnt (36,35,nsntl   )
      call setpnt (37,36,nmntl   )
      call setpnt (38,37,nsntl   )
      call setpnt (39,38,nmntl   )
      call setpnt (40,39,nsntl   )
      call setpnt (41,40,nmntl   )
      call setpnt (42,41,nsntl   )
      call setpnt (43,42,nmntl   )
      call setpnt (44,43,3*nsntl )
      call setpnt (45,44,3*nmntl )
      call setpnt (46,45,2*nsntl )
      call setpnt (47,46,2*nmntl )
      call setpnt (48,47,nsntl   )
      call setpnt (49,48,nmntl   )
      call setpnt (50,49,4*nslnmx)
      call setpnt (51,50,nsntl   )
      call setpnt (52,51,nmntl   )
      call setpnt (53,52,k47-k43 )
c
c     ensure sufficient storage for small eigenvalue problems
c
      ntlen=k53
      if (numfrq.ne.0) then
      nv=min(neq,numfrq+8)
      ne=3*nv*nv+5*nv+2*numnp+1+3*neq
      nc=ntlen-k17
      ntlen=k17+max(nc,ne)
      endif
      call expndm(ntlen)
c
!c.... in-core
      if(ioofc.eq.1)then
        mwspac=0
!c.... out-of-core
      else
        mwspac=ntlen-2*numnp
        ntpe2=lfna(3)
      endif
c
      if(nsl.ne.0)
     1call movvar(a(k23old),a(k23),k53-k23)
      if(irtin.eq.1)
     1call movvar(a(k20old),a(k81+1),numnp)
      if(irller.eq.1)
     1call movvar(a(k19old),a(k71),numnp)
      call movvar(a(k17old),a(k57),2*numnp2)
!!c
!!c     read and store loading information
!!c
      m6 =k16
      m7 =m6+nptm
      m8 =m7+nptm
      m9 =m8+4*nload
      m10=m9+nload
      m11=m10+nload
      m12=m11+nload
      m13=m12+nload
      m14=m13+nlcur+1
      m15=m14+2*nptm*nlcur
!!c
!!c     read in load curves
!!c
      call lodcvs (a(m6),a(m7),a(m13),a(m14),nlcur)
!!c
!!c     read in concentrated nodal loads
!!c
      call cnlods (a(k03),a(k04),a(m8),a(m9),a(m10),a(m11),a(m12))
!c
      if (nptst.eq.0) go to 90
      call blkcpy(a(m8),a(k58),nptst+nlcur+1+8*nload)
!c
   90 numnp1=numnp+1
!!c
!!c     generate pressure load data
!!c
  130 call loadp (a(k03),a(k04),a(k65),a(k66),a(k67),a(k68),a(k69),
     1 a(k70),nit)
!!c
!!c     generate displacement boundary condition cards
!!c
      call dbcint (a(k57),a(k72),a(k73),a(k74),a(k75),numdc)
!!c
!!c     generate body force load data
!!c
      call loadb
!!c
!!c     read in concentrated nodal masses and dampers
!!c
      call nodmas (a(k03),a(k04),a(k76),a(k77),a(k78),a(k79),
     1 a(k80))
!!c
!!c     initial conditions
!!c
      call inital(a(k03),a(k04),a(k18),a(k55),a(k56),a(k57),a(k16),icon)
!!c
!!c.... element body force
      if(nmbfi.ne.0)then
      call ebfin(a(k07))
      endif
!!c
!!c.... nodal temperatures
!!c
      if(ithopi.eq.2)then
      kbase=k82+1+(1+itopaz)*numnp
      kmult=kbase+numnp
      call tempin(a(kbase),a(kmult))
      endif
!c
!c.... j-integral
!c
!c     if(ijinti.ne.0)then
!c     call jinit(a(k10),a(k12),a(k03),a(k04),a(k83),a(k85),
!c    1 ityp2d,ijinti)
!c     endif
!c
      if(iconv.eq.1)call bye(1)
      nelpg=numelt
      call header
      write(lfnt(2),320) nelpg
!c*waeil*      write(*,320) nelpg                                          
!c
      
      return
!c
  150 call termin (txts,mssg,lcount,1)
!c
  160 format(a2)
  170 format(16i5)
  260 format(
     1/5x,'block ',i2
     2/8x,'first node of this block   .......................    =',i5
     3/8x,'last node of this block  .........................    =',i5)
  270 format(/5x,
     1/5x,'block ',i2
     2/8x,'first element of this block  .....................    =',i5
     3/8x,'last element of this block .......................    =',i5)
  300 format(///' n o d a l   p r i n t   b l o c k s ')
  310 format(///' e l e m e n t   p r i n t   b l o c k s ')
  320 format(///' number of elements per group =',i5)
      end
