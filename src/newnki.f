      subroutine newnki(nt)
!c     implicit double precision (a-h,o-z)                                    dp
!c
!c.... routine to read in control cards  to make sure input is correct
!c
      use mod_parameters
      use CN_Objects_manager
      
      real*8 hed                                                        
      character*80 txts,mesg
      character*2 nt
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
      common/cn0/iconv,lpb,nrcc,icon,iband,idirw,ftlst
      common/cn1/numati,numnpi,numeli,nblk1,nslidi,ntslvi,ntmsri,
     1           nnpbi,nepbi,ncnpi
      common/cn2/nlcuri,nlcmxi,ncnldi,npresi,ndbci,nbfri,nbfzi,
     1           nbfai,ncnmai,ncndpi,initci
      common/cn3/ibbari,intgi,nmbfi,ithopi,ithcri,ithini,iengri,
     1           ijinti
      common/cn4/itypei,ianali,ibwmni,icorei,ipcori,isolti,gamai,
     1           betai,raydmi
      common/cn5/delti,nstepi,ioprti,ioplti,irstri,irezi,mthsli,
     1 ditoli,entoli,nstrfi,nstiti,nitrfi,nrfsti
      common/cn6/nunldi,mthsui,idcndi,idcdri,iarcni,iadmpi,arcszi
      common/cn7/tollsi,tolsli,rftani,tolrzi,igeomi
      common/WMLthermal/thermalflag
      common/WMLBC/BCflag
!!!!!!	  common/couplinganalysis/ TDflag
      integer thermalflag,BCflag!!!!!!, TDflag
c


       WRITE(7777,*) '-- newnki.f'

!c.... read in problem and solution definition cards
      if(nt.eq.'91')then
      lcount=2
      mesg='error reading problem definition card 2'
      call gttxsg(txts,lcount)
!      read(unit=txts,fmt=10,err=100)numati,numnpi,numeli,nblk1,
      write(*,*)'before card2'
      read(txts,*,err=100)numati,numnpi,numeli,nblk1,
     1 nslidi,ntslvi,ntmsri,nnpbi,nepbi,ncnpi
     
!       write(*,*)'--------------'
!      write(*,*)numati,numnpi,numeli,nblk1,
!     1 nslidi,ntslvi,ntmsri,nnpbi,nepbi,ncnpi
      
      mesg='error reading problem definition card 3'
      call gttxsg(txts,lcount)
!      read(unit=txts,fmt=20,err=100)nlcuri,nlcmxi,ncnldi,npresi,
      read(txts,*,err=100)nlcuri,nlcmxi,ncnldi,npresi,
     1 ndbci,nbfri,nbfzi,nbfai,ncnmai,ncndpi,initci
!       write(*,*)'--------------'
!      write(*,*)nlcuri,nlcmxi,ncnldi,npresi,
!     1 ndbci,nbfri,nbfzi,nbfai,ncnmai,ncndpi,initci
     
      mesg='error reading problem definition card 4'
      call gttxsg(txts,lcount)
!      read(unit=txts,fmt=30,err=100)ibbari,intgi,nmbfi,ithopi,
      read(txts,*,err=100)ibbari,intgi,nmbfi,ithopi,
     1 ithcri,ithini,iengri,ijinti
!       write(*,*)'--------------'
!      write(*,*)ibbari,intgi,nmbfi,ithopi,
!     1 ithcri,ithini,iengri,ijinti
     
      mesg='error reading problem definition card 5'
      call gttxsg(txts,lcount)
!      read(unit=txts,fmt=40,err=100)itypei,ianali,ibwmni,icorei,
      read(txts,*,err=100)itypei,ianali,ibwmni,icorei,
     1 ipcori,isolti,gamai,betai,raydmi
!       write(*,*)'--------------'
!      write(*,*)itypei,ianali,ibwmni,icorei,
!     1 ipcori,isolti,gamai,betai,raydmi
     
      mesg='error reading solution definition card 1'
      call gttxsg(txts,lcount)
!      read(unit=txts,fmt=50,err=100)delti,nstepi,ioprti,ioplti,
      read(txts,*,err=100)delti,nstepi,ioprti,ioplti,
     1 irstri,irezi,mthsli,ditoli,entoli,nstrfi,nstiti,nitrfi,
     1 nrfsti
!       write(*,*)'--------------'
!      write(*,*)delti,nstepi,ioprti,ioplti,
!     1 irstri,irezi,mthsli,ditoli,entoli,nstrfi,nstiti,nitrfi,
!     1 nrfsti
     
      mesg='error reading solution definition card 2'
      call gttxsg(txts,lcount)
!      read(unit=txts,fmt=60,err=100)nunldi,mthsui,idcndi,idcdri,
      read(txts,*,err=100)nunldi,mthsui,idcndi,idcdri,
     1 iarcni,iadmpi,arcszi,thermalflag,BCflag!!!!!!!!!!!!!!!!,TDflag    !Thermal flag WML82210, BC flag WML 11811
!!!!      thermalflag   0=nothing  1= thermal  2=diff  3=both thermal and diff
     
!       write(*,*)'--------------'
!      write(*,*)nunldi,mthsui,idcndi,idcdri,
!     1 iarcni,iadmpi,arcszi,thermalflag,BCflag,TDflag    !Thermal flag WML82210, BC flag WML 11811
!!!!!!!!!!!            write(*,*)'TDflag',TDflag
      mesg='error reading solution definition card 3'
      call gttxsg(txts,lcount)
!      read(unit=txts,fmt=70,err=100)tollsi,tolsli,rftani,tolrzi,
      read(txts,*,err=100)tollsi,tolsli,rftani,tolrzi,
     1 igeomi
       write(*,*)'--------------'
      write(*,*)tollsi,tolsli,rftani,tolrzi,
     1 igeomi
      endif
c
c.... problem definition card 2
      nummat=numati
      
      numnp=numnpi
      numelt=numeli
      nsl=nslidi
      nsntl=ntslvi
      nmntl=ntmsri
      npb=nnpbi
      lpb=nepbi
      nrcc=ncnpi

c.... problem definition card 3
      nlcur=nlcuri
      nptm=nlcmxi
      nload=ncnldi
      numpc=npresi
      numdc=ndbci
      nthpy=nbfri
      nthpz=nbfzi
      nthps=nbfai
      imassn=ncnmai
      idampn=ncndpi
      icon=initci

!c.... problem definition card 4
      if(intgi.ne.1.or.ianali.gt.0)intgi=4
      maxint=intgi
      hgc=0.0001/16.*1.  !0.0025 for QS, hgc=alpha/16.*psi, suggested alpha=0.1, psi is fraction of critical damping for viscous, set to 1 for stiffness
      if(intgi.eq.4)hgc=0.
      if(ibbari.ne.1)ibbari=0
      if(intgi.eq.1.or.itypei.eq.2)ibbari=1
      ibar=ibbari+1
      if(ithopi.lt.-3.or.ithopi.gt.2)ithopi=0
      if(ithopi.lt.0)then
        itemp=ithopi
      else
        itemp=ithcri
      endif
      if(ithini.ne.1.or.ithopi.eq.0)ithini=0
      if(ithopi.eq.-3.and.ithini.eq.0)ithini=1
      irtin=ithini
      if(ianali.gt.0)then
        idump=0
      else
        idump=1
        if(iengri.lt.0.or.iengri.gt.2)iengri=0
      endif

!c.... problem definition card 5
      ityp2d=itypei
      if(ianali.lt.-2)ianali=0
      if(ianali.ge.-2.and.ianali.le.0)then
        imass=abs(ianali)
        numfrq=0
      else
        imass=1
        numfrq=ianali
      endif
      if(ibwmni.ne.1)ibwmni=0
      iband=ibwmni
      if(icorei.ne.1.and.icorei.ne.-1)icorei=0
      if(icorei.eq.1)then
        ioofc=0
      else
        ioofc=1
      endif
      if(ipcori.le.0.or.ipcori.gt.100)ipcori=25
      if(gamai.lt..5)gamai=.5
      cnwmk(1)=gamai
      cval=.25*((gamai+.5)**2)
      if(betai.lt.cval)betai=cval
      cnwmk(2)=betai

!c.... solution  definition card 1
      dt=delti
!      ------- set the DtCurrent in the module CN_Consts ---ismail2016-02-19
      DtCurrent=dt
!     -----------------------
      ntime=nstepi
      if(ioprti.le.0)ioprti=999999   !added 9 WML 3311
      iprint=ioprti
      if(ioplti.le.0)ioplti=1
      jprint=ioplti
      if(irstri.le.0)irstri=999999   !added 9 WML 3311
      irfreq=irstri
      if(irezi.le.0)irezi=999999     !added 9 WML 3311
      irzcnt=irezi
      if(mthsli.le.0.or.mthsli.gt.12)mthsli=1

!c.... acrlenth
      linsch=0
      termtm=dt*ntime
      mthsol=mthsli
      if(ditoli.le.0.)ditoli=1.e-3
      cvtl=ditoli
      if(entoli.le.0.)entoli=1.e-2
      ectl=entoli
      if(nstrfi.le.0)nstrfi=1
      isref=nstrfi
      if(nstiti.le.0)nstiti=1
      iequit=nstiti
      if(nitrfi.le.0)nitrfi=10
      if(ianali.gt.0)nitrfi=(ianali+9)/2
      ilimit=nitrfi
      if(nrfsti.le.0)nrfsti=15
      maxref=nrfsti
      maxref=15
!c.... solution definition card 2
      if(nunldi.lt.0)nunldi=0
      numspu=nunldi
      if(mthsui.le.0)mthsui=1
      mthunl=mthsui
      if(idcndi.le.0)idcndi=0
      idctrl=idcndi
      if((idcdri.gt.2.or.idcdri.le.1).and.idctrl.ne.0)then
        idctrl=0
      endif
      idirw=idcdri
      if(iarcni.ne.1)iarcni=0
      irco=iarcni
      if(iadmpi.ne.1)iadmpi=0
      idamp=iadmpi
      if(arcszi.lt.0.)arcszi=0.
      dsx=arcszi

!c.... solution definition card 3
      if(tollsi.le.0.)tollsi=.9
      tolls=tollsi
      if(tolsli.eq.0.)then
        tolsli=.001
      elseif(tolsli.lt.0.)then
        tolsli=0.
      endif
      sltol=tolsli
      if(rftani.eq.0.)then
        rftani=.01
      elseif(rftani.lt.0.)then
        rftani=0.
      endif
      slhrd=rftani
      tolrzi=ftlst
      if(igeomi.ne.1)igeomi=0
      igso=igeomi
c
      if(nt.eq.'91')then
      write(lfnt(2),109)
      write(lfnt(2),110)numati,numnpi,numeli,nblk1,
     1 nslidi,ntslvi,ntmsri,nnpbi,nepbi,ncnpi
      write(lfnt(2),120)nlcuri,nlcmxi,ncnldi,npresi,
     1 ndbci,nbfri,nbfzi
      write(lfnt(2),121)nbfai,ncnmai,ncndpi,initci
      write(lfnt(2),130)ibbari,intgi,nmbfi
      write(lfnt(2),131)ithopi,ithcri,ithini,iengri,ijinti
      write(lfnt(2),140)itypei,ianali
      write(lfnt(2),141)ibwmni,icorei,ipcori,isolti,gamai,betai,
     1 raydmi
      write(lfnt(2),149)
      write(lfnt(2),150)delti,nstepi,ioprti,ioplti,irstri,
     1 irezi,mthsli
      write(lfnt(2),151)ditoli,entoli,nstrfi,nstiti,nitrfi,
     1 nrfsti
      write(lfnt(2),160)nunldi,mthsui,idcndi
      write(lfnt(2),161)idcdri,iarcni,iadmpi,arcszi
      write(lfnt(2),170)tollsi,tolsli,rftani,tolrzi,igeomi
      endif
!----------call the CN manager to set up the necessary data ismail2016-02-17
!!!!      thermalflag   0=nothing  1= thermal  2=diff  3=both thermal and diff
      Call CNCreateObjects(int(numnpi),int(numeli)
     >                    ,int(thermalflag))
     
     
     
     
     
      return

  100 call termin(txts,mssg,lcount,1)
   10 format(i5,9i5)
   20 format(11i5)
   30 format(8i5)
   40 format(6i5,3e10.1)
   50 format(e10.3,i7,i7,4i5,2e10.1,4i5)
   60 format(6i5,e10.1,3i5)
   70 format(4e10.1,i5)
  109 format(//t25,'P R O B L E M     D E F I N I T I O N'//)
  110 format(
     1 ' number of materials ................................ ',i5,/,
     2 ' number of nodes .................................... ',i5,/,
     3 ' number of elements ................................. ',i5,/,
     4 ' not used ........................................... ',i5,/,
     5 ' number of slidelines ............................... ',i5,/,
     6 ' total number of slave nodes ........................ ',i5,/,
     7 ' total number of master nodes ....................... ',i5,/,
     8 ' number of nodal printout blocks .................... ',i5,/,
     9 ' number of element printout blocks .................. ',i5,/,
     & ' number of constrained nodal pairs .................. ',i5,/)
  120 format(
     1 ' number of load curves .............................. ',i5,/,
     2 ' maximum number of points defining any load curve ... ',i5,/,
     3 ' number of nodal loads or follower forces ........... ',i5,/,
     4 ' number of element sides with applied tractions ..... ',i5,/,
     5 ' number of displacement boundary conditions ......... ',i5,/,
     6 ' base acceleration in r-direction ................... ',i5,/,
     & '     eq.1: no                                         '   ,/,
     & '     eq.2: yes                                        '   ,/,
     7 ' base acceleration in z-direction ................... ',i5,/,
     & '     eq.1: no                                         '   ,/,
     & '     eq.2: yes                                        '     )
  121 format(
     8 ' body force loads due to angular velocity ........... ',i5,/,
     & '     lt.0: angular velocity about x-axis              '   ,/,
     & '     eq.0: no                                         '   ,/,
     & '     gt.0: angular velocity about z-axis              '   ,/,
     9 ' number of concentrated nodal masses ................ ',i5,/,
     & ' number of concentrated nodal dampers ............... ',i5,/,
     1 ' initial condition flag ............................. ',i5,/,
     & '     lt.0: initialize angular velocity                '   ,/,
     & '     eq.0: initialize nodal velocities to zero        '   ,/,
     & '     gt.0: initialize nodal velocities                '   ,/)
  130 format(
     1 ' element formulation flag ........................... ',i5,/,
     & '     eq.0: b-bar                                      '   ,/,
     & '     eq.1: old crystal2d formulation                  '   ,/,
     2 ' integration order flag ............................. ',i5,/,
     & '     eq.0: 2x2 gauss integration                      '   ,/,
     & '     eq.1: 1 point integration                        '   ,/,
     3 ' number of elements with body forces ................ ',i5)
  131 format(
     4 ' thermal option flag ................................ ',i5,/,
     & '     lt.0: temperatures read from disk file           '   ,/,
     & '     eq.0: no thermal effects                         '   ,/,
     & '     gt.0: temperatures are scaled by load curve      '   ,/,
     5 ' load curve number for temperature vs. time ......... ' i5,/,
     6 ' reference temperature flag ......................... ',i5,/,
     & '     eq.0: not specified on node cards                '   ,/,
     & '     eq.1: specified on node cards                    '   ,/,
     7 ' element dump flag .................................. ',i5,/,
     & '     eq.0: element energy is dumped in plot d.b.      '   ,/,
     & '     eq.1: element thickness is dumped in plot d.b.   '   ,/,
     & '     eq.2: element temperature is dumped in plot d.b. '   ,/,
     8 ' number of j-integral contours ...................... ',i5,/,
     & '     eq.0: no j-integral calculation                  '   ,/)
  140 format(
     1 ' geometry type ...................................... ',i5,/,
     & '     eq.0: axisymmetric                               '   ,/,
     & '     eq.1: plane strain                               '   ,/,
     2 ' analyis flag ....................................... ',i5,/,
     & '     eq.-2: dynamic analysis, statically initialized  '   ,/,
     & '     eq.-1: dynamic analysis                          '   ,/,
     & '     eq.0: quasistatic analysis                       '   ,/,
     & '     gt.0: number of eigenvalues to be extracted      '     )
  141 format(
     3 ' bandwidth minimization flag ........................ ',i5,/,
     & '     eq.0: no                                         '   ,/,
     & '     eq.1: yes                                        '   ,/,
     4 ' out-of-core solution flag .......................... ',i5,/,
     & '     eq.0: in-core solution                           '   ,/,
     & '     eq.1: out-of-core solution                       '   ,/,
     5 ' percent of computer memory to be used .............. ',i5,/,
     6 ' solution method .................................... ',i5,/,
     & '     eq.0: standard                                   '   ,/,
     & '     eq.1: island solution strategy                   '   ,/,
     7 ' newmark parameter gamma......................'      e12.5,/,
     8 ' newmark parameter beta ......................'      e12.5,/,
     9 ' mass proportional Rayleigh damping coeff ....'      e12.5,/)
  149 format(//t25,'S O L U T I O N     D E F I N I T I O N'//)
  150 format(
     1 ' time step size ..............................',     e12.5,/,
     2 ' number of time steps ............................... ',i6,/,
     3 ' step interval for printing ......................... ',i5,/,
     4 ' step interval for plotting ......................... ',i5,/,
     5 ' step interval for restart dumps .................... ',i5,/,
     6 ' step interval for automatic rezoning ............... ',i5,/,
     7 ' standard solution method flag ...................... ',i5,/,
     & '     eq.1: bfgs                                       '   ,/,
     & '     eq.2: broyden                                    '   ,/,
     & '     eq.3: davidon-fletcher-powell                    '   ,/,
     & '     eq.4: davidon symmetric                          '   ,/,
     & '     eq.5: modified newton                            '     ) 
  151 format(
     8 ' convergence tolerance on displacements ......',     e12.5,/,
     9 ' convergence tolerance on energy .............',     e12.5,/,
     & ' number of steps between reformations ............... ',i5,/,
     1 ' number of steps between iterations ................. ',i5,/,
     2 ' number of iterations between reformations .......... ',i5,/,
     3 ' number of reforms per step ......................... ',i5,/)
  160 format(
     1 ' number of arc length unloading steps ............... ',i5,/,
     2 ' arc length unloading method ........................ ',i5,/,
     & '     eq.1: bfgs                                       '   ,/,
     & '     eq.2: broyden                                    '   ,/,
     & '     eq.3: dfp                                        '   ,/,
     & '     eq.4: davidon                                    '   ,/,
     & '     eq.5: modified newton                            '   ,/,
     3 ' arc length displacement control flag ............... ',i5,/,
     & '     eq.0: displacement norm is used                  '   ,/,
     & '     gt.0: node number for displacement control       '     )
  161 format(
     4 ' d.o.f for nodal displacement control ............... ',i5,/,
     & '     eq.1: r - direction                              '   ,/,
     & '     eq.2: z - direction                              '   ,/,
     5 ' arc length constraint method ....................... ',i5,/,
     & '     eq.0: crisfields                                 '   ,/,
     & '     eq.1: ramms                                      '   ,/,
     6 ' arc length damping flag ............................ ',i5,/,
     & '     eq.0: no                                         '   ,/,
     & '     eq.1: yes                                        '   ,/,
     7 ' initial arc length size .....................',     e12.5,/)
  170 format(
     1 ' line search tolerance .......................',     e12.5,/,
     2 ' stiffness inserting tolerance ...............',     e12.5,/,
     3 ' reduction factor for friction ...............',     e12.5,/,
     4 ' rezoner least squares tolerance .............',     e12.5,/,
     5 ' geometric stiffness flag ........................... ',i5,/,
     & '     eq.0: no                                         '   ,/,
     & '     eq.1: yes                                        '   ,/)
      end
