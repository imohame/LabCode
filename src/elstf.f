      subroutine elstf (matp,fail)
c     implicit double precision (a-h,o-z)                                    dp
c
c.... compute stiffness and rhs
c
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
      common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
      common/bk08/kprint,nstep,ite,ilimit,newstf
      common/bk09/maxref,rhsn,rhsvn,cvtl,iteref,ectl,tolls
      common/bk12/ntlen
      common/bk14/lfna(15),lfnt(6)
      common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
      common/bk23/itemp,itherm,irtin
      common/bk29/numfrq,clengt
      common/bk32/nsref,nequit,time,timep,lprint,nprint
      common/bk34/aa(1)
      common/wdskfc/iesf,numelf,alm(8),astf(36)
      common/slar3/nsl,nsntl,nmntl,nslnmx,sltol,slhrd
      common/fissn0/maxneq,mwspac,ntpe0,ntpe1,nfissl(3)
      common/fissn1/melemt,nnns,ntpe2,n2g,llls
      common/fissn2/nwpblk,numblk,mwsusd,mxnepb,maxch,
     > matpr,mench,ifa(2)
      common/colht/icolht,iopta
      logical rezone
      common/rezone/rezone,nrzn,nctr,irzcnt,nrezon
      common/xcom0/imeth,interq,imess,istart,igraf
      common/cn3/ibbari,intgi,nmbfi,ithopi,ithcri,ithini,iengri,
     >           ijinti
      common/double/iprec,ncpw,unit
      common/main_block/ a(1)
	  common/mbsize/numelt2, numnp2, neq2
	  common /sigfrac/ sigmacrit0, sigmacrit1,sigmacrit2,sigmacrit3,
     >    DecayCount, f_decay, penalty,fractFlag
	  
      dimension matp(*),fail(*)
      logical icolht
      data icolht/.false./,mprntl/0/
      real sigmacrit0, sigmacrit1, sigmacrit2, sigmacrit3



!c        write(7777,*) '-- elstf.f'


      if (nstep.gt.2) then
      do 4 i=1,numelt
      if (fail(i).lt.0.24) then
      matp(i)=0
      endif
    4 continue
      endif
      call setary(fail,1.*unit,numelt2)
      iesf =0
      iopta=-4
      if(numfrq.ne.0.and.imass.ge.1)imass=0
!c
!c---- compute and assemble external load term ----------------------
!c
!c.... find value of load curves at current time
      if (.not.rezone)then
        call ldcset (a(k57+2*numnp2),a(k63),a(k64),time)
      endif
!c
!c.... add concentrated nodal load contribution to a(k19)
      call loadcn (a(k57),a(k18),a(k19),a(k63),a(k64),a(k58),a(k59),a
     1 (k60),a(k61),a(k62),time)
!c
!c.... add pressure load contribution to a(k19)
      call loadpr (a(k57),a(k18),a(k19),a(k63),a(k64),a(k65),a(k66),
     1 a(k67),a(k68),a(k69),a(k70),time)
!c
!c.... zero external load term a(k19) for prescribed displacement dof
      if(numdc.ne.0)then
        call zrorhs(a(k57),a(k72),a(k19),a(k73),numdc)
      endif
!c
!c---- compute temperatures --------------------------------------
!c
      if(ithopi.ge.-3.and.ithopi.ne.0)then
      call gtemp(a(k81),a(k81+1),a(k82),a(k82+1))
      endif
!c
!c-----------------------------------------------------------------
      call azero (a(k17),neq)
      melemt=0
!c
!c---- only for new stiffness form and factor----------------------
!c
      if (newstf.eq.0)then
!c.... for finite elements find column heights
      icolht=.true.
      call blkcpy (a(k16),a(k13),nwebuf)

!cw      call ovrlay ('crystal2d',3,0,'recall')
!!!        call ovrlay (3,0)
        call fem2dm()


!c
!c.... for slideline elements find column heights, lhs, and rhs,
!c     and assemble rhs
!ck      if(nsl.ne.0)call slidln
      icolht=.false.
!c
!c.... determine matrix profile / block structure
      call bsolvr (a(k17),a(ntlen),1,4)
!c.... print information on matrix profile / block structure
      mprnt=mwsusd+matpr+maxch+mench+numblk
      if(mprnt.ne.mprntl)then
      mprntl=mprnt
      if(imeth.eq.0)write(lfnt(2),71)nstep
      if(ioofc.eq.1)then
      write(lfnt(2),61)neq,matpr,maxch,mench
      else
      write(lfnt(2),62)neq,matpr,maxch,mench,mwspac,mwsusd,numblk
      endif
      if(timep.eq.0.)then
c*waeil*      if(imeth.eq.0) write(*,71)nstep
      if(ioofc.eq.1)then
c*waeil*      write(*,61)neq,matpr,maxch,mench
      else
c*waeil*      write(*,62)neq,matpr,maxch,mench,mwspac,mwsusd,numblk
      endif
      endif
      endif
      if (ioofc.eq.1)mwspac=-mwsusd
c
      endif
c
c---- in any case -----------------------------------------------
c
c.... compute and assemble finite element contributions
      lpar(2)=lpar(4)
c.... define 8 dof for finite elements
      nnns=8
c.... define 8+36=44 words in finite element lhs array
      llls=44
c.... compute finite element rhs and assemble and if necessary
c     compute lhs and assemble
      call blkcpy(a(k16),a(k13),nwebuf)

cw      call ovrlay('crystal2d',3,0,'recall')
!!!      call ovrlay(3,0)
        call fem2dm()


      if (.not.rezone)
     1call blkcpy(a(k13),a(k16),nwebuf)
c
c.... j-integral
c     if(ijinti.ne.0)call jmast
c
c.... compute and assemble concentrated nodal masses,
c     dampers, and rollers
      if(imassn+idampn+irller.ne.0)then
      call ldmas (a(k19),a(k18),a(k18),a(k55),a(k56),a(k57),a(k76),
     1 a(k77),a(k78),a(k57+2*numnp2),a(k79),a(k80),0,a(k71),
     2 aa(1),aa(1),iprec)
      endif
c
      if (lprint.eq.0) call header
c
c.... compute eigenvalues if necessary

      if(numfrq.ne.0) then

cw        call ovrlay ('crystal2d',4,0,'recall')
          call ovrlay(4,0)
      endif


c
c     calculate norm of rhs a(k19) (the rhs (residual) in this
c     case is the incremental external load)
      rhsn=dotprd(a(k19),a(k19))
c
      return
   61 format('     in-core solution of equations'/
     1       '     number of equations',t32,i7,/,
     2       '     words in stiffness matix',t32,i7,/,
     3       '     maximum column height',t32,i7,/,
     4       '     average column height',t32,i7)
   62 format('     out-of-core solution of equations'/
     1       '     number of equations',t32,i7,/,
     2       '     words in stiffness matix',t32,i7,/,
     3       '     maximum column height',t32,i7,/,
     4       '     average column height',t32,i7,/,
     5       '     working space available',t31,i8,/,
     6       '     working space used',t32,i7,/,
     7       '     number of blocks',t32,i7)
   71 format(//'   equation profile information at step',i5)
      end









