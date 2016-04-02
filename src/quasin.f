      subroutine quasin(ui,usi,r,tvc2,u,tvc1,tt,ierr)
c     implicit double precision (a-h,o-z)                                    dp
c
c     quasi-newton iteration methods
c
c     bfgs method by matthies and strang ---  published in
c     international journal for numerical methods in engineering
c     volume 14, pages 1613-1626, 1979, broyden's method
c     davidon-fletcher-powell method, and modified newton
c
      parameter (nume=40000)
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
      common/bk08/kprint,nstep,ite,ilimit,newstf
      common/bk09/maxref,rhsn,rhsvn,cvtl,iteref,ectl,tolls
      common/bk12/ntlen
      common/bk13/xnorm0(6),xnormc(6)!,xnk2d(20)
      common/bk14/lfna(15),lfnt(6)
      common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
      character*4 mess,mesg
      common/bk25/mess
      common/bk35/nbfgs,numupd
      common/wdskfc/iesf,numelf,alm(8),astf(36)
      common/rhsevl/numrhs
      common/total/itrlas,irflas,irhlas,itrtot,irftot,irhtot
      common/xcom0/imeth,interq,imess
      common/double/iprec,ncpw,unit
      dimension ui(*),usi(*),tvc2(*),r(*),u(*),tvc1(*)
      common/main_block/ a(1)
	  common /overlapping/ intersec(4, nume), area_coeff(nume), update_flag
	  
	  integer update_flag
c
c---- initialize for quasinewton procedure ------------------------
c
      ite=0
      itr=0
      nbfgs=0
      iteref=0
      numupd=0
      numrhs=1
      icount=0
      iphase=3
      stpper=0.05
      if (numdc.ne.0) call zrorhs (a(k57),a(k72),r,a(k73),numdc)
      if(mthsol.ne.5)then
      call blkcpi (0,r,a,k54,neq)
      endif
      g0=dotprd(r,ui)
      gc=g0
      do 20 i=1,neq
      tvc2(i)=ui(i)
      ui(i)=0.0
   20 continue
      if (numdc.ne.0) call swtche (a(k57),a(k72),tvc2,ui,a(k73),numdc)
c
c------------------- B E G I N   I T E R A T I O N ------------------
c
   30 continue
c.... increment # of equilibrium iterations for present lhs
      ite=ite+1
c.... increment # of equilibrium itermations for step
      itr=itr+1
c.... initialize line search step size
      step=1.0
      call intrup  (mesg,nwr,1)
      if (mess.eq.'    ' .and.mesg.ne.'sw2.') mess=mesg
c
c---- evaluate rhs -------------------------------------------------
c
      call assmrs (ui,usi,tvc2,r,u,tvc1,tt,step,0,g)
c
c---- line search --------------------------------------------------
c
cc.... perform line search if necessary
c      if (abs(g).gt.(tolls*abs(g0)).and.g0*g.le.0.)then
c      call lnesch (ui,usi,r,tvc2,u,tvc1,tt,step,g,g0,rlnew,riksf,alfa)
cc
cc.... when line search step is small set flags for reform
c      if(step.le.stpper)then
c      if(ite.eq.1)icount=icount+1
c      if(ite.ne.1.or.icount.le.7)then
c      if(step.lt.stpper)then
c        if(imeth.eq.0)then
c        ite=ilimit+1
c        else
c        iline=1
c        endif
c      endif
c      endif
c      endif
c      endif
c
c---- update displacement, compute norms ----------------------------
c
      do 50 i=1,neq
      ui(i)=ui(i)+step*tvc2(i)
      usi(i)=ui(i)
      tvc1(i)=u(i)+ui(i)
   50 continue
      dn1=dotprd(tvc1,tvc1)
      unm=dotprd(tvc2,tvc2)
      tol=max(dn1,dn2)*cvtl
      nstp1=nstep+1
      rhsvn=dotprd(r,r)
      dinorm=sqrt(unm)*step
      if(mesg.eq.'sw2.')then
        write (*,120) tt,nstp1,ite,numupd,numrhs,
     1  iteref,step,gc,g0,unm,tol
        mesg='    '
      elseif(mesg.eq.'sw8.')then
        imess=1
        mesg='    '
      endif
c
c.... place norms in array format
c     original or largest displacement increment norm
      xnorm0(1)=sqrt(max(dn1,dn2))
c.... current displacement increment norm
      xnormc(1)=sqrt(unm)
c.... original energy
      xnorm0(2)=gc
c.... current energy
      xnormc(2)=g0
c.... original or largest residual norm
      xnorm0(3)=sqrt(rhsn)
c.... current residual norm
      xnormc(3)=sqrt(rhsvn)
      write(*,*) xnorm0(1),xnormc(1),gc,g0,xnorm0(3),xnormc(3)
c
c---- double secret convergence check --------------------------------
c
c.... solution is way out
      if (abs(g0/(gc+1.e-06)).gt.100000.0)then
      icase=0
      else
c.... convergence check if new residual less than old
      if(rhsvn.le.rhsn)then
      if((unm.le.tol).and.
     1   (abs(g0).le.ectl*abs(gc).or.unm/10000.0.le.tol))then
      icase=3
      else
      if(ite.gt.ilimit)then
        write(lfnt(2),140)
        if(iteref+1.gt.maxref)then
          ierr=1
          return
        endif
        icase=1
      else
        icase=2
      endif
      endif
c.... residual larger than old one
      else
      write(lfnt(2),160)ite
      rhsn=rhsvn
      if(iteref+1.gt.maxref)then
        ierr=1
        return
      endif
      icase=1
      endif
      mthsll=mthsol
      endif
c
      rhsn=max(rhsvn,rhsn)
c
c---- stop execution of program (icase=0) ------------------------
c
      if(icase.eq.0)then
      ierr=1
      return
c
c---- reform stiffness and solve (icase=1) -----------------------
c
      elseif(icase.eq.1)then
      iteref=iteref+1
      if(imeth.eq.0)write(lfnt(2),150)
      ite=0
      iesf=0
      newstf=0
      nbfgs=0
      numupd=0
      iphase=4
      call azero (a(k17),neq)
      call assmrs (ui,usi,r,r,u,tvc1,tt,step,1,g)
      call blkcpi (0,r,a,k54,neq)
      call blkcpy (r,tvc2,neq)
      call bsolvr (tvc2,a(ntlen),2,8)
      g0=dotprd(tvc2,r)
      newstf=1
      iphase=3
      go to 30
c
c---- stiffness update and solve (icase=2) -------------------------
c
      elseif(icase.eq.2)then
      if (mthsll.eq.1) call bfgs (tvc2,r,tvc1,step,g,g0,dinorm)
      if (mthsll.eq.2) call broy (tvc2,r,tvc1,step)
      if (mthsll.eq.3) call dfp (a,tvc2,r,tvc1,step,neq,ntlen,nbfgs,
     1 numupd)
      if (mthsll.eq.4) call davd(a,tvc2,r,tvc1,step,neq,ntlen,nbfgs,
     1 numupd)
      if (mthsll.eq.5) call modn (tvc2,r)
      if(mthsll.ne.5)then
        call blkcpi (1,a,r,k54+nbfgs,neq)
      endif
      g0=dotprd(tvc2,r)
      go to 30
c
c------------------- E N D   I T E R A T I O N ------------------
c
c---- solution has converged (icase=3) -------------------------
c
      else
      iphase=2
      dn2=max(dn1,dn2)
      nstp1=nstep+1
ck      write(lfnt(2),130) nstp1,itr,iteref,numrhs,gc,g0,unm,tol
      itrtot=itrtot+itr
      irftot=irftot+iteref
      irhtot=irhtot+numrhs
      itrlas=itr
      irflas=irflas+iteref
      irhlas=numrhs
      iteref=0
      return
      endif
c
c----------------------------------------------------------------
c
  120 format('equilibrium iterations at time=',e12.4,' cycle=',i5
     1      /'   equilibrium iteration number           =',i12
     2      /'   number of stiffness updates            =',i12
     3      /'   number of right hand side evaluations  =',i12
     4      /'   stiffness matrix reformation number    =',i12
     5      /'   step size from line search             =',e12.4
     6      /'   initial energy norm                    =',e12.4
     7      /'   current energy norm                    =',e12.4
     8      /'   displacement norm                      =',e12.4
     9      /'   required tolerance                     =',e12.4)
  130 format(///
     1      /' equilibrium iteration in time step       =',i5
     2      /' number of iterations to converge         =',i5
     3      /' number of stiffness reformations         =',i5
     4      /' number of right hand side evaluations    =',i5
     5      /' initial energy norm                      =',e12.4
     6      /' final energy norm                        =',e12.4
     7      /' final displacement norm                  =',e12.4
     8      /' required tolerance                       =',e12.4)
  140 format(///' iteration limit reached  maximum permitted  ')
  150 format(///' stiffness matrix will now be reformed ')
  160 format(///' residual load vector larger than incremental loads',
     1' after iteration ',i5)
      end
