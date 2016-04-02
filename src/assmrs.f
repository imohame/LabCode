      subroutine assmrs(ui,usi,tvc2,r,u,tvc1,tt,step,iopt,g)
c     implicit double precision (a-h,o-z)                                    dp
c
c.... compute lhs and rhs for equilibrium iterations
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
      common/bk12/ntlen
      common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
      common/bk34/aa(1)
      common/slar3/nsl,nsntl,nmntl,nslnmx,sltol,slhrd
      common/rhsevl/numrhs
      common/fissn0/maxneq,mwspac,ntpe0,ntpe1,nfissl(3)
      common/fissn1/melemt,nnns,ntpe2,n2g,llls
      common/fissn2/nwpblk,numblk,mwsusd,mxnepb,maxch,matpr,mench,ifa(2)
      common/colht/icolht,iopta
      dimension ui(*),usi(*),tvc2(*),r(*),u(*),tvc1(*)
      common/main_block/ a(1)
      common/kelin/keli(3)
      common/double/iprec,ncpw,unit
	  common/mbsize/numelt2, numnp2, neq2
      logical icolht
      equivalence (lpar(7),mxn)
c


c      write(7777,*)  '--- assmrs.f'

      call setary(a(k08+numelt2),1.*unit,numelt2)
c.... increment number of right hand sides
      numrhs=numrhs+1
c.... set flag for zeroing stiffness on first call to fissln
      iopta=-4
c
c---- update displacement vectors ----------------------------------
c
      if(iopt.ne.1)then
      do 10 i=1,neq
c.... update step displacement with line search times incremental
c     displacement
      ui(i)=ui(i)+step*tvc2(i)
      usi(i)=ui(i)
c.... update total displacement with step displacement
      tvc1(i)=u(i)+ui(i)
   10 continue
      endif
c
c---- compute and assemble external load term ----------------------
c
c.... add concentrated nodal load contribution to r
      call loadcn(a(k57),tvc1,r,a(k63),a(k64),a(k58),a(k59),a(k60),a
     1 (k61),a(k62),tt)
c
c.... add pressure load contribution to r
      call loadpr(a(k57),tvc1,r,a(k63),a(k64),a(k65),a(k66),a(k67),
     1 a(k68),a(k69),a(k70),tt)
c
c.... zero external load term r for prescribed displacement dof
      if(numdc.ne.0)then
        call zrorhs(a(k57),a(k72),r,a(k73),numdc)
      endif
c
c--------------------------------------------------------------------
      melemt=0
      lpar(2)=lpar(4)
c
c---- only for new stifness form and factor -------------------------
c
      if(iphase.eq.4)then
c.... for finite elements find column heights
      icolht=.true.
      call blkcpy (a(k16),a(k13),nwebuf)

cw      call ovrlay ('crystal2d',3,0,'recall')
        call ovrlay(3,0)

c
c.... for slideline elements find column heights, lhs, rhs, and
c     assemble rhs
      iphase=4
ck      if(nsl.ne.0)call slidln
      icolht=.false.
c
c.... determine matrix profile / block strucure
      call bsolvr (a(k17),a(ntlen),1,4)
      if (ioofc.eq.1) mwspac=-mwsusd
      iphase=4
c
c---- no stiffness to be formed -----------------------------------

c      else
c      if(nsl.ne.0)call slidln
      endif
c
c---- in any case -------------------------------------------------
c
c.... define 8 dof for finite elements
      nnns=8
c.... define 8+36=44 words in finite element lhs array
      llls=44
c.... compute finite element rhs and assemble and if necessary
c     compute lhs and assemble
      call blkcpy (a(k16),a(k13),nwebuf)

cw      call ovrlay ('crystal2d',3,0,'recall')
        call ovrlay (3,0)

c
c.... compute and assemble concentrated nodal masses, dampers,
c     and rollers
      if(imassn+idampn+irller.ne.0)then
      call ldmas (r,tvc1,ui,a(k55),a(k56),a(k57),a(k76),a(k77),a(k78),a
     1 (k57+2*numnp),a(k79),a(k80),1,a(k71),aa(1),aa(1),iprec)
      endif
c
c---- reset displacement vectors -----------------------------------
c
      if(iopt.ne.1)then
      do 70 i=1,neq
c.... set step displacement back to value at previous iteration
      ui(i)=ui(i)-step*tvc2(i)
   70 usi(i)=ui(i)
c.... compute norm of new residual*increment displacement
      g=dotprd(tvc2,r)
      endif
c
      return
      end







