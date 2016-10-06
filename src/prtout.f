      subroutine prtout
c     implicit double precision (a-h,o-z)                                    dp
c
	  parameter (nume=40000)
	  parameter (nume2=20000)
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
      common/bk05/ifil,iadd,maxsiz,head(12)
      common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
      common/bk08/kprint,nstep,ite,ilimit,newstf
      common/bk14/lfna(15),lfnt(6)
      common/bk16/maxint,hgc
      common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
      common/bk18/nummat,ityp2d,ako(31)
      common/bk20/ntotal
      common/bk23/itemp,itherm,irtin
      common/bk32/nsref,nequit,time,timep,lprint,nprint
      common/bk34/bb(1)
      common/newz1/ibase,iadd5,iadd6,icnt6,locst6
      common/countstotal/numnptotal, numeltotal
      common/main_block/ a(1)
	  common /pcracktip/connect(4,nume2), node(2,nume2), penta(nume2),
     >	                 ndflag(2,nume2), numnpt, numeltu, ndc
	  common/wblock8/  abc(573,nume,4), his(573,nume,4)
      
      dimension tim(5)
	  
	  integer nummat, numnp, numelt, numnpt, numeltu

      if (mprint.le.0) go to 10
      lprint=0

      call hspwr (a(k18),a(k55),a(k56),a(k57),a(k81+1),a(k03)
     >,a(k04))
      
      !-- to write file 29, the elem connectivity and coordinaes
      CALL EC_WriteConnectivity(a(k02),a(k03),a(k04),a(k18)
     >,a(k57),a(k08))
!!     CALL EC_WriteConnectivity(ElemConnect,NodesCoordx, NodesCoordy,NodesDispl,DofIds,ElemMaterial)
!!!!c     calculate nodes position for current step
!!!      connect=0
!!!      penta=0
!!!      ndflag=0
!!!      node=0.0
!!!      call datapost(a(k02))
!!!	  
!!!      write(29,40) nummat, numnp, numelt, numnpt, numeltu
!!!      if (nstep.eq.0) then
!!!          numnptotal=0
!!!          numeltotal=0
!!!      endif
!!!
!!!!----------------------------ismail20150715
!!!      numnptotal=numnptotal+numnp
!!!      numeltotal=numeltotal+numelt
!!!      write(7017,*) numnptotal,numeltotal
!!!!      -- dump the entire history matrix for post-processing
!!!      write(7021)his(1:573,1:numelt,1)
!!!!----------------------------------------------       
!!!	  call outconnect(a(k08))
      

      go to 20
   10 nprint=0
      tim(1)=timep
      ifctr=iadd/maxsiz
      icnt=5+nummat+8*numnp+(21+idump)*numelt
      if (iadd+icnt.gt.(ifctr+1)*maxsiz) iadd=(ifctr+1)*maxsiz
      iaddsv=iadd
      iadd=iadd+5+nummat
!c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      if(ibase.ge.1)then
      ifctr6=iadd6/maxsiz
      icnt6=2+maxint*lpar(9)*numelt
      if(iadd6+icnt6.gt.(ifctr6+1)*maxsiz)iadd6=(ifctr6+1)*maxsiz
      endif
!c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

   20 nnel=0
      nprnt=0
      locstr=0
!c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      if(kprint.gt.0)then
!c.... write first 2 words of extended database
      if(ibase.ge.1)then
      call wrabsg(lfna(11),timep,1,iadd6,bb(ntotal+1),0)
!!!!!!!!!!!!!      call riosta(lfna(11))
      iadd6=iadd6+1
      maxlpr=maxint*lpar(9)
      call wrabsg(lfna(11),maxlpr,1,iadd6,bb(ntotal+1),1)
!!!!!!!!!!!!!      call riosta(lfna(11))
      iadd6=iadd6+1
      endif
      endif
!c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      call blkcpy (a(k16),a(k13),nwebuf)
      call prtstr (a(k05),a(k08),a(k10))
c
      if (kprint.le.0) go to 30
c
      lce=k16+lpar(9)-1
      if(ityp2d.le.1)then
        lct=lce
      else
        lct=lce-1
      endif
      call outdv4 (a(k18),a(k55),a(k56),a(k57),numnp,a(k13),2*numnp)         nk
!c     call outdv4 (a(k18),a(k55),a(k56),a(k57),numnp,a(k13),                 pl
!c    1 2*numnp,a(k81+1))                                                     pl
      if (idump.ne.0) call blkcpy (a(k16),a(k13),nwebuf)
      call outeng (a(k18),a(k55),a(k56),a(k57),numnp,a(lce),a(k09),
     1 a(k08),a(k02),a(k10),a(k11),a(k81+1),tim(2),tim(3),tim(5),
     2 a(k13),a(k13),a(k03),a(k04),a(lct),lpar(9))
!!!!!!!!!!ck      call wrabsg(lfna(7),tim,5,iaddsv,bb(ntotal+1),0)
!!!!!!!!!!ck      call riosta(lfna(7))
!!!!!!!!!!      if (idump.ne.0) then
!!!!!!!!!!ck      call wrabsg(lfna(7),a(k13),nummat,iaddsv+5,bb(ntotal+1),0)
!!!!!!!!!!ck      call riosta(lfna(7))
!!!!!!!!!!      endif
!!!!!!!!!!c
   30 if (kprint.le.0) kprint=1
      if (mprint.le.0) mprint=1

   40 format(5i5)
      return

      end
