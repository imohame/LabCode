      subroutine bsolvr(c,a,iopt,is)
c     implicit double precision (a-h,o-z)                                    dp
c
c.... interface routine to call fissln equation solver
c
      common/bk02/ioofc,iphase,imass,lpar(9)
      common/bk14/lfna(15),lfnt(6)
      common/bk15/cpuio(3,12),cpuip(3,12)
      common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
      common/bk34/s(1)
      common/wdskfc/iesf,numelf,alm(8),astf(36)
      common/fissn0/maxneq,mwspac,ntpe0,ntpe1,nfissl(3)
      common/fissn1/melemt,nnns,ntpe2,n2g,llls
      common/double/iprec,ncpw,unit
      common /main_block/ b(1)
      dimension c(*),a(*)
      real time1,time2,elpt1
      call CPU_TIME(time1)     

c
cw      call timin (cpuio,cpuip,is,2)
c
c---- in-core option ---------------------------------------------
      if(ioofc.ne.0)then
c
c.... determine matrix profile/block structure, c(*), and
c     set memory starting at a(1)
      if(iopt.eq.1)then
      call fissln (-1,c,a,s)
      mwspac=iabs(mwspac)
c
c.... factorize coefficient matrix a(*), and solve ax=c with
c     rhs c(*) and solution x overwriting c(*)
      elseif(iopt.eq.2)then
      melemt=0
      call fissln (-6,c,a,s)
      call fissln(-3,c,a,s)
c
c.... solve ax=c
      else
      call fissln (-3,c,a,s)
      endif
c
c---- out-of-core option ----------------------------------------
      else

c.... clear first mwspac+maxneq words in blank common
      call wrabsf(lfna(6),b,(mwspac+maxneq)*iprec,0)
      call riosta(lfna(6))
c.... copy rhs c into b(1)
      call blkcpy(c,b,maxneq)
c
c.... determine matrix profile block/ structure
      if (iopt.eq.1)then
      call fissln(1,b,b(1+maxneq),s)
c
c.... factorize coefficint matrix and solve equations
      elseif(iopt.eq.2)then
c.... number of element dof
      nnns=8
c.... number of words in stiffness + equation pointers
      llls=44
      melemt=numelf
      call fissln(2,b,b(1+maxneq),s)
      call fissln(3,b,b(1+maxneq),s)
c
c.... solve equations
      else
      call fissln(3,b,b(1+maxneq),s)
      endif
c
c.... read back blank common
      call rdabsf(lfna(6),b(1+maxneq),mwspac*iprec,
     1 maxneq*iprec,ioerr)
      call riosta(lfna(6))
c.... place output into c(*)
      call blkcpy(b,c,maxneq)
      call rdabsf(lfna(6),b,maxneq*iprec,0,ioerr)
      call riosta(lfna(6))
      endif
c-----------------------------------------------------------------
c
cw      call timin (cpuio,cpuip,is,3)
      iesf=0
      call CPU_TIME(time2)     
      elpt1=time2-time1     
      write(7016,*)'- time of bsolvr(c,a,iopt,is)',iopt,is,elpt1
      return
      end
