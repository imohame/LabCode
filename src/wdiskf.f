      subroutine wdiskf(a,numel)
c     implicit double precision (a-h,o-z)                                    dp
      common/bk14/lfna(15),lfnt(6)
      common/fissn0/maxneq,mwspac,ntpe0,ntpe1,nfissl(3)
      common/fissn1/melemt,nnns,ntpe2,n2g,llls
      common/fissn2/nwpblk,numblk,mwsusd,mxnepb,maxch,matpr,mench,ifa(2)
      common/fissn3/ifissl,kfissl(3)
      common/wdskfc/iesf,numelf,alm(8),astf(36)
      common/double/iprec,ncpw,unit
      dimension a(*)


      if (iesf.eq.0) numelf=0
      numelf=numelf+numel
      if (numel.eq.0) return
      if (llls.ne.44) go to 10
      lenn=44*numel*iprec
      call wrabsf(lfna(3),a,lenn,iesf)
      call riosta (lfna(3))
      iesf=iesf+lenn
      return
   10 nnns1=nnns+1
      lenn=44*iprec
      icoun=0
      do 20 j=1,44
   20 alm(j)=0
      do 50 i=1,numel
      do 30 j=1,nnns
   30 alm(j)=a(icoun+j)
      do 40 j=nnns1,llls
   40 astf(j-nnns)=a(icoun+j)
      call wrabsf(lfna(3),alm,lenn,iesf)
      call riosta (lfna(3))
      iesf=iesf+lenn
   50 icoun=icoun+llls
      return
      end
