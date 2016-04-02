       subroutine fisslb (iopt2,b,a,s)
c     implicit double precision (a-h,o-z)                                    dp

c....  in core fissle like solver based on actcol
c
       common/bk12/ntlen
       common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
       common/fissn0/maxneq,mwspac,ntpe0,ntpe1,nfissl(3)
       common/fissn1/melemt,nnns,ntpe2,n2g,llls
       common/fissn2/nwpblk,numblk,mwsusd,mxnepb,ma,nwk,mb,ifa(2)
       common/double/iprec,ncpw,unit
       dimension b(*),a(*),s(*)
       common /main_block/ aa(1)

       data nlc/0/,ilenm/1/,numblk/1/
c
       iopt=iopt2
       if (iopt.eq.-6) iopt=2
       if (iopt.lt.-3) go to 30
       if (abs(iopt).eq.3) go to 60
       if (abs(iopt).eq.2) go to 50
       ilenm=max(ntlen,ilenm)
       nn=ntlen+neq+1
       if(nn.gt.ilenm)then
       call expndm(nn)
       ilenm=nn
       endif
       neq1=neq+1
       do 10 i=1,neq1
   10 a(i)=0.
c
c      compute address of diagonal terms of stiffness matrix
c
       call addres (a,b,nwk,ma,mb)
c
c      compute storage
c
       nn=ntlen+neq+1+nwk
       if(nn.gt.ilenm)then
       call expndm(nn)
       ilenm=nn
       endif
       mwspac=ilenm-ntlen
       mwsusd=mwspac
       nwkn=nwk+neq+1
       neq2=neq+2
       do 20 i=neq2,nwkn
   20 a(i)=0.
c
       return
c
c      assemble stiffness matrix
c
c 30  call addeb(a(1),a(1+neq),s(1),s(1+nnns/2),nnns,llls-nnns/2,
c    1     2*llls-nnns,melemt)
   30 call addeb(a(1),a(1+neq),s(1),s(1+nnns),nnns,llls,
     1     llls*iprec,melemt)
c
       return
   50 call tridec (b,a,a(2+neq),neq)
c
       return
c
   60 call fwdbak (b,a,a(2+neq),neq)
c
       return
       end
