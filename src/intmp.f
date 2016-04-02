      subroutine intmp(told,temp1,tnew,temp2)
c     implicit double precision (a-h,o-z)                                    dp
c
c.... routine to get initialize temperature
c
      common/bk08/kprint,nstep,ite,ilimit,newstf
      common/bk32/nsref,nequit,time,timep,lprint,nprint
      common/cn1/numati,numnpi,numeli,nblk1,nslidi,ntslvi,ntmsri,
     1           nnpbi,nepbi,ncnpi
      common/cn3/ibbari,intgi,nmbfi,ithopi,ithcri,ithini,iengri,
     1           ijinti
      common/taux1/itopaz,ithadd
      common/main_block/ a(1)
      dimension temp1(*),temp2(*)
c
      tt=0.
      numnp3=(1+itopaz)*numnpi
      if(ithopi.gt.0)then
      knpc=igtpnt(63)
      kp=igtpnt(64)
      endif
c
c.... evaluate temperature from load curve (ithopi=1)
      if(ithopi.eq.1)then
c.... reference temperature used
      if(ithini.ge.1)then
        told=time
        tnew=time
        call blkcpy(temp1,temp2,numnpi)
      else
        told=tt
        call teval1(temp1,a(knpc),a(kp),ithcri,tt)
        tnew=tt
        call blkcpy(temp1,temp2,numnpi)
      endif
c
c.... evaluate temps from load curve and vectors (ithopi=2)
      elseif(ithopi.eq.2)then
      if(ithini.ge.1)then
        told=time
        tnew=time
        call blkcpy(temp1,temp2,numnpi)
      else
        told=tt
        kbase=igtpnt(82)+1+numnp3
        kmult=kbase+numnpi
        call teval2(temp1,a(knpc),a(kp),a(kbase),a(kmult),
     1   ithcri,tt)
        tnew=tt
        call blkcpy(temp1,temp2,numnpi)
      endif
c
c.... read from topaz plotfiles (ithopi<0)
      elseif(ithopi.lt.0)then
      ithadd=64+2*numnpi+5*numeli
      if(ithopi.gt.-3)then
        if(ithini.ge.1)then
          told=time
          tnew=time
          call blkcpy(temp1,temp2,numnpi)
        else
          call teval3(told,temp1)
          tnew=told
          call blkcpy(temp1,temp2,numnp3)
        endif
      else
        call teval3(tnew,temp2)
        told=tt
        tnew=tt
        call blkcpy(temp1,temp2,numnpi)
      endif
      endif
c
      return
      end
