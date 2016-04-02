      subroutine gtemp(told,temp1,tnew,temp2)
c     implicit double precision (a-h,o-z)                                    dp
c
c.... routine to get new temperatures
c
      common/bk08/kprint,nstep,ite,ilimit,newstf
      common/bk14/lfna(15),lfnt(6)
      common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
      common/bk32/nsref,nequit,time,timep,lprint,nprint
      common/cn1/numati,numnpi,numeli,nblk1,nslidi,ntslvi,ntmsri,
     1           nnpbi,nepbi,ncnpi
      common/cn3/ibbari,intgi,nmbfi,ithopi,ithcri,ithini,iengri,
     1           ijinti
      common/taux1/itopaz,ithadd
      common/main_block/ a(1)
      dimension temp1(*),temp2(*)
c
c.... update old temperature
      told=tnew
      numnp3=(1+itopaz)*numnpi
      do 10 j=1,numnp3
      temp1(j)=temp2(j)
   10 continue
c.... update for last step
      if(nstep.eq.ntime)return
      if(ithopi.gt.0)then
      knpc=igtpnt(63)
      kp=igtpnt(64)
      endif
c
c.... evaluate temperaure from load curve (ithopi=1)
      if(ithopi.eq.1)then
      tnew=time
      call teval1(temp2,a(knpc),a(kp),ithcri,time)
c
c.... evaluate temps from load curve and vectors (ithopi=2)
      elseif(ithopi.eq.2)then
      kbase=igtpnt(82)+1+numnp3
      kmult=kbase+numnpi
      tnew=time
      call teval2(temp2,a(knpc),a(kp),a(kbase),a(kmult),
     1 ithcri,time)
c
c.... read from topaz plotfiles (ithopi<0)
      elseif(ithopi.lt.0)then
   20 continue
      isave=ithadd
c.... read next state
      call teval3(tnew,temp2)
c.... interpolation (ithopi=-2)
      if(ithopi.lt.-1)then
      if (abs(tnew-time).lt.1.e-12) return
      if(time.gt.tnew)goto 20
      scale=(time-told)/(tnew-told)
      tnew=time
      do 40 j=1,numnp3
      dtemp=temp2(j)-temp1(j)
   40 temp2(j)=temp1(j)+scale*dtemp
      ithadd=isave
c.... interpolation (ithopi=-3)
      elseif(ithopi.eq.-3)then
      scale=float(nstep+1)/float(ntime)
      tnew=time
      do 50 j=1,numnp3
      dtemp=temp2(j)-temp1(j)
   50 temp2(j)=temp1(j)+scale*dtemp
      ithadd=isave
      endif
      endif
c
      return
      end
