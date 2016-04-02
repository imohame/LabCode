      subroutine teval3(tnew,temp)
c     implicit double precision (a-h,o-z)                                    dp
c
c.... routine to read next temperature state from topaz plotfile
c
      common/bk14/lfna(15),lfnt(6)
      common/bk34/bb(1)
      common/cn1/numati,numnpi,numeli,nblk1,nslidi,ntslvi,ntmsri,
     1           nnpbi,nepbi,ncnpi
      common/taux1/itopaz,ithadd
      dimension temp(*)
c
      call rdabsg(lfna(4),tnew,1,ithadd,ioerr,bb,0)
      call riosta(lfna(4))
      if(tnew.gt.1.e20)then
      ifactr=(ithadd-1)/262144+1
      ithadd=262144*ifactr
      call rdabsg(lfna(4),tnew,1,ithadd,ioerr,bb,0)
      call riosta(lfna(4))
      endif
      ithadd=ithadd+5+numati
      call rdabsg(lfna(4),temp,numnpi,ithadd,ioerr,bb,0)
      call riosta(lfna(4))
      ithadd=ithadd+2*numnpi
      if(itopaz.eq.0)then
      ithadd=ithadd+2*numnpi
      else
      call rdabsg(lfna(4),temp(1+numnpi),2*numnpi,ithadd,ioerr,
     1 bb,0)
      call riosta(lfna(4))
      ithadd=ithadd+2*numnpi
      endif
      ifactr=(ithadd-1)/262144+1
      if (ithadd+5+numati+4*numnpi.gt.ifactr*262144) then
      ithadd=262144*ifactr
      endif
c
      return
      end
