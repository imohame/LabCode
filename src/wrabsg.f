      subroutine wrabsg(ibuff,a,len,kloc,b,iop)
c     double precision a                                                     dp
c
c.... routine to write single precision plot files when either single
c     or double precision real variables are used
c
      common/double/iprec,ncpw,unit
      dimension a(*),b(*)

c
c.... write data single precision
      if(iprec.eq.1)then
      call wrabsf(ibuff,a,len,kloc)
c
c.... write data for double precision
      else
      if(iop.eq.1)then
      call wrabsf(ibuff,a,len,kloc)
      else
      ngroup=(len-1)/1000+1
      len1=1000
      ipoint=1
      do 10 i=1,ngroup
      jloc=kloc+(i-1)*1000
      if(i.eq.ngroup)len1=len-(ngroup-1)*1000
      do 20 j=1,len1
      b(j)=real(a(ipoint))
      ipoint=ipoint+1
   20 continue
      call wrabsf(ibuff,b,len1,jloc)
   10 continue
      endif
      endif
c
      return
      end
