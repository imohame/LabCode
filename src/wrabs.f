      subroutine wrabs  (fit,w,nw,da)                                   vax750


c
c.... entry for random write
c
c     calling sequence: call wrabsf(fit,w,nw,da)
c
c     input arguments
c            fit      the file information table
c            w        data to be written to disk
c            nw       number of words to write to disk
c            da       zero base disk address
c
      implicit integer (a-z)                                            vax750
      dimension fit(8),w(nw)                                            vax750
c
      common/frfcm1/mxfrf,ifrf,buflen,fcsize,dskloc ,curlen,kop,ier     vax750
      character*8 frfn,frn,kfn                                          vax750
      common/frfcm2/frfn(2,16),frn,kfn                                  vax750
c
c.... get family size, family root name, and name of open family member
      ifrf=fit(2)                                                       vax750
      fcsize=fit(4)                                                     vax750
      frn=frfn(1,ifrf)                                                  vax750
      kfn=frfn(2,ifrf)                                                  vax750
c.... get buffer pointers
      buflen=abs(fit(3))                                                vax750
      dskloc =fit(5)                                                    vax750
      curlen=fit(6)                                                     vax750
      kop=1                                                             vax750
      l=nw                                                              vax750
      m=0                                                               vax750
c.... set up access to correct family member
   30 kd=da+m                                                           vax750
      call asgrfm (kd,fit)                                              vax750
      ll=min(l,fit(4)-kd)                                               vax750
c.... move data into buffer
   40 i=kd-dskloc                                                       vax750
      if (i.lt.0) go to 50                                              vax750
      bloc=i                                                            vax750
      blen=min(ll,buflen-i)                                             vax750
      i=0                                                               vax750
      go to 60                                                          vax750
   50 bloc=0                                                            vax750
      blen=ll+i                                                         vax750
      if (blen.gt.buflen) go to 80                                      vax750
c.... branch if no overlap
   60 if (blen.le.0) go to 80                                           vax750
      ll=blen                                                           vax750
      do 70 k=1,ll                                                      vax750
      fit(k+bloc+7)=w(k+m-i)                                            vax750
   70 continue                                                          vax750
      fit(3)=-buflen                                                    vax750
      if (i.lt.0) m=m-ll                                                vax750
      l=l-ll                                                            vax750
      m=m+ll                                                            vax750
c.... loop if all requested data has not been transferred
      if (l.ne.0) go to 30                                              vax750
      return                                                            vax750
c.... flush the buffer if data is present which is not on disk
   80 if (fit(3).lt.0) then                                             vax750
      fit(3)=-fit(3)                                                    vax750
      call wdiska (fit(1),fit(8),buflen,dskloc )                        vax750
      fit(7)=max(fit(7),dskloc +buflen)                                 vax750
      endif                                                             vax750
c.... for blocks larger than the buffer, write the data directly
   90 if (ll.lt.buflen) go to 110                                       vax750
      if (mod(kd,buflen).ne.0) go to 110                                vax750
      nr=ll/buflen                                                      vax750
      do 100 n=1,nr                                                     vax750
      call wdiska (fit(1),w(m+1),buflen,kd)                             vax750
      m=m+buflen                                                        vax750
      kd=kd+buflen                                                      vax750
  100 continue                                                          vax750
      fit(7)=max(fit(7),kd)                                             vax750
      l=l-nr*buflen                                                     vax750
      if (l.eq.0) return                                                vax750
      ll=ll-nr*buflen                                                   vax750
      if (ll.eq.0) go to 30                                             vax750
c.... initialize the buffer if required
  110 continue                                                          vax750
      dskloc =kd-mod(kd,buflen)                                         vax750
      fit(5)=dskloc                                                     vax750
      curlen=min(dim(fit(7),dskloc ),buflen)                            vax750
      if (curlen.lt.1) go to 120                                        vax750
      if ((kd.ne.dskloc ).or.(ll.lt.curlen)) then                       vax750
      call rdiska (fit(1),fit(8),curlen,dskloc )                        vax750
      endif                                                             vax750
  120 curlen=buflen                                                     vax750
      fit(6)=curlen                                                     vax750
      go to 40                                                          vax750
c
      end                                                               vax750
