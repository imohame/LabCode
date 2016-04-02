c     subroutine memory(a,n)                                            vms
c     implicit double precision (a-h,o-z)                               vms  dp
c
c     this routine will check dimension of a array
c
c     common/array/maxa,maxadd,ifield                                   vms
c     common/bk14/lfna(15),lfnt(6)                                      vms
c     dimension a(1)                                                    vms
c     data length/1/                                                    vms
c
c     length=length+n                                                   vms
c     if (length.lt.maxadd.and.n.gt.0) then                             vms
c     do 10 i=1,n                                                       vms
c  10 a(i)=0.                                                           vms
c     endif                                                             vms
c     if (length.lt.maxadd) return                                      vms
c     nshort=length-maxadd                                              vms
c     write(lfnt(2),20)nshort                                           vms
c     write(*,20)nshort                                                 vms
c  20 format(//' --- error --- common block is ',i5,' words short')     vms
c     call bye(2)                                                       vms
c     end                                                               vms
c     subroutine mmemry(a,n)                                            unics
c     dimension a(*)                                                    unics
c     data length/1/                                                    unics
c     length=length+n                                                   unics
c     call memory ('UC',n)                                              unics
c     end                                                               unics
c     subroutine memory(a,n)                                            nersc
c
c     this routine will check dimension of a array
c
c     common/array/maxa,maxadd,ifield                                   nersc
c     common/bk14/lfna(15),lfnt(6)                                      nersc
c     dimension a(1)                                                    nersc
c     data length/1/                                                    nersc
c
c     length=length+n                                                   nersc
c     if (length.lt.maxadd.and.n.gt.0) then                             nersc
c     do 10 i=1,n                                                       nersc
c  10 a(i)=0.                                                           nersc
c     endif                                                             nersc
c     if (length.lt.maxadd) return                                      nersc
c     nshort=length-maxadd                                              nersc
c     write(lfnt(2),20)nshort                                           nersc
c     write(*,20)nshort                                                 nersc
c  20 format(//' --- error --- common block is ',i5,' words short')     nersc
c     call adios(2)                                                     nersc
c     end                                                               nersc
      subroutine rdabsg(ibuff,a,len,kloc,ioerr,b,iop)
c     double precision a                                                     dp
c
c.... routine to read single precision plot files when either single
c     or double precision real variables are used
c
      common/double/iprec,ncpw,unit
      dimension a(*),b(*)
c
c.... read data single precision
      if(iprec.eq.1)then
      call rdabsf(ibuff,a,len,kloc,ioerr)
c
c.... read data for double precision
      else
      if(iop.eq.1)then
      call rdabsf(ibuff,a,len,kloc,ioerr)
      else
      ngroup=(len-1)/1000+1
      len1=1000
      ipoint=1
      do 10 i=1,ngroup
      jloc=kloc+(i-1)*1000
      if(i.eq.ngroup)len1=len-(ngroup-1)*1000
      call rdabsf(ibuff,b,len1,jloc,ioerr)
      do 20 j=1,len1
      a(ipoint)=dble(b(j))
      ipoint=ipoint+1
   20 continue
   10 continue
      endif
      endif
c
      return
      end
