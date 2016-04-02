      subroutine memory(a,n)                                            wkstn
c     implicit double precision (a-h,o-z)                               wkstndp
c
c     this routine will check dimension of a array
c
      common/array/maxa,maxadd,ifield                                   wkstn
      common/bk14/lfna(15),lfnt(6)                                      wkstn
      dimension a(1)                                                    wkstn
      data length/1/                                                    wkstn
c
      length=length+n                                                   wkstn
      if (length.lt.maxadd.and.n.gt.0) then                             wkstn
      do 10 i=1,n                                                       wkstn
   10 a(i)=0.                                                           wkstn
      endif                                                             wkstn
      if (length.lt.maxadd) return                                      wkstn
      nshort=length-maxadd                                              wkstn
      write(lfnt(2),20)nshort                                           wkstn
      write(*,20)nshort                                                 wkstn
   20 format(//' --- error --- common block is ',i10,' words short')     wkstn
      call bye(2)                                                       wkstn
      end                                                               wkstn
