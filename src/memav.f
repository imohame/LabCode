      integer function memav(idum)                                      vax750
c
c.... this routine will return available memory
c
      common/array/maxa,maxadd,ifield                                   vax750
      common/main_block/ a(1)                                           vax750

      memav=maxa-ifield                                                 vax750
      return                                                            vax750
      end                                                               vax750
