      subroutine getnam (n,name)                                        vax750
c     implicit double precision (a-h,o-z)                                    dp
      character*8  names, name                                          vax750
      common/filen/names(25)                                            vax750
      name=names(n)                                                     vax750
      return                                                            vax750
      end                                                               vax750
