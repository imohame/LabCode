      subroutine blkcpy(a,b,n)                                          vax750
c     implicit double precision (a-h,o-z)                                    dp
      dimension a(1),b(1)                                               vax750
      do 10 i=1,n                                                       vax750
  10  b(i)=a(i)                                                         vax750
      return                                                            vax750
      end                                                               vax750
