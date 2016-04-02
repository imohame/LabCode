      subroutine exita(n)                                               vax750
c     implicit double precision (a-h,o-z)                                    dp
      common/bk14/lfna(15),lfnt(6)                                      vax750
	  integer::n
	  write(*,*) n
      write(lfnt(2),100) n                                              vax750
 100  format(2x,'program stop',i5)                                      vax750
c     call exit(n)                                                      nersc
      stop                                                              vax750
      end                                                               vax750
