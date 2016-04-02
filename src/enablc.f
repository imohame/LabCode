      subroutine enablc                                                 unix
!!!!!      external ctrlco                                                   unix
!!!!!      integer signum,flag                                               unix
!!!!!      integer signal                                                    unix
!!!!!
!!!!!      write(7777,*) '-- enablc.f'
!!!!!
!!!!!      signum=2                                                          unix
!!!!!      flag  =-1                                                         unix
!!!!!c      inum  =signal(signum,ctrlco,flag)                                 wkstn
!!!!!c     call fsigctl('REGISTER','SIGINT',ctrlco)                          unics
!!!!!      return                                                            unix
      end                                                               unix
