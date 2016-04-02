      subroutine rdabsf (ibuff,iw,len,kloc,ioerr)                       vax750
      common/iobufs/iobuf(700,13)                                       vax750
      dimension iw(1)                                                   vax750
      call rdabs (iobuf(1,ibuff),iw,len,kloc,ioerr)                     vax750
      return                                                            vax750
      end                                                               vax750
