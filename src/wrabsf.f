      subroutine wrabsf (ibuff,iw,len,kloc)                             vax750
      common/iobufs/iobuf(700,13)                                       vax750
      dimension iw(1)                                                   vax750

 
      call wrabs (iobuf(1,ibuff),iw,len,kloc)                           vax750

      return                                                            vax750
      end                                                               vax750
