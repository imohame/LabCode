c     subroutine enable_ctrlc                                           vms
c     integer*4 tt_chan,sys$qiow,dvi$_devclass,item_list(3),            vms
c    1 sys$getdvi,dc$k_term                                             vms
c     common/ttchn/tt_chan                                              vms
c     external ctrlc_rout                                               vms
c     parameter (dvi$_devclass=4)                                       vms
c     parameter (dc$k_term=66)                                          vms
c     include '($iodef)'                                                vms
c     ibuf=0                                                            vms
c     ilen=0                                                            vms
c     item_list(1)=ior(ishft(dvi$_devclass,16),4)                       vms
c     item_list(2)=%loc(ibuf)                                           vms
c     item_list(3)=%loc(ilen)                                           vms
c     istat=sys$getdvi(,,'sys$input',item_list,,,,)                     vms
c     if(.not. istat) call lib$stop(%val(istat))                        vms
c     if(ibuf.ne.dc$k_term) return                                      vms
c     istat=sys$qiow(,%val(tt_chan),%val(io$_setmode.or.io$m_ctrlcast), vms
c    1 ,,,ctrlc_rout,,%val(3),,,)                                       vms
c     if(.not. istat) call lib$stop(%val(istat))                        vms
c     return                                                            vms
c     end                                                               vms
c     subroutine ctrlc_rout                                             vms
c     character*4 msg                                                   vms
c     common/msg/msg                                                    vms
c     character*1 letter                                                vms
c     dimension letter(80)                                              vms
c     data msg/'    '/                                                  vms
c     call enable_ctrlc                                                 vms
c     write(*,40)                                                       vms
c     read(*,30) (letter(i),i=1,80)                                     vms
c     do 10 i=1,80                                                      vms
c  10 call lwrcas(letter(i),letter(i))                                  vms
c     do 20 i=1,80                                                      vms
c     if(letter(i).eq.' ') go to 20                                     vms
c     letter(i+3)='.'                                                   vms
c     msg=letter(i)//letter(i+1)//letter(i+2)//letter(i+3)              vms
c     return                                                            vms
c  20 continue                                                          vms
c     return                                                            vms
c  30 format(80a1)                                                      vms
c  40 format($,1x,'enter sense switch:')                                vms
c     end                                                               vms
      subroutine intrup (mess,nwr,i)                                    vax750
      character*4 mess,msg                                              vax750
      common/msg/msg                                                    vax750
c     call mfeint                                                       nersc
      mess=msg                                                          vax750
      msg='    '                                                        vax750
      return                                                            vax750
      end                                                               vax750
