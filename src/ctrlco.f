      subroutine ctrlco(signum)                                         unix
!!!!!!!!      integer signum                                                    unix
!!!!!!!!      character*4 msg                                                   unix
!!!!!!!!      common/msg/msg                                                    unix
!!!!!!!!      character*1 letter                                                unix
!!!!!!!!      dimension letter(80)                                              unix
!!!!!!!!      save iopen                                                        unix
!!!!!!!!      data msg/'    '/                                                  unix
!!!!!!!!      data iopen/0/                                                     unix
!!!!!!!!      write( *,41)                                                      unix
!!!!!!!!      write( *,40)                                                      unix
!!!!!!!!      read ( *,30)(letter(i),i=1,50)                                    unix
!!!!!!!!      call enablc                                                       unix
!!!!!!!!      do 20 i=1,46                                                      unix
!!!!!!!!      if (letter(i).eq.' ') go to 20                                    unix
!!!!!!!!      letter(i+3)='.'                                                   unix
!!!!!!!!      msg=letter(i)//letter(i+1)//letter(i+2)//letter(i+3)              unix
!!!!!!!!      return                                                            unix
!!!!!!!!   20 continue                                                          unix
!!!!!!!!      return                                                            unix
!!!!!!!!   30 format(80a1)                                                      unix
!!!!!!!!   40 format(1x,'.enter sense switch:',$)                               unix
!!!!!!!!   41 format(2x,/)                                                      unix
      end                                                               unix
