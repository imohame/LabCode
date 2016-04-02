      subroutine lwrcas(dest,sour)                                      vax750
      character*1 dest,sour                                             vax750
      data ia/65/,iz/90/,idist/32/                                      vax750
      is=ichar(sour)                                                    vax750
      if (is.ge.ia.and.is.le.iz) then                                   vax750
      id=is+idist                                                       vax750
      dest=char(id)                                                     vax750
      else                                                              vax750
      dest=sour                                                         vax750
      endif                                                             vax750
      return                                                            vax750
      end                                                               vax750
