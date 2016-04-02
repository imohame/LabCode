      subroutine ovrlay(number,iparm)

cw      subroutine ovrlay(name,number,iparm,ikeep)
c     implicit double precision (a-h,o-z)                                    dp
c
c
*      write(7777,*) '*** ovrlay.f'

      go to (10,20,30,40),number

   10 call nik2di
      return

   20 call initil
      return

   30 call fem2dm
      return

   40 call eigen

      return
      end
