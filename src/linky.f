      subroutine linky(names)                                            vax750
c     implicit double precision (a-h,o-z)                                    dp
c
c     reads and cracks execute line
c
      common/iobufs/iobuf(700,13)                                       vax750
      character*1 name(8),letter(80)                                    vax750
      character*8 names                                                 vax750
      character*80 lnkarg                                               unix
      common/args/lnkarg,numargs                                        wkstn
c     common/args/lnkarg,numargs                                        unics
      dimension names(*)                                                vax750
c
      write(7777,*) '-- link.f '
      
      if(numargs.eq.0)then                                              wkstn
c     if(numargs.eq.0)then                                              unics
      write(*,110)                                                      vax750
      read(*,100)letter                                                 vax750
      do 10 i=1,80                                                      vax750
   10 call lwrcas(letter(i),letter(i))                                  vax750
      else                                                              wkstn
      do 12 i=1,80                                                      wkstn
      call lwrcas(letter(i),lnkarg(i:i))                                wkstn
   12 continue                                                          wkstn
      endif                                                             wkstn
c     else                                                              unics
c     do 12 i=1,80                                                      unics
c     call lwrcas(letter(i),lnkarg(i:i))                                unics
c  12 continue                                                          unics
c     endif                                                             unics
c
      i=0                                                               vax750
   20 i=i+1                                                             vax750
      if (i.gt.80) go to 90                                             vax750
      if (letter(i).ne.'=' ) go to 20                                   vax750
      j=i                                                               vax750
   30 j=j-1                                                             vax750
      if (letter(j).eq.' ' ) go to 30                                   vax750
      l=i                                                               vax750
   40 l=l+1                                                             vax750
      if (letter(l).eq.' ' ) go to 40                                   vax750
      do 50 ilet=1,8                                                    vax750
      name(ilet)=' '                                                    vax750
   50 continue                                                          vax750
      name(1)= letter(l)                                                vax750
      m=l                                                               vax750
   60 m=m+1                                                             vax750
      if (letter(m).eq.' ' ) go to 70                                   vax750
      if (letter(m).eq.',' ) go to 70                                   vax750
      if (m-l+1.gt.8) go to 70                                          vax750
      name(m-l+1)= letter(m)                                            vax750
      go to 60                                                          vax750
   70 n=0                                                               vax750
      if (letter(j).eq.'i' ) n=18                                       vax750
      if(n.eq.0) go to 80                                               vax750
      names(n)=name(1)//name(2)//name(3)//name(4)//name(5)//name(6)     vax750
     1//name(7)//name(8)                                                vax750
      write(*,*)' in linky.f ',names(n)
   80 i=m                                                               vax750
      go to 20                                                          vax750
c
   90 if(names(18).ne.'xyz123')                                         vax750
     1open(unit=18,file=names(18),status='old',                         vax750
     1               form='formatted')                                  vax750
       ierr=setvbuf3f_local(18,1,100)
       write(*,*)' in linky.f ',names(18)
!      This is the result.out file
      open(unit=17,file=names(17),status='unknown',                     vax750
     1               form='formatted')                                  vax750
       ierr=setvbuf3f_local(17,1,100)
c
      return                                                            vax750
  100 format (80a1)                                                     vax750
  110 format(////,2x,'default file names :'//                           vax750
     1'  o=n2hsp    g=n2hsp    d=n2dump'//                              vax750
     2'  please define input file names or change defaults:')           vax750
      end                                                               vax750
