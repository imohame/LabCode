      subroutine uertst (ier,name)
c                                  specifications for arguments
      integer            ier
      integer            name(1)
c                                  specifications for local variables
      integer            i,ieq,ieqdf,iounit,level,levold,nameq(6),
     *                   namset(6),namupk(6),nin,nmtb
      data               namset/1hu,1he,1hr,1hs,1he,1ht/
      data               nameq/6*1h /
      data               level/4/,ieqdf/0/,ieq/1h=/
c                                  unpack name into namupk
c                                  first executable statement
co      call uspkd (name,6,namupk,nmtb)
c                                  get output unit number
      call ugetio(1,nin,iounit)
c                                  check ier
      if (ier.gt.999) go to 25
      if (ier.lt.-32) go to 55
      if (ier.le.128) go to 5
      if (level.lt.1) go to 30
c                                  print terminal message
      if (ieqdf.eq.1) write(iounit,35) ier,nameq,ieq,namupk
      if (ieqdf.eq.0) write(iounit,35) ier,namupk
      go to 30
    5 if (ier.le.64) go to 10
      if (level.lt.2) go to 30
c                                  print warning with fix message
      if (ieqdf.eq.1) write(iounit,40) ier,nameq,ieq,namupk
      if (ieqdf.eq.0) write(iounit,40) ier,namupk
      go to 30
   10 if (ier.le.32) go to 15
c                                  print warning message
      if (level.lt.3) go to 30
      if (ieqdf.eq.1) write(iounit,45) ier,nameq,ieq,namupk
      if (ieqdf.eq.0) write(iounit,45) ier,namupk
      go to 30
   15 continue
c                                  check for uerset call
      do 20 i=1,6
         if (namupk(i).ne.namset(i)) go to 25
   20 continue
      levold = level
      level = ier
      ier =levold
      if (level.lt.0) level = 4
      if (level.gt.4) level = 4
      go to 30
   25 continue
      if (level.lt.4) go to 30
c                                  print non-defined message
      if (ieqdf.eq.1) write(iounit,50) ier,nameq,ieq,namupk
      if (ieqdf.eq.0) write(iounit,50) ier,namupk
   30 ieqdf = 0
      return
   35 format(19h *** terminal error,10x,7h(ier = ,i3,
     1       20h) from imsl routine ,6a1,a1,6a1)
   40 format(27h *** warning with fix error,2x,7h(ier = ,i3,
     1       20h) from imsl routine ,6a1,a1,6a1)
   45 format(18h *** warning error,11x,7h(ier = ,i3,
     1       20h) from imsl routine ,6a1,a1,6a1)
   50 format(20h *** undefined error,9x,7h(ier = ,i5,
     1       20h) from imsl routine ,6a1,a1,6a1)
c
c                                  save p for p = r case
c                                    p is the page namupk
c                                    r is the routine namupk
   55 ieqdf = 1
      do 60 i=1,6
   60 nameq(i) = namupk(i)
   65 return
      end
