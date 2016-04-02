      subroutine ugetio(iopt,nin,nout)
c                                  specifications for arguments
      integer            iopt,nin,nout
c                                  specifications for local variables
      integer            nind,noutd
      data               nind/1/,noutd/2/
c                                  first executable statement
      if (iopt.eq.3) go to 10
      if (iopt.eq.2) go to 5
      if (iopt.ne.1) go to 9005
      nin = nind
      nout = noutd
      go to 9005
    5 nind = nin
      go to 9005
   10 noutd = nout
 9005 return
      end
