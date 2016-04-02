      subroutine rdiska (lus,w,nw,da)                                   vax750
c
c     interface to direct access i/o for rwabsf random i/o routines
c
c.... entry to transfer one record from disk to a buffer
c
c     input arguments
c            lus      the file unit specifier (logical unit no.)
c            nw       number of words to read from disk
c                     (must be a multiple of 512)
c            da       zero base disk word.address
c                     (must be on a sector boundary)
c
c     output arguments
c            w        data read from disk
c
      implicit integer(a-z)                                             vax750
      dimension w(nw)                                                   vax750
c
      common/frfcm1/mxfrf,ifrf,buflen,fcsize,dskloc ,curlen,kop,ier     vax750
c
      lda=da/buflen+1                                                   vax750
      read (lus,rec=lda,iostat=ios) w                                   vax750
      return                                                            vax750
c
c.... entry to transfer one record from a buffer to disk
c
c     input arguments
c            lus      the file unit specifier (logical unit no.)
c            w        data to be written to disk
c            nw       number of words to write to disk
c                     (must be a multiple of 512)
c            da       zero base disk word.address
c                     (must be on a sector boundary)
c
      entry wdiska(lus,w,nw,da)                                         vax750
      lda=da/buflen+1                                                   vax750
      write (lus,rec=lda) w                                             vax750
      return                                                            vax750
c
      end                                                               vax750
