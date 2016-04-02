      subroutine fissla(ist,b,a,s)
!c     implicit double precision (a-h,o-z)                                    dp
!c
!c     a fast implicit solver for systems of linear equations
!c
!c     steven j. sackett
!c     lawrence livermore laboratory
!c     livermore, california 94550
!c
!c     robert l. taylor
!c     dept. of civil engineering
!c     university of california
!c     berkeley, california 94720
!c
!c     february 11, 1983
!c
!c-----argument usage----------------------------------------------------
      dimension b(1),a(1),s(1)
!c
!c     input arguments
!c           ist       stage flag
!c                         = 1    determine coefficient matrix profile
!c                                and block structure and set memory
!c                                pointers
!c                         = 2    assemble and factorize coefficient
!c                                matrix
!c                         = 3    solve equations using previously
!c                                factorized coefficient matrix
!c                         = 4    zero coefficient matrix and assemble
!c                                contributions from nel elements
!c                         = 5    assemble coefficient matrix contri-
!c                                butions from nel elements
!c                         = 6    factorize assembled coefficient matrix
!c            b        coefficient matrix column heights (excluding
!c                     diagonal), if ist = 1; the right hand side
!c                     (load vector), if ist = 3
!c            a        working space for blocked coefficient matrix
!c            s        array to hold a group of element records during
!c                     equation (coefficient matrix) assembly, if the
!c                     assembly is done by fissle, or containing
!c                     information for a user supplied version of
!c                     subroutine asblk, if assembly is done by the user
!c
!c     if ist = 1, b is assumed to be an integer array; if ist > 1, it is
!c     assumed to be real. if ist = 2, 4, or 5, however, and the equation
!c     assembly is done by fissle, b is not used. the last argument, s,
!c     is required only if ist = 2, 4, or 5. if equation assembly is done
!c     by fissle, this array must be dimensioned for at least ng*ls
!c     words, where ng is the number of elements per group and ls is the
!c     element record length (see common block /fissl1/).
!c
!c     output arguments
!c            b        coefficient matrix diagonal pointers, if ist = 1;
!c                     the right hand side (load vector), if ist = 2, 4,
!c                     or 5 and the user supplies a version of asblk
!c                     which includes assembly of the load vector;
!c                     the solution vector, if ist = 3
!c            a        inverse of d from u'du factorization of the
!c                     coefficient matrix, if ist = 2, 3, or 6
!c
!c     the first neq locations of the working space array, a, may be used
!c     as scratch space between calls to fissle. the remaining nwu-neq
!c     locations may be preserved or destroyed at the option of the user.
!c     if they are preserved, negative values of ist should be used to
!c     inform fissle of this fact and avoid unnecessary disk transfers.
!c     for problems that fit in one block, use of this option can
!c     significantly reduce solution cost.
!c
!c-----common blocks set by the user-------------------------------------
      common/bk12/ntlen
      common/fissn0/neq,mxw,no,n1,nfissl(3)
!c           neq       number of equations
!c           mxw       maximum dimension for the working space
!c           no        unit specifier for error messages
!c           n1        unit specifier for the coefficient matrix file
!c           nfissl(i) maximum number of warning messages to print
!c                     (default is nfissl(i)=20)
!c                       i = 1    zero pivot warnings
!c                       i = 2    equation sign change warnings
!c                       i = 3    loss of accuracy warnings
!c
!c     if mxw = 0, the working space array, a, is assumed to be located
!c     at the end of the users memory space. in this case fissle will
!c     expand memory to yield a working space size which will minimize
!c     the number of blocks required but not exceed the total memory
!c     available. if mxw < 0, fissle proceeds as if mxw = 0 except that
!c     abs(mxw) words of memory are assumed to be already assigned to a.
!c     in all cases the actual number of words of working space used is
!c     stored in the third word of common block /fissl2/.
!c
      common/fissn1/nel,ns,n2,ng,ls
!c           nel       number of elements
!c           ns        maximum number of degrees of freedom per element
!c           n2        unit specifier for file containing element records
!c           ng        number of elements per group (maximum=128)
!c           ls        element record length (minimum=ns+ns*(ns+1)/2)
!c
!c     the file assigned to unit n2 is assumed to contain nel records,
!c     one for each element. the first ns words of each record give
!c     the equation numbers for that element's degrees of freedom (if
!c     an element has less than ns degrees of freedom, zeroes must be
!c     used for the excess equation numbers). the next ns*(ns+1)/2
!c     words contain the upper triangle of the element coefficient
!c     (stiffness) matrix, stored columnwise. for vectorization these
!c     records are processed ng at a time. if ng > 0 they are assumed
!c     to be random-absolute, i.e. written with wrabsf calls. if ng < 0
!c     they are assumed to be sequential-binary, i.e. written with
!c     write(n2) calls (in this case a single call must, except for
!c     possibly a short last record, transmit a logical record
!c     containing ng element records). if ng is not specified by the
!c     user it is set to its maximum value, 128. if ls is not specified
!c     it is set to ns+ns*(ns+1)/2. if n2 = 0, the s array is assumed
!c     to already contain data for nel elements and an element file
!c     is not explicitly required.
!c
!c     this information is used by fissle to assemble the system
!c     coefficient matrix. if assembly is done with a user supplied
!c     version of subroutine asblk, /fissl1/ and the element file
!c     described above are not explicitly required.
!c
!c-----common blocks set by fissle---------------------------------------
      common/fissn2/nwb,nblk,nwu,meb,mxc,mpr,mch,ia1,ia2
!c           nwb       number of words per block
!c           nblk      number of blocks
!c           nwu       total working space used
!c           meb       maximum number of equations in any one block
!c           mxc       maximum column height
!c           mpr       matrix profile (sum of column heights)
!c           mch       mean (average) column height
!c           ia1       offsets to access data in the working space
!c           ia2
!c
      common/fissn3/ifissl,kfissl(3)
!c           ifissl    fissle fatal error flag
!c                         = 0    if there are no fatal errors
!c                         =-1    if a fatal error is detected
!c           kfissl(i) number of warning errors in last call to fissle
!c                       i = 1    zero pivot warnings
!c                       i = 2    equation sign change warnings
!c                       i = 3    loss of accuracy warnings
!c
      common/fissl4/je,jc
c           je        eqn. number of the first eqn. in the current
c                     block -1
c           jc        number of columns in the current equation block
c
      common/fissl5/dt(2)
!c           dt(1)     set negative if the sign of a diagonal term
!c                     changes during factorization
!c           dt(2)     set to zero if more than six significant digits
!c                     are lost in a diagonal term during factorization
!c-----------------------------------------------------------------------
!c
!c     note- this package contains coding for several different computer
!c     systems. statements flagged with "7600" in columns 74-77 apply
!c     to the cdc-7600 only. statements flagged with "cray" in columns
!c     74-77 apply to the cray-1 only. statements flagged with "ansi" in
!c     columns 74-77 are ansi standard fortran equivalents to vectorized
!c     cray coding.
!c
!c-----------------------------------------------------------------------
!!!      common/kkk/km,kkjj
!c
!c     parameter giving the minimum block size
!c     (this should be a multiple of the disk sector size)
      parameter (mbs0=64)
      parameter (mbs1=mbs0-1,mbs2=2*mbs0,mbs3=mbs2-1)
c
      data iai/0/
c
      if(abs(ist).gt.3) go to 201
      if(abs(ist)-2) 101,201,301
c
  101 continue
      nblk   = 1
      ifissl = 0

c.... compute coefficient matrix profile and set diagonal pointers
      mpr    = mapad(b,neq,mxc)
      mch    = mpr/neq

c.... estimate optimal block size
      nwp = mxw
      nwa = nwp
      if(nwp.gt.0) go to 115

      nwa = memav(0) - nwp
      nwp = -nwp
  115 nwu = mpr + neq + 3
      nwa = mbs2*(min(nwu+mbs3,nwa-neq)/mbs2)
      nwb = nwu

      if(nwb.le.nwa) go to 121
      nwb = nwa/2
      nwc = mch+1
      ncb = (nwb-3)/nwc
      mbs = max(ncb*nwc,mxc+1) + 3
      mbs = mbs0*((mbs+mbs1)/mbs0)
      nwb = min(nwb,mbs)
      nwu = 2*nwb

      if(nwb.eq.mbs) go to 121
      ifissl = -1
      go to 901


c.... set memory pointers and expand field length
  121 ia1 = neq
      nwu = nwu + neq
      ia2 = nwu - nwb
      nwr = nwu - nwp
      if(nwr.gt.0) call expndm((ntlen+nwu))

c.... set up block structure

      call mkblk(a(ia1+1),a(1),b(1),ist)
      if(ifissl.lt.0) go to 901
      return
c
  201 continue
      ifissl=0
c.... assemble equations (coefficient matrix) in blocks
      if(abs(ist).eq.6) go to 251
      call asblk(a(ia1+1),a(ia2+1),b(1),s(1),ist)
ck-------------------------------------------------
ck-------------------------------------------------
      ltta=(neq*(neq-1))/2+neq
ck      write(25,202)
  202 format(////5x,'Global stiffness matrix :'/)
      do 204 ii=1,ltta 
c      write(*,203) ii,a(11+ii)
  203 format(5x,'K(',i3,') =',e14.6)
  204 continue 
ck-------------------------------------------------
ck-------------------------------------------------
      if(abs(ist).ne.2) return
c.... factorize the coefficient matrix
  251 kfissl(1)=0
      kfissl(2)=0
      kfissl(3)=0
      call trifac(a(ia1+1),a(ia2+1),a(1),ist)
      return
c
  301 continue
ck---------------------------------------------------------
ck---------------------------------------------------------
!!!      if (km.eq.1) then
!!!ck      write(25,302)
!!!  302 format(//5x,'Load vector :'/) 
!!!      do 307 i=1,neq
!!!c      write(*,304) i,b(i)
!!!  304 format(5x,'F(',i2,') = ',f10.4)
!!!  307 continue
!!!      endif  
ck---------------------------------------------------------
ck--------------------------------------------------------- 
c.... reduce the load vector and backsubstitute
      call forred(a(ia1+1),a(ia2+1),a(1),b(1),ist)
c     if(and(nblk,1).ne.0) go to 311                                    cray1
      if(and(nblk,1).ne.0) go to 311                                    unix
c     if(iand(nblk,1).ne.0) go to 311                                   vms
      call bacsub(a(ia1+1),a(ia2+1),b(1))
      return
  311 call bacsub(a(ia2+1),a(ia1+1),b(1))
!!!!!ck----------------------------------------------
!!!!!ck----------------------------------------------
!!!!!      if (km.eq.1) then
!!!!!ck      write(25,332)
!!!!!  332 format(//5x,'Displacements :'/)
!!!!!      do 334 i=1,neq  
!!!!!ck      write(25,333) i,b(i)
!!!!!  333 format(5x,'U(',i2,') = ',f10.4)
!!!!!  334 continue
!!!!!      endif
!!!!!ck----------------------------------------------
!!!!!ck----------------------------------------------  
      return
c
c.... error in memory allocation
  901 mnw=2*max(nblk,mbs)+neq
      write(no,1001)neq,nblk,nwp,nwa,mnw
      return
 1001 format(5x,'**fatal error** insufficient memory to block equations'
     1/i15,' equations'/i15,' blocks'
     2/i15,' words of working space provided'
     3/i15,' words available'
     4/i15,' words required')
      end
