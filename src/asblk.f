      subroutine asblk(ia,ja,b,s,ist)

       use mod_parameters

c     implicit double precision (a-h,o-z)                                    dp
c
c     assemble equations (coefficient matrix) in blocks
c
c     this subroutine does a finite-element type equation assembly. the
c     coefficient matrix for the system is formed by adding together
c     contributions from individual element matrices stored in compact
c     form. the equation numbers given for each element's degrees of
c     freedom map the element matrices from the compact storage form
c     to the appropriate place in the system coefficient matrix (see
c     comments in subroutine addem for further details).
c
c     a user supplied version of this subroutine is required if direct
c     control over equation assembly is desired.
c
      dimension ia(4),ja(4),b(*),s(*)
c
c     input arguments
c           ia        working space for a block of equations
c           ja        working space for a second block of equations
c            b        load vector (not used in this version)
c            s        array to hold a group of element records
c
c     ia and ja are both nwb words in length. space for two blocks of
c     equations is provided because it is available and the assembly of
c     more than one block at a time reduces the number of element file
c     searches required. the load vector, b, is included in the argument
c     list for flexibility in the replacement of this version of asblk
c     by a user supplied version. this allows for the possibility of
c     load vector assembly at equation assembly time.
c
c     each block begins with three control words- (1) the number of
c     the first block required to factorize the block, (2) the global
c     equation number of the last equation in the previous block,
c     and (3) the number of columns in the block. this is followed
c     by pointers to the diagonal terms of the coefficient matrix
c     columns in the block and, finally, the actual column data.
c
      common/fissn0/neq,mxw,no,n1,nfissl(3)
      common/fissn1/nel,ns,n2,ng,ls
      common/fissn2/nwb,nblk,nwu,meb,mxc,mpr,mch,ia1,ia2
      common/fissn3/ifissl,kfissl(3)
      common/double/iprec,ncpw,unit
c=======================================================================
c Addition by KHALiL - 04/11/07
c      common/bk48/stress(4),ft,p1,p2,ag,d(4,4),ipt,NEL1,nstate
c=======================================================================

c
      if(ng.eq.0) ng=nelemg
      if(ls.eq.0) ls=ns+ns*(ns+1)/2
      nepg=abs(ng)
      if((nepg.le.nelemg).and.(ls.ge.(ns+ns*(ns+1)/2))) go to 11
      ifissl=-1
      write(no,1001)nelemg
      return
c
   11 nw=2*nwb*iprec
      n1da=0
      n=1
      if((n.eq.nblk).and.(ist.lt.0)) go to 21
c.... retreive two blocks of equations from disk
   20 continue
      if(n.eq.nblk) nw=nwb*iprec
      call rdabsf(n1,ia,nw,n1da,ioerr)
      call riosta (n1)
c.... initialize equation counters and zero coefficients for block n
   21 je=ia(2)
      ke=ia(3)+je
      le=ke
      if(abs(ist).eq.5) go to 31
      call fzero(ia(4),ia(3)+1,nwb-3)
   31 if(n.eq.nblk) go to 51
c.... initialize equation counters and zero coefficients for block n+1
      le=ja(3)+ke
      if(abs(ist).eq.5) go to 51
      call fzero(ja(4),ja(3)+1,nwb-3)
c
c.... loop through element file
   51 if(ng.lt.0) rewind n2
      keda=0
      ner=nel
	  	  
      do 100 l=1,nel,nepg
c.... read in coefficient matrix data for a group of elements
      neg=min(nepg,ner)
      lr=neg*ls*iprec
      ner=ner-neg
      if(ng.ge.0) go to 55
      read(n2)(s(i),i=1,lr)
      go to 61
   55 if(n2.eq.0) go to 61
      call rdabsf(n2,s,lr,keda,ioerr)
      call riosta (n2)
      keda=keda+lr
c.... add contributions to equations in blocks n and n+1
   61 call addem(je,ke,le,ia(4),ia(4),s(1),s(ns+1),ns,ls,
     1 ls*iprec,neg)
c
  100 continue
c
c.... write updated blocks to disk
      if((nblk.eq.1).and.(ist.lt.0)) go to 201
      call wrabsf(n1,ia,nw,n1da)
      call riosta (n1)
      n1da=n1da+nw
      n=n+2
      if(n.le.nblk) go to 20
c
  201 continue
      return
 1001 format(5x,'**fatal error** ng > %i or ls < ns+ns*(ns+1)/2')
      end
