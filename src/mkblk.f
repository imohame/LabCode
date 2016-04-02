      subroutine mkblk(ia,ipt,idiag,ist)
c     implicit double precision (a-h,o-z)                                    dp
c
c     set up block structure for equations
c
c     each block begins with three control words- (1) the number of
c     the first block required to factorize the block, (2) the global
c     equation number of the last equation in the previous block,
c     and (3) the number of columns in the block. this is followed
c     by pointers to the diagonal terms of the coefficient matrix
c     columns in the block and, finally, the actual column data.
c
      dimension ia(4),ipt(*),idiag(*)
      common/fissn0/neq,mxw,no,n1,nfissl(3)
      common/fissn2/nwb,nblk,nwu,meb,mxc,mpr,mch,ia1,ia2
      common/fissn3/ifissl,kfissl(3)
      common/double/iprec,ncpw,unit

c
      if(ia1.ne.ia2) go to 11
c
c.... establish pointers for single block problem
      ia(1)=1
      ia(2)=0
      ia(3)=neq
      meb=neq
      do 4 i=1,neq
      ia(i+3)=idiag(i)+neq
    4 continue

      if(ist.lt.0) return

c.... write block information to file n1
      n1da=0
      nw=nwb*iprec

      call wrabsf(n1,ia,nw,n1da,ioerr)
!!!!!!!!!!      call riosta (n1)
      return
c
c.... establish pointers for multi-block problem
   11 meb=0
      ie=0
      ic=0
      la=0
      ir=0
      nw=nwb*iprec
      ipt(1)=0
      n1da=-nw
      do 200 n=1,neq
      id=idiag(n)
      if((id-la+ic+3).lt.nwb) go to 101
      if(ic.ne.0) go to 21
      write(no,1001)n
   15 ifissl=-1
      return
   21 ipt(nblk)=ipt(nblk)+ic
      do 50 i=1,nblk
      ib=i
      if(ie.lt.ipt(ib)) go to 61
   50 continue
   61 ia(1)=ib
      ia(2)=ipt(nblk)-ic
      ia(3)=ic
      meb=max(meb,ic)
      do 70 k=1,ic
      ia(k+3)=ia(k+3)+ic
   70 continue

c.... write block information to file n1
      n1da=n1da+nw

      call wrabsf(n1,ia,nw,n1da,ioerr)
!!!!!!!!!      call riosta (n1)

c.... reestablish pointers
      nblk=nblk+1
      if(nblk.gt.nwb) go to 15
      ic=0
      la=ir
      ie=n-1
      ipt(nblk)=ie
  101 ic=ic+1
      ia(ic+3)=id-la
      ie=min(ie,n-id+ir)
      ir=id
  200 continue
c
      ipt(nblk)=ipt(nblk)+ic
      do 250 i=1,nblk
      ib=i
      if(ie.lt.ipt(ib)) go to 261
  250 continue
  261 ia(1)=ib
      ia(2)=ipt(nblk)-ic
      ia(3)=ic
      meb=max(meb,ic)
      do 270 k=1,ic
      ia(k+3)=ia(k+3)+ic
  270 continue
c.... write last block to file n1
      n1da=n1da+nw

      call wrabsf(n1,ia,nw,n1da,ioerr)
!!!!!!!!!!      call riosta (n1)

      return

 1001 format(5x,'**fatal error** insufficient memory to store column '
     1,i6)
      end
