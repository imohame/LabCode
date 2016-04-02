      subroutine forred(ia,ja,d,b,ist)
c     implicit double precision (a-h,o-z)                                    dp
c
c     reduce load vector block by block
c
      dimension ia(4),ja(4),d(*),b(*)
      common/fissn0/neq,mxw,no,n1,nfissl(3)
      common/fissn2/nwb,nblk,nwu,meb,mxc,mpr,mch,ia1,ia2
      common/fissl4/je,jc
      common/double/iprec,ncpw,unit
c
c     note- this routine has been vectorized for the cray-1. statements
c     flagged with "cray" in columns 74-77 apply to the cray only.
c     statements flagged with "ansi" in columns 74-77 are ansi standard
c     fortran equivalent to the cray coding.
c
      n1da=0
      nw=nwb*iprec
      n=1
      if((n.eq.nblk).and.(ist.lt.0)) go to 101
c.... retrieve block 1 from disk
      call rdabsf(n1,ia,nw,n1da,ioerr)
   20 continue
      call riosta (n1)
      if(n.eq.nblk) go to 101
c.... retrieve block n+1 from disk
      n1da=n1da+nw
      call rdabsf(n1,ja,nw,n1da,ioerr)
c.... reduce b using block n
  101 je=ia(2)
      jc=ia(3)
      call recol(-1,1,jc,ia(4),ia(4),b(je+1),d(je+1))
      if(n.eq.nblk) go to 201
      call riosta (n1)
      if((n+1).eq.nblk) go to 121
c.... retrieve block n+2 from disk
      n1da=n1da+nw
      call rdabsf(n1,ia,nw,n1da,ioerr)
c.... reduce b using block n+1
  121 je=ja(2)
      jc=ja(3)
      call recol(-1,1,jc,ja(4),ja(4),b(je+1),d(je+1))
      n=n+2
      if(n.le.nblk) go to 20
c.... divide reduced load vector by diagonal (inverses are stored)
  201 continue
      do 300 j=1,neq
c     b(j)=cvmgz(b(j),d(j)*b(j),d(j))                                   cray1
      if(d(j).ne.0)b(j)=d(j)*b(j)                                       vax75
  300 continue
      return
      end
