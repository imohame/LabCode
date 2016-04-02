      subroutine trifac(ia,ja,d,ist)
c     implicit double precision (a-h,o-z)                                    dp
c
c     factorize coefficient matrix using block by block column reduction
c
      dimension ia(4),ja(4),d(*)
      common/fissn0/neq,mxw,no,n1,nfissl(3)
      common/fissn2/nwb,nblk,nwu,meb,mxc,mpr,mch,ia1,ia2
      common/fissl4/je,jc
      common/double/iprec,ncpw,unit
c
      je=0
      jc=neq
      n1da=0
      nw=nwb*iprec
      n=1
      if(n.eq.nblk) go to 101
c.... retrieve block n from disk
   20 continue
      call rdabsf(n1,ja,nw,n1da,ioerr)
      call riosta (n1)
      jb=ja(1)
      je=ja(2)
      jc=ja(3)
      nm=n-jb
      if(nm.eq.0) go to 101
      n2da=(jb-2)*nw
      do 100 m=1,nm
c.... retrieve block n-jb+m from disk
      n2da=n2da+nw
      call rdabsf(n1,ia,nw,n2da,ioerr)
      call riosta (n1)
      ie=ia(2)
      ic=ia(3)
c.... reduce block n using previously reduced block n-jb+m
      call actred(ia(4),ia(4),ja(4),ja(4),ie,ic,d(je+1))
  100 continue
c.... reduce block n using block n
  101 call actred(ja(4),ja(4),ja(4),ja(4),je,jc,d(je+1))
c.... write reduced block onto disk
      if((nblk.eq.1).and.(ist.lt.0)) return
      call wrabsf(n1,ja,nw,n1da)
      call riosta (n1)
      n1da=n1da+nw
      n=n+1
      if(n.le.nblk) go to 20
      return
      end




