      subroutine bacsub(ia,ja,b)
c     implicit double precision (a-h,o-z)                                    dp
c
c     backsubstitute equations block by block
c
      dimension ia(4),ja(4),b(*)
      common/fissn0/neq,mxw,no,n1,nfissl(3)
      common/fissn2/nwb,nblk,nwu,meb,mxc,mpr,mch,ia1,ia2
      common/fissl4/je,jc
      common/double/iprec,ncpw,unit
c
      nw=nwb*iprec
c.... backsubstitute from block nblk
      call colbac(ja(4),ja(4),b(je+1))
      if(nblk.eq.1) return
c.... retrieve block nblk-2 from disk
      n1da=(nblk-3)*nw
      call rdabsf(n1,ja,nw,n1da,ioerr)
c.... backsubstitute from block nblk-1
      je=ia(2)
      jc=ia(3)
      call colbac(ia(4),ia(4),b(je+1))
      do 200 n=3,nblk,2
      call riosta (n1)
      if(n.eq.nblk) go to 101
c.... retrieve block nblk-n from disk
      n1da=n1da-nw
      call rdabsf(n1,ia,nw,n1da,ioerr)
c.... backsubstitute from block nblk-n+1
  101 je=ja(2)
      jc=ja(3)
      call colbac(ja(4),ja(4),b(je+1))
      if(n.eq.nblk) go to 200
      call riosta (n1)
      if((n+1).eq.nblk) go to 121
c.... retrieve block nblk-n-1 from disk
      n1da=n1da-nw
      call rdabsf(n1,ja,nw,n1da,ioerr)
c.... backsubstitute from block nblk-n
  121 je=ia(2)
      jc=ia(3)
      call colbac(ia(4),ia(4),b(je+1))
  200 continue
      return
      end
