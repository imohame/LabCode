      subroutine modn(d,w)
c     implicit double precision (a-h,o-z)                                    dp
      common/bk12/ntlen
      common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
      dimension d(*),w(*)
      common/main_block/ a(1)

      call blkcpy (w,d,neq)

      call bsolvr (d,a(ntlen),3,8)

      return
      end
