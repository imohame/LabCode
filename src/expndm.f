      subroutine expndm(n)
c     implicit double precision (a-h,o-z)                                    dp
c
      common/array/maxa,maxadd,ifield
      common /main_block/ b(1)
c
c     expand large core memory
c
      write(*,*) ifield, n
      lmin=min(ifield,n)
c     if (n-ifield.ne.0)call memory(b(lmin),n-ifield)                   cray1
c     if (n-ifield.ne.0)call memory(b(lmin),n-ifield)                   vms
      if (n-ifield.ne.0)call memory(b(lmin),n-ifield)                   wkstn
c     if (n-ifield.ne.0)call mmemry(b(lmin),n-ifield)                   unics
c     if (n-ifield.ne.0)call memory(b(lmin),n-ifield)                   nersc
      ifield=n
c
      return
      end
