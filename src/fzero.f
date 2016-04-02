      subroutine fzero(a,ns,ne)
c     implicit double precision (a-h,o-z)                                    dp
      dimension a(*)
      do 10 i=ns,ne
      a(i)=0.
   10 continue
      return
      end
