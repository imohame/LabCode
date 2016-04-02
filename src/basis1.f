      subroutine basis1 (r,s,h,p)
c     implicit double precision (a-h,o-z)                                    dp
c
      dimension h(*),p(2,*)
c
      rp=1.0+r
      sp=1.0+s
      rm=1.0-r
      sm=1.0-s
      h(1)=0.25*rp*sp
      h(2)=0.25*rm*sp
      h(3)=0.25*rm*sm
      h(4)=0.25*rp*sm
      p(1,1)=0.25*sp
      p(1,2)=-p(1,1)
      p(1,3)=-0.25*sm
      p(1,4)=-p(1,3)
      p(2,1)=0.25*rp
      p(2,2)=0.25*rm
      p(2,3)=-p(2,2)
      p(2,4)=-p(2,1)
      return
      end
