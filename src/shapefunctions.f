      subroutine shapefunctions(psi,eta,N,ShpFcnDeriv,maxint)
      integer maxint
      real*8 eta(maxint),psi(maxint)
      real*8 N(maxint,4), ShpFcnDeriv(maxint,2,4)

      do i=1,maxint
      N(i,1)=1./4.*(1.-psi(i))*(1.-eta(i))
      N(i,2)=1./4.*(1.+psi(i))*(1.-eta(i))
      N(i,3)=1./4.*(1.+psi(i))*(1.+eta(i))
      N(i,4)=1./4.*(1.-psi(i))*(1.+eta(i))
      ShpFcnDeriv(i,1,1)=-(1-eta(i))/4.
      ShpFcnDeriv(i,1,2)=(1-eta(i))/4.
      ShpFcnDeriv(i,1,3)=(1+eta(i))/4.
      ShpFcnDeriv(i,1,4)=-(1+eta(i))/4.
      ShpFcnDeriv(i,2,1)=-(1-psi(i))/4.
      ShpFcnDeriv(i,2,2)=-(1+psi(i))/4.
      ShpFcnDeriv(i,2,3)=(1+psi(i))/4.
      ShpFcnDeriv(i,2,4)=(1-psi(i))/4.
      enddo

      return
      end