      subroutine m2dsolve(dt,time)

      common /myvar/ a,b,c,e,taur,twomu,xmhigh1,xmhigh2
      common /myvar1/ gdotr(2),tau(2)
      dimension ystart(2)
      external RKQC
      external derivs

      nvar      = 2
      ystart(1) = tau(1)
      ystart(2) = tau(2)
      x1        = time
      x2        = time+dt
      eps       = 1.00e-04
      h1        = dt
      hmin      = h1/1.0e+3

      call M2ODEINT(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs
     >              ,RKQC)

      tau(1)    = ystart(1)
      tau(2)    = ystart(2)

      return
      end
