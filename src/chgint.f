      subroutine chgint
c     implicit double precision (a-h,o-z)                                    dp
      common/bk11/delt,alfa,iequit,iprint,isref
      common/bk26/dt,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
      a0=1./(alfa*dt*dt)
      a1=(delt/alfa)/dt
      a2=1./(alfa*dt)
      a3=.5/alfa-1.
      a4=(delt/alfa)-1.
      a5=dt*(.5*(delt/alfa)-1.)
      a6=a0
      a7=-a2
      a8=-a3
      a9=dt*(1.-delt)
      a10=delt*dt
      return
      end
