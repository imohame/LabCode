      subroutine updat (u,udt,udtt,ui)
!!!!  call       updat (a(k18),a(k55),a(k56),a(k19))    
c     implicit double precision (a-h,o-z)                                    dp
c
c.... update displacement, velocity, and acceleraton (newmark-beta)
c
      common/bk02/ioofc,iphase,imass,lpar(9)
      common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
      common/bk26/dt,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
      dimension u(*),udt(*),udtt(*),ui(*)
c
c.... update total displacement u(*) with step displacemnt ui(*)
      do 10 i=1,neq
   10 u(i)=u(i)+ui(i)
c
c.... update velocity ut(*) and acceleration utt(*) if dynamics
      if(imass.eq.1)then
      do 20 i=1,neq
      udttld=udtt(i)
      udtt(i)=a6*ui(i)+a7*udt(i)+a8*udttld
   20 udt(i)=udt(i)+a9*udttld+a10*udtt(i)
c
c.... change to dynamic problem if imass=2
      elseif(imass.eq.2)then
      imass=1
      endif
c
      return
      end
