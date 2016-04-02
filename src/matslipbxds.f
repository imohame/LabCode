       SUBROUTINE  matslipbxds(sig,prop,ln)

c   Bi-Crystal Double-Slip System
C
C-----------------------------------------------------------------------
C
C     DESCRIPTION:
C
C
       use mod_parameters
      REAL phi1,phi2

      common/myvar/a,b,c,e,taur,twomu,xmhigh1,xmhigh2
      common/myvar1/gdotr(2),tau(2)
      common/myvar2/xmhigh3,xmhigh4
      common/myvar3/a01,b1,c1,e1,a02,b2,c2,e2,a03,b3,c3,e3
      common/intgrt/nintg
      common/bk16/maxint,hgc
      common/hokao/lst,nnn2(nume,4)
      common/berg/dels11(nelemg),dels22(nelemg),dels33(nelemg),
     > dels12(nelemg)
      common/vect21/sig11(nelemg),sig22(nelemg),sig33(nelemg),
     > sig12(nelemg)

      common/vect3/
     1 dgi1(nelemg,4),dgi2(nelemg,4),dgi3(nelemg,4),dgi4(nelemg,4),
     > dgi5(nelemg,4),
     2 dgt1(nelemg,4),dgt2(nelemg,4),dgt3(nelemg,4),dgt4(nelemg,4),
     3 f11v(nelemg),f22v(nelemg),f12v(nelemg),f21v(nelemg),
     > dsd5(nelemg),
     4 sig11s(nelemg),sig22s(nelemg),sig33s(nelemg),sig12s(nelemg),
     5 ddp1(nelemg,4),ddp2(nelemg,4),ddp3(nelemg,4),ddp4(nelemg,4),
     > ddp5(nelemg,4)

      common/bk08/kprint,nstep,ite,ilimit,newstf
      common/bk32/nsref,nequit,time,timep,lprint,nprint
      common/meng/
     1 f11dd(nelemg),f22dd(nelemg),f12dd(nelemg),f21dd(nelemg),
     > f33dd(nelemg),
     2 f11dv(nelemg),f22dv(nelemg),f12dv(nelemg),f21dv(nelemg),
     > f33dv(nelemg),
     3 vg11(nelemg),vg22(nelemg),vg12(nelemg),vg21(nelemg),
     > vg33(nelemg),fdet(nelemg),
     4 d1(nelemg),d2(nelemg),d3(nelemg),d4(nelemg),spin(nelemg)
      
      common/history/his1(nume,4),his2(nume,4),his3(nume,4),
     > his4(nume,4)
     1 ,his5(nume,4),his6(nume,4),his7(nume,4),his8(nume,4),
     > his9(nume,4)
     2 ,his10(nume,4),his11(nume,4),his12(nume,4),his13(nume,4)
     3 ,his14(nume,4),his15(nume,4),his16(nume,4),his17(nume,4)
     4 ,his18(nume,4),his19(nume,4)

      common/hist/abc1(nume,4),abc2(nume,4),abc3(nume,4),
     > abc4(nume,4)
     1 ,abc5(nume,4),abc6(nume,4),abc7(nume,4),abc8(nume,4),
     > abc9(nume,4)
     2 ,abc10(nume,4),abc11(nume,4),abc12(nume,4),abc13(nume,4)
     3 ,abc14(nume,4),abc15(nume,4),abc16(nume,4),abc17(nume,4)
     4 ,abc18(nume,4),abc19(nume,4)

      common/custr/sign1(nelemg),sign2(nelemg),sign3(nelemg),
     > sign4(nelemg)
      common/range/mft,mlt,lft,llt,nftm1
      common/bk26/dt,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
      common/vect8/dsave(4,4,nelemg)
      common/cccc/ntt
      common/density/denim(4),denm(4),gdot(4)
      common/poutt/nst1,nst2,nst3,nst4
      common/nslipsys/nnns
      common/alamda/alamda1(nume,4),alamda2(nume,4),dnim(4),dnm(4)

      dimension sig(ln,*),p(2,2),p1(2,2),prop(*)
c
      open(30,file='rotation.f',status='unknown')
      open(31,file='sherslip.f',status='unknown')

c.....resolved shear stress(thisis TWO slip systems)
      open(32,file='restrss1.f',status='unknown')
      open(33,file='restrss2.f',status='unknown')
      open(36,file='restress.f',status='unknown')
      open(37,file='sliprat1.f',status='unknown')
      open(38,file='sliprat2.f',status='unknown')
      open(39,file='tempreture.out',status='unknown')

      open(222,file='interface.out',status='unknown')

!      data pi     /3.141592654/
!      data eta    /68.170/
!      data tempr  /1.00/
!      data xi     /0.5/
c
c
ck---------------------------
ck....timep : total time
ck....dt    : time step size
ck....nstep : no. of step
ck....nintg : intergr. point
ck....ntt   : iteration
ck---------------------------

      nnn1=nstep

ck---------------------------------
ck....define total current stresses
ck---------------------------------
      DO 5 i=mft,mlt
        sign1(i) = sig(1,i)
        sign2(i) = sig(2,i)
        sign3(i) = sig(3,i)
        sign4(i) = sig(4,i)
    5 CONTINUE

ck---------------------
ck....define properties
ck---------------------
      tauy    = prop(3)
      gdot01  = prop(5)
      gdot02  = prop(5)
      gdot03  = prop(5)
      gdot04  = prop(5)
      xmhigh1 = prop(4)
      xmhigh2 = prop(4)
      xmhigh3 = prop(4)
      xmhigh4 = prop(4)
      gdotcr  = prop(6)
      ielem   = prop(8)
      nst1    = prop(7)
      nst2    = 2*nst1
      nst3    = 3*nst1
      nst4    = 4*nst1


      IF (nstep.EQ.nst1.OR.nstep.EQ.nst2.OR.nstep.EQ.nst3
     $    .OR.nstep.EQ.nst4) THEN

c          WRITE(30,*) '------------> nstep = ',nstep
c          WRITE(31,*) '------------> nstep = ',nstep
c          WRITE(32,*) '------------> nstep = ',nstep
c          WRITE(33,*) '------------> nstep = ',nstep
c          WRITE(36,*) '------------> nstep = ',nstep
c          WRITE(37,*) '------------> nstep = ',nstep
c          WRITE(38,*) '------------> nstep = ',nstep
c          WRITE(39,*) '------------> nstep = ',nstep

      ENDIF        

ck
      DO 10 i=mft,mlt 

        ink=i+nftm1


ck---------------------------------
ck....define deviatoric deformation
ck---------------------------------
        traced=d1(i)+d2(i)+d3(i)
        dd1=d1(i)-traced/3.00
        dd2=d2(i)-traced/3.00
        dd3=d3(i)-traced/3.00
        dd4=d4(i)

ck------------------------------
ck....define deviatoric stresses
ck------------------------------
        press=(sign1(i)+sign2(i)+sign3(I))/3.00
        sdev1=sign1(i)-press
        sdev2=sign2(i)-press
        sdev3=sign3(i)-press
        sdev4=sign4(i)

ck-------------------------------------------------------
ck....update internal variables for multiple slip systems
ck-------------------------------------------------------
        IF(nnn1.eq.nnn2(ink,nintg)) THEN

           chi    = abc1(ink,nintg)*pi/180.0
           gamma  = abc2(ink,nintg)
           tau(1) = abc3(ink,nintg)
           tau(2) = abc4(ink,nintg)
           taur   = abc7(ink,nintg)
           gdot(1)= abc8(ink,nintg)
           gdot(2)= abc9(ink,nintg)
           temp   = (abc10(ink,nintg)+273.0)/293.0
           temp   = max(temp,tempr)

        ELSE

           chi    = abc1(ink,nintg)*pi/180.0
           gamma  = abc2(ink,nintg)
           tau(1) = abc3(ink,nintg)
           tau(2) = abc4(ink,nintg)
           taur   = abc7(ink,nintg)
           gdot(1)= abc8(ink,nintg)
           gdot(2)= abc9(ink,nintg)
           temp   = (abc10(ink,nintg)+273.0)/293.0
           temp   = max(temp,tempr)


ck  
           abc1(ink,nintg)  = his1(ink,nintg)
           abc2(ink,nintg)  = his2(ink,nintg)
           abc3(ink,nintg)  = his3(ink,nintg)
           abc4(ink,nintg)  = his4(ink,nintg)
           abc7(ink,nintg)  = his7(ink,nintg)
           abc8(ink,nintg)  = his8(ink,nintg)
           abc9(ink,nintg)  = his9(ink,nintg)
           abc10(ink,nintg) = his10(ink,nintg)

           abc11(ink,nintg) = abc10(ink,nintg)/20.


ck---------------------------------------------------
ck....printout internal variables for current contour
ck---------------------------------------------------
       IF (nstep.eq.nst1.OR.nstep.eq.nst2.OR.nstep.eq.nst3
     $.OR.nstep.eq.nst4) THEN

c//---
c// Only to print out the angles to the file <interface.out>
c//
           IF(ink.EQ.2) THEN
              write(222,10020) (180.*phi1/pi),(180.*phi2/pi)
           ELSEIF(ink.EQ.(ielem+2)) THEN
              write(222,10021) (180.*phi1/pi),(180.*phi2/pi)
           ENDIF


           WRITE(30,1001) abc1(ink,nintg)
           WRITE(31,1001) abc2(ink,nintg)
           WRITE(32,1001) abc3(ink,nintg)
           WRITE(33,1001) abc4(ink,nintg)
           WRITE(36,1001) abc7(ink,nintg)
           WRITE(37,1001) abc8(ink,nintg)
           WRITE(38,1001) abc9(ink,nintg)
           WRITE(39,1001) abc11(ink,nintg)


           IF((ink.GE.71).AND.(ink.LE.80)) THEN
             WRITE(222,*) '    '
             WRITE(222,*) '==> Element No. :',ink
             WRITE(222,*) '-------------------------'
             WRITE(222,10011) abc1(ink,nintg)
             WRITE(222,10012) abc2(ink,nintg)
             WRITE(222,10013) abc3(ink,nintg)
             WRITE(222,10014) abc4(ink,nintg)
             WRITE(222,10015) abc7(ink,nintg)
             WRITE(222,10016) abc8(ink,nintg)
             WRITE(222,10017) abc9(ink,nintg)
             WRITE(222,10018) abc11(ink,nintg)
             WRITE(222,*) '    '
             WRITE(222,*) '    '
          ENDIF

 1001    FORMAT(5x,e15.8)
10011    FORMAT('Rotation ............................ = ',e15.8)
10012    FORMAT('Accumulated plastic shear strain .....= ',e15.8)
10013    FORMAT('Resolved shear stress for slip (1) .. = ',e15.8)
10014    FORMAT('Resolved shear stress for slip (2) .. = ',e15.8)
10015    FORMAT('Reference shear stress .............. = ',e15.8)
10016    FORMAT('Shear strain (slip) rate for slip (1) = ',e15.8)
10017    FORMAT('Shear strain (slip) rate for slip (2) = ',e15.8)
10018    FORMAT('Current temp. of solid (T) / Tr ..... = ',e15.8)

10020    FORMAT('Crystal (A): phi(1) = ',f6.2,' ; phi(2) = ',f6.2)
10021    FORMAT('Crystal (B): phi(1) = ',f6.2,' ; phi(2) = ',f6.2)


       ENDIF

      ENDIF

ck-------------------------------------------
ck....define initial angles from loading axis
ck-------------------------------------------
ck....single crystal case :
ck      if(ink.eq.1) then
ck....bicrystal case :

      IF(ink.le.ielem) THEN
         phi1   = (pi/180.0)*35.0 
         phi2   = (pi/180.0)*35.0

         ym     = prop(1)
         g      = prop(1)/(2.0*(1.0+prop(2)))      
         twomu  = g*2.0
         proppp = prop(1)/(3.*(1.-2.*prop(2)))
         alamdt = proppp*dt
      
      ELSE

         phi1   = (pi/180.0)*35.0
         phi2   = (pi/180.0)*35.0

         ym     = prop(49)
         g      = prop(49)/(2.0*(1.0+prop(2)))
         twomu  = g*2.0
         proppp = prop(49)/(3.*(1.-2.*prop(2)))
         alamdt = proppp*dt

      ENDIF


ck---------------------------------------
ck....compute tangent constitutive matrix
ck---------------------------------------
ck....single crystal case :
ck      if(ink.eq.1) then
ck....bicrystal case :

      IF(ink.le.ielem) THEN
         dsave(1,1,i) = prop(10)
         dsave(1,2,i) = prop(11)
         dsave(1,3,i) = prop(11)
         dsave(1,4,i) = 0.
         dsave(2,1,i) = prop(11)
         dsave(2,2,i) = prop(10)
         dsave(2,3,1) = prop(11)
         dsave(2,4,i) = 0.
         dsave(3,1,i) = prop(11)
         dsave(3,2,i) = prop(11)
         dsave(3,3,i) = prop(10)
         dsave(3,4,i) = 0.
         dsave(4,1,i) = 0.
         dsave(4,2,i) = 0.
         dsave(4,3,i) = 0.
         dsave(4,4,i) = prop(25)
      ELSE
         dsave(1,1,i) = prop(58)
         dsave(1,2,i) = prop(59)
         dsave(1,3,i) = prop(59)
         dsave(1,4,i) = 0.
         dsave(2,1,i) = prop(59)
         dsave(2,2,i) = prop(58)
         dsave(2,3,1) = prop(59)
         dsave(2,4,i) = 0.
         dsave(3,1,i) = prop(59)
         dsave(3,2,i) = prop(59)
         dsave(3,3,i) = prop(58)
         dsave(3,4,i) = 0.
         dsave(4,1,i) = 0.
         dsave(4,2,i) = 0.
         dsave(4,3,i) = 0.
         dsave(4,4,i) = prop(73)
      ENDIF

ck--------------------------
ck....check if gdot > gdotcr
ck--------------------------
      IF(gdot(1).ge.gdotcr) THEN
         xmhigh1 = 1.00
         gdot01  = gdotcr
      ENDIF

      IF(gdot(2).ge.gdotcr) THEN
         xmhigh2 = 1.00
         gdot02  = gdotcr
      ENDIF

ck--------------------------------------
ck....calculate Pij for each slip system
ck--------------------------------------
ck....single crystal case :
ck      if(ink.eq.1) then
ck....bicrystal case :

      IF(ink.le.ielem) THEN
         angl1  = (phi1-chi)
         angl2  = (phi2+chi)

         c2t  = cos(2.0*angl1)
         s2t  = sin(2.0*angl1)
         c2t1 = cos(2.0*angl2)
         s2t1 = sin(2.0*angl2)

ck....slip system 1 :
         p(1,1)  = -s2t/2.0
         p(1,2)  = -c2t/2.0
         p(2,2)  = s2t/2.0
         wp      = 0.5

ck....slip system 2 :
         p1(1,1) = -s2t1/2.0
         p1(1,2) = c2t1/2.0
         p1(2,2) = s2t1/2.0
         wp1     = -0.5

      ELSE

         angl1  = (phi1-chi)
         angl2  = (phi2+chi)

         c2t  = cos(2.0*angl1)
         s2t  = sin(2.0*angl1)
         c2t1 = cos(2.0*angl2)
         s2t1 = sin(2.0*angl2)

ck....slip system 1 :
         p(1,1) = -s2t/2.0
         p(1,2) = -c2t/2.0
         p(2,2) =  s2t/2.0
         wp     =  0.5

ck....slip system 2 :
        p1(1,1) = -s2t1/2.0
        p1(1,2) =  c2t1/2.0
        p1(2,2) =  s2t1/2.0
        wp1     = -0.5

      ENDIF

ck----------------------------------------
ck....perform multiplication of DDK*P(I,J)
ck----------------------------------------

      gdotr(1) = dd1*p(1,1)+dd2*p(2,2)+2.*dd4*p(1,2)
      gdotr(2) = dd1*p1(1,1)+dd2*p1(2,2)+2.*dd4*p1(1,2)

ck-------------------------
ck....define hardening laws
ck-------------------------
C:------>>
c.... corection factor for the tempr
       corfctr = (tempr/temp)**xi

ck....power laws :
       taur=tauy*((100.0*gamma+1.0)**0.1)


ck....hyberbolic hardening laws :
c       taur=tauy+(0.8*tauy*tanh(11.125*gamma))

ck....hyberbolic hardening (including thermal effects) :
ck      x=(tempr/temp)**v
ck      taur=tauy+(0.8*tauy*tanh(11.125*gamma))*x

ck....hardening from dislocation evolution :
C      taur=0.976923e-03+0.5*g*bv*(denim(1)**0.5+denim(2)**0.5
C     $+denim(3)**0.5+denim(4)**0.5)


ck------------------------------------------------------------
ck....define coefficients for nonlinear differential equations
ck....(obtain tau(i), i=1,12 for f.c.c.)
ck------------------------------------------------------------

      a = gdot01*(p(1,1)*p(1,1)+p(2,2)*p(2,2)+
     #    2.0*p(1,2)*p(1,2))
      b = gdot02*(p(1,1)*p1(1,1)+p(2,2)*p1(2,2)
     #    +2.0*p(1,2)*p1(1,2))
      c = gdot01*(p1(1,1)*p(1,1)+p1(2,2)*p(2,2)
     #    +2.0*p1(1,2)*p(1,2))
      e = gdot02*(p1(1,1)*p1(1,1)
     #    +p1(2,2)*p1(2,2)+2.0*p1(1,2)*p1(1,2))


ck--------------------------------------
ck
      IF (nstep.ne.0.and. ntt.eq.1) THEN
         timexx = timep-dt
      ELSE
         timexx = timep
      ENDIF

ck--------------------------------------

      CALL m2dsolve(dt,timexx)

ck--------------------------------------------------------
ck....calculate shear strain (shear strain should have the
ck....same sign as the resolved shear stress)
ck--------------------------------------------------------
      sat    = 1.0/1.00
      gdot(1)=gdot01*((abs(tau(1)/taur)*sat)**(xmhigh1-1.))
     1        *(tau(1)/taur)*sat
      gdot(2)=gdot02*((abs(tau(2)/taur)*sat)**(xmhigh2-1.))
     1        *(tau(2)/taur)*sat


c      write(7777,*) '- gdot(1): ', gdot(1)
c      write(7777,*) '- gdot(2): ', gdot(2)


ck-----------------------------------
ck....define plastic deformation rate
ck-----------------------------------
      dep1 = p(1,1)*gdot(1)+p1(1,1)*gdot(2)
      dep2 = p(2,2)*gdot(1)+p1(2,2)*gdot(2)
      dep4 = p(1,2)*gdot(1)+p1(1,2)*gdot(2)

ck-----------------------------
ck....calculate OMEGAP & OMEGAE
ck-----------------------------
      omegap = gdot(1)*wp+gdot(2)*wp1
      omegae = (spin(i)+omegap)*dt
ck
      o12    = dt*omegap*(sdev1-sdev2)
      o11    = -2.00*dt*omegap*(sdev4)
      o22    = -o11

ck------------------------------
ck....update deviatoric stresses
ck------------------------------
      ssdev4 = sdev4+twomu*dt*(dd4-dep4)+o12
      ssdev1 = sdev1+twomu*dt*(dd1-dep1)+o11
      ssdev2 = sdev2+twomu*dt*(dd2-dep2)+o22
      ssdev3 = sdev3+twomu*dt*dd3

      ss=abs(ssdev1*dep1)+abs(ssdev2*dep2)+abs(2.00*dep4*ssdev4)
      temp   = temp+dt*eta*ss

ck-------------------------------------------------------------------
ck....update pressure and add deviatoric part to get current stresses
ck-------------------------------------------------------------------
      press    = press+alamdt*traced
ck
      sign1(i) = ssdev1+press
      sign2(i) = ssdev2+press
      sign3(i) = ssdev3+press
      sign4(i) = ssdev4
c
      sig(1,i) = sign1(i)
      sig(2,i) = sign2(i)
      sig(3,i) = sign3(i)
      sig(4,i) = sign4(i)

ck
      chi             = chi+omegae
      his1(ink,nintg) = chi*(180.0/pi)
      gamma=dt*(abs(gdot(1))+abs(gdot(2)))+gamma
      his2(ink,nintg) = gamma
      his3(ink,nintg) = tau(1)
      his4(ink,nintg) = tau(2)
      his7(ink,nintg) = taur
      his8(ink,nintg) = gdot(1)
      his9(ink,nintg) = gdot(2)
      his10(ink,nintg)=(temp*293.0)-273.0

      nnn2(ink,nintg) = nnn1

  10  CONTINUE

c------------------------------
c.....obtain stresses increment
c------------------------------
      DO 440 i=mft,mlt
         dels11(i) = sign1(i)-sig(1,i)
         dels22(i) = sign2(i)-sig(2,i)
         dels33(i) = sign3(i)-sig(3,i)
         dels12(i) = sign4(i)-sig(4,i)

c----------------------------
c.....update rotated stresses 
c----------------------------
        sig11(i)   = sig(1,i)
        sig22(i)   = sig(2,i)
        sig33(i)   = sig(3,i)
        sig12(i)   = sig(4,i)
ck
        sig11s(i)  = sig(1,i)
        sig22s(i)  = sig(2,i)
        sig33s(i)  = sig(3,i)
        sig12s(i)  = sig(4,i)
  440 CONTINUE


      RETURN
      END
