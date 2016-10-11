      subroutine stfnf (yz,matp,hgs,nmel,prop,den)
!c.....this subroutine calculations for B, Bbar, K, Dij, Wij, Fij
!c     implicit double precision (a-h,o-z)                                    dp
      use mod_parameters
      use EC_Objects_manager
      
      common/bk02/ioofc,iphase,imass,lpar(9)
      common/bk03/numdc,imassn,idampn,irller,penstf
      common/bk08/kprint,nstep,ite,ilimit,newstf
      common/bk13/xnorm0(6),xnormc(6)!,xnk2d(20)
      common/bk16/maxint,hgc
      common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
      common/bk34/s(44,1)
      common/intgp/c(16),ipt,nel,nelsub
      common/range/mft,mlt,lft,llt,nftm1
      common/intgrt/nintg
      common/riksw1/sp1,ds0,cs01,cs02,rrsp,linsch,igso,irco,idamp
      common/xcom0/imeth,interq,imess,istart,igraf
      common/vect0/
     1 q11(nelemg),q22(nelemg),q12(nelemg),q21(nelemg),
     2 r11(nelemg),r22(nelemg),r12(nelemg),r21(nelemg),
     3 s11(nelemg),s22(nelemg),s12(nelemg),s21(nelemg)
      common/vect1/r1(nelemg),r2(nelemg),r3(nelemg),r4(nelemg),
     > r5(nelemg),r6(nelemg),
     1 r7(nelemg),r8(nelemg),mtype(nelemg),mte(nelemg)
      common/vect2/
     1 ed11(nelemg),ed21(nelemg),ed12(nelemg),ed22(nelemg),
     2 ed13(nelemg),ed23(nelemg),ed14(nelemg),ed24(nelemg),
     3 fd11(nelemg),fd21(nelemg),fd12(nelemg),fd22(nelemg),
     4 fd13(nelemg),fd23(nelemg),fd14(nelemg),fd24(nelemg)
      common/vect3/
     1 dgi1(nelemg,4),dgi2(nelemg,4),dgi3(nelemg,4),dgi4(nelemg,4),
     > dgi5(nelemg,4),
     2 dgt1(nelemg,4),dgt2(nelemg,4),dgt3(nelemg,4),dgt4(nelemg,4),
     3 f11v(nelemg),f22v(nelemg),f12v(nelemg),f21v(nelemg),
     > dsd5(nelemg),
     4 sig11s(nelemg),sig22s(nelemg),sig33s(nelemg),sig12s(nelemg),
     5 ddp1(nelemg,4),ddp2(nelemg,4),ddp3(nelemg,4),ddp4(nelemg,4),
     > ddp5(nelemg,4)
      common/vect4/
     1 py1(nelemg),py2(nelemg),py3(nelemg),py4(nelemg),
     2 pz1(nelemg),pz2(nelemg),pz3(nelemg),pz4(nelemg),
     3 ph1(nelemg),ph2(nelemg),ph3(nelemg),ph4(nelemg)
      common/vect6/
     1 cg11(nelemg),cg21(nelemg),cg12(nelemg),cg22(nelemg),
     2 cg13(nelemg),cg23(nelemg),cg14(nelemg),cg24(nelemg),
     3 ch11(nelemg),ch21(nelemg),ch12(nelemg),ch22(nelemg),
     4 ch13(nelemg),ch23(nelemg),ch14(nelemg),ch24(nelemg),
     5 ci11(nelemg),ci21(nelemg),ci12(nelemg),ci22(nelemg),
     6 ci13(nelemg),ci23(nelemg),ci14(nelemg),ci24(nelemg)
      common/vect7/
     1 diavg(nelemg),dfavg(nelemg),volnd(nelemg),
     > volcd(nelemg),vlinc(nelemg),volgp(nelemg),
     2 dimod(nelemg),sig11(nelemg),sig22(nelemg),
     > sig33(nelemg),sig12(nelemg)
      common/vect8/
     1 d(16,nelemg)
      common/vect10/volume(nelemg),
     1 pavgy1(nelemg),pavgy2(nelemg),pavgy3(nelemg),pavgy4(nelemg),
     2 pavgz1(nelemg),pavgz2(nelemg),pavgz3(nelemg),pavgz4(nelemg),
     3 pavgh1(nelemg),pavgh2(nelemg),pavgh3(nelemg),pavgh4(nelemg)
      common/vect11/row1(nelemg),row2(nelemg),row3(nelemg),
     > row4(nelemg),
     1 bmy1(nelemg),bmy2(nelemg),bmy3(nelemg),bmy4(nelemg),
     2 bmz1(nelemg),bmz2(nelemg),bmz3(nelemg),bmz4(nelemg),
     3 bmh1(nelemg),bmh2(nelemg),bmh3(nelemg),bmh4(nelemg)
      common/vect12/svvol(nelemg,4),
     1 spy1(nelemg,4),spy2(nelemg,4),spy3(nelemg,4),spy4(nelemg,4),
     2 spz1(nelemg,4),spz2(nelemg,4),spz3(nelemg,4),spz4(nelemg,4),
     3 sph1(nelemg,4),sph2(nelemg,4),sph3(nelemg,4),sph4(nelemg,4)
      common/vect13/
     1 x1112(nelemg),x1314(nelemg),x1114(nelemg),x1213(nelemg),
     2 x2122(nelemg),x2324(nelemg),x2124(nelemg),x2223(nelemg),
     3 hgmod1(nelemg),hgmod2(nelemg),hsv1(nelemg),hsv2(nelemg),
     > hsv3(nelemg),hsv4(nelemg),
     4 hg1(nelemg),hg2(nelemg)
      common/vect15/
     1 d1(nelemg),d2(nelemg),d3(nelemg),d4(nelemg)
      common/vect18/
     1 sclo(nelemg),volmd(nelemg),volmd4(nelemg,4)
      common/vect19/
     1 sdy1(nelemg),sdy2(nelemg),sdy3(nelemg),sdy4(nelemg),
     2 sdz1(nelemg),sdz2(nelemg),sdz3(nelemg),sdz4(nelemg)
      common/vect20/
     1 sdgt1(nelemg,4),sdgt2(nelemg,4),sdgt3(nelemg,4),
     2 sdgt4(nelemg,4),sdgt5(nelemg,4),dfdet(nelemg,4)
      common/vect99/streff(nelemg)
!c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      common/berg/dels11(nelemg),dels22(nelemg),dels33(nelemg),
     > dels12(nelemg),rb(nelemg,8)
!c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      common/knn/k1
      common/hokao/lst
      common/cccc/ntt
      common/effmod/effmod(nelemg)

      common/main_block/ a(1)

!c----------------------------------------------
!c.....add up common block meng and bk26 by m.k.
!c----------------------------------------------
      common/vect13b/
     1 hgmod1c(nelemg),hgmod2c(nelemg),hsv1c(nelemg),
     > hsv2c(nelemg),hsv3c(nelemg),
     2 hsv4c(nelemg),hg1c(nelemg),hg2c(nelemg)
      common/vect4b/
     1py1c(nelemg),py2c(nelemg),py3c(nelemg),py4c(nelemg),
     2pz1c(nelemg),pz2c(nelemg),pz3c(nelemg),pz4c(nelemg),
     3ph1c(nelemg),ph2c(nelemg),ph3c(nelemg),ph4c(nelemg)
      common/meng/
     1 f11dd(nelemg),f22dd(nelemg),f12dd(nelemg),f21dd(nelemg),
     > f33dd(nelemg),
     2 f11dv(nelemg),f22dv(nelemg),f12dv(nelemg),f21dv(nelemg),
     > f33dv(nelemg),
     3 vg11(nelemg),vg22(nelemg),vg12(nelemg),vg21(nelemg),
     > vg33(nelemg),fdet(nelemg),    
     4 dg11(nelemg),dg22(nelemg),dg33(nelemg),dg12(nelemg),
     > spin(nelemg)
      common/bk26/dt,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
      common/hourglass/fhg(40000,8),fhghis(40000,8),fhg1(nelemg),
     > fhg2(nelemg),
     1 fhg3(nelemg),fhg4(nelemg),fhg5(nelemg),fhg6(nelemg),
     > fhg7(nelemg),fhg8(nelemg)
      common/hgenergypass/hgener(nelemg),totener(nelemg),
     > inertia(nelemg)
      common/intener/intener(nelemg)
      common/hourglassstress/hgstress1c(nelemg),hgstress2c(nelemg)
	  common /stressflag/ ElemFractCode(nume),ElemDecayCount(nume)
	  common /overlapping/ intersec(4, nume), area_coeff(nume),
     > update_flag
	  common /strvect/ stra(nelemg,4)
!c.....END    
      integer ink,llt,lft,nftm1
      dimension matp(*),yz(8,nelemg),hgs(*),prop(48,*),den(*)
	  dimension fd11c(nelemg),fd21c(nelemg),fd12c(nelemg),
     > fd22c(nelemg),
     1  fd13c(nelemg),fd23c(nelemg),fd14c(nelemg),fd24c(nelemg)
      dimension cg11c(nelemg),cg21c(nelemg),cg12c(nelemg),
     > cg22c(nelemg),
     1 cg13c(nelemg),cg23c(nelemg),cg14c(nelemg),cg24c(nelemg)
      equivalence (lpar(1),model),(lpar(5),ityp2d)
      real intener, area_coeff, stra
      real Aratio
	  integer ElemFractCode, update_flag
!c.....Parameter for hourglass control type
!c     3=stiffness hourglass control
!c     2=viscous hourglass control
      data ihgtyp/3/

!      do i=1,1000
!      enddo
!c        write(7777,*) '-- stfnf.f'
!      write(*,*) 'in stfnf, lft=,llt=',lft,llt
!c      write(*,*) ' '
!c.....initial nodal coordinates 
      ntt=ntt+1
!$OMP PARALLEL DO       
      do 10 i=lft,llt
      ci11(i)=yz(1,i)
      ci21(i)=yz(2,i)
      ci12(i)=yz(3,i)
      ci22(i)=yz(4,i)
      ci13(i)=yz(5,i)
      ci23(i)=yz(6,i)
      ci14(i)=yz(7,i)
      ci24(i)=yz(8,i)
!   10 continue
!c.....updating nodal coordinates at step n+1
!      do 20 i=lft,llt
      cg11(i)=ci11(i)+ed11(i)
      cg21(i)=ci21(i)+ed21(i)
      cg12(i)=ci12(i)+ed12(i)
      cg22(i)=ci22(i)+ed22(i)
      cg13(i)=ci13(i)+ed13(i)
      cg23(i)=ci23(i)+ed23(i)
      cg14(i)=ci14(i)+ed14(i)
      cg24(i)=ci24(i)+ed24(i)
!   20 continue
!c.....updating nodal coordinates at step n+1/2
!      do 30 i=lft,llt
      ch11(i)=cg11(i)-.50*fd11(i)
      ch21(i)=cg21(i)-.50*fd21(i)
      ch12(i)=cg12(i)-.50*fd12(i)
      ch22(i)=cg22(i)-.50*fd22(i)
      ch13(i)=cg13(i)-.50*fd13(i)
      ch23(i)=cg23(i)-.50*fd23(i)
      ch14(i)=cg14(i)-.50*fd14(i)
      ch24(i)=cg24(i)-.50*fd24(i)
!   30 continue
!c
!      do 40 i=lft,llt
      diavg(i)=0.
      dfavg(i)=0.
      volnd(i)=0.
      volcd(i)=0.
      x1112(i)=ci11(i)-ci12(i)
      x2122(i)=ci21(i)-ci22(i)
      x1114(i)=ci11(i)-ci14(i)
      x2124(i)=ci21(i)-ci24(i)
      x1314(i)=ci13(i)-ci14(i)
      x2324(i)=ci23(i)-ci24(i)
      x1213(i)=ci12(i)-ci13(i)
      x2223(i)=ci22(i)-ci23(i)
!   40 continue
   10 continue
!c
!c    ityp2d: Axisym=0, pl strain=1, pl stress=2
      scale=float(min(1,ityp2d))

      do 70 lst=1,maxint
      k1=lst
      ipt=lst
      nintg=ipt
      if (maxint.eq.1) nintg=5
!c     obtain derivatives of shape functions
!c     ciij are the initial coords., vlinc is the Jacobian
!c     this gives the derivatives of the shape functions with respect to the
!c     reference configuration X, pyi and pzi where i=1..4
      call strdsp (ci11,ci12,ci13,ci14,vlinc)

c
!$OMP PARALLEL DO       
      do 50 i=lft,llt
	  if (vlinc(i).lt.0.0) write(*,*) 'Negative Jacobian',i
!c     obtain du/dX through sum(d(Ni)/dX*ui)
      dgi1(i,ipt)=py1(i)*fd11(i)+py2(i)*fd12(i)+py3(i)*fd13(i)+py4(i)
     1 *fd14(i)
      dgi2(i,ipt)=pz1(i)*fd21(i)+pz2(i)*fd22(i)+pz3(i)*fd23(i)+pz4(i)
     1 *fd24(i)
      dgi3(i,ipt)=pz1(i)*fd11(i)+pz2(i)*fd12(i)+pz3(i)*fd13(i)+pz4(i)
     1 *fd14(i)
      dgi4(i,ipt)=py1(i)*fd21(i)+py2(i)*fd22(i)+py3(i)*fd23(i)+py4(i)
     1 *fd24(i)
!c     dgi5 is an axisymmetric term
      dgi5(i,ipt)=ph1(i)*fd11(i)+ph2(i)*fd12(i)+ph3(i)*fd13(i)+ph4(i)
     1 *fd14(i)
!c     volume calculation used later with half-step
      volmd4(i,ipt)=1./vlinc(i)
!   50 continue
!c     obtain d(x+u)/dX through sum(d(Ni)/dX*(xi+ui)) where xi is at step n
!      do 60 i=lft,llt
      sdgt1(i,ipt)=py1(i)*cg11(i)+py2(i)*cg12(i)+py3(i)*cg13(i)+py4(i)
     1 *cg14(i)
      sdgt2(i,ipt)=pz1(i)*cg21(i)+pz2(i)*cg22(i)+pz3(i)*cg23(i)+pz4(i)
     1 *cg24(i)
      sdgt3(i,ipt)=pz1(i)*cg11(i)+pz2(i)*cg12(i)+pz3(i)*cg13(i)+pz4(i)
     1 *cg14(i)
      sdgt4(i,ipt)=py1(i)*cg21(i)+py2(i)*cg22(i)+py3(i)*cg23(i)+py4(i)
     1 *cg24(i)
!c     sdgt5 is an axisymmetric term, scale=1 for plane strain
      sdgt5(i,ipt)=ph1(i)*cg11(i)+ph2(i)*cg12(i)+ph3(i)*cg13(i)+ph4(i)
     1 *cg14(i)+scale
!c     determinant of the deformation gradient at step n+1
      dfdet(i,ipt)=(sdgt1(i,ipt)*sdgt2(i,ipt)-sdgt3(i,ipt)
     1 *sdgt4(i,ipt))*sdgt5(i,ipt)
!c     volume increment to or from x to psi to deformed
      dfavg(i)=dfavg(i)+vlinc(i)*dfdet(i,ipt)
!c     Fij at step n obtained by d(x+u)/dX-du/dX
      dgt1(i,ipt)=sdgt1(i,ipt)-dgi1(i,ipt)
      dgt2(i,ipt)=sdgt2(i,ipt)-dgi2(i,ipt)
      dgt3(i,ipt)=sdgt3(i,ipt)-dgi3(i,ipt)
      dgt4(i,ipt)=sdgt4(i,ipt)-dgi4(i,ipt)
!c     volume change to or from x to psi
      volnd(i)=volnd(i)+vlinc(i)
!   60 continue
!      do i=lft,llt
!c     obtain strain component du/dy, dv/dz, du/dz, dv/dy
      stra(i,1)=py1(i)*ed11(i)+py2(i)*ed12(i)+py3(i)*ed13(i)+py4(i)   !du/dy
     1 *ed14(i)
      stra(i,2)=pz1(i)*ed21(i)+pz2(i)*ed22(i)+pz3(i)*ed23(i)+pz4(i)   !dv/dz
     1 *ed24(i)
      stra(i,3)=pz1(i)*ed11(i)+pz2(i)*ed12(i)+pz3(i)*ed13(i)+pz4(i)   !du/dz
     1 *ed14(i)
      stra(i,4)=py1(i)*ed21(i)+py2(i)*ed22(i)+py3(i)*ed23(i)+py4(i)   !dv/dy
     1 *ed24(i)
!      end do
   50 continue
   70 continue
!c     changing basis for strain-displacement calculations to half-step
!$OMP PARALLEL DO       
      do 80 i=lft,llt
      x1112(i)=ch11(i)-ch12(i)
      x2122(i)=ch21(i)-ch22(i)
      x1114(i)=ch11(i)-ch14(i)
      x2124(i)=ch21(i)-ch24(i)
      x1314(i)=ch13(i)-ch14(i)
      x2324(i)=ch23(i)-ch24(i)
      x1213(i)=ch12(i)-ch13(i)
      x2223(i)=ch22(i)-ch23(i)
   80 continue

      do 110 lst=1,maxint
!c.....maxint = max. integeration points
      k1=lst
      ipt=lst
      nintg=ipt
      if (maxint.eq.1) nintg=5
!c     
!c     calling strain-displacement to set up strain increment at 
!c     half-step, used by Nike2D subroutines
      call strdsp (ch11,ch12,ch13,ch14,vlinc)

!c     calculation of deformation gradient increment in terms of the
!c     half-step configuration
!!!!!$OMP PARALLEL DO       
      do 90 i=lft,llt
      ddp1(i,ipt)=py1(i)*fd11(i)+py2(i)*fd12(i)+py3(i)*fd13(i)+py4(i)
     1 *fd14(i)
      ddp2(i,ipt)=pz1(i)*fd21(i)+pz2(i)*fd22(i)+pz3(i)*fd23(i)+pz4(i)
     1 *fd24(i)
      ddp3(i,ipt)=pz1(i)*fd11(i)+pz2(i)*fd12(i)+pz3(i)*fd13(i)+pz4(i)
     1 *fd14(i)
      ddp4(i,ipt)=py1(i)*fd21(i)+py2(i)*fd22(i)+py3(i)*fd23(i)+py4(i)
     1 *fd24(i)
      ddp5(i,ipt)=ph1(i)*fd11(i)+ph2(i)*fd12(i)+ph3(i)*fd13(i)+ph4(i)
     1 *fd14(i)
!c     used in conjunction with volmd4 calculated from full step
      volmd4(i,ipt)=volmd4(i,ipt)*vlinc(i)
!   90 continue
!      do 100 i=lft,llt
      diavg(i)=diavg(i)+vlinc(i)*(ddp1(i,ipt)+ddp2(i,ipt)+ddp5(i,ipt))
      volcd(i)=volcd(i)+vlinc(i)
!  100 continue
   90 continue

  110 continue

!c     changing basis for strain-displacement calculations to full step
!$OMP PARALLEL DO       
      do 120 i=lft,llt
      x1112(i)=cg11(i)-cg12(i)
      x2122(i)=cg21(i)-cg22(i)
      x1114(i)=cg11(i)-cg14(i)
      x2124(i)=cg21(i)-cg24(i)
      x1314(i)=cg13(i)-cg14(i)
      x2324(i)=cg23(i)-cg24(i)
      x1213(i)=cg12(i)-cg13(i)
      x2223(i)=cg22(i)-cg23(i)
      diavg(i)=diavg(i)/(volcd(i)+1.e-30)
      dfavg(i)=dfavg(i)/(volnd(i)+1.e-30)
  120 continue

!c     calculation of deviatoric for B-bar calculations?
      if (maxint.eq.1) go to 150
      do 140 ipt=1,4
      do 130 i=lft,llt
      dimod(i)=sclo(i)*(diavg(i)-ddp1(i,ipt)-ddp2(i,ipt)-
     > ddp5(i,ipt))/3.
      ddp1(i,ipt)=ddp1(i,ipt)+dimod(i)
      ddp2(i,ipt)=ddp2(i,ipt)+dimod(i)
  130 ddp5(i,ipt)=ddp5(i,ipt)+dimod(i)
  140 continue
   
!c     element volume need for B-Bar
  150 call scm0 (volume,13*nmel)
      do 190 lst=1,maxint
      k1=lst
      ipt=lst
      nintg=ipt
      if (maxint.eq.1) nintg=5
c
!c     calling strain-displacement to set up strain increment at 
!c     full-step, used for calculations of B, Bbar
      call strdsp (cg11,cg12,cg13,cg14,volgp)
!c
!c     ibar equals 2 for the single integration point case.  
!c     it might be worth switching the 0 to 2 here, because
!c     these calculations are only used for the stffns.f
!c     subroutine, which is only called when ibar=1 
!c     ibar=ibbari+1, ibbari input as 0 for Bbar element formulation
!c     and 1 for Nike2D element formulation.  For the cases of 
!c     single point integration and plane stress, ibbari is changed
!c     to 1, and therefore ibar=2

!c     B-Bar Calculation 
      if (ibar.eq.0) go to 170
!$OMP PARALLEL DO       
      do 160 i=lft,llt
      volume(i)=volume(i)+volgp(i)
      pavgy1(i)=pavgy1(i)+volgp(i)*py1(i)
      pavgy2(i)=pavgy2(i)+volgp(i)*py2(i)
      pavgy3(i)=pavgy3(i)+volgp(i)*py3(i)
      pavgy4(i)=pavgy4(i)+volgp(i)*py4(i)
      pavgz1(i)=pavgz1(i)+volgp(i)*pz1(i)
      pavgz2(i)=pavgz2(i)+volgp(i)*pz2(i)
      pavgz3(i)=pavgz3(i)+volgp(i)*pz3(i)
      pavgz4(i)=pavgz4(i)+volgp(i)*pz4(i)
      pavgh1(i)=pavgh1(i)+volgp(i)*ph1(i)
      pavgh2(i)=pavgh2(i)+volgp(i)*ph2(i)
      pavgh3(i)=pavgh3(i)+volgp(i)*ph3(i)
      pavgh4(i)=pavgh4(i)+volgp(i)*ph4(i)
  160 continue
!$OMP PARALLEL DO       
  170 do 180 i=lft,llt
      spy1(i,lst)=py1(i)
      spy2(i,lst)=py2(i)
      spy3(i,lst)=py3(i)
      spy4(i,lst)=py4(i)
      spz1(i,lst)=pz1(i)
      spz2(i,lst)=pz2(i)
      spz3(i,lst)=pz3(i)
      spz4(i,lst)=pz4(i)
      sph1(i,lst)=ph1(i)
      sph2(i,lst)=ph2(i)
      sph3(i,lst)=ph3(i)
      sph4(i,lst)=ph4(i)
      svvol(i,lst)=volgp(i)
  180 continue
  190 continue
      if (ibar.eq.0) go to 210
      third=1./3.
!c     1/3*Bbari
!$OMP PARALLEL DO       
      do 200 i=lft,llt
      pavgy1(i)=third*(pavgy1(i)+pavgh1(i))/volume(i)
      pavgz1(i)=third*pavgz1(i)/volume(i)
      pavgy2(i)=third*(pavgy2(i)+pavgh2(i))/volume(i)
      pavgz2(i)=third*pavgz2(i)/volume(i)
      pavgy3(i)=third*(pavgy3(i)+pavgh3(i))/volume(i)
      pavgz3(i)=third*pavgz3(i)/volume(i)
      pavgy4(i)=third*(pavgy4(i)+pavgh4(i))/volume(i)
      pavgz4(i)=third*pavgz4(i)/volume(i)
  200 continue

  210 do 430 lst=1,maxint
      ipt=lst
      nintg=ipt

!$OMP PARALLEL DO       
      do 220 i=lft,llt
      py1(i)=spy1(i,lst)
      py2(i)=spy2(i,lst)
      py3(i)=spy3(i,lst)
      py4(i)=spy4(i,lst)
      pz1(i)=spz1(i,lst)
      pz2(i)=spz2(i,lst)
      pz3(i)=spz3(i,lst)
      pz4(i)=spz4(i,lst)
      ph1(i)=sph1(i,lst)
      ph2(i)=sph2(i,lst)
      ph3(i)=sph3(i,lst)
      ph4(i)=sph4(i,lst)
      volgp(i)=svvol(i,lst)
      volmd(i)=volmd4(i,lst)
  220 continue

!c     deformation gradient dx/dX at step n
      do 230 i=lft,llt
      f11v(i)=dgt1(i,ipt)
      f22v(i)=dgt2(i,ipt)
      f12v(i)=dgt3(i,ipt)
      f21v(i)=dgt4(i,ipt)
  230 continue

!c     get rotation tensor at step n, stored as r's
!c     rotation calculated in terms of Fij at step n
      call getvrt (r11,r22,r12,r21)

!c     deformation gradient dx/dX at step n+1/2      
!$OMP PARALLEL DO       
	  do 240 i=lft,llt
      f11v(i)=dgt1(i,ipt)+.50*dgi1(i,ipt)
      f22v(i)=dgt2(i,ipt)+.50*dgi2(i,ipt)
      f12v(i)=dgt3(i,ipt)+.50*dgi3(i,ipt)
      f21v(i)=dgt4(i,ipt)+.50*dgi4(i,ipt)
  240 continue
!c
!c     get rotation tensor at step n+1/2, stored as q's
!c     rotation calculated in terms of Fij at step n+1/2
      call getvrt (q11,q22,q12,q21)
!c
!c     deformation gradient dx/dX at step n+1
!$OMP PARALLEL DO       
      do 250 i=lft,llt
      f11v(i)=dgt1(i,ipt)+dgi1(i,ipt)
      f22v(i)=dgt2(i,ipt)+dgi2(i,ipt)
      f12v(i)=dgt3(i,ipt)+dgi3(i,ipt)
      f21v(i)=dgt4(i,ipt)+dgi4(i,ipt)
  250 continue
!c
!c--------------------------------------------------------
!c.....obtain velocity gradient from F(n) & F(n+1) 
!c--------------------------------------------------------
!c     this uses n+1 deformation gradient and the difference
!c     between n+1 and n deformation gradient.  The values
!c     dg11, dg22, dg33, dg12, and spin are sent to matpoly
!$OMP PARALLEL DO       
      do 255 i=lft,llt
      f11dd(i)=sdgt1(i,ipt)
      f22dd(i)=sdgt2(i,ipt)
      f12dd(i)=sdgt3(i,ipt)
      f21dd(i)=sdgt4(i,ipt)
      f33dd(i)=sdgt5(i,ipt)
      f11dv(i)=dgi1(i,ipt)/dt
      f22dv(i)=dgi2(i,ipt)/dt
      f12dv(i)=dgi3(i,ipt)/dt
      f21dv(i)=dgi4(i,ipt)/dt
      f33dv(i)=dgi5(i,ipt)/dt
     
      fdet(i)=f11dd(i)*f22dd(i)-f12dd(i)*f21dd(i)
      vg11(i)=(f11dv(i)*f22dd(i)-f12dv(i)*f21dd(i))/fdet(i)
      vg22(i)=(f11dd(i)*f22dv(i)-f21dv(i)*f12dd(i))/fdet(i)
      vg12(i)=(f11dd(i)*f12dv(i)-f11dv(i)*f12dd(i))/fdet(i)
      vg21(i)=(f21dv(i)*f22dd(i)-f22dv(i)*f21dd(i))/fdet(i)
      vg33(i)=f33dv(i)/f33dd(i)
!c------------------------------
!c.....green - st. venant strain
!c------------------------------
!ck      gstrain(1)=0.5*(f11dd(i)*f11dd(i)+f21dd(i)*f21dd(i)-1.0)
!ck      gstrain(2)=0.5*(f12dd(i)*f12dd(i)+f22dd(i)*f22dd(i)-1.0)
!ck      gstrain(3)=0.5*(f11dd(i)*f12dd(i)+f21dd(i)*f22dd(i))
!c
!c------------------------------------------------------------------
!c.....obtain the rate of deformation gradient & spin tensor by m.k.
!c------------------------------------------------------------------
      dg11(i)=vg11(i)
      dg22(i)=vg22(i)
      dg33(i)=vg33(i)
      dg12(i)=0.5*(vg12(i)+vg21(i))
c
      spin(i)=0.5*(vg21(i)-vg12(i))      
      sp11=0.

c      write(*,*) dg11(i),dg22(i),dg33(i),dg12(i)
!c
!ck      if (lst.eq.1) then
!ck      write(27,248) sdgt1(i,ipt),sdgt3(i,ipt),sdgt4(i,ipt),
!ck     .sdgt2(i,ipt),dfdet(i,ipt)
!ck  248 format(5x,'F : f11=',e14.8,2x,'f12=',e14.8,2x,'f21=',e14.8,
!ck     .       2x,'f22=',e14.8,2x,'fdet=',e14.8/)
!ck      write(27,251) r11(i),r12(i),r21(i),r22(i)
!ck  251 format(5x,'R : r11=',e14.8,2x,'r12=',e14.8,2x,'r21=',e14.8,
!ck     .       2x,'r22=',f7.4/)
!ck      write(27,252) vg11(i),vg12(i),vg21(i),vg22(i)
!ck  252 format(5x,'L : v11=',e14.8,2x,'v12=',e14.8,2x,'v21=',e14.8,
!ck     .       2x,'v22=',e14.8/)
!ck      write(27,253) dg11(i),dg12(i),dg12(i),dg22(i),dg33(i)
!ck  253 format(5x,'D : d11=',e14.8,2x,'d12=',e14.8,2x,'d21=',e14.8,
!ck     .       2x,'d22=',e14.8,2x,'d33=',e14.8/)
!ck      write(27,254) sp11,spin(i),spin(i),sp11
!ck  254 format(5x,'W : w11=',e14.8,2x,'w12=',e14.8,2x,'w21= -',e14.8,
!ck     .       2x,'w22=',e14.8/)
!ck      write(27,267) gstrain(3)
!ck  267 format(5x,'strain-xy =',f8.4/)
!ck      endif
!c  
  255 continue 
!c.....END
!c
!c     get rotation tensor at step n+1, stored as s's
!c     rotation calculated in terms of Fij at step n+1
!      call getvrt (s11,s22,s12,s21)
!c
!c     This rotates the half-step deformation gradient increment used by 
!c     Nike2D subroutines back to the original configuration by qTdq
!c     where q is the half-step rotation tensor
      call rotstr (ddp1(1,ipt),ddp2(1,ipt),ddp3(1,ipt),ddp4(1,ipt)
     > ,ddp5(1,ipt))

      mft=1
      mlt=1
      llt1=llt-1
      nel=nftm1+1
      if (lft.eq.llt) go to 280
      do 270 i=lft,llt1

      if (matp(i).eq.matp(i+1)) go to 270
      model=mtype(mft)
      nelsub=nftm1+mft
!c
!c  this is where we call our material subroutines to update any necessary parts of Kij
!c     and get the stress values
      call sslcs (a,model)
!c.... island material group gather
!c     imeth is set to 0, so the next few if statements do not apply
      matpc=-matp(mft)
      if(imeth.gt.0)call xpegat(matpc)

      mft=mlt+1
  270 mlt=mlt+1
  280 model=mtype(mft)
      nelsub=nftm1+mft

      call sslcs (a,model)
!c.... island material group gather
      matpc=-matp(mft)
      if(imeth.gt.0)call xpegat(matpc)
!c.... island element group gather
      if(imeth.gt.0)call xpegat(3)

!c     this does not apply, appears to be some kind of incremental BT*sigma
      if(imeth.gt.0.and.imass.ne.0)then
!$OMP PARALLEL DO       
      do 301 i=lft,llt
      dels11(i)=dels11(i)*volgp(i)
      dels22(i)=dels22(i)*volgp(i)
      dels33(i)=dels33(i)*volgp(i)
      dels12(i)=dels12(i)*volgp(i)
!  301 continue
!      do 302 i=lft,llt
      sdy1(i)=py1(i)*dels11(i)+pz1(i)*dels12(i)
      sdz1(i)=pz1(i)*dels22(i)+py1(i)*dels12(i)
      sdy2(i)=py2(i)*dels11(i)+pz2(i)*dels12(i)
      sdz2(i)=pz2(i)*dels22(i)+py2(i)*dels12(i)
      sdy3(i)=py3(i)*dels11(i)+pz3(i)*dels12(i)
      sdz3(i)=pz3(i)*dels22(i)+py3(i)*dels12(i)
      sdy4(i)=py4(i)*dels11(i)+pz4(i)*dels12(i)
      sdz4(i)=pz4(i)*dels22(i)+py4(i)*dels12(i)
!  302 continue
!      do 303 i=lft,llt
      rb(i,1)=rb(i,1)+sdy1(i)+ph1(i)*dels33(i)
      rb(i,2)=rb(i,2)+sdz1(i)
      rb(i,3)=rb(i,3)+sdy2(i)+ph2(i)*dels33(i)
      rb(i,4)=rb(i,4)+sdz2(i)
      rb(i,5)=rb(i,5)+sdy3(i)+ph3(i)*dels33(i)
      rb(i,6)=rb(i,6)+sdz3(i)
      rb(i,7)=rb(i,7)+sdy4(i)+ph4(i)*dels33(i)
      rb(i,8)=rb(i,8)+sdz4(i)
!  303 continue
!!!!!!!!!!      do 304 i=lft,llt
!!!!!!!!!      xnk2d(1)=xnk2d(1)+rb(i,1)*fd11(i)+rb(i,2)*fd21(i)
!!!!!!!!!     1 +rb(i,3)*fd12(i)+rb(i,4)*fd22(i)
!!!!!!!!!     2 +rb(i,5)*fd13(i)+rb(i,6)*fd23(i)
!!!!!!!!!     3 +rb(i,7)*fd14(i)+rb(i,8)*fd24(i)
!!!!!!!!!!  304 continue
  301 continue
      endif
!c
!c
!c     Beginning of one point integration of BT*sigma
      if (maxint.ne.1) go to 310
      do 300 i=lft,llt
      volgp(i)=4.*volgp(i)
  300 continue

  310 do 320 i=lft,llt
      sig11(i)=sig11s(i)*volgp(i)
      sig22(i)=sig22s(i)*volgp(i)
      sig33(i)=sig33s(i)*volgp(i)
  320 sig12(i)=sig12s(i)*volgp(i)
!c
!c     This is modifying the stiffness values d that are 
!c     passed from matpoly
      if(iphase.ne.3)then
      do 330 i=lft,llt
      do 330 j=1,16
      d(j,i)=d(j,i)*volgp(i)
  330 continue
      endif
!c
!c     once again, this appears to be part of the Bbar formulation
!c     which only applies to the stffns.f subroutine
      if (ibar.eq.0) go to 350
      third=-1./3.
!$OMP PARALLEL DO       
      do 340 i=lft,llt
      row1(i)=sig11(i)+sig22(i)+sig33(i)
      bmy1(i)=sclo(i)*(pavgy1(i)+third*(py1(i)+ph1(i)))
      bmz1(i)=sclo(i)*(pavgz1(i)+third*pz1(i))
      bmy2(i)=sclo(i)*(pavgy2(i)+third*(py2(i)+ph2(i)))
      bmz2(i)=sclo(i)*(pavgz2(i)+third*pz2(i))
      bmy3(i)=sclo(i)*(pavgy3(i)+third*(py3(i)+ph3(i)))
      bmz3(i)=sclo(i)*(pavgz3(i)+third*pz3(i))
      bmy4(i)=sclo(i)*(pavgy4(i)+third*(py4(i)+ph4(i)))
      bmz4(i)=sclo(i)*(pavgz4(i)+third*pz4(i))
  340 continue
c
c     This is the final result of the stress update from BT*sigma
c     Is this incremental??
!!!!!$OMP PARALLEL DO       
  350 do 360 i=lft,llt
      sdy1(i)=py1(i)*sig11(i)+pz1(i)*sig12(i)
      sdz1(i)=pz1(i)*sig22(i)+py1(i)*sig12(i)
      sdy2(i)=py2(i)*sig11(i)+pz2(i)*sig12(i)
      sdz2(i)=pz2(i)*sig22(i)+py2(i)*sig12(i)
      sdy3(i)=py3(i)*sig11(i)+pz3(i)*sig12(i)
      sdz3(i)=pz3(i)*sig22(i)+py3(i)*sig12(i)
      sdy4(i)=py4(i)*sig11(i)+pz4(i)*sig12(i)
      sdz4(i)=pz4(i)*sig22(i)+py4(i)*sig12(i)
  360 continue
!c
!c     nodal forces
!$OMP PARALLEL DO       
      do 370 i=lft,llt
      r1(i)=r1(i)+sdy1(i)+ph1(i)*sig33(i)
      r2(i)=r2(i)+sdz1(i)
      r3(i)=r3(i)+sdy2(i)+ph2(i)*sig33(i)
      r4(i)=r4(i)+sdz2(i)
      r5(i)=r5(i)+sdy3(i)+ph3(i)*sig33(i)
      r6(i)=r6(i)+sdz3(i)
      r7(i)=r7(i)+sdy4(i)+ph4(i)*sig33(i)
      r8(i)=r8(i)+sdz4(i)
!c      write(*,*) r1(i),r2(i),r3(i),r4(i)
!c      write(*,*) r5(i),r6(i),r7(i),r8(i),i
  370 continue		      

      if (maxint.ne.1) go to 396
!c
!c     Hourglass Stiffness Algorithms
!c     Our ihgtyp is set to 1, go to 381
      if (ihgtyp.eq.2) go to 388
	  if (ihgtyp.eq.3) go to 392
      if (ihgtyp.eq.1) go to 381
!$OMP PARALLEL DO       
      do 380 i=lft,llt
      hgmod1(i)=hgs(i)*(ed11(i)-ed12(i)+ed13(i)-ed14(i))
      hgmod2(i)=hgs(i)*(ed21(i)-ed22(i)+ed23(i)-ed24(i))
      r1(i)=r1(i)+hgmod1(i)
      r2(i)=r2(i)+hgmod2(i)
      r3(i)=r3(i)-hgmod1(i)
      r4(i)=r4(i)-hgmod2(i)
      r5(i)=r5(i)+hgmod1(i)
      r6(i)=r6(i)+hgmod2(i)
      r7(i)=r7(i)-hgmod1(i)
      r8(i)=r8(i)-hgmod2(i)
  380 continue

      go to 396
	  
!c     Hourglass stiffness algorithm 2
!!!!!$OMP PARALLEL DO       
  381 do 383 i=lft,llt
      hgmod1(i)=cg11(i)-cg12(i)+cg13(i)-cg14(i)
      hgmod2(i)=cg21(i)-cg22(i)+cg23(i)-cg24(i)
      hsv1(i)  =1.-hgmod1(i)*py1(i)-hgmod2(i)*pz1(i)
      hsv2(i)  =-1.-hgmod1(i)*py2(i)-hgmod2(i)*pz2(i)
      hsv3(i)  =1.-hgmod1(i)*py3(i)-hgmod2(i)*pz3(i)
      hsv4(i)  =-1.-hgmod1(i)*py4(i)-hgmod2(i)*pz4(i)
  383 continue
!!!!$OMP PARALLEL DO       
      do 385 i=lft,llt
      hg1(i)=hgs(i)*(hsv1(i)*ed11(i)+hsv2(i)*ed12(i)
     1             + hsv3(i)*ed13(i)+hsv4(i)*ed14(i))
      hg2(i)=hgs(i)*(hsv1(i)*ed21(i)+hsv2(i)*ed22(i)
     2             + hsv3(i)*ed23(i)+hsv4(i)*ed24(i))
  385 continue
!c
!c     Contribution of hourglass stiffness to nodal forces
!$OMP PARALLEL DO       
      do 387 i=lft,llt
      r1(i)=r1(i)+hg1(i)*hsv1(i)
      r2(i)=r2(i)+hg2(i)*hsv1(i)
      r3(i)=r3(i)+hg1(i)*hsv2(i)
      r4(i)=r4(i)+hg2(i)*hsv2(i)
      r5(i)=r5(i)+hg1(i)*hsv3(i)
      r6(i)=r6(i)+hg2(i)*hsv3(i)
      r7(i)=r7(i)+hg1(i)*hsv4(i)
      r8(i)=r8(i)+hg2(i)*hsv4(i)
  387 continue
      go to 396
!c     hourglass stiffnes algorithm 3
!!!!!$OMP PARALLEL DO       
  388 do 389 i=lft,llt
!c     Hourglass viscosity setup WML 12610
!c     add in material props (stiffness, density, etc)
!c     put hourglass stiffnesses into common block (vect13b)
      ink=i+nftm1
!c     Hourglass Coefficient CQ*psi (hgc=alpha/16*psi)
        hgs(i)=1.*volgp(i)*(py1(i)*py1(i)+py2(i)
     1 *py2(i)+py3(i)*py3(i)+py4(i)*py4(i)+pz1(i)*pz1(i)
     2  +pz2(i)*pz2(i)+pz3(i)*pz3(i)+pz4(i)*pz4(i))/2.*
     3  effmod(i)*hgc
      if (hgs(i).lt.0.) hgs(i)=0.
!c      write(*,*) den(matp(ink)),matp(ink),ink
!c        write(*,*) hgs(i)
!c      write(*,*) hgs(i)

!c     rotate x vector to corotational
      cg11c(i)=cg11(i)*r11(i)+cg21(i)*r21(i)
	  cg21c(i)=cg11(i)*r12(i)+cg21(i)*r22(i)
	  cg12c(i)=cg12(i)*r11(i)+cg22(i)*r21(i)
	  cg22c(i)=cg12(i)*r12(i)+cg22(i)*r22(i)
      cg13c(i)=cg13(i)*r11(i)+cg23(i)*r21(i)
	  cg23c(i)=cg13(i)*r12(i)+cg23(i)*r22(i)
	  cg14c(i)=cg14(i)*r11(i)+cg24(i)*r21(i)
	  cg24c(i)=cg14(i)*r12(i)+cg24(i)*r22(i)
!c     rotate increment to corotational
      fd11c(i)=(fd11(i)*r11(i)+fd21(i)*r21(i))
	  fd21c(i)=(fd11(i)*r12(i)+fd21(i)*r22(i))
	  fd12c(i)=(fd12(i)*r11(i)+fd22(i)*r21(i))
	  fd22c(i)=(fd12(i)*r12(i)+fd22(i)*r22(i))
      fd13c(i)=(fd13(i)*r11(i)+fd23(i)*r21(i))
	  fd23c(i)=(fd13(i)*r12(i)+fd23(i)*r22(i))
	  fd14c(i)=(fd14(i)*r11(i)+fd24(i)*r21(i))
	  fd24c(i)=(fd14(i)*r12(i)+fd24(i)*r22(i))
	  
      hgmod1c(i)=cg11c(i)-cg12c(i)+cg13c(i)-cg14c(i)
      hgmod2c(i)=cg21c(i)-cg22c(i)+cg23c(i)-cg24c(i)

	  
!c     get shape function derivatives with respect to corotational

      x1112(i)=cg11c(i)-cg12c(i)
      x2122(i)=cg21c(i)-cg22c(i)
      x1114(i)=cg11c(i)-cg14c(i)
      x2124(i)=cg21c(i)-cg24c(i)
      x1314(i)=cg13c(i)-cg14c(i)
      x2324(i)=cg23c(i)-cg24c(i)
      x1213(i)=cg12c(i)-cg13c(i)
      x2223(i)=cg22c(i)-cg23c(i)



!c.....maxint = max. integeration points
 
      nintg=5

      call strdsp2 (cg11c,cg12c,cg13c,cg14c)

!c     hourglass vectors gamma
      hsv1c(i)  =1.-hgmod1c(i)*py1c(i)-hgmod2c(i)*pz1c(i)
      hsv2c(i)  =-1.-hgmod1c(i)*py2c(i)-hgmod2c(i)*pz2c(i)
      hsv3c(i)  =1.-hgmod1c(i)*py3c(i)-hgmod2c(i)*pz3c(i)
      hsv4c(i)  =-1.-hgmod1c(i)*py4c(i)-hgmod2c(i)*pz4c(i)
  389 continue
!!!!!$OMP PARALLEL DO       
      do 390 i=lft,llt
!c     Hourglass Stresses Qi=CQ*psi*qdoti=CQ*psi*gammaT*vi
      hg1c(i)=hgs(i)*(hsv1c(i)*fd11c(i)/dt+hsv2c(i)*fd12c(i)/dt
     1             + hsv3c(i)*fd13c(i)/dt+hsv4c(i)*fd14c(i)/dt)
      hg2c(i)=hgs(i)*(hsv1c(i)*fd21c(i)/dt+hsv2c(i)*fd22c(i)/dt
     2             + hsv3c(i)*fd23c(i)/dt+hsv4c(i)*fd24c(i)/dt)
  390 continue
!c
!c     Contribution of hourglass stiffness to nodal forces
!!!!!!!$OMP PARALLEL DO       
      do 391 i=lft,llt
      ink=i+nftm1	  
      fhg1(i)=hg1c(i)*hsv1c(i)*volgp(i)
      fhg2(i)=hg2c(i)*hsv1c(i)*volgp(i)
      fhg3(i)=hg1c(i)*hsv2c(i)*volgp(i)
      fhg4(i)=hg2c(i)*hsv2c(i)*volgp(i)
      fhg5(i)=hg1c(i)*hsv3c(i)*volgp(i)
      fhg6(i)=hg2c(i)*hsv3c(i)*volgp(i)
      fhg7(i)=hg1c(i)*hsv4c(i)*volgp(i)
      fhg8(i)=hg2c(i)*hsv4c(i)*volgp(i)
	  
      r1(i)=r1(i)+r11(i)*(fhg1(i))+r12(i)*(fhg2(i))
      r2(i)=r2(i)+r21(i)*(fhg1(i))+r22(i)*(fhg2(i))
      r3(i)=r3(i)+r11(i)*(fhg3(i))+r12(i)*(fhg4(i))
      r4(i)=r4(i)+r21(i)*(fhg3(i))+r22(i)*(fhg4(i))
      r5(i)=r5(i)+r11(i)*(fhg5(i))+r12(i)*(fhg6(i))
      r6(i)=r6(i)+r21(i)*(fhg5(i))+r22(i)*(fhg6(i))
      r7(i)=r7(i)+r11(i)*(fhg7(i))+r12(i)*(fhg8(i))
      r8(i)=r8(i)+r21(i)*(fhg7(i))+r22(i)*(fhg8(i))
	  
	  


      hgener(i)=hg1c(i)*hg1c(i)/hgs(i)*dt*volgp(i)+
     1          hg2c(i)*hg2c(i)/hgs(i)*dt*volgp(i)+hgener(i)
      inertia(i)=1./8.*volgp(i)*den(matp(ink))*(fd11(i)*
     3     fd11(i)+fd21(i)*fd21(i)+fd12(i)*fd12(i)+fd22(i)*fd22(i)
     4     +fd13(i)*fd13(i)+fd23(i)*fd23(i)+fd14(i)*fd14(i)+
     5     fd24(i)*fd24(i))/dt/dt
      totener(i)=intener(i)*volgp(i)*dt+totener(i)
!c      write(*,*) r1(i),r2(i),r3(i),r4(i)
!c      write(*,*) r5(i),r6(i),r7(i),r8(i),i,'after hg'
  391 continue
      go to 396
  
!!!!!$OMP PARALLEL DO       
  392 do 393 i=lft,llt
  
!c      write(*,*) ink,effmod(i),prop(10,matp(ink)),hgs(i)
!
!c     Hourglass stiffness setup WML 21511
!c     add in material props (stiffness, density, etc)
!c     put hourglass stiffnesses into common block (vect13b)
      ink=i+nftm1
!c     hourglass coefficient CQ, hgc=alpha/16
        hgs(i)=1.*volgp(i)*(py1(i)*py1(i)+py2(i)
     1 *py2(i)+py3(i)*py3(i)+py4(i)*py4(i)+pz1(i)*pz1(i)
     2  +pz2(i)*pz2(i)+pz3(i)*pz3(i)+pz4(i)*pz4(i))/2.*
     3  effmod(i)*hgc
      if (hgs(i).lt.0.) hgs(i)=0.
!c      hgs(i)=hgs(i)*effmod(i)/prop(10,matp(ink))
!c        write(*,*) hgs(i)
!c      write(*,*) hgs(i)
! 
!c     rotate x vector to corotational
      cg11c(i)=cg11(i)*r11(i)+cg21(i)*r21(i) 
      cg21c(i)=cg11(i)*r12(i)+cg21(i)*r22(i)
      cg12c(i)=cg12(i)*r11(i)+cg22(i)*r21(i)
      cg22c(i)=cg12(i)*r12(i)+cg22(i)*r22(i)
      cg13c(i)=cg13(i)*r11(i)+cg23(i)*r21(i) 
      cg23c(i)=cg13(i)*r12(i)+cg23(i)*r22(i)
      cg14c(i)=cg14(i)*r11(i)+cg24(i)*r21(i)
      cg24c(i)=cg14(i)*r12(i)+cg24(i)*r22(i)
!c     rotate DISPLACEMENT to corotational 
      fd11c(i)=(fd11(i)*r11(i)+fd21(i)*r21(i))
      fd21c(i)=(fd11(i)*r12(i)+fd21(i)*r22(i))
      fd12c(i)=(fd12(i)*r11(i)+fd22(i)*r21(i))
      fd22c(i)=(fd12(i)*r12(i)+fd22(i)*r22(i))
      fd13c(i)=(fd13(i)*r11(i)+fd23(i)*r21(i))
      fd23c(i)=(fd13(i)*r12(i)+fd23(i)*r22(i))
      fd14c(i)=(fd14(i)*r11(i)+fd24(i)*r21(i))
      fd24c(i)=(fd14(i)*r12(i)+fd24(i)*r22(i))
	  
      hgmod1c(i)=cg11c(i)-cg12c(i)+cg13c(i)-cg14c(i)
      hgmod2c(i)=cg21c(i)-cg22c(i)+cg23c(i)-cg24c(i)

	  
!c     get shape function derivatives with respect to corotational

      x1112(i)=cg11c(i)-cg12c(i)
      x2122(i)=cg21c(i)-cg22c(i)
      x1114(i)=cg11c(i)-cg14c(i)
      x2124(i)=cg21c(i)-cg24c(i)
      x1314(i)=cg13c(i)-cg14c(i)
      x2324(i)=cg23c(i)-cg24c(i)
      x1213(i)=cg12c(i)-cg13c(i)
      x2223(i)=cg22c(i)-cg23c(i)



!c.....maxint = max. integeration points
 
      nintg=5

      call strdsp2 (cg11c,cg12c,cg13c,cg14c)

!c     hourglass vectors gamma
      hsv1c(i)  =1.-hgmod1c(i)*py1c(i)-hgmod2c(i)*pz1c(i)
      hsv2c(i)  =-1.-hgmod1c(i)*py2c(i)-hgmod2c(i)*pz2c(i)
      hsv3c(i)  =1.-hgmod1c(i)*py3c(i)-hgmod2c(i)*pz3c(i)
      hsv4c(i)  =-1.-hgmod1c(i)*py4c(i)-hgmod2c(i)*pz4c(i)
  393 continue
!$OMP PARALLEL DO       
      do 394 i=lft,llt
!c     hourglass strains times CQ, or Qidot
      hg1c(i)=hgs(i)*(hsv1c(i)*fd11c(i)+hsv2c(i)*fd12c(i)
     1             + hsv3c(i)*fd13c(i)+hsv4c(i)*fd14c(i))/dt
      hg2c(i)=hgs(i)*(hsv1c(i)*fd21c(i)+hsv2c(i)*fd22c(i)
     2             + hsv3c(i)*fd23c(i)+hsv4c(i)*fd24c(i))/dt
!c     hourglass stresses Qi
      hgstress1c(i)=hgstress1c(i)+hg1c(i)*dt
      hgstress2c(i)=hgstress2c(i)+hg2c(i)*dt
  394 continue
      
!c	  do i=lft, llt
      fhg1(i)=hgstress1c(i)*hsv1c(i)*volgp(i)
      fhg2(i)=hgstress2c(i)*hsv1c(i)*volgp(i)
      fhg3(i)=hgstress1c(i)*hsv2c(i)*volgp(i)
      fhg4(i)=hgstress2c(i)*hsv2c(i)*volgp(i)
      fhg5(i)=hgstress1c(i)*hsv3c(i)*volgp(i) 
      fhg6(i)=hgstress2c(i)*hsv3c(i)*volgp(i)
      fhg7(i)=hgstress1c(i)*hsv4c(i)*volgp(i)
      fhg8(i)=hgstress2c(i)*hsv4c(i)*volgp(i)
!c      end do
!c      fhg1(i)=fhg1(i)+hg1c(i)*hsv1c(i)*volgp(i)
!c      fhg2(i)=fhg2(i)+hg2c(i)*hsv1c(i)*volgp(i)
!c      fhg3(i)=fhg3(i)+hg1c(i)*hsv2c(i)*volgp(i)
!c      fhg4(i)=fhg4(i)+hg2c(i)*hsv2c(i)*volgp(i)
!c      fhg5(i)=fhg5(i)+hg1c(i)*hsv3c(i)*volgp(i) 
!c      fhg6(i)=fhg6(i)+hg2c(i)*hsv3c(i)*volgp(i)
!c      fhg7(i)=fhg7(i)+hg1c(i)*hsv4c(i)*volgp(i)
!c      fhg8(i)=fhg8(i)+hg2c(i)*hsv4c(i)*volgp(i)
!c
!c     Contribution of hourglass stiffness to nodal forces
!$OMP PARALLEL DO       
      do 395 i=lft,llt
      r1(i)=r1(i)+r11(i)*(fhg1(i))+r12(i)*(fhg2(i))
      r2(i)=r2(i)+r21(i)*(fhg1(i))+r22(i)*(fhg2(i))
      r3(i)=r3(i)+r11(i)*(fhg3(i))+r12(i)*(fhg4(i))
      r4(i)=r4(i)+r21(i)*(fhg3(i))+r22(i)*(fhg4(i))
      r5(i)=r5(i)+r11(i)*(fhg5(i))+r12(i)*(fhg6(i))
      r6(i)=r6(i)+r21(i)*(fhg5(i))+r22(i)*(fhg6(i))
      r7(i)=r7(i)+r11(i)*(fhg7(i))+r12(i)*(fhg8(i))
      r8(i)=r8(i)+r21(i)*(fhg7(i))+r22(i)*(fhg8(i))
      inertia(i)=0.
      hgener(i)=hgstress1c(i)*hg1c(i)/hgs(i)*dt*volgp(i)+
     1    hgstress2c(i)*hg2c(i)/hgs(i)*dt*volgp(i)+hgener(i)
      totener(i)=intener(i)*volgp(i)*dt+totener(i)

!c      write(*,*) r1(i),r2(i),r3(i),r4(i)
!c      write(*,*) r5(i),r6(i),r7(i),r8(i),i,'after hg'
  395 continue

!c      write(*,*)'out of loop'

  396 if(iphase.ne.3.and.numdc.ne.0)goto 410
      if(iphase.eq.3)goto 430
  400 if (newstf) 430,410,430
  410 continue
!c
!c     igso controls geometric stiffness, set to 0
      if (igso.eq.1) call nonlin(s(9,1),ityp2d)
!c
!c     stffns is where Bbar method is implemented for 4 node
!c     integration
      if (ibar.eq.1) call stffns (s(9,1),d(1,1))
      if (ibar.eq.1) go to 430
!c     ityp2d is 0 for axi, 1 for pl. str., and 2 for pl. stress
      if (ityp2d.ne.0)
!c     calculates BTDB for one point integration or pl. stress
     1call pnbtcb (s(9,1),d(1,1),lst)
      if (ityp2d.eq.0)
!c     axisymmetric only
     1call axibtc (s(9,1),d(1,1))
      if (maxint.gt.1) go to 430
!c
!c     Hourglass stiffness control subroutines
!c     setcoe sets hgs as a function of the max s(i,j)
!c     hgstfb incorporates hourglass stiffness in stiffness matrix
      if (hgc.ne.0.00.and.ihgtyp.lt.2) call setcoe (hgs,s(9,1))
!c      if (ihgtyp.eq.3) then
!c      do 420 i=lft,llt
!c      ink=i+nftm1
!c      hgs(i)=hgs(i)*effmod(i)/prop(10,matp(ink))

!c 420  continue
!c      endif
      if (ihgtyp.eq.0) call hgstd2 (hgs,s(9,1))
      if (ihgtyp.eq.1) call hgstfb (hgs,s(9,1))
      if (ihgtyp.eq.2) call hgstfb2 (hgs,s(9,1))
      if (ihgtyp.eq.3) call hgstfb3 (hgs,s(9,1))

!c
  430 continue
!c     update element internal force and stiffness for cracked element
      do i=lft,llt
        ink=i+nftm1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		if(ElemFractCode(ink)==2) then !-- already decayed
		if(EC_GetElemSplit(ink)>0) then !-- already decayed
          Aratio=EC_GetElemAreaRatio (ink)
!!!        if((ink==200).or.(ink==401))then
!!!!            write(*,*)ink,Aratio
!!        endif
		  r1(i)=Aratio*r1(i)
		  r2(i)=Aratio*r2(i)
		  r3(i)=Aratio*r3(i)
		  r4(i)=Aratio*r4(i)
		  r5(i)=Aratio*r5(i)
		  r6(i)=Aratio*r6(i)
		  r7(i)=Aratio*r7(i)
		  r8(i)=Aratio*r8(i) 
		  do j=1,36
		    s(8+j,i)=Aratio*s(8+j,i)
	      end do
		end if
	  end do
      return
      end




