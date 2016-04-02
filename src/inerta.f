      subroutine inerta(xm,udt,udtt,lm,s,iequit,iprec)
c     implicit double precision (a-h,o-z)                                    dp
c
c     modify load vector and stiffness matrix to account for
c     inertial effects
c
       use mod_parameters
      common/bk02/ioofc,iphase,imass,lpar(9)
      common/bk13/xnorm0(6),xnormc(6)!,xnk2d(20)
      common/bk26/dt,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
      common/cn4/itypei,ianali,ibwmni,icorei,ipcori,isolti,gamai,
     1           betai,raydmi
      common/range/mft,mlt,lft,llt,nftm1
      
      common/vect1/r1(nelemg),r2(nelemg),r3(nelemg),
     > r4(nelemg),r5(nelemg),r6(nelemg),
     > r7(nelemg),r8(nelemg),mtype(nelemg),mte(nelemg)
      
      common/vect2/
     1 ed1(nelemg),ed2(nelemg),ed3(nelemg),ed4(nelemg),
     2 ed5(nelemg),ed6(nelemg),ed7(nelemg),ed8(nelemg),
     3 fd1(nelemg),fd2(nelemg),fd3(nelemg),fd4(nelemg),
     4 fd5(nelemg),fd6(nelemg),fd7(nelemg),fd8(nelemg)
      
      common/vect5/
     1 wv1(nelemg),wv2(nelemg),wv3(nelemg),wv4(nelemg),
     2 wv5(nelemg),wv6(nelemg),wv7(nelemg),wv8(nelemg),
     3 vv1(nelemg),vv2(nelemg),vv3(nelemg),vv4(nelemg),
     4 vv5(nelemg),vv6(nelemg),vv7(nelemg),vv8(nelemg),
     5 av1(nelemg),av2(nelemg),av3(nelemg),av4(nelemg),
     6 av5(nelemg),av6(nelemg),av7(nelemg),av8(nelemg)
      common/vect9/scl(8,nelemg),yz(8,nelemg)
      dimension xm(4,*),udt(*),udtt(*),lm(44*iprec,*),s(44,*)
      if (imass.ne.1) return
      b0=0.
      if(iequit.eq.1) b0=a0
      do 10 i=lft,llt
      vv1(i)  = udt(lm(1,i))
      vv2(i)  = udt(lm(2,i))
      vv3(i)  = udt(lm(3,i))
      vv4(i)  = udt(lm(4,i))
      vv5(i)  = udt(lm(5,i))
      vv6(i)  = udt(lm(6,i))
      vv7(i)  = udt(lm(7,i))
      vv8(i)  = udt(lm(8,i))
      av1(i)  =udtt(lm(1,i))
      av2(i)  =udtt(lm(2,i))
      av3(i)  =udtt(lm(3,i))
      av4(i)  =udtt(lm(4,i))
      av5(i)  =udtt(lm(5,i))
      av6(i)  =udtt(lm(6,i))
      av7(i)  =udtt(lm(7,i))
      av8(i)  =udtt(lm(8,i))
   10 continue
      do 20 i=lft,llt
      wv1(i)=scl(1,i)*(-a2*vv1(i)-a3*av1(i)+b0*fd1(i))
      wv2(i)=scl(2,i)*(-a2*vv2(i)-a3*av2(i)+b0*fd2(i))
      wv3(i)=scl(3,i)*(-a2*vv3(i)-a3*av3(i)+b0*fd3(i))
      wv4(i)=scl(4,i)*(-a2*vv4(i)-a3*av4(i)+b0*fd4(i))
      wv5(i)=scl(5,i)*(-a2*vv5(i)-a3*av5(i)+b0*fd5(i))
      wv6(i)=scl(6,i)*(-a2*vv6(i)-a3*av6(i)+b0*fd6(i))
      wv7(i)=scl(7,i)*(-a2*vv7(i)-a3*av7(i)+b0*fd7(i))
      wv8(i)=scl(8,i)*(-a2*vv8(i)-a3*av8(i)+b0*fd8(i))
   20 continue
      do 30 i=lft,llt
      r1(i)=r1(i)+xm(1,i)*wv1(i)
      r2(i)=r2(i)+xm(1,i)*wv2(i)
      r3(i)=r3(i)+xm(2,i)*wv3(i)
      r4(i)=r4(i)+xm(2,i)*wv4(i)
      r5(i)=r5(i)+xm(3,i)*wv5(i)
      r6(i)=r6(i)+xm(3,i)*wv6(i)
      r7(i)=r7(i)+xm(4,i)*wv7(i)
      r8(i)=r8(i)+xm(4,i)*wv8(i)
   30 continue
      if(raydmi.ne.0.and.iequit.eq.1)then
      fac1=raydmi/dt
      do 35 i=lft,llt
      r1(i)=r1(i)+fac1*xm(1,i)*fd1(i)
      r2(i)=r2(i)+fac1*xm(1,i)*fd2(i)
      r3(i)=r3(i)+fac1*xm(2,i)*fd3(i)
      r4(i)=r4(i)+fac1*xm(2,i)*fd4(i)
      r5(i)=r5(i)+fac1*xm(3,i)*fd5(i)
      r6(i)=r6(i)+fac1*xm(3,i)*fd6(i)
      r7(i)=r7(i)+fac1*xm(4,i)*fd7(i)
      r8(i)=r8(i)+fac1*xm(4,i)*fd8(i)
   35 continue
      endif
      fac2=a0+raydmi/dt
      do 40 i=lft,llt
      s(1,i) = s(1,i)+fac2*xm(1,i)
      s(3,i) = s(3,i)+fac2*xm(1,i)
      s(6,i) = s(6,i)+fac2*xm(2,i)
      s(10,i)=s(10,i)+fac2*xm(2,i)
      s(15,i)=s(15,i)+fac2*xm(3,i)
      s(21,i)=s(21,i)+fac2*xm(3,i)
      s(28,i)=s(28,i)+fac2*xm(4,i)
      s(36,i)=s(36,i)+fac2*xm(4,i)
   40 continue
!c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!      do 50 i=lft,llt
!      xnk2d(2)=xnk2d(2)+
!     > fd1(i)**2*xm(1,i)+
!     > fd2(i)**2*xm(1,i)+
!     > fd3(i)**2*xm(2,i)+
!     > fd4(i)**2*xm(2,i)+
!     > fd5(i)**2*xm(3,i)+
!     > fd6(i)**2*xm(3,i)+
!     > fd7(i)**2*xm(4,i)+
!     > fd8(i)**2*xm(4,i)
!   50 continue
!c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      return
      end
