      subroutine s3mn (prop,propc,sig,epx,ener,ln)
c  this is a J@ Plasticity based on a radial return solution
c     implicit double precision (a-h,o-z)                                    dp
c
       use mod_parameters
      common/bk02/ioofc,iphase,imass,lpar(9)
      common/intgp/
     1 c11,c21,c31,c41,c12,c22,c32,c42,c13,c23,c33,c43,
     > c14,c24,c34,c44
      common/range/mft,mlt,lft,llt,nftm1
      common/double/iprec,ncpw,unit
      common/vect3/
     1 dgi1(nelemg,4),dgi2(nelemg,4),dgi3(nelemg,4),dgi4(nelemg,4),
     > dgi5(nelemg,4),
     2 dgt1(nelemg,4),dgt2(nelemg,4),dgt3(nelemg,4),dgt4(nelemg,4),
     3 f11v(nelemg),f22v(nelemg),f12v(nelemg),f21v(nelemg),
     > dsd5(nelemg),
     4 sig11s(nelemg),sig22s(nelemg),sig33s(nelemg),sig12s(nelemg),
     5 ddp1(nelemg,4),ddp2(nelemg,4),ddp3(nelemg,4),ddp4(nelemg,4),
     > ddp5(nelemg,4)
      common/vect8/
     1 dsave(4,4,nelemg)
      common/vect14/
     1 aj2(nelemg),ak2(nelemg),ak(nelemg),scle(nelemg),depi(nelemg),
     > deps(nelemg),
     2 p(nelemg),davg(nelemg),aj1(nelemg),t1(nelemg),t2(nelemg),
     > t3(nelemg),t4(nelemg),
     3 da1(nelemg),da2(nelemg),da3(nelemg),da4(nelemg),qhs(nelemg),
     > depn(nelemg),
     4 scale1(nelemg),scale0(nelemg)
      common/vect15/
     1 d1(nelemg),d2(nelemg),d3(nelemg),d4(nelemg)
c     common/plmvc1/plwk(128)                                                pl
c
      dimension prop(*),propc(4,*),sig(ln,*),epx(ln,*),
     > ener(ln,*),d(4,4)
c    1 ,seffo(128)                                                           pl
      equivalence (c11,d)
      data third/-.333333333333333/
c
c     properties defined in MAZE
      qb  =prop(5)
      qh  =prop(1)*prop(4)/(prop(1)-prop(4))
      qs  =prop(3)
      ym  =prop(1)
      pr  =prop(2)
      g   =ym/(1.+pr)
      qa  =1.-qb
      gd2 =.5*g
      blk =-ym/(1.-2.*pr)
      qbqh=qb*qh
c
      call rotat1 (sig,ener,ln)
c
c     compute trial stress
c
c     do 45 i=mft,mlt                                                        pl
c     seffo(i)=1.5*(sig(1,i)**2+sig(2,i)**2+sig(3,i)**2+2.*sig(4,i)**2)      pl
c    1       -.5*(sig(1,i)+sig(2,i)+sig(3,i))**2                             pl
c  45 continue                                                               pl
      do 10 i=mft,mlt
      qhs(i) =qh
      davg(i)=third*(d1(i)+d2(i)+d3(i))
   10 p(i)=blk*davg(i)
      do 20 i=mft,mlt
      da1(i)=sig(1,i)+p(i)+g*(d1(i)+davg(i))
      da2(i)=sig(2,i)+p(i)+g*(d2(i)+davg(i))
      da3(i)=sig(3,i)+p(i)+g*(d3(i)+davg(i))
   20 da4(i)=sig(4,i)+gd2*d4(i)
c
      if (prop(7).eq.0.0) go to 30
c
c     variable hardening moduli
c
      call vyield (prop(6),prop(14),epx(3,1),ln,ak,qhs)
c
      qb=1.
      qa=0.
      go to 50
c
   30 do 40 i=mft,mlt
   40 ak(i)=qs+qbqh*epx(3,i)
c
   50 do 60 i=mft,mlt
      p(i) =third*(da1(i)+da2(i)+da3(i))
      t1(i)=p(i)+da1(i)-epx(1,i)
      t2(i)=p(i)+da2(i)-epx(2,i)
      t3(i)=p(i)+da3(i)+epx(1,i)+epx(2,i)
   60 t4(i)=da4(i)-epx(4,i)
      do 70 i=mft,mlt
      aj2(i)=1.5*(t1(i)**2+t2(i)**2+t3(i)**2)+3.*t4(i)**2
   70 ak2(i)=aj2(i)-ak(i)*ak(i)
      do 80 i=mft,mlt
   80 scle(i)=.50+sign(.5*unit,ak2(i))
      if (iphase.eq.3) go to 110
      do 100 i=mft,mlt
      if (scle(i).eq.0.) go to 100
      scale=-1.0
      b1=da1(i)-sig(1,i)
      b2=da2(i)-sig(2,i)
      b3=da3(i)-sig(3,i)
      c4=da4(i)-sig(4,i)
      av=third*(b1+b2+b3)
      c1=b1+av
      c2=b2+av
      c3=b3+av
      ax=1.5/(ak(i)*ak(i))
      a=ax*(c1*c1+c2*c2+c3*c3+2.*c4*c4)
      b=ax*(c1*t1(i)+c2*t2(i)+c3*t3(i)+2.*c4*t4(i))
      phi=ak2(i)/(ak(i)*ak(i))
      q11=b*b-a*phi
      if (min(q11,a)-1.e-20.lt.0.0) go to 90
      scale=(-b+sqrt(q11))/a
   90 scale0(i)=min(-scale,1.*unit)
      scale1(i)=1.-scale0(i)
  100 continue
  110 fac1=1.5*g
      fac2=qa*qh/fac1
      do 120 i=mft,mlt
  120 aj1(i)=sqrt(aj2(i))+1.0-scle(i)
      do 130 i=mft,mlt
      depi(i) =scle(i)*(aj1(i)-ak(i))/(fac1+qhs(i))
      epx(3,i)=epx(3,i)+depi(i)
      deps(i) =fac1*depi(i)/aj1(i)
  130 depn(i) =fac2*deps(i)
      do 140 i=mft,mlt
      sig(1,i)=da1(i)  -deps(i)*t1(i)
      sig(2,i)=da2(i)  -deps(i)*t2(i)
      sig(3,i)=da3(i)  -deps(i)*t3(i)
      sig(4,i)=da4(i)  -deps(i)*t4(i)
      epx(1,i)=epx(1,i)+depn(i)*t1(i)
      epx(2,i)=epx(2,i)+depn(i)*t2(i)
  140 epx(4,i)=epx(4,i)+depn(i)*t4(i)
c     do 145 i=mft,mlt                                                       pl
c     seff=1.5*(sig(1,i)**2+sig(2,i)**2+sig(3,i)**2+2.*sig(4,i)**2)          pl
c    1       -.5*(sig(1,i)+sig(2,i)+sig(3,i))**2                             pl
c     plwk(i)=plwk(i)+.125*(sqrt(seff)+sqrt(seffo(i)))*depi(i)               pl
c 145 continue                                                               pl
c
      if (prop(7).eq.0.0.and.iphase.eq.3) go to 200
c
c     compute tangent constitutive matrix
c
      do 150 i=1,4
      do 150 j=1,4
  150 d(i,j)=propc(i,j)
c
      do 160 i=mft,mlt
      ak(i)=qs+qbqh*epx(3,i)
      dsave(1,1,i)=c11
      dsave(1,2,i)=c12
      dsave(1,3,i)=c13
      dsave(1,4,i)=c14
      dsave(2,1,i)=c21
      dsave(2,2,i)=c22
      dsave(2,3,i)=c23
      dsave(2,4,i)=c24
      dsave(3,1,i)=c31
      dsave(3,2,i)=c32
      dsave(3,3,i)=c33
      dsave(3,4,i)=c34
      dsave(4,1,i)=c41
      dsave(4,2,i)=c42
      dsave(4,3,i)=c43
  160 dsave(4,4,i)=c44
c
      do 190 i=mft,mlt
      if (scle(i).eq.0.) go to 190
      if (prop(7).eq.0.) go to 180
      call yield (prop(6),prop(14),epx(3,i),ak(i),qh)
      pres=third*(sig(1,i)+sig(2,i)+sig(3,i))
      q1=sig(1,i)+pres
      q2=sig(2,i)+pres
      q3=sig(3,i)+pres
      q4=sig(4,i)
      phi=1.5*(q1*q1+q2*q2+q3*q3+2.0*q4*q4)/(ak(i)*ak(i))-1.0
      if (phi.le.0.0) go to 170
      scale=1.0/sqrt(1.0+phi)
      sig(1,i)=scale*q1-pres
      sig(2,i)=scale*q2-pres
      sig(3,i)=scale*q3-pres
      sig(4,i)=scale*q4
  170 if (iphase.eq.3) go to 190
      qh=max(qh,.01*g)
  180 av=third*(sig(1,i)+sig(2,i)+sig(3,i))
      p1=sig(1,i)+av-epx(1,i)
      p2=sig(2,i)+av-epx(2,i)
      p3=sig(3,i)+av+epx(1,i)+epx(2,i)
      p4=sig(4,i)-epx(4,i)
      call rott2s(p1,p2,p4,i)
      p4=2.*p4
      q1=c11*p1+c12*(p2+p3)
      q2=c22*p2+c12*(p1+p3)
      q3=c33*p3+c12*(p1+p2)
      q4=c44*p4
      q11=p1*q1+p2*q2+p3*q3+p4*q4
      q11=q11+ak(i)*ak(i)*qh/2.25
      q11=1.0/q11
      q1q11=q1*q11
      q2q11=q2*q11
      q3q11=q3*q11
      d11=c11-q1q11*q1
      d12=c12-q1q11*q2
      d13=c13-q1q11*q3
      d14=   -q1q11*q4
      d22=c22-q2q11*q2
      d23=c23-q2q11*q3
      d24=   -q2q11*q4
      d33=c33-q3q11*q3
      d34=   -q3q11*q4
      d44=c44-q4*q11*q4
      scal1=c11*scale1(i)
      scal2=c12*scale1(i)
      dsave(1,1,i)=d11*scale0(i)+scal1
      dsave(1,2,i)=d12*scale0(i)+scal2
      dsave(1,3,i)=d13*scale0(i)+scal2
      dsave(1,4,i)=d14*scale0(i)
      dsave(2,1,i)=dsave(1,2,i)
      dsave(2,2,i)=d22*scale0(i)+scal1
      dsave(2,3,i)=d23*scale0(i)+scal2
      dsave(2,4,i)=d24*scale0(i)
      dsave(3,1,i)=dsave(1,3,i)
      dsave(3,2,i)=dsave(2,3,i)
      dsave(3,3,i)=d33*scale0(i)+scal1
      dsave(3,4,i)=d34*scale0(i)
      dsave(4,1,i)=dsave(1,4,i)
      dsave(4,2,i)=dsave(2,4,i)
      dsave(4,3,i)=dsave(3,4,i)
      dsave(4,4,i)=d44*scale0(i)+c44*scale1(i)
c
  190 continue
c
  200 call rotat2 (sig,ener,ln)
c
      do 210 i=mft,mlt
      sig11s(i)=sig(1,i)
      sig22s(i)=sig(2,i)
      sig33s(i)=sig(3,i)
  210 sig12s(i)=sig(4,i)
c
      return
      end
