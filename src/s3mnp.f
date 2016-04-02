      subroutine s3mnp(prop,propc,sig,epx,ener,thick,ln)
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
     > deps(nelemg), p(nelemg),davg(nelemg),aj1(nelemg),t1(nelemg),
     > t2(nelemg),t3(nelemg),t4(nelemg),
     3 da1(nelemg),da2(nelemg),da3(nelemg),da4(nelemg),qhs(nelemg),
     > depn(nelemg),
     4 scale1(nelemg),scale0(nelemg)
      common/vect15/
     1 d1(nelemg),d2(nelemg),d3(nelemg),d4(nelemg)
      common/vect16/
     1 t456(nelemg),fc1qhs(nelemg)
c
      dimension prop(*),propc(4,*),sig(ln,*),epx(ln,*),ener(ln,*),
     1 d(4,4),thick(ln,*)
      equivalence (c11,d)
      data third/-.333333333333333/
c
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
      q1=ym*pr/((1.0+pr)*(1.0-2.0*pr))
      q2=.5*g
      q3=q1+2.0*q2
      q3=1./q3
c
      call rotat1 (sig,ener,ln)
c
c     compute trial stress
c
      do 20 i=mft,mlt
      qhs(i)=qh
      da4(i)=sig(4,i)+gd2*d4(i)
      d3(i) =-(sig(3,i)+q1*(d1(i)+d2(i)))*q3
      t4(i) =da4(i)-epx(4,i)
   20 continue
c
      if (prop(7).eq.0.0) then
      do 30 i=mft,mlt
   30 ak(i)=qs+qbqh*epx(3,i)
      else
c
c     variable hardening moduli
c
      call vyield (prop(6),prop(14),epx(3,1),ln,ak,qhs)
c
      qb=1.
      qa=0.
      endif
c
c
      fac1=1.5*g
      do 44 i=mft,mlt
      t456(i)=3.*t4(i)*t4(i)
      davg(i)=third*(d1(i)+d2(i)+d3(i))
      p(i)=blk*davg(i)
      da1(i)=sig(1,i)+p(i)+g*(d1(i)+davg(i))
      da2(i)=sig(2,i)+p(i)+g*(d2(i)+davg(i))
      da3(i)=sig(3,i)+p(i)+g*(d3(i)+davg(i))
      p(i) =third*(da1(i)+da2(i)+da3(i))
      t1(i)=p(i)+da1(i)-epx(1,i)
      t2(i)=p(i)+da2(i)-epx(2,i)
      t3(i)=p(i)+da3(i)+epx(1,i)+epx(2,i)
      aj2(i)=1.5*(t1(i)**2+t2(i)**2+t3(i)**2)+t456(i)
      ak2(i)=aj2(i)-ak(i)*ak(i)
      scle(i)=.50+sign(.5*unit,ak2(i))
      aj1(i) =0.0
      depi(i)=0.0
      deps(i)=0.0
      fc1qhs(i)=1./(fac1+qhs(i))
   44 continue
      scl=0.
      do 45 i=mft,mlt
   45 scl=scl+scle(i)
      if (nint(scl).eq.0.) then
      do 46 i=mft,mlt
      sig(1,i)=da1(i)
      sig(2,i)=da2(i)
      sig(3,i)=0.
      sig(4,i)=da4(i)
      thick(1,i)=thick(1,i)+d3(i)*thick(1,i)
   46 continue
      go to 75
      endif
      do 60 i=mft,mlt
      if (scle(i).eq.0.)  go to 60
      sg3old=da3(i)  -deps(i)*t3(i)
      sg3lst=sg3new
      sg3new=sg3old
      epslst=epsnew
      epsnew=d3(i)
      d3(i)=-d1(i)-d2(i)
      do 50 iter=2,20
      davg(i)=third*(d1(i)+d2(i)+d3(i))
      p(i)=blk*davg(i)
      da1(i)=sig(1,i)+p(i)+g*(d1(i)+davg(i))
      da2(i)=sig(2,i)+p(i)+g*(d2(i)+davg(i))
      da3(i)=sig(3,i)+p(i)+g*(d3(i)+davg(i))
      p(i) =third*(da1(i)+da2(i)+da3(i))
      t1(i)=p(i)+da1(i)-epx(1,i)
      t2(i)=p(i)+da2(i)-epx(2,i)
      t3(i)=p(i)+da3(i)+epx(1,i)+epx(2,i)
      aj2(i)=1.5*(t1(i)**2+t2(i)**2+t3(i)**2)+t456(i)
      ak2(i)=aj2(i)-ak(i)*ak(i)
      scle(i)=.50+sign(.5*unit,ak2(i))
      aj1(i) =sqrt(aj2(i))+1.0-scle(i)
      depi(i)=scle(i)*(aj1(i)-ak(i))*fc1qhs(i)
      deps(i)=fac1*depi(i)/aj1(i)
      sg3old=da3(i)  -deps(i)*t3(i)
      sg3lst=sg3new
      sg3new=sg3old
      epslst=epsnew
      epsnew=d3(i)
      demn=1.e-20+sg3new-sg3lst
      if (abs(demn).lt.1.e-14) go to 60
      d3(i)=epslst-sg3lst*(epsnew-epslst)/demn
      if(abs(epsnew-epslst)/abs(d3(i)).lt.00.00001) go to 60
   50 continue
      write (*,240)
   60 continue
      if (iphase.eq.3) go to 67
      do 66 i=mft,mlt
      if (scle(i).eq.0.) go to 66
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
      a=ax*(c1*c1+c2*c2+c3*c3+2.*c4**2)
      b=ax*(c1*t1(i)+c2*t2(i)+c3*t3(i)+2.*c4*t4(i))
      phi=ak2(i)/(ak(i)*ak(i))
      q11=b*b-a*phi
      if (q11.lt.0.0) go to 64
      scale=(-b+sqrt(q11))/a
   64 scale0(i)=min(-scale,1.*unit)
      scale1(i)=1.-scale0(i)
   66 continue
   67 fac2=qa/fac1
      do 70 i=mft,mlt
      depn(i)=fac2*qhs(i)*deps(i)
      sig(1,i)=da1(i) -deps(i)*t1(i)
      sig(2,i)=da2(i) -deps(i)*t2(i)
      sig(3,i)=0.0
      sig(4,i)=da4(i) -deps(i)*t4(i)
      epx(1,i)=epx(1,i)+depn(i)*t1(i)
      epx(2,i)=epx(2,i)+depn(i)*t2(i)
      epx(3,i)=epx(3,i)+depi(i)
      epx(4,i)=epx(4,i)+depn(i)*t4(i)
      thick(1,i)=thick(1,i)+d3(i)*thick(1,i)
   70 continue
   75 continue
      if (prop(7).ne.0.0) then
      call vyield (prop(6),prop(14),epx(3,1),ln,ak,qhs)
      do 130 i=mft,mlt
      t456(i)=3.*sig(4,i)*sig(4,i)
      p(i) =third*(sig(1,i)+sig(2,i)+sig(3,i))
      t1(i)=p(i)+sig(1,i)
      t2(i)=p(i)+sig(2,i)
      t3(i)=p(i)+sig(3,i)
      aj2(i)=1.5*(t1(i)**2+t2(i)**2+t3(i)**2)+t456(i)
      ak2(i)=aj2(i)-ak(i)*ak(i)
      scle(i)=.50+sign(.5*unit,ak2(i))
      deps(i)=1.-scle(i)+
     > scle(i)*sqrt((ak(i)*ak(i))/(aj2(i)+1.-scle(i)))
      sig(1,i)=deps(i)*sig(1,i)
      sig(2,i)=deps(i)*sig(2,i)
      sig(4,i)=deps(i)*sig(4,i)
  130 continue
      endif
c
      if (iphase.eq.3) go to 200
c
c     compute tangent constitutive matrix
c
      do 140 i=1,4
      do 140 j=1,4
  140 d(i,j)=propc(i,j)
      p1=-c31/c33
      p2=-c32/c33
      p4=-c34/c33
      gp11 = c11 + p1*c13
      gp12 = c12 + p2*c13
      gp13 = p4*c13 + c14
      gp21 = c21 + p1*c23
      gp22 = c22 + p2*c23
      gp23 = p4*c23 + c24
      gp31 = c31 + p1*c33
      gp32 = c32 + p2*c33
      gp33 = p4*c33 + c34
      gp41 = c41 + p1*c43
      gp42 = c42 + p2*c43
      gp43 = p4*c43 + c44
      c11  = gp11 + p1*gp31
      c21  = gp21 + p2*gp31
      c31  = 0.
      c41  = p4*gp31 + gp41
      c12  = gp12 + p1*gp32
      c22  = gp22 + p2*gp32
      c32  = 0.
      c42  = p4*gp32 + gp42
      c13  = 0.
      c23  = 0.
      c33  = 0.
      c43  = 0.
      c14  = gp13 + p1*gp33
      c24  = gp23 + p2*gp33
      c34  = 0.
      c44  = p4*gp33 + gp43
c
      do 150 i=mft,mlt
      ak(i)=qs+qbqh*epx(3,i)
      dsave(1,1,i)=thick(1,i)*c11
      dsave(1,2,i)=thick(1,i)*c12
      dsave(1,3,i)=0.
      dsave(1,4,i)=thick(1,i)*c14
      dsave(2,1,i)=thick(1,i)*c21
      dsave(2,2,i)=thick(1,i)*c22
      dsave(2,3,i)=0.
      dsave(2,4,i)=thick(1,i)*c24
      dsave(3,1,i)=0.
      dsave(3,2,i)=0.
      dsave(3,3,i)=0.
      dsave(3,4,i)=0.
      dsave(4,1,i)=thick(1,i)*c41
      dsave(4,2,i)=thick(1,i)*c42
      dsave(4,3,i)=0.
  150 dsave(4,4,i)=thick(1,i)*c44
c
      do 155 i=1,4
      do 155 j=1,4
  155 d(i,j)=propc(i,j)
c
      do 190 i=mft,mlt
      if (scle(i).eq.0.) go to 190
      if (prop(7).eq.0.) go to 180
      qh=max(qhs(i),.01*g)
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
      p1=-dsave(3,1,i)/dsave(3,3,i)
      p2=-dsave(3,2,i)/dsave(3,3,i)
      p4=-dsave(3,4,i)/dsave(3,3,i)
      gp11 = dsave(1,1,i) + p1*dsave(1,3,i)
      gp12 = dsave(1,2,i) + p2*dsave(1,3,i)
      gp13 = dsave(1,3,i)*p4 + dsave(1,4,i)
      gp21 = dsave(2,1,i) + p1*dsave(2,3,i)
      gp22 = dsave(2,2,i) + p2*dsave(2,3,i)
      gp23 = dsave(2,3,i)*p4 + dsave(2,4,i)
      gp31 = dsave(3,1,i) + p1*dsave(3,3,i)
      gp32 = dsave(3,2,i) + p2*dsave(3,3,i)
      gp33 = dsave(3,3,i)*p4 + dsave(3,4,i)
      gp41 = dsave(4,1,i) + p1*dsave(4,3,i)
      gp42 = dsave(4,2,i) + p2*dsave(4,3,i)
      gp43 = dsave(4,3,i)*p4 + dsave(4,4,i)
      dsave(1,1,i)=thick(1,i)*(gp11 + p1*gp31)
      dsave(2,1,i)=thick(1,i)*(gp21 + p2*gp31)
      dsave(3,1,i)=0.
      dsave(4,1,i)=thick(1,i)*(p4*gp31 + gp41)
      dsave(1,2,i)=thick(1,i)*(gp12 + p1*gp32)
      dsave(2,2,i)=thick(1,i)*(gp22 + p2*gp32)
      dsave(3,2,i)=0.
      dsave(4,2,i)=thick(1,i)*(p4*gp32 + gp42)
      dsave(1,3,i)=0.
      dsave(2,3,i)=0.
      dsave(3,3,i)=0.
      dsave(4,3,i)=0.
      dsave(1,4,i)=thick(1,i)*(gp13 + p1*gp33)
      dsave(2,4,i)=thick(1,i)*(gp23 + p2*gp33)
      dsave(3,4,i)=0.
      dsave(4,4,i)=thick(1,i)*(p4*gp33 + gp43)
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
  240 format(' convergence failure in plane stress plasticity subr.')
      end
