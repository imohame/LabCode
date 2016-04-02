      subroutine gpsbw (nr,nd,lv,ls,lu,kc,nc,ka)
c     implicit double precision (a-h,o-z)                                    dp
c
c***********************************************************************
c     this routine uses the gibbs-poole-stockmeyer algorithm to        *
c     determine a nodal numbering scheme which results in minimal      *
c     bandwidth/profile                                                *
c     (see siam j. num. anal., vol.13, no.2, april 1976, pp236-250)    *
c***********************************************************************
      dimension nr(*),nd(2,*),lv(*), ls(*), lu(*), kc(*), nc(*), ka(*)
      common /sabr1/ mind, maxd, ibw, nbw, ipr, npr, nnc, nzn, nnp
      common /sabr2/ kbw, kpr
      common /gps1/ lw, ld
      common /gps2/ inbr, lnbr
c.....part i: finding endpoints of a pseudo-diameter.
      iv=0
   10 iv=iv+1
      if (nd(1,iv).ne.mind) go to 10
      call lvstr (iv,nd,lv,ka)
      do 20 j=1,nnp
      if (lv(j).eq.0) go to 80
   20 continue
   30 i=0
      do 50 j=1,nnp
      if (lv(j).ne.ld) go to 50
      k=i
      i=i+1
      nr(i)=j
      jd=nd(1,j)
   40 if (k.eq.0) go to 50
      m=nr(k)
      if (nd(1,m).le.jd) go to 50
      nr(k+1)=m
      nr(k)=j
      k=k-1
      go to 40
   50 continue
      lvw=lw
      lvd=ld
      luw=nnp
      lud=0
      do 70 j=1,i
      is=nr(j)
      call lvstr (is,nd,ls,ka)
      if (ld.le.lvd) go to 60
      iv=is
c     call blkcpy(ls,lv,nnp)                                            cray1
      call blicpy(ls,lv,nnp)                                            vax750
      go to 30
   60 if (luw.le.lw) go to 70
      iu=is
      luw=lw
      lud=ld
c     call blkcpy(ls,lu,nnp)                                            cray1
      call blicpy(ls,lu,nnp)                                            vax750
   70 continue
      if (lud.eq.lvd) go to 90
   80 write (*,400)
      call bye (2)
c.....part ii: minimizing level width.
   90 continue
      do 100 j=1,nnp
      lu(j)=lvd+1-lu(j)
      ls(j)=0
  100 continue
      ii=0
      do 110 j=1,nnp
      if (lv(j).ne.lu(j)) go to 110
      ls(j)=-1
      ii=ii+1
  110 continue
      nif=0
      if (ii.eq.nnp) go to 290
      ii=0
      kk=0
      do 120 j=1,nnp
      if (ls(j).ne.0) go to 120
      kk=kk+1
      nr(kk)=ii
      call disco (j,nd,ls,ka,kc(ii+1),kk)
      nc(kk)=lw
      ii=ii+lw
  120 continue
      do 130 k=2,kk
      nk=nc(k-1)
      do 130 l=k,kk
      nl=nc(l)
      if (nl.le.nk) go to 130
      nc(k-1)=nl
      nc(l)=nk
      nk=nl
      nn=nr(k-1)
      nr(k-1)=nr(l)
      nr(l)=nn
  130 continue
      i=0
      do 140 k=1,kk
      n=nc(k)
c     call blkcpy (kc(nr(k)+1),ls(i+1),n)                               cray1
      call blicpy (kc(nr(k)+1),ls(i+1),n)                               vax750
      i=i+n
      ls(i)=-ls(i)
  140 continue
      do 150 j=1,lvd
      nr(j)=0
  150 continue
      do 160 j=1,nnp
      if (lv(j).ne.lu(j)) go to 160
      i=lv(j)
      nr(i)=nr(i)+1
  160 continue
      do 170 i=1,ii
      nc(i)=0
      kc(i)=0
  170 continue
      m=1
      do 280 i=1,ii
      j=ls(i)
      k=abs(j)
      l=lv(k)
      kc(l)=kc(l)+1
      l=lu(k)
      nc(l)=nc(l)+1
      if (j.gt.0) go to 280
      ls(i)=k
      n1=0
      n2=0
      do 190 l=1,lvd
      if (kc(l).eq.0) go to 180
      n1=max(n1,nr(l)+kc(l))
  180 if (nc(l).eq.0) go to 190
      n2=max(n2,nr(l)+nc(l))
  190 continue
      if (n1-n2) 210,200,230
  200 if (luw.lt.lvw) go to 230
  210 continue
      do 220 l=1,lvd
      nr(l)=nr(l)+kc(l)
  220 continue
      go to 260
  230 continue
      if (m.eq.1) nif=m
      do 240 l=1,lvd
      nr(l)=nr(l)+nc(l)
  240 continue
      do 250 l=m,i
      j=ls(l)
      lv(j)=lu(j)
  250 continue
  260 m=i+1
      do 270 l=1,lvd
      kc(l)=0
      nc(l)=0
  270 continue
  280 continue
c.....part iii: numbering
  290 if (nd(1,iv).le.nd(1,iu)) go to 310
      nif=nif+1
      is=iv
      iv=iu
      iu=is
      do 300 j=1,nnp
      lv(j)=lvd+1-lv(j)
  300 continue
  310 nif=nif-1
      i0=0
      i1=0
      inbr=1
      lnbr=1
      nr(inbr)=iv
      lv(iv)=0
      go to 330
  320 i2=inbr
      lnbr=lnbr+1
      call numrl (i1,i2,nr,nd,lv,ka)
      i0=i2
      i1=i2
  330 continue
      call numrl (i0,nnp,nr,nd,lv,ka)
      id=maxd+1
      is=0
      do 340 j=1,nnp
      if (lv(j).ne.lnbr) go to 340
      jd=nd(1,j)
      if (jd.ge.id) go to 340
      id=jd
      is=j
  340 continue
      if (is.eq.0) go to 350
      i0=inbr
      inbr=inbr+1
      nr(inbr)=is
      lv(is)=0
      go to 330
  350 if (lnbr.lt.lvd) go to 320
      if (nif.eq.0) go to 370
      nh=nnp/2
      do 360 i=1,nh
      nrt=nr(i)
      nr(i)=nr(nnp-i+1)
      nr(nnp-i+1)=nrt
  360 continue
  370 continue
      call comb1 (nd,nr,lv,ka)
      if (kpr.ge.npr) go to 390
      npr=kpr
      nbw=kbw
      if (nzn.eq.0) go to 390
      do 380 j=1,nzn
      do 380 i=1,nnp
      if (nr(i).lt.nr(j-nzn)) go to 380
      nr(i)=nr(i)+1
  380 continue
  390 return
  400 format (/' error termination in gps --- model has disjoint parts')
      end
