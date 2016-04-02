      subroutine eqnum(id,y,z,b,matp,a,numnp,numel,neq,iband,ncn,nrcc,
     1 roller,idctrl,idirw)
c     implicit double precision (a-h,o-z)                                    dp
c

      common/WMLthermal/thermalflag   !!!!,thermalconstraint(40000),Tinit(40000),Rqold(40000)
      integer thermalflag,numnp
      common/WMLthermalreordr/idth(40000)
      integer idth,neq2
      common/bk14/lfna(15),lfnt(6)
      dimension y(*),z(*),id(2,*),b(*),matp(*),a(*),ncn(3,*),roller(*)
c
      call header
      nprnt=0
      neq=0
      ier=0
      if (iband.eq.0) go to 10  ! bandwidth minimization flag
c
      call reordr (a,b,matp,id,numnp,neq,numel,ncn,nrcc,ier,iband)
c
      write(lfnt(2),130)
      nprnt=87
      if (ier.eq.0) go to 80
c
   10 do 70 n=1,numnp
      il=1
      iu=2
      if (nrcc.eq.0) go to 30
      do 20 j=1,nrcc
      if (n.eq.ncn(2,j)) then
      if (ncn(3,j).eq.3) go to 70
      if (ncn(3,j).eq.1) il=2
      if (ncn(3,j).eq.2) iu=1
      endif
   20 continue
   30 do 60 i=il,iu
      if (id(i,n)) 50,40,50
   40 neq=neq+1
      id(i,n)=neq
      go to 60
   50 id(i,n)=0
   60 continue
   70 continue
   
      if (thermalflag > 0) then
      neq2=0
   71 do 78 n=1,numnp
      il=1
      iu=1
      do 77 i=il,iu
      if (idth(n)) 73,72,73
   72 neq2=neq2+1
      idth(n)=neq2
      go to 77
   73 idth(n)=0
   77 continue
   78 continue
      endif
c
c------------------------------------------
      write(*,79) neq
   79 format(5x,'Number of equations =',i3)
c------------------------------------------      
c
   80 if (nrcc.eq.0) go to 100
      do 90 j=1,nrcc
      n=ncn(1,j)
      m=ncn(2,j)
      if (ncn(3,j).eq.3) then
      id(1,m)=id(1,n)
      id(2,m)=id(2,n)
      elseif (ncn(3,j).eq.2) then
      id(2,m)=id(2,n)
      elseif (ncn(3,j).eq.1) then
      id(1,m)=id(1,n)
      endif
   90 continue
c
  100 continue
      do 120 n=1,numnp
      bcc=0.0
      if (id(1,n).eq.0.and.id(2,n).ne.0) bcc=1.0
      if (id(1,n).ne.0.and.id(2,n).eq.0) bcc=2.0
      if (id(1,n).eq.0.and.id(2,n).eq.0) bcc=3.0
      if (nprnt.gt.0) go to 110
      nprnt=50
      call header
      write(lfnt(2),130)
  110 nprnt=nprnt-1
      if(roller(n).ne.0.) bcc=roller(n)
      write(lfnt(2),140) n,bcc,y(n),z(n),id(1,n),id(2,n)
  120 continue
c
      if (idctrl.ne.0) idctrl=id(idirw,idctrl)
      return
c
c
  130 format(///' g e n e r a t e d  n o d a l  d a t a '/
     1/' node',5x,' b.c.',10x,'y',20x,'z',20x,'eq-y',5x,'eq-z')
  140 format(i5,5x,f5.0,5x,e12.4,8x,e12.4,13x,i5,5x,i5)
      end
