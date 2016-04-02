      subroutine reordr(nrv,b,matp,id,numnp,neq,numel,ncn,nrcc,ier,
     1 iband)
c     implicit double precision (a-h,o-z)                                    dp
c
c     obtain nodal reorder vector and number equations
c
      dimension nrv(*),b(*),matp(*),id(2,*),ncn(3,*)

      common/main_block/ a(1)

      common/bk00/
     1k01,k02,k03,k04,k05,k06,k07,k08,k09,k10,k11,k12,
     2k13,k14,k15,k16,k17,k18,k19,k20,k21,k22,k23,k24,
     3k25,k26,k27,k28,k29,k30,k31,k32,k33,k34,k35,k36,
     4k37,k38,k39,k40,k41,k42,k43,k44,k45,k46,k47,k48,
     5k49,k50,k51,k52,k53,k54,k55,k56,k57,k58,k59,k60,
     6k61,k62,k63,k64,k65,k66,k67,k68,k69,k70,k71,k72,
     7k73,k74,k75,k76,k77,k78,k79,k80,k81,k82,k83,k84,
     8k85,k86,k87,k88,k89,k90,k91,k92,k93,k94,k95,k96
      common/bk14/lfna(15),lfnt(6)
      common/bk15/cpuio(36),cpuip(36)
      common/slar3/nsl,nsntl,nmntl,nslnmx,sltol,slhrd
      common/sabr0/maxnd
      common/sabr1/mind,maxd,ibw,nbw,ipr,npr,nnc,nzn,nnp
      common/WMLthermal/thermalflag !!!!!,thermalconstraint(40000),Tinit(40000),Rqold(40000)
      integer thermalflag,numnp
      common/WMLthermalreordr/idth(40000),nbwth
      integer idth,neq2,nbwth
      character*80 txts,mssg
c
      maxnd=256
cw      call timin (cpuio,cpuip,2,2)
      if (iband.eq.2) go to 20
      write(lfnt(2),130)
      write (*,130)
      nv1=1
      nv2=nv1+2*numnp
      nv3=nv2+2*numnp
      nv4=nv3+numnp
      nv5=nv4+numnp
      nv6=nv5+numnp
      nv7=nv6+numnp
      nv8=nv7+numnp
      lk=maxnd*numnp
      la=nv8+lk
c
c     expand memory
c
      nn=k53+la
      call expndm(nn)
      numnp2=2*numnp
      do 10 i=1,numnp2
   10 nrv(i+nv2-1)=0
c
c     construct connectivity and nodal degree arrays for nike mesh
c
      call nsetka (numel,nrv(nv2),nrv(nv8),b,matp,ncn,nrcc,a(k23)
     1 ,a(k24),a(k25),a(k27),a(k28),a(k29),a(k30),nsl)
      call rezzap (numnp,nrv,nrv(nv1),nrv(nv2),nrv(nv8))
      nv1=nv1+nzn
c
c     contract memory to minimum required
c
      nnn=k53+nv8+nnc
      if(nnn.lt.nn)then
      call expndm(nnn)
      endif
c
c     minimize bandwidth/profile using gps algorithm
c
      call gpsbw (nrv(nv1),nrv(nv2),nrv(nv3),nrv(nv4),
     1          nrv(nv5),nrv(nv6),nrv(nv7),nrv(nv8))
c
      nbw1=2*ibw+1
      nbw2=2*nbw+1
      nbwth=nbw+1
      write(*,*) nbwth
      ipr=ipr+nnp
      npr=npr+nnp
      write(lfnt(2),110) nbw1,ipr,nbw2,npr
      write (*,110) nbw1,ipr,nbw2,npr
c    

      if (npr.lt.ipr) go to 30
      ier=1
      go to 100
c
c     read the reorder vector from the input file
c
   20 nblks=(numnp-1)/16+1
      lfst=1
      llst=min(16,numnp)
      mssg=' error encountered reading the nodal reorder vector'
      do 25 j=1,nblks
      call gttxsg (txts,lcount)
      read (unit=txts,fmt=120,err=102)(nrv(i),i=lfst,llst)
      lfst=llst+1
      llst=min(llst+16,numnp)
   25 continue
c
c     number equations
c

   30 do 90 i=1,numnp

      n=nrv(i)
      il=1
      iu=2
      if (n.eq.0) n=i
      if (nrcc.eq.0) go to 50
      do 40 j=1,nrcc
      if (n.eq.ncn(2,j)) then
      if (ncn(3,j).eq.3) go to 90
      if (ncn(3,j).eq.1) il=2
      if (ncn(3,j).eq.2) iu=1
      endif
   40 continue
   50 do 80 k=il,iu
      if (id(k,n)) 60,60,70
   60 neq=neq+1
      id(k,n)=neq
      go to 80
   70 id(k,n)=0
   80 continue
   90 continue
      write(*,*) 'before thermalnumbering, thermalflag',thermalflag
c     number equations for thermal FEM
      if (thermalflag > 0) then   
      neq2=0
   91 do 99 i=1,numnp
      n=nrv(i)
      il=1
      iu=1
      if (n.eq.0) n=i

   92 do 98 k=il,iu
      if (idth(n)) 93,93,94
   93 neq2=neq2+1
      idth(n)=neq2
      go to 98
   94 idth(n)=0
   98 continue
   99 continue
      endif
c
c      do i=1,numnp
c     write(*,*) idth(i)
c     enddo
  100 continue
cw      call timin (cpuio,cpuip,2,3)
      return
c
  102 call termin (txts,mssg,lcount,1)
c
  110 format(/' bandwidth (profile) before minimization =',i5,' ('
     1,i7,')'/' bandwidth (profile) after minimization  =',i5,' ('
     2,i7,')')
  120 format(16i5)
  130 format(///' gps bandwidth/profile minimization attempted')
      end
