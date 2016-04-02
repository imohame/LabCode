      subroutine hspwr (u,udt,udtt,id,temp,y,z)
c     implicit double precision (a-h,o-z)                                    dp
c
c     write nodal data into hsp file
c
      parameter (nume=40000)
	  parameter (nume2=20000)
      common/bk02/ioofc,iphase,imass,lpar(9)
      common/bk08/kprint,nstep,ite,ilimit,newstf
      common/bk10/npb,nodep(2,8)
      common/bk14/lfna(15),lfnt(6)
      common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
      common/bk23/itemp,itherm,irtin
      common/bk32/nsref,nequit,time,timep,lprint,nprint
	  common /pcracktip/ connect(4,nume2),node(2,nume2),penta(nume2),
     >	                 ndflag(2,nume2), numnpt, numeltu, ndc
c
      dimension u(*),udt(*),udtt(*),id(2,*),temp(*),d(7)
      dimension y(*),z(*),d1(2)

	  real node
c	  
	  do i=1,nume2
	      do j=1,2
		      node(j,i)=0.0
		  end do
	  end do
c
c     print displacements, velocities, and accelerations
c
      ic=0
      do 50 ib=1,npb
      node1=nodep(1,ib)
      if (node1.eq.0) go to 50
      node2=nodep(2,ib)
      do 40 ii=node1,node2
      do 10 i=1,7
   10 d(i)=0.
      do 20 i=1,2
      kk=id(i,ii)
      if (kk.eq.0) go to 20
      d(i)=u(kk)
      if (imass.eq.0) go to 20
      d(i+2)=udt(kk)
      d(i+4)=udtt(kk)
   20 continue
      if (itemp.ne.0) d(7)=temp(ii)
c      if (ic.gt.0) go to 30
c      ic=50
c      call header
c      write(lfnt(2),60) nstep,timep
c   30 ic=ic-1
c      write(lfnt(2),70) ii,(d(l),l=1,7)
ck----for history plot by m.k.
          d1(1)=d(1)+y(ii)
          d1(2)=d(2)+z(ii)
c      write(29,100) ii,d1(1),d1(2)
          node(1,ii)=d1(1)
		  node(2,ii)=d1(2)
ck----end
   40 continue
	  
   50 continue
      return
c
c
c   60 format(///' n o d a l   p r i n t   o u t   f o r   t i m e   s t
c     1e p ',i5,31x,' ( at time ',1pe10.4,' )'/
c     2/' nodal point y-displacement  z-displacement    y-velocity      z
c     3-velocity  y-acceleration  z-acceleration  temperature')
c   70 format (i9,7e15.4)
  100 format(i5,5x,e15.6,5x,e15.6)
      end
