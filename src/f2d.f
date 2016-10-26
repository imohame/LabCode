      subroutine f2d (cht,hgs,d,a,rhs,usi,ix,y,z,beta,freep,matp,
!c     calculate nodal forces for the nominal stress-strain curve
!c     this is only ofcourse for uniaxial loading
     1 matype,den,ym,nbck,udt,udtt,fval,iequit,gs,id)                       
!c    1matype,den,ym,nbck,udt,udtt,fval,iequit,gs,id,plwork)                pl
!c     implicit double precision (a-h,o-z)                                  dp
!c
      use mod_parameters
      
      common/bk00/
     1k01,k02,k03,k04,k05,k06,k07,k08,k09,k10,k11,k12,
     2k13,k14,k15,k16,k17,k18,k19,k20,k21,k22,k23,k24,
     3k25,k26,k27,k28,k29,k30,k31,k32,k33,k34,k35,k36,
     4k37,k38,k39,k40,k41,k42,k43,k44,k45,k46,k47,k48,
     5k49,k50,k51,k52,k53,k54,k55,k56,k57,k58,k59,k60,
     6k61,k62,k63,k64,k65,k66,k67,k68,k69,k70,k71,k72,
     7k73,k74,k75,k76,k77,k78,k79,k80,k81,k82,k83,k84,
     8k85,k86,k87,k88,k89,k90,k91,k92,k93,k94,k95,k96
      common/bk02/ioofc,iphase,imass,lpar(9)
      common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
      common/bk08/kprint,nstep,ite,ilimit,newstf
      common/bk13/xnorm0(6),xnormc(6) !,xnk2d(20)
      common/bk16/maxint,hgc
      common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
      common/bk18/nummat,ityp2d,ako(31)
      common/bk32/nsref,nequit,time,timep,lprint,nprint
      common/bk34/lmm(1)
      common/fissn1/melemt,nnns,ntpe2,n2g,llls
      common/range/mft,mlt,lft,llt,nftm1
      
      common/vect1/r1(nelemg),r2(nelemg),r3(nelemg),
     > r4(nelemg),r5(nelemg),r6(nelemg),
     > r7(nelemg),r8(nelemg),mtype(nelemg),mte(nelemg)
      common/vect2/
     1 ed1(nelemg),ed2(nelemg),ed3(nelemg),ed4(nelemg),
     2 ed5(nelemg),ed6(nelemg),ed7(nelemg),ed8(nelemg),
     3 fd1(nelemg),fd2(nelemg),fd3(nelemg),fd4(nelemg),
     4 fd5(nelemg),fd6(nelemg),fd7(nelemg),fd8(nelemg)
     
      common/vect9/scale(8,nelemg),yz(8,nelemg)
      common/vect17/betap(nelemg)
      common/vect18/sclo(nelemg)
!c     common/plmvc1/plwk(nelemg)                                                pl
      common/colht/icolht,iopta
      common/double/iprec,ncpw,unit
      common/xcom0/imeth,interq,imess,istart,igraf
!!!      common/kkk/km,kkjj
      common/forces/ fff(8,40000)
      common/WMLBC/BCflag
      common/hourglass/fhg(40000,8),fhghis(40000,8),fhg1(nelemg),
     > fhg2(nelemg),fhg3(nelemg),fhg4(nelemg),fhg5(nelemg),
     > fhg6(nelemg),fhg7(nelemg),fhg8(nelemg)
     
      common/hourglass2/hgsstore(40000),hgshis(40000)
      
      common/hourglassk/khg11(nelemg),khg12(nelemg),khg13(nelemg),
     > khg14(nelemg),khg22(nelemg),khg23(nelemg),khg24(nelemg),
     > khg33(nelemg),khg34(nelemg), khg44(nelemg)
     
      common/hourglassks/khg11s(40000),khg12s(40000),
     1khg13s(40000),khg14s(40000),khg22s(40000),khg23s(40000),
     2khg24s(40000),khg33s(40000),khg34s(40000),khg44s(40000)
      common/hourglasskh/khg11h(40000),khg12h(40000),
     1khg13h(40000),khg14h(40000),khg22h(40000),khg23h(40000),
     2khg24h(40000),khg33h(40000),khg34h(40000),khg44h(40000)
      common/hgenergy/hgenerstore(40000),hgenerhis(40000)
      common/totalenergy/totenerstore(40000),totenerhis(40000),
     > inertener(40000)
      common/hgenergypass/hgener(nelemg),
     > totener(nelemg),inertia(nelemg)
      common/hourglassstress/hgstress1c(nelemg),hgstress2c(nelemg)
      common/hgstress/hgstress1store(40000),hgstress2store(40000),
     >   hgstress1his(40000),hgstress2his(40000)
!!	  common /overlapping/ intersec(4, nume), area_coeff(nume), 
!!     > update_flag
!!	  common /stressflag/ ElemFractCode(nume),ElemDecayCount(nume)
!!	  common /sigfrac/ sigmacrit0, sigmacrit1, sigmacrit2,sigmacrit3,
!!     >   DecayCount, f_decay, penalty,fractFlag
      
      
  
      integer BCflag
      dimension cht(*),hgs(*),d(*),a(*),rhs(*),ix(4,*),y(*),z(*),
     1  udt(*),usi(*),beta(*),freep(5,*),matp(*),den(*),nbck(*),
     2  udtt(*),matype(*),ym(4,*),fval(*),gs(*),id(2,*)
!c    3  ,plwork(*)                                                           pl
      logical icolht
!c	  real rp1(128), rp2(128), rp3(128), rp4(128), rp5(128)
!c	 1     rp6(128), rp7(128), rp8(128)
      real rp(8,nelemg)
!!!      real  area_coeff
	  real xo, xp, yo, yp, dx, dy, penalty
	  integer ElemFractCode, ElemDecayCount, DecayCount, q,fractFlag
	  real sigmacrit0, sigmacrit1, sigmacrit2, sigmacrit3
!c
      equivalence (lpar(2),numel)
!c
!c      penalty=2.281E+5


     
      
      numel=numelt
!!c&&&&&&&&&&&&&&&&&&
!      if(imeth.gt.0)then
!      xnk2d(1)=0.
!      xnk2d(2)=0.
!      endif
!!c&&&&&&&&&&&&&&&&
      nelg=(numel-1)/nelemg+1
      nmel=nelemg
      
      if(imeth.gt.0.and.(.not.icolht))call xpegat(1)
      do 110 ng=1,nelg
      lft=1
      llt=min(nelemg,numel-nelemg*(ng-1))
      nft=1+(ng-1)*nelemg
      nlt=min(numel,ng*nelemg)
      nftm1=nft-1
!!!!      write(*,*)'numel numelt lft llt nft nftm1 nelemg ng nelg', 
!!!!     > 'ntime numnp neq'
!!!!      write(*,*)numel,numelt,lft,llt,nft,nftm1,nelemg,ng,nelg, 
!!!!     > ntime,numnp,neq
!!!!!      write(*,*)nftm1,lft,llt
      
      if (icolht) go to 90

      call scm0 (lmm,44*nmel)
      call scm0 (r1 ,10*nmel)
      call scm0 (ed1,16*nmel)
!c     call azero(plwk,128)                                                   pl

!!!!!$OMP PARALLEL DO       
      do 10 i=lft,llt
      lprec=(i-1)*44*iprec
      khg11(i)=khg11s(i+nftm1)
      khg12(i)=khg12s(i+nftm1)
      khg13(i)=khg13s(i+nftm1)
      khg14(i)=khg14s(i+nftm1)
      khg22(i)=khg22s(i+nftm1)
      khg23(i)=khg23s(i+nftm1)
      khg24(i)=khg24s(i+nftm1)
      khg33(i)=khg33s(i+nftm1)
      khg34(i)=khg34s(i+nftm1)
      khg44(i)=khg44s(i+nftm1)
      fhg1(i)=fhg(i+nftm1,1)
      fhg2(i)=fhg(i+nftm1,2)
      fhg3(i)=fhg(i+nftm1,3)
      fhg4(i)=fhg(i+nftm1,4)
      fhg5(i)=fhg(i+nftm1,5)
      fhg6(i)=fhg(i+nftm1,6)
      fhg7(i)=fhg(i+nftm1,7)
      fhg8(i)=fhg(i+nftm1,8)
      hgs(i)=hgsstore(i+nftm1)
      totener(i)=totenerstore(i+nftm1)
      hgener(i)=hgenerstore(i+nftm1)
      hgstress1c(i)=hgstress1store(i+nftm1)
      hgstress2c(i)=hgstress2store(i+nftm1)
!c	  write(*,*) i+nftm1, fhg2(i), 'fhg2'
!c     equation numbers
      lmm(lprec+1)=id(1,ix(1,i+nftm1))	
      lmm(lprec+2)=id(2,ix(1,i+nftm1))
      lmm(lprec+3)=id(1,ix(2,i+nftm1))
      lmm(lprec+4)=id(2,ix(2,i+nftm1))
      lmm(lprec+5)=id(1,ix(3,i+nftm1))
      lmm(lprec+6)=id(2,ix(3,i+nftm1))
      lmm(lprec+7)=id(1,ix(4,i+nftm1))
      lmm(lprec+8)=id(2,ix(4,i+nftm1))
!c     initial coordinates
      yz(1,i) =y(ix(1,i+nftm1))
      yz(2,i) =z(ix(1,i+nftm1))
      yz(3,i) =y(ix(2,i+nftm1))
      yz(4,i) =z(ix(2,i+nftm1))
      yz(5,i) =y(ix(3,i+nftm1))
      yz(6,i) =z(ix(3,i+nftm1))
      yz(7,i) =y(ix(4,i+nftm1))
      yz(8,i) =z(ix(4,i+nftm1))
!c     displacement from reference to step n+1
      ed1(i)  =d(lmm(lprec+1))
      ed2(i)  =d(lmm(lprec+2))
      ed3(i)  =d(lmm(lprec+3))
      ed4(i)  =d(lmm(lprec+4))
      ed5(i)  =d(lmm(lprec+5))
      ed6(i)  =d(lmm(lprec+6))
      ed7(i)  =d(lmm(lprec+7))
      ed8(i)  =d(lmm(lprec+8))
!c      write(*,*) ed1(i),ed2(i),ed3(i),ed4(i)
!c      write(*,*) ed5(i),ed6(i),ed7(i),ed8(i),i
!c     incremental displacement from step n to n+1
      fd1(i)  =usi(lmm(lprec+1))
      fd2(i)  =usi(lmm(lprec+2))
      fd3(i)  =usi(lmm(lprec+3))
      fd4(i)  =usi(lmm(lprec+4))
      fd5(i)  =usi(lmm(lprec+5))
      fd6(i)  =usi(lmm(lprec+6))
      fd7(i)  =usi(lmm(lprec+7))
      fd8(i)  =usi(lmm(lprec+8))
   10 continue
      do 20 i=lft,llt
      mtype(i)=0
      if (matp(nftm1+i).ne.0) mtype(i)=matype(matp(nftm1+i))
   20 sclo(i) =ako(mtype(i)+1)

!!!!!$OMP PARALLEL DO       
      do 30 i=lft,llt
      lprec=(i-1)*44*iprec
      scale(1,i)=float((1+sign(1,lmm(lprec+1)-1))/2)
      scale(2,i)=float((1+sign(1,lmm(lprec+2)-1))/2)
      scale(3,i)=float((1+sign(1,lmm(lprec+3)-1))/2)
      scale(4,i)=float((1+sign(1,lmm(lprec+4)-1))/2)
      scale(5,i)=float((1+sign(1,lmm(lprec+5)-1))/2)
      scale(6,i)=float((1+sign(1,lmm(lprec+6)-1))/2)
      scale(7,i)=float((1+sign(1,lmm(lprec+7)-1))/2)
      scale(8,i)=float((1+sign(1,lmm(lprec+8)-1))/2)
      betap(i)=beta(nftm1+i)
   30 continue
!$OMP PARALLEL DO       
      do 40 i=lft,llt
!c     this is where displacements are zeroed out due to BC's      
      ed1(i)=scale(1,i)*ed1(i)
      ed2(i)=scale(2,i)*ed2(i)
      ed3(i)=scale(3,i)*ed3(i)
      ed4(i)=scale(4,i)*ed4(i)
      ed5(i)=scale(5,i)*ed5(i)
      ed6(i)=scale(6,i)*ed6(i)
      ed7(i)=scale(7,i)*ed7(i)
      ed8(i)=scale(8,i)*ed8(i)
      fd1(i)=scale(1,i)*fd1(i)
      fd2(i)=scale(2,i)*fd2(i)
      fd3(i)=scale(3,i)*fd3(i)
      fd4(i)=scale(4,i)*fd4(i)
      fd5(i)=scale(5,i)*fd5(i)
      fd6(i)=scale(6,i)*fd6(i)
      fd7(i)=scale(7,i)*fd7(i)
      fd8(i)=scale(8,i)*fd8(i)
   40 continue
!c
!c
!c     update stress, stress divergence, and compute element stiffnesses
!c
!c      write(*,*) 'calling stfnf'
      call stfnf (yz,matp,hgs(nft),nmel,a(k12),den)
!c----------------------------------------------------------------
!c----------------------------------------------------------------
!c      if (km.eq.1) then
!c      do 43 i=lft,llt
!c      iicc=(i-1)*44*iprec
!c      write(25,42) i,(lmm(iicc+j),j=9,17),(lmm(iicc+jj),jj=18,26),
!c     .             (lmm(iicc+k),k=27,35),(lmm(iicc+kk),kk=36,44)
!c   42 format(//5x,'Local stiffness matrix :   Element (',i3,')'//
!c     .       5x,9(f8.4)/5x,9(f8.4)/5x,9(f8.4)/5x,9(f8.4))
!c   43 continue
!c      endif
!c
!c      do 48 i=lft,llt
!c      write(25,44) i
!c   44 format(//5x,'Element (',i3,')'/5x,'Iforce-y (1)',
!c     .2x,'Iforce-z (1)',2x,'Iforce-y (2)',2x,'Iforce-z (2)',
!c     .2x,'Iforce-y (3)',2x,'Iforce-z (3)',2x,'Iforce-y (4)',
!c     .2x,'Iforce-z (4)')   
!c      write(25,47) r1(i),r2(i),r3(i),r4(i),r5(i),r6(i),r7(i),r8(i)
!c   47 format(3x,8(e14.4)/)
!c   48 continue
!c
!$OMP PARALLEL DO       
      do 45 i=lft,llt
      fff(1,i+nftm1)=r1(i)
      fff(2,i+nftm1)=r2(i)
      fff(3,i+nftm1)=r3(i)
      fff(4,i+nftm1)=r4(i)
      fff(5,i+nftm1)=r5(i)
      fff(6,i+nftm1)=r6(i)
      fff(7,i+nftm1)=r7(i)
      fff(8,i+nftm1)=r8(i)
   45 continue
   
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c                    Penalty force for convergence
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!$OMP PARALLEL DO       
      rp(1:8,lft:llt)=0.0
!!      do i=lft, llt
!!	      do j=1, 8
!!		      rp(j,i)=0.0
!!		  end do
!!	  end do
	  
!c	  do i=lft,llt
!c
!c		  ink = nftm1+i
!c		  if(ink == 2) then
!c		      if(ElemFractCode(ink)==2) then
!c			      do k=1,4
!c				  if (k == 2) then
!c				  xp = y(ix(k,ink)) + d(id(1,ix(k,ink)))
!c	              xo = y(ix(k,4+1)) + d(id(1,ix(k,4+1)))
!c	              yp = z(ix(k,ink)) + d(id(2,ix(k,ink)))
!c	              yo = z(ix(k,4+1)) + d(id(2,ix(k,4+1)))
!c				  dx = xp - xo
!c	              dy = yp - yo
!c		          rp(2*k-1,i)=-penalty*dx*(0.40000**(ctr_flag(ink)-DecayCount))
!c				  rp(2*k,i)=-penalty*dy*(0.40000**(ctr_flag(ink)-DecayCount))
!c				  else if(k==3) then
!c		          xo = y(ix(k,ink)) + d(id(1,ix(k,ink)))
!c	              xp = y(ix(k,4+1)) + d(id(1,ix(k,4+1)))
!c	              yo = z(ix(k,ink)) + d(id(2,ix(k,ink)))
!c	              yp = z(ix(k,4+1)) + d(id(2,ix(k,4+1)))
!c	              dx = xp - xo
!c	              dy = yp - yo
!c		          rp(2*k-1,i)=penalty*dx*(0.40000**(ctr_flag(ink)-DecayCount))
!c				  rp(2*k,i)=penalty*dy*(0.40000**(ctr_flag(ink)-DecayCount))
!c				  end if
!c				  end do
!c				  do q=1,8
!c	                 rp(q,4+1)=-rp(q,2)
!c	              end do
!c			  end if
!c		  else if(ink==1) then
!c		      if(ElemFractCode(ink)==2) then
!c			      do k=1,4
!c				  if (k .le. 2) then
!c				  xp = y(ix(k,ink)) + d(id(1,ix(k,ink)))
!c	              xo = y(ix(k,4+2)) + d(id(1,ix(k,4+2)))
!c	              yp = z(ix(k,ink)) + d(id(2,ix(k,ink)))
!c	              yo = z(ix(k,4+2)) + d(id(2,ix(k,4+2)))
!c				  dx = xp - xo
!c	              dy = yp - yo
!c		          rp(2*k-1,i)=-penalty*dx*(0.40000**(ctr_flag(ink)-DecayCount))
!c				  rp(2*k,i)=-penalty*dy*(0.40000**(ctr_flag(ink)-DecayCount))
!c				  else 
!c		          xo = y(ix(k,ink)) + d(id(1,ix(k,ink)))
!c	              xp = y(ix(k,4+2)) + d(id(1,ix(k,4+2)))
!c	              yo = z(ix(k,ink)) + d(id(2,ix(k,ink)))
!c	              yp = z(ix(k,4+2)) + d(id(2,ix(k,4+2)))
!c	              dx = xp - xo
!c	              dy = yp - yo
!c		          rp(2*k-1,i)=penalty*dx*(0.40000**(ctr_flag(ink)-DecayCount))
!c				  rp(2*k,i)=penalty*dy*(0.40000**(ctr_flag(ink)-DecayCount))
!c				  end if
!c				  end do
!c				  do q=1,8
!c	                 rp(q,4+2)=-rp(q,1)
!c	              end do
!c			   end if
!c		  end if
!c      end do
	  
	  
!$OMP PARALLEL DO       
	  do i=lft, llt
	      r1(i)=r1(i)+rp(1,i)
		  r2(i)=r2(i)+rp(2,i)
		  r3(i)=r3(i)+rp(3,i)
	      r4(i)=r4(i)+rp(4,i)
	      r5(i)=r5(i)+rp(5,i)
	      r6(i)=r6(i)+rp(6,i)
	      r7(i)=r7(i)+rp(7,i)
	      r8(i)=r8(i)+rp(8,i)
	  end do  
!c-----------------------------------------------------------------
!c-----------------------------------------------------------------   
!c
!c
!$OMP PARALLEL DO       
      do 50 i=lft,llt
      fhghis(i+nftm1,1)=fhg1(i)
      fhghis(i+nftm1,2)=fhg2(i)
      fhghis(i+nftm1,3)=fhg3(i)
      fhghis(i+nftm1,4)=fhg4(i)
      fhghis(i+nftm1,5)=fhg5(i)
      fhghis(i+nftm1,6)=fhg6(i)
      fhghis(i+nftm1,7)=fhg7(i)
      fhghis(i+nftm1,8)=fhg8(i)
      hgshis(i+nftm1)=hgs(i)
      totenerhis(i+nftm1)=totener(i)
      hgenerhis(i+nftm1)=hgener(i)
      hgstress1his(i+nftm1)=hgstress1c(i)
      hgstress2his(i+nftm1)=hgstress2c(i)
      inertener(i+nftm1)=inertia(i)
   50 continue
!c     account for inertial effects
!c
      call inerta (ym(1,nft),udt,udtt,lmm,lmm(8*iprec+1),iequit,
     1 iprec)
!c
!c
!c     magnetic body force loads
!c
      call lodmbf (fval,freep(1,nft),ym(1,nft),den,matp(nft))
!c
!c
!c     body force loads
!c
      call lodbdy (fval,ym(1,nft))
!c
!c
!c     displacement boundary conditions
!c
!c      write(*,*) 'BCflag', BCflag
      if (BCflag.eq.0) then
      call ulbc(iequit,lmm,lmm(8*iprec+1),fval,a(k73),a(k74),a(k75),
     1 nbck(nft),ix(1,nft),id,a(k72),iprec)
!c      write(*,*) 'in first if', BCflag
      endif
      if (BCflag.eq.1) then
      call ulbc2(iequit,lmm,lmm(8*iprec+1),fval,a(k73),a(k74),a(k75),
     1 nbck(nft),ix(1,nft),id,a(k72),iprec)	
!c      write(*,*) 'in second if', BCflag
      endif  
!c
!c
!c     element birth and death options
!c
      call brthdt(freep(1,nft),mtype,beta(nft),matp(nft),lmm,
     1 lmm(8*iprec+1),iprec)
!c
!c.... assemble plastic work contributions into global element vector
!c     do 50 i=lft,llt                                                        pl
!c     plwork(nftm1+i)=plwk(i)                                                pl
!c  50 continue                                                               pl
!c  			  
!c      write(*,*) 'in f2d before rhs update rhs(12*44+2)', rhs(12*44+2)
!c      write(*,*) 'in f2d before rhs update r2(13)', r2(13)

!!!!$OMP PARALLEL DO       
      do 60 i=lft,llt
      lprec=(i-1)*44*iprec
!c      write(*,*) r1(i),r2(i),r3(i),r4(i)
!c      write(*,*) r5(i),r6(i),r7(i),r8(i)
      rhs(lmm(lprec+1))=rhs(lmm(lprec+1))-scale(1,i)*r1(i)
      rhs(lmm(lprec+2))=rhs(lmm(lprec+2))-scale(2,i)*r2(i)
      rhs(lmm(lprec+3))=rhs(lmm(lprec+3))-scale(3,i)*r3(i)
      rhs(lmm(lprec+4))=rhs(lmm(lprec+4))-scale(4,i)*r4(i)
      rhs(lmm(lprec+5))=rhs(lmm(lprec+5))-scale(5,i)*r5(i)
      rhs(lmm(lprec+6))=rhs(lmm(lprec+6))-scale(6,i)*r6(i)
      rhs(lmm(lprec+7))=rhs(lmm(lprec+7))-scale(7,i)*r7(i)
      rhs(lmm(lprec+8))=rhs(lmm(lprec+8))-scale(8,i)*r8(i)
   60 continue
!c
!c      write(*,*) 'in f2d after rhs update rhs(12*44+2)', rhs(12*44+2)
      if(iphase.eq.3)goto 110
   70 if (newstf) 110,80,110
   80 continue
      melemt=llt-lft+1
      if (ioofc.eq.1) call fissln (iopta,gs,gs,lmm)
      if (ioofc.eq.0) call wdiskf (lmm,melemt)
      iopta=-5
      go to 110
   90 continue
   
!!!!!$OMP PARALLEL DO       
      do 100 i=lft,llt
      lprec=(i-1)*44*iprec
      lmm(lprec+1)=id(1,ix(1,i+nftm1))
      lmm(lprec+2)=id(2,ix(1,i+nftm1))
      lmm(lprec+3)=id(1,ix(2,i+nftm1))
      lmm(lprec+4)=id(2,ix(2,i+nftm1))
      lmm(lprec+5)=id(1,ix(3,i+nftm1))
      lmm(lprec+6)=id(2,ix(3,i+nftm1))
      lmm(lprec+7)=id(1,ix(4,i+nftm1))
      lmm(lprec+8)=id(2,ix(4,i+nftm1))
  100 continue
!c     icolht is true for new stiffness formation
      if (icolht) call clht (cht,lmm,nlt-nft+1,iprec)
  110 continue
      if (icolht) return
      hgc=0.0001/16.*1.  !0.0025 works for QS, hgc=alpha/16.*psi, suggested alpha=0.1, psi is fraction of critical damping
      if(imeth.gt.0)call xpegat(4)
      return

      end












