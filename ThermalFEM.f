!!!      subroutine ThermalFEM(y,z,ele,nod,maxint)
!!!!c     This program is for thermal finite element calculation
!!!!c     Variable Rqglobal, Rqold, Rqolder, correspond to variable at 
!!!!c     current, last, last two time step.
!!!!c     Variable Tinit, Tinitd, Tinitdold, correspond to variable at 
!!!!c     current, last, last two time step. 
!!!          use  mod_dtimeSpecs
!!!          use  mod_file_units
!!!          use  mod_parameters
!!!          
!!!!!!!      parameter(nume = 40000)    
!!!      integer ele,nod,maxint
!!!      dimension y(nod), z(nod)
!!!      
!!!      integer kk,ii,bb,INFO_flag
!!!
!!!      real*8 rhocp(ele),eta(maxint),psi(maxint)
!!!      real*8 N(maxint,4),ShpFcnDeriv(maxint,2,4),coords(maxint,2)
!!!      real*8 Jmatrix(2,2),Jacob(ele,maxint),Gamma(2,2),B(maxint,2,4)
!!!      
!!!      real*8 kappa(2,2),B1(2,4),B2(2,4),B3(2,4),B4(2,4)
!!!      real*8 B1T(4,2),B2T(4,2),B3T(4,2),B4T(4,2)
!!!      real*8 N1T(4,1),N2T(4,1),N3T(4,1),N4T(4,1)
!!!      real*8 Nh(1,4),NhT(4,1),N1(1,4),N2(1,4),N3(1,4),N4(1,4)
!!!      real*8 ci1(4,4),ci2(4,4),ci3(4,4),ci4(4,4),c(ele,4,4)
!!!      real*8 ke1a(2,4),ke2a(2,4),ke3a(2,4),ke4a(2,4)
!!!      real*8 ke1b(4,4),ke2b(4,4),ke3b(4,4),ke4b(4,4)
!!!      real*8 ke1(4,4),ke2(4,4),ke3(4,4),ke4(4,4),ke(ele,4,4)
!!!	  real*8 ks1a(4,4),ks2a(4,4),ks3a(4,4),ks4a(4,4)     ! for stress assisted diffusion
!!!	  real*8 ks1b(4,1),ks2b(4,1),ks3b(4,1),ks4b(4,1)
!!!	  real*8 ks1c(4,4),ks2c(4,4),ks3c(4,4),ks4c(4,4)
!!!	  real*8 ks1(4,4),ks2(4,4),ks3(4,4),ks4(4,4)
!!!      real*8 Rq1(4,1),Rq2(4,1),Rq3(4,1),Rq4(4,1)
!!!      
!!!      real*8 Cglobal(nod,nod), Kglobal(nod,nod), Rqglobal(nod,1)
!!!      real*8 RHS(nod,1), Rq(ele,4,1),DtCondDiff !,TkelvinDelta(nod,1),TkelvinCurrent(nod,1)
!!!      real*8 he(ele,4,4),rhe(ele,4,1), h12e(4,4),h23e(4,4),
!!!     > h34e(4,4)
!!!      real*8 h41e(4,4), rh12e(4,1),rh23e(4,1),rh34e(4,1),
!!!     > rh41e(4,1)
!!!      real*8 rh12egp(4,1),rh23egp(4,1),rh34egp(4,1),rh41egp(4,1)
!!!      real*8 h12i(4,4),h23i(4,4),h34i(4,4),h41i(4,4),W
!!!      real*8 thermalkd(ele),Tfld(ele),DijSijd(ele),Tinitd(nod,1)
!!!      !,hd(ele)
!!!      real*8 beta,Tcoeff1(nod,nod),Rnew(nod,1), sigp(4,1)
!!!      real*8 Rold(nod,1),J12,J23,J34,J41,Rolder(nod,1)
!!!      real*8 RHS2(nod,1),LHS(nod,nod),TdSol(nod,1),alpha
!!!      real*8 Tcoeff2(nod,nod)
!!!     
!!!	  real*8 gradp1a(2,1), gradp2a(2,1), gradp3a(2,1),gradp4a(2,1)
!!!	  real*8 gradp(2,1), gradpabs
!!!      real*8 gradslip1a(2,1), gradslip2a(2,1), gradslip3a(2,1)
!!!	  real*8 gradslip4a(2,1), gradslipa(2,1)	  
!!!      integer  IPIV(nod)
!!!      
!!!      common/bk00/
!!!     1k01,k02,k03,k04,k05,k06,k07,k08,k09,k10,k11,k12,
!!!     2k13,k14,k15,k16,k17,k18,k19,k20,k21,k22,k23,k24,
!!!     3k25,k26,k27,k28,k29,k30,k31,k32,k33,k34,k35,k36,
!!!     4k37,k38,k39,k40,k41,k42,k43,k44,k45,k46,k47,k48,
!!!     5k49,k50,k51,k52,k53,k54,k55,k56,k57,k58,k59,k60,
!!!     6k61,k62,k63,k64,k65,k66,k67,k68,k69,k70,k71,k72,
!!!     7k73,k74,k75,k76,k77,k78,k79,k80,k81,k82,k83,k84,
!!!     8k85,k86,k87,k88,k89,k90,k91,k92,k93,k94,k95,k96
!!!      common /main_block/ a(1)
!!!      common/bk08/kprint,nstep,ite,ilimit,newstf
!!!      common/bk11/cnwmk(2),iequit,iprint,isref,iPrintOutputFlag 
!!!      integer iequit,iprint,isref,iPrintOutputFlag 
!!!      common/bk26/dt,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
!!!!!!!      real*8 dt,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
!!!      common/WMLthermali/connect(nume,4)
!!!      integer connect
!!!      common /WMLthermalBC/const,constraint(nume)
!!!      integer const,constraint
!!!      common /WMLthermal2/thermalk(nume),thermalh(nume),
!!!     >      thermalRo(nume),thermalcp(nume),thermalx(nume)
!!!     >      ,thermalD(1000)
!!!!!!!!!      real*8 thermalk,thermalh,thermalRo,thermalcp,thermalKa     
!!!      common /WMLthermal/ thermalflag,thermalconstraint(nume),   
!!!     > Tinit(nume),Rqold(nume)
!!!      integer thermalflag,thermalconstraint
!!!!!!!!!      real*8 Tinit,Rqold
!!!      common/WMLthermalSolve/Rqolder(nume,1),Tinitdold(nume) !,dummy(40000,1)
!!!!!!!!!      real*8 Rqolder,Tinitdold
!!!      common/WMLthermalmatpoly/DijSije(nume),Tele(nume)
!!!!!!!!!      real*8 DijSije,Tele
!!!	  common/couplinganalysis/ TDflag
!!!	  integer TDflag      ! 0 for thermal analysis, 1 for diffusion analysis
!!!	  common/hydrodiffusion/ hycon(nume)
!!!!!!!!      real*8 hycon
!!!	  common /propag/ sigalt(4,nume)
!!!!!!      real*8 sigalt
!!!	  common /grad_pressure/ gradpdata(3,nume)
!!!!!!      real*8 gradpdata
!!!	  common/intgrt/   nintg
!!!      integer nintg
!!!      
!!!      
!!!      integer nstep,intmeth
!!!	  real*8  Vh, Temp, press(nume), pressnode(nume)
!!!      real*8 rDiffCoeffElem
!!!	  integer nelenode(nume)
!!!!!!!!!!         write(*,*) '----------- thermalkd(i)',thermalkd(1:10)
!!!!!!!!!!            write(*,*) '----------- y(i,1)',y(1),y(2),y(12),y(13)
!!!!!!!!!!            write(*,*) '----------- y(i,1)',z(1),z(2),z(12),z(13)
!!!!!!!!!      common/fissn1/melemt,nnns,ntpe2,n2g,llls
!!!!!!!!!!!!	  common/wblock2/  g_source, g_immob, g_minter, g_recov, b_v,
!!!!!!!!!!!!     1                             b_vvec(87),nmo,nim
!!!!!!!!!!!!      common/wblock8/  abc(573,nume,4), his(573,nume,4)
!!!
!!!      
!!!!!      common/WMLthermalperiod/period(40000)
!!!!!      common/WMLthermalreordr/idth(nume),nbwth
!!!      
!!!!!!!!      integer idth,nbwth,period
!!!!c	  common /GND_loop/ gradslip(2,24,40000), rho_gnd(2,24,40000)
!!!!!!!	  real  slipnode(nume,24), slipp(4,1)
!!!!!!!!	  integer nelenode2(nume), nsi
!!!!!!!	  real gradslip, rho_gnd, ly, lz
!!!!!!!!!!      write(*,*) 'in ThermalFEM'
!!!
!!!      DtCondDiff=dt
!!!      DtCondDiff=dt/idtimeStepsSol(1)
!!!! ##################################################
!!!      do idt=1,idtimeStepsSol(1) !---- this is the refinement loop, chops dt by n-sections
!!!        Cglobal=0.0
!!!        Kglobal=0.0
!!!        Rqglobal=0.0
!!!        Tcoeff1=0.0
!!!        Tcoeff2=0.0
!!!        RHS2=0.0
!!!        LHS=0.0
!!!        TdSol=0.0
!!!        IPIV=0.0
!!!        INFO_flag=0.0
!!!        DijSijd=0
!!!        
!!!        
!!!!!!!      do idt=1,1 !---- this is the refinement loop, chops dt by n-sections
!!!!!!!        DtCondDiff=dt/1
!!!!!!!! ##################################################
!!!
!!!      Vh=2.0E-6           ! partial molar volume of hydrogen
!!!      Temp=293.0
!!!	  
!!!  	  if(TDflag==0) then
!!!          ks1=0.0
!!!          ks2=0.0
!!!          ks3=0.0
!!!          ks4=0.0
!!!!!!	      do i=1,4
!!!!!!		      do j=1,4
!!!!!!			      ks1(i,j)=0.0
!!!!!!				  ks2(i,j)=0.0
!!!!!!				  ks3(i,j)=0.0
!!!!!!				  ks4(i,j)=0.0
!!!!!!			  end do
!!!!!!		  end do
!!!      else if(TDflag==1) then !!!----diffusion
!!!!!	      do i=1,nod
!!!		      pressnode=0.0
!!!			  nelenode=0          ! number of elements connecting node i
!!!!!!		  end do
!!!	      do i=1,ele
!!!		      press(i)=(sigalt(1,i)+sigalt(2,i)+sigalt(3,i))/3.0
!!!			  do j=1,4
!!!			      pressnode(connect(i,j))=pressnode(connect(i,j))+
!!!     >              press(i)
!!!				  nelenode(connect(i,j))=nelenode(connect(i,j))+1
!!!			  end do
!!!		  end do
!!!		  do i=1,nod
!!!		      if(nelenode(i)==0) then
!!!			      write(*,*) 'warning a/0, check ThermalFEM.f'
!!!			  end if
!!!		      pressnode(i)=pressnode(i)/nelenode(i)   
!!!		  end do
!!!	  end if
!!!	  
!!!!c     Initializing temperatures based on the equation numbers
!!!      if (nstep.eq.0) then
!!!!!!!c     Setting initial vectors for intmeth 2 and 3      
!!!          Tinitdold(1:nod,1)=0
!!!          Rqolder(1:nod,1)=0
!!!          Rqold(1:nod,1)=0
!!!      endif
!!!
!!!!c     Setting Element Material Parameters, removing "!"s would make the
!!!!c     non-dimensional terms dimensional
!!!!      ---use the mat diff. coeff.--- this has to be done before calling DiffCoeffTableGetD
!!!      thermalkd(1:ele)=thermalk(1:ele) 
!!!      do i=1,ele
!!!	     if (TDflag==0) then    ! thermal analysis
!!!             rhocp(i)=(thermalRo(i)*thermalcp(i)) !/thermalKa(i)     !1./etae(i)!*E/To   !!!!             rhocp(i)=1./etae(i)!*E/To
!!!		 else if (TDflag==1) then ! diffusion analysis
!!!		     rhocp(i)=1.0
!!!!             --- find the D from the table or equation
!!!         rDiffCoeffElem=0.0
!!!         call DiffCoeffTableGetD(i,rDiffCoeffElem,Tele(i,1),a(k08))
!!!         thermalkd(i)=rDiffCoeffElem
!!!		 end if
!!!      end do
!!!      if (TDflag==1) then ! diffusion analysis
!!!        write(iFU_DiffCoeff_out,26) thermalkd(1:ele)
!!!      endif
!!!!      
!!!      DijSijd(1:ele)=DijSije(1:ele)!    ! Q rate of heat generation, textbook by Cook
!!!      
!!!!!!!!!   ----------------------------------testing
!!!!!!!!      write(*,*)'Tinit',Tinit(1:nod,1)
!!!!!!!!      write(*,*)const,constraint(1:const)
!!!!!!!!      write(*,*)const,constraint(1:nod)
!!!!!!!!!   ----------------------------------testing
!!!!!!  we can update this here b/c Tinit is saved         
!!!      Tinitd(1:nod,1)=Tinit(1:nod,1)
!!!      
!!!!c      Defining integration points, corresponding with node 1
!!!!c      at psi1,eta1=sqrt(3)/3, w=1 and maxint=4 and going ccw 
!!!      psi(1)=-.57735
!!!      eta(1)=-.57735
!!!      psi(2)=.57735
!!!      eta(2)=-.57735
!!!      psi(3)=.57735
!!!      eta(3)=.57735
!!!      psi(4)=-.57735
!!!      eta(4)=.57735
!!!      
!!!          
!!!!c     Setting the Weight Factors for Gauss integration, hardcoded below
!!!!c     to maxint=4      
!!!      if (maxint.eq.1) then
!!!           W=4
!!!      elseif (maxint.eq.4) then
!!!           W=1
!!!      endif 
!!!      
!!!      call shapefunctions(psi,eta,N,ShpFcnDeriv,maxint)
!!!! ##################################################
!!!      
!!!!!!!c     Zeroing out Global C, K matrix and Rq vector
!!!!!!      do i=1,nod
!!!!!!          do j=1,nod
!!!!!!              Cglobal(i,j)=0.000
!!!!!!              Kglobal(i,j)=0.
!!!!!!          end do
!!!!!!          Rqglobal(i,1)=0.
!!!!!!      end do
!!!!!!!!      write(*,*) '----------- Tele'
!!!      
!!!!c     starting the big element loop
!!!      do i=1,ele
!!!!c           thermal conductivity      
!!!            kappa(1,1)=thermalkd(i) 
!!!            kappa(1,2)=0.
!!!            kappa(2,1)=0.
!!!            kappa(2,2)=thermalkd(i)
!!!!c      reading nodes associated with element i 
!!!            node1=connect(i,1)
!!!            node2=connect(i,2)
!!!       	    node3=connect(i,3)
!!!            node4=connect(i,4)
!!!            Coords(1,1)=y(node1)!*bv          ! current configuration
!!!            Coords(1,2)=z(node1)!*bv
!!!            Coords(2,1)=y(node2)!*bv
!!!            Coords(2,2)=z(node2)!*bv
!!!            Coords(3,1)=y(node3)!*bv
!!!            Coords(3,2)=z(node3)!*bv
!!!            Coords(4,1)=y(node4)!*bv
!!!            Coords(4,2)=z(node4)!*bv
!!!!!!!!!!            write(*,*) '----------- kappa',kappa
!!!!!!!!!!            write(*,*) '----------- Coords',Coords
!!!!c     given shape function derivatives wrt isoparametric coords
!!!!c     and coordinate matrix for the 4 Gauss points, find the B-matrix  
!!!       do j=1,maxint
!!!          Jmatrix=matmul(ShpFcnDeriv(j,1:2,1:4),Coords)
!!!          Jacob(i,j)=Jmatrix(1,1)*Jmatrix(2,2)-Jmatrix(2,1)*Jmatrix(1,2)
!!!          Gamma(1,1)=1/Jacob(i,j)*Jmatrix(2,2)
!!!          Gamma(1,2)=-1/Jacob(i,j)*Jmatrix(1,2)
!!!          Gamma(2,1)=-1/Jacob(i,j)*Jmatrix(2,1)
!!!          Gamma(2,2)=1/Jacob(i,j)*Jmatrix(2,2)
!!!          B(j,1:2,1:4)=matmul(Gamma,ShpFcnDeriv(j,1:2,1:4))
!!!       enddo
!!!        B1(1:2,1:4)=B(1,1:2,1:4)
!!!        B2(1:2,1:4)=B(2,1:2,1:4)
!!!        B3(1:2,1:4)=B(3,1:2,1:4)
!!!        B4(1:2,1:4)=B(4,1:2,1:4)
!!!        
!!!        call matrixtrans(B1,B1T,2,4)
!!!        call matrixtrans(B2,B2T,2,4)
!!!        call matrixtrans(B3,B3T,2,4)
!!!        call matrixtrans(B4,B4T,2,4)
!!!        
!!!        call matrixtrans(N(1,1:4),N1T,1,4)
!!!        call matrixtrans(N(2,1:4),N2T,1,4)
!!!        call matrixtrans(N(3,1:4),N3T,1,4)
!!!        call matrixtrans(N(4,1:4),N4T,1,4)  
!!!
!!!        N1(1,1:4)=N(1,1:4)
!!!        N2(1,1:4)=N(2,1:4)
!!!        N3(1,1:4)=N(3,1:4)
!!!        N4(1,1:4)=N(4,1:4)      
!!!
!!! 
!!!        ci1=matmul(N1T,N1)
!!!        ci2=matmul(N2T,N2)
!!!        ci3=matmul(N3T,N3)
!!!        ci4=matmul(N4T,N4)        
!!!
!!!!c     Element c matrix for the 4 Gauss points        
!!!        do ii=1,4
!!!        do jj=1,4
!!!        ci1(ii,jj)=ci1(ii,jj)*rhocp(i)*Jacob(i,1)*W
!!!        ci2(ii,jj)=ci2(ii,jj)*rhocp(i)*Jacob(i,2)*W
!!!        ci3(ii,jj)=ci3(ii,jj)*rhocp(i)*Jacob(i,3)*W
!!!        ci4(ii,jj)=ci4(ii,jj)*rhocp(i)*Jacob(i,4)*W
!!!        enddo
!!!        enddo
!!!        
!!!        ke1a=matmul(kappa,B1)
!!!        ke2a=matmul(kappa,B2)
!!!        ke3a=matmul(kappa,B3)
!!!        ke4a=matmul(kappa,B4)
!!!        
!!!        ke1b=matmul(B1T,ke1a)
!!!        ke2b=matmul(B2T,ke2a)        
!!!        ke3b=matmul(B3T,ke3a)        
!!!        ke4b=matmul(B4T,ke4a)  
!!!!c      Element k matrix for the 4 Gauss points            
!!!        ke1(1:4,1:4)=ke1b(1:4,1:4)*Jacob(i,1)*W
!!!        ke2(1:4,1:4)=ke2b(1:4,1:4)*Jacob(i,2)*W
!!!        ke3(1:4,1:4)=ke3b(1:4,1:4)*Jacob(i,3)*W
!!!        ke4(1:4,1:4)=ke4b(1:4,1:4)*Jacob(i,4)*W
!!!
!!!!c      Ks matrix caused by stress assisted diffusion		
!!!		if(TDflag==1) then  !!!-- diffusion
!!!		    ks1a(1:4,1:4)=ke1b(1:4,1:4)*Vh/(8.3142*Temp)
!!!			ks2a(1:4,1:4)=ke2b(1:4,1:4)*Vh/(8.3142*Temp)
!!!			ks3a(1:4,1:4)=ke3b(1:4,1:4)*Vh/(8.3142*Temp)
!!!			ks4a(1:4,1:4)=ke4b(1:4,1:4)*Vh/(8.3142*Temp)
!!!			
!!!			sigp(1,1)=pressnode(connect(i,1))
!!!			sigp(2,1)=pressnode(connect(i,2))
!!!			sigp(3,1)=pressnode(connect(i,3))
!!!			sigp(4,1)=pressnode(connect(i,4))
!!!			
!!!			ks1b=matmul(ks1a,sigp)
!!!			ks2b=matmul(ks2a,sigp)
!!!			ks3b=matmul(ks3a,sigp)
!!!			ks4b=matmul(ks4a,sigp)
!!!			
!!!			ks1c=matmul(ks1b,N1)
!!!			ks2c=matmul(ks2b,N2)
!!!			ks3c=matmul(ks3b,N3)
!!!			ks4c=matmul(ks4b,N4)
!!!			
!!!			ks1(1:4,1:4)=ks1c(1:4,1:4)*Jacob(i,1)*W
!!!			ks2(1:4,1:4)=ks2c(1:4,1:4)*Jacob(i,2)*W
!!!			ks3(1:4,1:4)=ks3c(1:4,1:4)*Jacob(i,3)*W
!!!			ks4(1:4,1:4)=ks4c(1:4,1:4)*Jacob(i,4)*W
!!!			
!!!!c       calculate pressure gradient
!!!            gradp1a=matmul(B1, sigp)
!!!			gradp2a=matmul(B2, sigp)
!!!			gradp3a=matmul(B3, sigp)
!!!			gradp4a=matmul(B4, sigp)
!!!			gradp(1:2,1)=gradp1a(1:2,1)+gradp2a(1:2,1)
!!!     >                   +gradp3a(1:2,1)+gradp4a(1:2,1)
!!!            gradp(1,1)=gradp(1,1)/4.0
!!!            gradp(2,1)=gradp(2,1)/4.0
!!!            gradpabs=sqrt((gradp(1,1))**2+(gradp(2,1))**2)
!!!			
!!!            gradpdata(1,i)=gradp(1,1)
!!!            gradpdata(2,i)=gradp(2,1)
!!!            gradpdata(3,i)=gradpabs	
!!!        end if					
!!!		
!!!!c      Element Rq matrix given plastic work for each element 
!!!        Rq1(1:4,1)=N1T(1:4,1)*DijSijd(i)*W*Jacob(i,1)
!!!        Rq2(1:4,1)=N2T(1:4,1)*DijSijd(i)*W*Jacob(i,2)
!!!        Rq3(1:4,1)=N3T(1:4,1)*DijSijd(i)*W*Jacob(i,3)
!!!        Rq4(1:4,1)=N4T(1:4,1)*DijSijd(i)*W*Jacob(i,4)
!!!       
!!!        ke(i,1:4,1:4)=ke1(1:4,1:4)+ke2(1:4,1:4)+ke3(1:4,1:4)+
!!!     >            ke4(1:4,1:4)-ks1(1:4,1:4)-ks2(1:4,1:4)
!!!     >            -ks3(1:4,1:4)-ks4(1:4,1:4)
!!!        c(i,1:4,1:4)=ci1(1:4,1:4)+ci2(1:4,1:4)+ci3(1:4,1:4)+ 
!!!     >            ci4(1:4,1:4)
!!!        Rq(i,1:4,1)=Rq1(1:4,1)+Rq2(1:4,1)+Rq3(1:4,1)+Rq4(1:4,1)
!!!      enddo
!!!
!!!
!!!!!!!!c     Adding convection terms to the K matrix and R vector
!!!!!!!      do i=1,ele
!!!!!!!        do j=1,4
!!!!!!!          do kk=1,4
!!!!!!!             he(i,j,kk)=0.
!!!!!!!          enddo
!!!!!!!          rhe(i,j,1)=0.
!!!!!!!        enddo
!!!!!!!      enddo
!!!      he=0.0
!!!      rhe(1:ele,1:4,1)=0.
!!!
!!!!c     Assembling global matrices
!!!      do i=1,ele
!!!        node1=connect(i,1)
!!!        node2=connect(i,2)
!!!        node3=connect(i,3)
!!!        node4=connect(i,4)
!!!    
!!!        Kglobal((node1),(node1))=Kglobal((node1),
!!!     1       (node1))+ke(i,1,1)+he(i,1,1)
!!!
!!!        Kglobal((node1),(node2))=Kglobal((node1),
!!!     1       (node2))+ke(i,1,2)+he(i,1,2)
!!!
!!!        Kglobal((node1),(node3))=Kglobal((node1),
!!!     1       (node3))+ke(i,1,3)+he(i,1,3)
!!!
!!!        Kglobal((node1),(node4))=Kglobal((node1),
!!!     1       (node4))+ke(i,1,4)+he(i,1,4)
!!!
!!!        Kglobal((node2),(node1))=Kglobal((node2),
!!!     1       (node1))+ke(i,2,1)+he(i,2,1)
!!!
!!!        Kglobal((node2),(node2))=Kglobal((node2),
!!!     1       (node2))+ke(i,2,2)+he(i,2,2)
!!!
!!!        Kglobal((node2),(node3))=Kglobal((node2),
!!!     1       (node3))+ke(i,2,3)+he(i,2,3)
!!!
!!!        Kglobal((node2),(node4))=Kglobal((node2),
!!!     1       (node4))+ke(i,2,4)+he(i,2,4)
!!!
!!!        Kglobal((node3),(node1))=Kglobal((node3),
!!!     1       (node1))+ke(i,3,1)+he(i,3,1)
!!!
!!!        Kglobal((node3),(node2))=Kglobal((node3),
!!!     1       (node2))+ke(i,3,2)+he(i,3,2)
!!!
!!!        Kglobal((node3),(node3))=Kglobal((node3),
!!!     1       (node3))+ke(i,3,3)+he(i,3,3)
!!!
!!!        Kglobal((node3),(node4))=Kglobal((node3),
!!!     1       (node4))+ke(i,3,4)+he(i,3,4)
!!!
!!!        Kglobal((node4),(node1))=Kglobal((node4),
!!!     1       (node1))+ke(i,4,1)+he(i,4,1)
!!!
!!!        Kglobal((node4),(node2))=Kglobal((node4),
!!!     1       (node2))+ke(i,4,2)+he(i,4,2)
!!!
!!!        Kglobal((node4),(node3))=Kglobal((node4),
!!!     1       (node3))+ke(i,4,3)+he(i,4,3)
!!!
!!!        Kglobal((node4),(node4))=Kglobal((node4),
!!!     1       (node4))+ke(i,4,4)+he(i,4,4)
!!!
!!!    
!!!        Cglobal((node1),(node1))=Cglobal((node1),
!!!     1       (node1))+c(i,1,1)
!!!        Cglobal((node1),(node2))=Cglobal((node1),
!!!     1       (node2))+c(i,1,2)
!!!        Cglobal((node1),(node3))=Cglobal((node1),
!!!     1       (node3))+c(i,1,3)
!!!        Cglobal((node1),(node4))=Cglobal((node1),
!!!     1       (node4))+c(i,1,4)
!!!        Cglobal((node2),(node1))=Cglobal((node2),
!!!     1       (node1))+c(i,2,1)
!!!        Cglobal((node2),(node2))=Cglobal((node2),
!!!     1       (node2))+c(i,2,2)
!!!        Cglobal((node2),(node3))=Cglobal((node2),
!!!     1       (node3))+c(i,2,3)
!!!        Cglobal((node2),(node4))=Cglobal((node2),
!!!     1       (node4))+c(i,2,4)
!!!        Cglobal((node3),(node1))=Cglobal((node3),
!!!     1       (node1))+c(i,3,1)
!!!        Cglobal((node3),(node2))=Cglobal((node3),
!!!     1       (node2))+c(i,3,2)
!!!        Cglobal((node3),(node3))=Cglobal((node3),
!!!     1       (node3))+c(i,3,3)
!!!        Cglobal((node3),(node4))=Cglobal((node3),
!!!     1       (node4))+c(i,3,4)
!!!        Cglobal((node4),(node1))=Cglobal((node4),
!!!     1       (node1))+c(i,4,1)
!!!        Cglobal((node4),(node2))=Cglobal((node4),
!!!     1       (node2))+c(i,4,2)
!!!        Cglobal((node4),(node3))=Cglobal((node4),
!!!     1       (node3))+c(i,4,3)
!!!        Cglobal((node4),(node4))=Cglobal((node4),
!!!     1       (node4))+c(i,4,4)
!!! 
!!!        Rqglobal((node1),1)=Rqglobal((node1),1)+Rq(i,1,1)+
!!!     1       rhe(i,1,1)
!!!
!!!        Rqglobal((node2),1)=Rqglobal((node2),1)+Rq(i,2,1)+
!!!     1       rhe(i,2,1)
!!!
!!!        Rqglobal((node3),1)=Rqglobal((node3),1)+Rq(i,3,1)+
!!!     1       rhe(i,3,1)
!!!
!!!        Rqglobal((node4),1)=Rqglobal((node4),1)+Rq(i,4,1)+
!!!     1       rhe(i,4,1)
!!!
!!!
!!!    
!!!      enddo
!!!!     Crank-Nicolson
!!!      intmeth=2
!!!!!!!!      write(*,*) '---------intmeth=2'
!!!
!!!      if (intmeth == 2) then
!!!      beta=0.75	
!!!      alpha=-0.25
!!!      
!!!!!!!!       write(*,*) '---------do i=1,nod=2'
!!!
!!!      do i=1,nod
!!!        do j=1,nod
!!!        LHS(i,j)=1./DtCondDiff*Cglobal(i,j)+(beta)*Kglobal(i,j)
!!!     1            *(1.+alpha)
!!!       
!!!        Tcoeff1(i,j)=1./DtCondDiff*Cglobal(i,j)+(beta)*Kglobal(i,j)*
!!!     1               (alpha)-(1.-beta)*(1.+alpha)*Kglobal(i,j)
!!!        Tcoeff2(i,j)=(1.-beta)*Kglobal(i,j)*alpha
!!!        enddo
!!!        Rold(i,1)=-(beta)*(alpha)*Rqold(i,1)+(1.-beta)*(1.+alpha)
!!!     1           *Rqold(i,1)
!!!        Rnew(i,1)=(beta)*Rqglobal(i,1)*(1.+alpha)
!!!        Rolder(i,1)=-(1.-beta)*(alpha)*Rqolder(i,1)
!!!      enddo
!!!       RHS2(1:nod,1)=matmul(Tcoeff1(1:nod,1:nod),Tinitd(1:nod,1))
!!!     1             +matmul(Tcoeff2(1:nod,1:nod),Tinitdold(1:nod,1))
!!!     2             +Rold(1:nod,1)+Rnew(1:nod,1)+Rolder(1:nod,1)
!!!!c     modifying C and K matrices for prescribed temperatures
!!!!!!!!!      write(*,*) nod,const,constraint(1:nod)
!!!!!      do i=1,nod
!!!        do j=1,const
!!!            nodeID=constraint(j)
!!!!!        if (i.eq.int(constraint(j))) then
!!!!!!!           write(*,*) i
!!!!!            do kk=1,nod
!!!                LHS(nodeID,1:nod)=0.
!!!!!!            enddo
!!!            LHS(nodeID,nodeID)=1.
!!!            RHS2(nodeID,1)=Tinitd(nodeID,1)
!!!!!        endif
!!!      enddo
!!!!!      enddo
!!!!      to apply the effect of thermal loading
!!!      call ThermalSetLoading(nod,RHS2,LHS,Tinitd)
!!!
!!!!c     modifications for periodic BC's
!!!!      do i=1,nod
!!!!        if (period(i).gt.0) then
!!!!            do kk=1,nod
!!!!                LHS((i),(kk))=0.
!!!!            enddo
!!!!            LHS((i),(i))=1.
!!!!            RHS2((i),1)=Tinitd((period(i)),1)
!!!!c            write(*,*) i,period(i)
!!!!        endif
!!!!      enddo
!!!
!!!!c      call gespd(LHS,RHS2,TdSol,nod,1)
!!!!c      call CG(LHS,RHS2,TdSol,Tinitd,nod)
!!!!c      call LU(LHS,RHS2,TdSol,nod,nbwth)
!!!
!!!      
!!!!!!!       write(*,*) '---------DGESV=2'
!!!       call DGESV(nod,1,LHS,nod,IPIV,RHS2,nod,INFO_flag)
!!!       
!!!       if(INFO_flag.eq.0)then
!!!             TdSol(1:nod,1)=RHS2(1:nod,1)     
!!!	     write(*,*)'solve successful',DtCondDiff
!!!       endif
!!!       
!!!       if(INFO_flag.ne.0)then    
!!!	     write(*,*)'solve unsuccessful',DtCondDiff
!!!	     call bye(2)
!!!       endif
!!!      endif
!!!
!!!
!!!!c------------------------------------------------     
!!!!c this section takes any nodal temperatures that
!!!!c are less than the initial temp of 293K and 
!!!!c reinitializes them to 293K. The temperatures
!!!!c possibly go below 293K because of numerical 
!!!!c instabilties.
!!!      if(nstep.gt.0)then
!!!	      if(TDflag==1) then   ! H2 difffusion
!!!              do i=1,nod
!!!                  if(TdSol(i,1).lt.0.0)then
!!!	                  TdSol(i,1)=0.0
!!!	              end if
!!!	              if(Tinitd(i,1).lt.0.0)then
!!!	                  Tinitd(i,1)=0.0
!!!	              end if
!!!              end do
!!!          elseif(TDflag==0) then     ! thermal conduction
!!!		      do i=1,nod
!!!                  if(TdSol(i,1).lt.293)then
!!!	                  TdSol(i,1)=293
!!!	              end if
!!!	              if(Tinitd(i,1).lt.293)then
!!!	                  Tinitd(i,1)=293
!!!	              end if
!!!              end do
!!!		  end if
!!!      end if
!!!!------------------------------------------------	 
!!!!      update only the nodes that has ICs not BCs
!!!!   ----------------------------------testing
!!!
!!!        do j=1,const
!!!            NodeID=int(constraint(j))
!!!            TdSol(NodeID,1) =Tinit(NodeID,1)
!!!        enddo
!!!      
!!!        call ThermalLdSetToSol(nod,TdSol)
!!!! !!!  ----------------------------------testing
!!!!!      write(*,*) 'TdSol(1:nod,1)'
!!!!!      write(*,*) Tinit(1:nod,1)
!!!!!!!------------------------------------------------	 
!!!        
!!!      Tinit(1:nod,1)=TdSol(1:nod,1)!/To
!!!      Rqold(1:nod,1)=Rqglobal(1:nod,1)
!!!      Rqolder(1:nod,1)=Rqold(1:nod,1)       ! save value at n step
!!!      Tinitdold(1:nod,1)=Tinitd(1:nod,1)
!!!      
!!!!!!  we can't update this here b/c Tinitd is not saved in any block but Tinit is saved   
!!!!!!!!!!      Tinitd(1:nod,1)=Tinit(1:nod,1)       ! this is the current step solution 
!!!      
!!!      if(TDflag==0) then !if thermal
!!!          do i=1,ele
!!!              Tele(i,1)=0.25*Tinit((connect(i,1)),1)+0.25*
!!!     1        Tinit((connect(i,2)),1)+0.25*Tinit((connect(i,3)),1)
!!!     2        +0.25*Tinit((connect(i,4)),1)
!!!          end do
!!!      elseif(TDflag==1) then !if diffusion
!!!	      do i=1,ele
!!!              hycon(i,1)=0.25*Tinit((connect(i,1)),1)+0.25*
!!!     1        Tinit((connect(i,2)),1)+0.25*Tinit((connect(i,3)),1)
!!!     2        +0.25*Tinit((connect(i,4)),1)
!!!          end do
!!!	  end if
!!!!!!!      write(*,*)'Tinit(1:nod,1)',Tinit(1:nod,1)
!!!!!!!      write(*,*)'hycon(1:nod,1)',hycon(1:ele,1)
!!!!!!!! ##################################################
!!!      if(iPrintOutputFlag .eq. 1) then
!!!        idtPrintOutput=mod(idt,idtimeStepsOutput(1))
!!!!!!        idtPrintOutput=mod(idt,1)
!!!        if(idtPrintOutput == 0) then ! there is a reminder  ex. 2200/1000 = 200 but 3000/1000 =0
!!!          ! if nstep is not divisible with iprint then do not print
!!!            if(TDflag==0) then !if thermal
!!!                !            write the temperature output for all elements
!!!                write(iFU_temperElem_out,26) Tele(1:ele,1)
!!!                !            write the temperature output for all nodes
!!!                write(iFU_temperNodes_out,26) TdSol(1:nod,1)
!!!            elseif(TDflag==1) then !if diffusion
!!!                !            write the temperature output for all nodes
!!!                write(iFU_DiffNodes_out,26) TdSol(1:nod,1)
!!!                !            write the temperature output for all elements
!!!                write(iFU_DiffElem_out,26) hycon(1:ele,1)
!!!            endif
!!!        endif
!!!        
!!!      endif
!!!      
!!!!!!!!!                --print output for node 32 for testing
!!!!!!!!      write(iFU_temperNodeID_out,26) TdSol(40,1) ! nodes ID
!!!!!!!!      write(iFU_temperElemID_out,26) TdSol(26,1) !elem ID
!!!!!!!!! ##################################################
!!!      enddo !do idt=1,idtimeStepsSol(1) !---- this is the refinement loop, chops dt by n-sections
!!!! ##################################################
!!!      if(iPrintOutputFlag .eq. 1) then
!!!          if(TDflag==0) then !if thermal
!!!            do i=1,ele
!!!                DijSijd(i)=DijSije(i)/rhocp(i)
!!!            enddo
!!!!            write the DijSije output for all elements
!!!             write(iFU_DijSijElem_out,26) DijSijd(1:ele)
!!!          endif
!!!          
!!!      endif
!!!      iPrintOutputFlag=0
!!!   26 format(e15.6)
!!!      
!!!      end
!!!      
!!!      
!!!    
