      subroutine ThermalFEM(y,z,ele,nod,maxint)
!c     This program is for thermal finite element calculation
!c     Variable Rqglobal, Rqold, Rqolder, correspond to variable at 
!c     current, last, last two time step.
!c     Variable Tinit, Tinitd, Tinitdold, correspond to variable at 
!c     current, last, last two time step. 
          use mod_dtimeSpecs
          
      parameter(nume = 40000)    
      integer ele,nod
      dimension y(nod), z(nod)
      integer maxint,kk,ii,bb,INFO_flag

      real*8 rhocp(ele),eta(maxint),psi(maxint)
      real*8 N(maxint,4), ShpFcnDeriv(maxint,2,4), coords(maxint,2)
      real*8 Jmatrix(2,2), Jacob(ele,maxint), Gamma(2,2),B(maxint,2,4)
      real*8 kappa(2,2),B1(2,4),B2(2,4),B3(2,4),B4(2,4)
      real*8 B1T(4,2),B2T(4,2),B3T(4,2),B4T(4,2)
      real*8 N1T(4,1),N2T(4,1),N3T(4,1),N4T(4,1)
      real*8 Nh(1,4),NhT(4,1),N1(1,4),N2(1,4),N3(1,4),N4(1,4)
      real*8 ci1(4,4),ci2(4,4),ci3(4,4),ci4(4,4),c(ele,4,4)
      real*8 ke1a(2,4),ke2a(2,4),ke3a(2,4),ke4a(2,4)
      real*8 ke1b(4,4),ke2b(4,4),ke3b(4,4),ke4b(4,4)
      real*8 ke1(4,4),ke2(4,4),ke3(4,4),ke4(4,4),ke(ele,4,4)
	  real*8 ks1a(4,4),ks2a(4,4),ks3a(4,4),ks4a(4,4)     ! for stress assisted diffusion
	  real*8 ks1b(4,1),ks2b(4,1),ks3b(4,1),ks4b(4,1)
	  real*8 ks1c(4,4),ks2c(4,4),ks3c(4,4),ks4c(4,4)
	  real*8 ks1(4,4),ks2(4,4),ks3(4,4),ks4(4,4)
      real*8 Rq1(4,1),Rq2(4,1),Rq3(4,1),Rq4(4,1), Rq(ele,4,1)
      real*8 Cglobal(nod,nod), Kglobal(nod,nod), Rqglobal(nod,1)
      real*8 RHS(nod,1),DtCondDiff !,TkelvinDelta(nod,1),TkelvinCurrent(nod,1)
      real*8 he(ele,4,4),rhe(ele,4,1), h12e(4,4),h23e(4,4),
     > h34e(4,4)
      real*8 h41e(4,4), rh12e(4,1),rh23e(4,1),rh34e(4,1),
     > rh41e(4,1)
      real*8 rh12egp(4,1),rh23egp(4,1),rh34egp(4,1),rh41egp(4,1)
      real*8 h12i(4,4),h23i(4,4),h34i(4,4),h41i(4,4),W
      real*8 thermalkd(ele),Tfld(ele),DijSijd(ele),
     > Tinitd(nod,1) !,hd(ele)
      real*8 beta,Tcoeff1(nod,nod),Rnew(nod,1), sigp(4,1)
      real*8 Rold(nod,1),J12,J23,J34,J41,Rolder(nod,1)
      real*8 RHS2(nod,1),LHS(nod,nod),TdSol(nod,1),alpha,
     > Tcoeff2(nod,nod)
	  real*8 gradp1a(2,1), gradp2a(2,1), gradp3a(2,1),gradp4a(2,1)
	  real*8 gradp(2,1), gradpabs
      real*8 gradslip1a(2,1), gradslip2a(2,1), gradslip3a(2,1)
	  real*8 gradslip4a(2,1), gradslipa(2,1)
	  
      integer  IPIV(nod)

      common/bk08/kprint,nstep,ite,ilimit,newstf
      common/bk11/cnwmk(2),iequit,iprint,isref,PrintOutputFlag
      
      common/bk26/dt,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
	  common/wblock2/  g_source, g_immob, g_minter, g_recov, b_v,
     1                             b_vvec(87),nmo,nim
      common/wblock8/  abc(573,nume,4), his(573,nume,4)
      common/WMLthermali/connect(nume,4)
      integer connect,nstep,intmeth
      
      common /WMLthermalBC/const,constraint(nume)
      integer const,constraint
      common /WMLthermal2/thermalk(nume),thermalh(nume),
     >      thermalRo(nume),thermalcp(nume),thermalKa(nume)
      common/WMLthermal/thermalflag,thermalconstraint(nume),
     1    Tinit(nume,1),Rqold(nume,1)

!!      common/WMLthermalperiod/period(40000)
      common/WMLthermalmatpoly/DijSij(nume),Tele(nume,1)
      real DijSij
!!      common/WMLthermalreordr/idth(nume),nbwth
      common/WMLthermalSolve/Rqolder(nume,1),
     1       Tinitdold(nume,1) !,dummy(40000,1)
      integer idth,nbwth,period
      common/fissn1/melemt,nnns,ntpe2,n2g,llls
	  common/couplinganalysis/ TDflag
	  common/hydrodiffusion/ hycon(nume,1)
	  common /propag/ sigalt(4,nume)
	  common /grad_pressure/ gradpdata(3,nume)
!c	  common /GND_loop/ gradslip(2,24,40000), rho_gnd(2,24,40000)
	  common/intgrt/   nintg
	  integer TDflag      ! 0 for thermal analysis, 1 for diffusion analysis
	  real hycon, sigalt, Vh, Temp, press(nume), pressnode(nume)
	  integer stressflag, nelenode(nume)
	  real gradpdata, slipnode(nume,24), slipp(4,1)
	  integer nelenode2(nume), nsi
	  real gradslip, rho_gnd, ly, lz
      write(*,*) 'in ThermalFEM'

      if(TDflag==1) then      ! Diffusion analysis
          stressflag=1        ! 0 for no stress assisted, 1 for yes
	  elseif(TDflag==0) then  ! thermal analysis
	      stressflag=0
	  end if
      Vh=2.0E-6           ! partial molar volume of hydrogen
      Temp=293.0
	  
  	  if(stressflag==0) then    ! 0 for no stress assisted
	      do i=1,4
		      do j=1,4
			      ks1(i,j)=0.0
				  ks2(i,j)=0.0
				  ks3(i,j)=0.0
				  ks4(i,j)=0.0
			  end do
		  end do
	  else if(stressflag==1) then   ! 0 for stress assisted
	      do i=1,nod
		      pressnode(i)=0.0
			  nelenode(i)=0          ! number of elements connecting node i
		  end do
	      do i=1,ele
		      press(i)=(sigalt(1,i)+sigalt(2,i)+sigalt(3,i))/3.0
			  do j=1,4
			      pressnode(connect(i,j))=pressnode(connect(i,j))+
     >              press(i)
				  nelenode(connect(i,j))=nelenode(connect(i,j))+1
			  end do
		  end do
		  do i=1,nod
		      if(nelenode(i)==0) then
			      write(*,*) 'warning a/0, check ThermalFEM.f'
			  end if
		      pressnode(i)=pressnode(i)/nelenode(i)   
		  end do
	  end if
	  
!c     Initializing temperatures based on the equation numbers
      if (nstep.eq.0) then
!!!!c     Setting initial vectors for intmeth 2 and 3      
          Tinitdold(1:nod,1)=0
          Rqolder(1:nod,1)=0
          Rqold(1:nod,1)=0
      end if

      
!c     Setting Element Material Parameters, removing "!"s would make the
!c     non-dimensional terms dimensional
      do i=1,ele
	     if (TDflag==0) then    ! thermal analysis
             rhocp(i)=(thermalRo(i)*thermalcp(i)) !/thermalKa(i)     !1./etae(i)!*E/To
!!!!             rhocp(i)=1./etae(i)!*E/To
		 else if (TDflag==1) then ! diffusion analysis
		     rhocp(i)=1.0
		 end if
         thermalkd(i)=thermalk(i)               !*E*cspeed*bv/To
!!!!!!!         hd(i)=thermalh(i)                      !*E*cspeed/To
!!!         Tfld(i)=Tfl(i)!*To
         DijSijd(i)=DijSij(i)!/bv*cspeed*E    ! Q rate of heat generation, textbook by Cook
      end do
!!!!!!   ----------------------------------testing
!!!!!      write(*,*)'Tinit',Tinit(1:nod,1)
!!!!!      write(*,*)'TkelvinCurrent',TkelvinCurrent(1:nod,1)
!!!!!      write(*,*)'Tinit(1:const)'
!!!!!      write(*,*)Tinit(1:nod)
!!!!!      write(*,*)'constraint(1:const)'
!!!!!      write(*,*)const,constraint(1:const)
!!!!!      write(*,*)const,constraint(1:nod)
!!!!!!   ----------------------------------testing
        do i=1,nod
          Tinitd(i,1)=Tinit(i,1)!*To
!!!!!!          TdSol(i,1)=TkelvinCurrent(i,1)!*To
         end do
      
!!!!!!!!c     Change in Temperature for intmeth=3 
!!!!!!!      do i=1,nod
!!!!!!!          TkelvinDelta(i,1)=Tinitd(i,1)-Tinitdold(i,1)
!!!!!!!      end do      
!!!!!!!      if(nstep.eq.0) then
!!!!!!!          do i=1,nod
!!!!!!!              TkelvinDelta(i,1)=0.
!!!!!!!          end do
!!!!!!!      end if
                  
      
     
!c      Defining integration points, corresponding with node 1
!c      at psi1,eta1 and going ccw
     
      psi(1)=-.57735
      eta(1)=-.57735
      psi(2)=.57735
      eta(2)=-.57735
      psi(3)=.57735
      eta(3)=.57735
      psi(4)=-.57735
      eta(4)=.57735
      
      
      call shapefunctions(psi,eta,N,ShpFcnDeriv,maxint)
      
!c     Setting the Weight Factors for Gauss integration, hardcoded below
!c     to maxint=4      
      if (maxint.eq.1) then
           W=4
      elseif (maxint.eq.4) then
           W=1
      endif 
! ##################################################
! ##################################################
      do idt=1,idtimeStepsSol(1) !---- thia is the refinement loop, chops dt by n-sections
        DtCondDiff=dt/idtimeStepsSol(i)
! ##################################################
!c     Zeroing out Global C, K matrix and Rq vector
      do i=1,nod
          do j=1,nod
              Cglobal(i,j)=0.000
              Kglobal(i,j)=0.
          end do
          Rqglobal(i,1)=0.
      end do
      
      
!c     starting the big element loop
      do i=1,ele
!c           thermal conductivity      
            kappa(1,1)=thermalkd(i) 
            kappa(1,2)=0.
            kappa(2,1)=0.
            kappa(2,2)=thermalkd(i)

!c      reading nodes associated with element i 
            node1=connect(i,1)
            node2=connect(i,2)
       	    node3=connect(i,3)
            node4=connect(i,4)
            Coords(1,1)=y(node1)!*bv          ! current configuration
            Coords(1,2)=z(node1)!*bv
            Coords(2,1)=y(node2)!*bv
            Coords(2,2)=z(node2)!*bv
            Coords(3,1)=y(node3)!*bv
            Coords(3,2)=z(node3)!*bv
            Coords(4,1)=y(node4)!*bv
            Coords(4,2)=z(node4)!*bv
!c     given shape function derivatives wrt isoparametric coords
!c     and coordinate matrix for the 4 Gauss points, find the B-matrix  
       do j=1,maxint
          Jmatrix=matmul(ShpFcnDeriv(j,1:2,1:4),Coords)
          Jacob(i,j)=Jmatrix(1,1)*Jmatrix(2,2)-Jmatrix(2,1)*Jmatrix(1,2)
          Gamma(1,1)=1/Jacob(i,j)*Jmatrix(2,2)
          Gamma(1,2)=-1/Jacob(i,j)*Jmatrix(1,2)
          Gamma(2,1)=-1/Jacob(i,j)*Jmatrix(2,1)
          Gamma(2,2)=1/Jacob(i,j)*Jmatrix(2,2)
          B(j,1:2,1:4)=matmul(Gamma,ShpFcnDeriv(j,1:2,1:4))
       enddo
        B1(1:2,1:4)=B(1,1:2,1:4)
        B2(1:2,1:4)=B(2,1:2,1:4)
        B3(1:2,1:4)=B(3,1:2,1:4)
        B4(1:2,1:4)=B(4,1:2,1:4)
        
        call matrixtrans(B1,B1T,2,4)
        call matrixtrans(B2,B2T,2,4)
        call matrixtrans(B3,B3T,2,4)
        call matrixtrans(B4,B4T,2,4)
        
        call matrixtrans(N(1,1:4),N1T,1,4)
        call matrixtrans(N(2,1:4),N2T,1,4)
        call matrixtrans(N(3,1:4),N3T,1,4)
        call matrixtrans(N(4,1:4),N4T,1,4)  

        N1(1,1:4)=N(1,1:4)
        N2(1,1:4)=N(2,1:4)
        N3(1,1:4)=N(3,1:4)
        N4(1,1:4)=N(4,1:4)      

 
        ci1=matmul(N1T,N1)
        ci2=matmul(N2T,N2)
        ci3=matmul(N3T,N3)
        ci4=matmul(N4T,N4)        

!c     Element c matrix for the 4 Gauss points        
        do ii=1,4
        do jj=1,4
        ci1(ii,jj)=ci1(ii,jj)*rhocp(i)*Jacob(i,1)*W
        ci2(ii,jj)=ci2(ii,jj)*rhocp(i)*Jacob(i,2)*W
        ci3(ii,jj)=ci3(ii,jj)*rhocp(i)*Jacob(i,3)*W
        ci4(ii,jj)=ci4(ii,jj)*rhocp(i)*Jacob(i,4)*W
        enddo
        enddo
        
        ke1a=matmul(kappa,B1)
        ke2a=matmul(kappa,B2)
        ke3a=matmul(kappa,B3)
        ke4a=matmul(kappa,B4)
        
        ke1b=matmul(B1T,ke1a)
        ke2b=matmul(B2T,ke2a)        
        ke3b=matmul(B3T,ke3a)        
        ke4b=matmul(B4T,ke4a)  
!c      Element k matrix for the 4 Gauss points            
        ke1(1:4,1:4)=ke1b(1:4,1:4)*Jacob(i,1)*W
        ke2(1:4,1:4)=ke2b(1:4,1:4)*Jacob(i,2)*W
        ke3(1:4,1:4)=ke3b(1:4,1:4)*Jacob(i,3)*W
        ke4(1:4,1:4)=ke4b(1:4,1:4)*Jacob(i,4)*W

!c      Ks matrix caused by stress assisted diffusion		
		if(stressflag==1) then
		    ks1a(1:4,1:4)=ke1b(1:4,1:4)*Vh/(8.3142*Temp)
			ks2a(1:4,1:4)=ke2b(1:4,1:4)*Vh/(8.3142*Temp)
			ks3a(1:4,1:4)=ke3b(1:4,1:4)*Vh/(8.3142*Temp)
			ks4a(1:4,1:4)=ke4b(1:4,1:4)*Vh/(8.3142*Temp)
			
			sigp(1,1)=pressnode(connect(i,1))
			sigp(2,1)=pressnode(connect(i,2))
			sigp(3,1)=pressnode(connect(i,3))
			sigp(4,1)=pressnode(connect(i,4))
			
			ks1b=matmul(ks1a,sigp)
			ks2b=matmul(ks2a,sigp)
			ks3b=matmul(ks3a,sigp)
			ks4b=matmul(ks4a,sigp)
			
			ks1c=matmul(ks1b,N1)
			ks2c=matmul(ks2b,N2)
			ks3c=matmul(ks3b,N3)
			ks4c=matmul(ks4b,N4)
			
			ks1(1:4,1:4)=ks1c(1:4,1:4)*Jacob(i,1)*W
			ks2(1:4,1:4)=ks2c(1:4,1:4)*Jacob(i,2)*W
			ks3(1:4,1:4)=ks3c(1:4,1:4)*Jacob(i,3)*W
			ks4(1:4,1:4)=ks4c(1:4,1:4)*Jacob(i,4)*W
			
!c       calculate pressure gradient
            gradp1a=matmul(B1, sigp)
			gradp2a=matmul(B2, sigp)
			gradp3a=matmul(B3, sigp)
			gradp4a=matmul(B4, sigp)
			gradp(1:2,1)=gradp1a(1:2,1)+gradp2a(1:2,1)
     >                   +gradp3a(1:2,1)+gradp4a(1:2,1)
            gradp(1,1)=gradp(1,1)/4.0
            gradp(2,1)=gradp(2,1)/4.0
            gradpabs=sqrt((gradp(1,1))**2+(gradp(2,1))**2)
			
            gradpdata(1,i)=gradp(1,1)
            gradpdata(2,i)=gradp(2,1)
            gradpdata(3,i)=gradpabs	
        end if					
		
!c      Element Rq matrix given plastic work for each element 
        Rq1(1:4,1)=N1T(1:4,1)*DijSijd(i)*W*Jacob(i,1)
        Rq2(1:4,1)=N2T(1:4,1)*DijSijd(i)*W*Jacob(i,2)
        Rq3(1:4,1)=N3T(1:4,1)*DijSijd(i)*W*Jacob(i,3)
        Rq4(1:4,1)=N4T(1:4,1)*DijSijd(i)*W*Jacob(i,4)
       
        ke(i,1:4,1:4)=ke1(1:4,1:4)+ke2(1:4,1:4)+ke3(1:4,1:4)+
     >            ke4(1:4,1:4)-ks1(1:4,1:4)-ks2(1:4,1:4)
     >            -ks3(1:4,1:4)-ks4(1:4,1:4)
        c(i,1:4,1:4)=ci1(1:4,1:4)+ci2(1:4,1:4)+ci3(1:4,1:4)+ 
     >            ci4(1:4,1:4)
        Rq(i,1:4,1)=Rq1(1:4,1)+Rq2(1:4,1)+Rq3(1:4,1)+Rq4(1:4,1)
      enddo


!c     Adding convection terms to the K matrix and R vector
      do i=1,ele
        do j=1,4
          do kk=1,4
             he(i,j,kk)=0.
          enddo
          rhe(i,j,1)=0.
        enddo
      enddo


!c     Assembling global matrices
      do i=1,ele
        node1=connect(i,1)
        node2=connect(i,2)
        node3=connect(i,3)
        node4=connect(i,4)
    
        Kglobal((node1),(node1))=Kglobal((node1),
     1       (node1))+ke(i,1,1)+he(i,1,1)

        Kglobal((node1),(node2))=Kglobal((node1),
     1       (node2))+ke(i,1,2)+he(i,1,2)

        Kglobal((node1),(node3))=Kglobal((node1),
     1       (node3))+ke(i,1,3)+he(i,1,3)

        Kglobal((node1),(node4))=Kglobal((node1),
     1       (node4))+ke(i,1,4)+he(i,1,4)

        Kglobal((node2),(node1))=Kglobal((node2),
     1       (node1))+ke(i,2,1)+he(i,2,1)

        Kglobal((node2),(node2))=Kglobal((node2),
     1       (node2))+ke(i,2,2)+he(i,2,2)

        Kglobal((node2),(node3))=Kglobal((node2),
     1       (node3))+ke(i,2,3)+he(i,2,3)

        Kglobal((node2),(node4))=Kglobal((node2),
     1       (node4))+ke(i,2,4)+he(i,2,4)

        Kglobal((node3),(node1))=Kglobal((node3),
     1       (node1))+ke(i,3,1)+he(i,3,1)

        Kglobal((node3),(node2))=Kglobal((node3),
     1       (node2))+ke(i,3,2)+he(i,3,2)

        Kglobal((node3),(node3))=Kglobal((node3),
     1       (node3))+ke(i,3,3)+he(i,3,3)

        Kglobal((node3),(node4))=Kglobal((node3),
     1       (node4))+ke(i,3,4)+he(i,3,4)

        Kglobal((node4),(node1))=Kglobal((node4),
     1       (node1))+ke(i,4,1)+he(i,4,1)

        Kglobal((node4),(node2))=Kglobal((node4),
     1       (node2))+ke(i,4,2)+he(i,4,2)

        Kglobal((node4),(node3))=Kglobal((node4),
     1       (node3))+ke(i,4,3)+he(i,4,3)

        Kglobal((node4),(node4))=Kglobal((node4),
     1       (node4))+ke(i,4,4)+he(i,4,4)

    
        Cglobal((node1),(node1))=Cglobal((node1),
     1       (node1))+c(i,1,1)
        Cglobal((node1),(node2))=Cglobal((node1),
     1       (node2))+c(i,1,2)
        Cglobal((node1),(node3))=Cglobal((node1),
     1       (node3))+c(i,1,3)
        Cglobal((node1),(node4))=Cglobal((node1),
     1       (node4))+c(i,1,4)
        Cglobal((node2),(node1))=Cglobal((node2),
     1       (node1))+c(i,2,1)
        Cglobal((node2),(node2))=Cglobal((node2),
     1       (node2))+c(i,2,2)
        Cglobal((node2),(node3))=Cglobal((node2),
     1       (node3))+c(i,2,3)
        Cglobal((node2),(node4))=Cglobal((node2),
     1       (node4))+c(i,2,4)
        Cglobal((node3),(node1))=Cglobal((node3),
     1       (node1))+c(i,3,1)
        Cglobal((node3),(node2))=Cglobal((node3),
     1       (node2))+c(i,3,2)
        Cglobal((node3),(node3))=Cglobal((node3),
     1       (node3))+c(i,3,3)
        Cglobal((node3),(node4))=Cglobal((node3),
     1       (node4))+c(i,3,4)
        Cglobal((node4),(node1))=Cglobal((node4),
     1       (node1))+c(i,4,1)
        Cglobal((node4),(node2))=Cglobal((node4),
     1       (node2))+c(i,4,2)
        Cglobal((node4),(node3))=Cglobal((node4),
     1       (node3))+c(i,4,3)
        Cglobal((node4),(node4))=Cglobal((node4),
     1       (node4))+c(i,4,4)
 
        Rqglobal((node1),1)=Rqglobal((node1),1)+Rq(i,1,1)+
     1       rhe(i,1,1)

        Rqglobal((node2),1)=Rqglobal((node2),1)+Rq(i,2,1)+
     1       rhe(i,2,1)

        Rqglobal((node3),1)=Rqglobal((node3),1)+Rq(i,3,1)+
     1       rhe(i,3,1)

        Rqglobal((node4),1)=Rqglobal((node4),1)+Rq(i,4,1)+
     1       rhe(i,4,1)


    
      enddo
            
      intmeth=2
      
!!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!!!!!c     Direct Integration from Cook, et. al.  
!!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!!!!      if (intmeth.eq.1) then
!!!!!!c     beta=1 Backward Difference
!!!!!!c     beta=2./3. Galerkin
!!!!!!c     beta=0.5 Crank-Nicholson
!!!!!!c     beta=0 Forward Difference
!!!!!      beta=0.5
!!!!!      
!!!!!
!!!!!      do i=1,nod
!!!!!        do j=1,nod
!!!!!        LHS(i,j)=1./DtCondDiff*Cglobal(i,j)+beta*Kglobal(i,j)
!!!!!        Tcoeff1(i,j)=1./DtCondDiff*Cglobal(i,j)-(1.-beta)*Kglobal(i,j)
!!!!!        enddo
!!!!!        Rold(i,1)=(1.-beta)*Rqold(i,1)
!!!!!        Rnew(i,1)=beta*Rqglobal(i,1)
!!!!!      enddo
!!!!!      RHS2=matmul(Tcoeff1(1:nod,1:nod),Tinitd)+Rold+Rnew
!!!!!
!!!!!!c     modifying C and K matrices for prescribed temperatures
!!!!!      do i=1,nod
!!!!!        do j=1,const
!!!!!        if (i.eq.int(constraint(j))) then
!!!!!!c        write(*,*) i
!!!!!            do kk=1,nod
!!!!!                LHS((i),(kk))=0.
!!!!!            enddo
!!!!!            LHS((i),(i))=1.
!!!!!            RHS2((i),1)=Tinitd((i),1)           
!!!!!        endif
!!!!!      enddo
!!!!!      enddo
!!!!!
!!!!!!c     modifications for periodic BC's
!!!!!!      do i=1,nod
!!!!!!        if (period(i).gt.0) then
!!!!!!            do kk=1,nod
!!!!!!                LHS((i),(kk))=0.
!!!!!!            enddo
!!!!!!            LHS((i),(i))=1.
!!!!!!            RHS2((i),1)=Tinitd((period(i)),1)
!!!!!!c            write(*,*) i,period(i)
!!!!!!        endif
!!!!!!      enddo
!!!!!
!!!!!!c      call gespd(LHS,RHS2,TdSol,nod,1)
!!!!!!c      call CG(LHS,RHS2,TdSol,Tinitd,nod)
!!!!!      call LU(LHS,RHS2,TdSol,nod,nbwth)
!!!!!
!!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!!!!!c     Method suggested in Cook, et. al. 12.4-14            
!!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!!!!      elseif (intmeth.eq.2) then
      if(intmeth == 2)then
        beta=0.75	
        alpha=-0.25
      

      do i=1,nod
        do j=1,nod
        LHS(i,j)=1./DtCondDiff*Cglobal(i,j)+(beta)*Kglobal(i,j)
     1            *(1.+alpha)
       
        Tcoeff1(i,j)=1./DtCondDiff*Cglobal(i,j)+(beta)*Kglobal(i,j)*
     1               (alpha)-(1.-beta)*(1.+alpha)*Kglobal(i,j)
        Tcoeff2(i,j)=(1.-beta)*Kglobal(i,j)*alpha
        enddo
        Rold(i,1)=-(beta)*(alpha)*Rqold(i,1)+(1.-beta)*(1.+alpha)
     1           *Rqold(i,1)
        Rnew(i,1)=(beta)*Rqglobal(i,1)*(1.+alpha)
        Rolder(i,1)=-(1.-beta)*(alpha)*Rqolder(i,1)
      enddo
       RHS2(1:nod,1)=matmul(Tcoeff1(1:nod,1:nod),Tinitd(1:nod,1))
     1             +matmul(Tcoeff2(1:nod,1:nod),Tinitdold(1:nod,1))
     2             +Rold(1:nod,1)+Rnew(1:nod,1)+Rolder(1:nod,1)
!c     modifying C and K matrices for prescribed temperatures
!!!!!!      write(*,*) nod,const,constraint(1:nod)
!!      do i=1,nod
        do j=1,const
            nodeID=constraint(j)
!!        if (i.eq.int(constraint(j))) then
!!!!           write(*,*) i
!!            do kk=1,nod
                LHS(nodeID,1:nod)=0.
!!!            enddo
            LHS(nodeID,nodeID)=1.
            RHS2(nodeID,1)=Tinitd(nodeID,1)
!!        endif
      enddo
!!      enddo
!      to apply the effect of thermal loading
      call ThermalSetLoading(nod,RHS2,LHS,Tinitd)

!!!!!c     modifications for periodic BC's
!!!!!      do i=1,nod
!!!!!        if (period(i).gt.0) then
!!!!!            do kk=1,nod
!!!!!                LHS((i),(kk))=0.
!!!!!            enddo
!!!!!            LHS((i),(i))=1.
!!!!!            RHS2((i),1)=Tinitd((period(i)),1)
!!!!!c            write(*,*) i,period(i)
!!!!!        endif
!!!!!      enddo
!!!!
!!!!!c      call gespd(LHS,RHS2,TdSol,nod,1)
!!!!!c      call CG(LHS,RHS2,TdSol,Tinitd,nod)
!!!!!c      call LU(LHS,RHS2,TdSol,nod,nbwth)

      
       call DGESV(nod,1,LHS,nod,IPIV,RHS2,nod,INFO_flag)
       
       if(INFO_flag == 0)then
             TdSol(1:nod,1)=RHS2(1:nod,1)     
	     write(*,*)'solve successful'
       endif
       
       if(INFO_flag /= 0)then    
	     write(*,*)'solve unsuccessful'
	     call bye(2)
       endif
      endif

!!!!!!!!!!!!	   
!!!!!!!!!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!!!!!!!!!!!!c     Direct integration as in Ansys manual
!!!!!!!!!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!!!!!!!!!!!      elseif (intmeth.eq.3) then
!!!!!!!!!!!!      beta=0.5
!!!!!!!!!!!!
!!!!!!!!!!!!
!!!!!!!!!!!!      do i=1,nod
!!!!!!!!!!!!        do j=1,nod
!!!!!!!!!!!!        LHS(i,j)=1./DtCondDiff/beta*Cglobal(i,j)+Kglobal(i,j)
!!!!!!!!!!!!        Tcoeff1(i,j)=1./DtCondDiff/beta*Cglobal(i,j)
!!!!!!!!!!!!        Tcoeff2(i,j)=(1.-beta)/beta*Cglobal(i,j)/DtCondDiff
!!!!!!!!!!!!        enddo
!!!!!!!!!!!!        Rnew(i,1)=Rqglobal(i,1)
!!!!!!!!!!!!      enddo
!!!!!!!!!!!!      RHS2=matmul(Tcoeff1,Tinitd)+
!!!!!!!!!!!!     1       matmul(Tcoeff2,TkelvinDelta)+Rnew
!!!!!!!!!!!!
!!!!!!!!!!!!!c     modifying C and K matrices for prescribed temperatures
!!!!!!!!!!!!      do i=1,nod
!!!!!!!!!!!!        do j=1,const
!!!!!!!!!!!!        if (i.eq.int(constraint(j))) then
!!!!!!!!!!!!!!!        write(*,*) i
!!!!!!!!!!!!            do kk=1,nod
!!!!!!!!!!!!                LHS((i),(kk))=0.
!!!!!!!!!!!!            enddo
!!!!!!!!!!!!            LHS((i),(i))=1.
!!!!!!!!!!!!            RHS2((i),1)=Tinitd((i),1)
!!!!!!!!!!!!            
!!!!!!!!!!!!        endif
!!!!!!!!!!!!      enddo
!!!!!!!!!!!!      enddo
!!!!!!!!!!!!
!!!!!!!!!!!!!c     modifications for periodic BC's
!!!!!!!!!!!!!      do i=1,nod
!!!!!!!!!!!!!        if (period(i).gt.0) then
!!!!!!!!!!!!!            do kk=1,nod
!!!!!!!!!!!!!                LHS((i),(kk))=0.
!!!!!!!!!!!!!            enddo
!!!!!!!!!!!!!            LHS((i),(i))=1.
!!!!!!!!!!!!!            RHS2((i),1)=Tinitd((period(i)),1)
!!!!!!!!!!!!!c            write(*,*) i,period(i)
!!!!!!!!!!!!!        endif
!!!!!!!!!!!!!      enddo
!!!!!!!!!!!!
!!!!!!!!!!!!!c      call gespd(LHS,RHS2,TdSol,nod,1)
!!!!!!!!!!!!!c      call CG(LHS,RHS2,TdSol,Tinitd,nod)
!!!!!!!!!!!!      call LU(LHS,RHS2,TdSol,nod,nbwth)
!!!!!!!!!!!!      do i=1,nod
!!!!!!!!!!!!      Tinitdold(i,1)=Tinitd(i,1)
!!!!!!!!!!!!      enddo
!!!!!!!!!!!!
!!!!!!!!!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!!!!!!!!!!!!     Explicit         
!!!!!!!!!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!!!!!!!!!!!      elseif (intmeth.eq.0) then
!!!!!!!!!!!!      
!!!!!!!!!!!!      do i=1,nod
!!!!!!!!!!!!      RHS(i,1)=0.
!!!!!!!!!!!!      enddo
!!!!!!!!!!!!
!!!!!!!!!!!!      RHS=matmul(Kglobal,Tinitd)
!!!!!!!!!!!!
!!!!!!!!!!!!      RHS(1:nod,1)=-RHS(1:nod,1)+Rqglobal(1:nod,1)
!!!!!!!!!!!!
!!!!!!!!!!!!      
!!!!!!!!!!!!!     modifying C and K matrices for prescribed temperatures
!!!!!!!!!!!!      do i=1,nod
!!!!!!!!!!!!        do j=1,const
!!!!!!!!!!!!        if (i.eq.int(constraint(j))) then
!!!!!!!!!!!!            do kk=1,nod
!!!!!!!!!!!!                Cglobal((i),(kk))=0.
!!!!!!!!!!!!                Cglobal((kk),(i))=0.
!!!!!!!!!!!!                Kglobal((i),(kk))=0.
!!!!!!!!!!!!            enddo
!!!!!!!!!!!!            Cglobal((i),(i))=1.
!!!!!!!!!!!!            Rqglobal((i),1)=0.
!!!!!!!!!!!!            RHS((i),1)=0.
!!!!!!!!!!!!        endif
!!!!!!!!!!!!      enddo
!!!!!!!!!!!!      enddo
!!!!!!!!!!!!      do i=1,nod
!!!!!!!!!!!!      RHS(i,1)=RHS(i,1)*DtCondDiff
!!!!!!!!!!!!      do j=1,nod
!!!!!!!!!!!!      Cglobal(i,j)=Cglobal(i,j)*DtCondDiff
!!!!!!!!!!!!      enddo
!!!!!!!!!!!!!c      RHS(1:nod,1)=-matmul(Cglobal(1:nod,1:nod),TkelvinDelta(1:nod,1))
!!!!!!!!!!!!!c     1          +RHS(1:nod,1)
!!!!!!!!!!!!      enddo
!!!!!!!!!!!!
!!!!!!!!!!!!      call gespd(Cglobal,RHS,TkelvinDelta,nod,1)
!!!!!!!!!!!!
!!!!!!!!!!!!      TdSol(1:nod,1)=Tinitd(1:nod,1)+TkelvinDelta(1:nod,1)
!!!!!!!!!!!!     1        *DtCondDiff
!!!!!!!!!!!!      endif

!c------------------------------------------------     
!c this section takes any nodal temperatures that
!c are less than the initial temp of 293K and 
!c reinitializes them to 293K. The temperatures
!c possibly go below 293K because of numerical 
!c instabilties.
      if(nstep > 0)then
	      if(TDflag==1) then   ! H2 difffusion
              do i=1,nod
                  if(TdSol(i,1) < 0.0)then
	                  TdSol(i,1)=0.0
	              end if
	              if(Tinitd(i,1) < 0.0)then
	                  Tinitd(i,1)=0.0
	              end if
              end do
          elseif(TDflag==0) then     ! thermal conduction
		      do i=1,nod
                  if(TdSol(i,1) < 293)then
	                  TdSol(i,1)=293
	              end if
	              if(Tinitd(i,1) < 293)then
	                  Tinitd(i,1)=293
	              end if
              end do
		  end if
      end if
!------------------------------------------------	 
!      update only the nodes that has ICs not BCs
!   ----------------------------------testing
!!!      write(*,*) '------------------------------------------'
!!!      write(*,*) 'TdSol(1:nod,1)'
!!!      write(*,*) TdSol(1:nod,1)
!!!      write(*,*) 'Tinit(1:nod,1)'
!!!      write(*,*) Tinit(1:nod,1)
!!!      write(*,*) 'const,constraint(1:const)'
!!!      write(*,*) const,constraint(1:const)
      
        do j=1,const
            NodeID=int(constraint(j))
            TdSol(NodeID,1) =Tinit(NodeID,1)
        enddo
      
        call ThermalLdSetToSol(nod,TdSol)
! !!!  ----------------------------------testing
!!      write(*,*) 'TdSol(1:nod,1)'
!!      write(*,*) TdSol(1:nod,1)
!!      write(*,*) 'Tinit(1:nod,1)'
!!      write(*,*) Tinit(1:nod,1)
!!!!------------------------------------------------	 
        
      Tinit(1:nod,1)=TdSol(1:nod,1)!/To
      
      Rqold(1:nod,1)=Rqglobal(1:nod,1)
      Rqolder(1:nod,1)=Rqold(1:nod,1)       ! save value at n step
      Tinitdold(1:nod,1)=Tinitd(1:nod,1)   ! this is the solution of the previous step not the current
      ! Tinitd(1:nod,1)=Tinit(1:nod,1)       ! this is the current step solution 
!      write the temp sol file
!c      do i=1,nod
!c      write(*,*) 'Tinit',Tinit(i,1),'Rqold',Rqold(i,1),'Tinitd',
!c     1     Tinitd(i,1),'Tinitdold',Tinitdold(i,1)
!c      enddo
      if(TDflag==0) then !if thermal
          do i=1,ele
              Tele(i,1)=0.25*Tinit((connect(i,1)),1)+0.25*
     1        Tinit((connect(i,2)),1)+0.25*Tinit((connect(i,3)),1)
     2        +0.25*Tinit((connect(i,4)),1)
          end do
      elseif(TDflag==1) then !if diffusion
	      do i=1,ele
              hycon(i,1)=0.25*Tinit((connect(i,1)),1)+0.25*
     1        Tinit((connect(i,2)),1)+0.25*Tinit((connect(i,3)),1)
     2        +0.25*Tinit((connect(i,4)),1)
          end do
	  end if
      
      if(PrintOutputFlag .eq. 1) then
        idtPrintOutput=mod(idtimeStepsSol(1),idtimeStepsOutput(1))
        if(idtPrintOutput == 0) then ! there is a reminder  ex. 2200/1000 = 200 but 3000/1000 =0
          ! if nstep is not divisible with iprint then do not print
          
!            write the temperature output for all nodes
            write(iFU_temperNodes_out,26) TdSol(1:nod,1)
!            write the temperature output for all elements
            write(iFU_temperElem_out,26) Tele(1:ele,1)
        endif
        PrintOutputFlag=0
      endif
      
   26 format(e15.6)
! ##################################################
      enddo
!!!!!c      open(1,file='temp.dat',status='unknown')
!!!!!c      write(1,*) ' TITLE="Values from the file: ggg.f"'
!!!!!c      write(1,*) ' VARIABLES="x","y","ggg.f"'
!!!!!c      write(1,*) 'ZONE T="Scalar Field",N=',nod,',E=',ele,
!!!!!c     1 ',ET=QUADRILATERAL,'
!!!!!c      write(1,*) ' F=FEPOINT'
!!!!!c      do i=1,nod
!!!!!c      write(1,*) y(i),z(i),TdSol(idth(i),1)!*293.
!!!!!c      enddo
!!!!!c      do i=1,ele
!!!!!c      write(1,*) connect(i,1),connect(i,2),connect(i,3),connect(i,4)
!!!!!c      enddo
!!!!!c      write(*,*) ele12(1:ele12size)
!!!!!c      write(*,*) ele23(1:ele23size)
!!!!!c      write(*,*) ele34(1:ele34size)
!!!!!c      write(*,*) ele41(1:ele41size)
!!!!!c      do i=1,ele
!!!!!c      write(*,*) thermalkd(i),hd(i),etae(i),DijSijd(i),TkelvinDelta(i,1),Jacob
!!!!!c     1   (i,1),B1(1,1),DtCondDiff
!!!!!c      enddo
      return
      end
      
      
    
