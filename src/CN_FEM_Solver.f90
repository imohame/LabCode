
!c     This program is for thermal finite element calculation
!c     Variable rRn, rRn_1, rRn_2, correspond to variable at 
!c     current, last, last two time step.
!c     Variable rTn, rTn_1, rTn_2, correspond to variable at 
!c     current, last, last two time step. 


subroutine CN_FEM_Solver(y, z, id, u,matp)!, ElemCount,NodeCount)
        use  mod_dtimeSpecs
        use CN_Objects_manager
        
!        implicit none
!        real,     intent(in) :: y(*), z(*),  u(*)!connect(4,*),
!        integer,  intent(in) :: id(2,*)!,ElemCount,NodeCount
!        real , intent(in) ::matp(*)
        dimension y(*), z(*), matp(*),  u(*),id(2,*)
      
            
      real*8 N(4,4),ShpFcnDeriv(4,2,4),coords(4,2)
      real*8 Jmatrix(2,2),Gamma(2,2),B(4,2,4)
      real*8 rhocp(ElemCountAct),Jacob(ElemCountAct,4)
      
      real*8 kappa(2,2)
      real*8 B1(2,4),B2(2,4),B3(2,4),B4(2,4)
      real*8 B1T(4,2),B2T(4,2),B3T(4,2),B4T(4,2)
      real*8 ke1a(2,4),ke2a(2,4),ke3a(2,4),ke4a(2,4)
      real*8 ke1b(4,4),ke2b(4,4),ke3b(4,4),ke4b(4,4)
      real*8 ke1(4,4),ke2(4,4),ke3(4,4),ke4(4,4)
	  real*8 ks1(4,4),ks2(4,4),ks3(4,4),ks4(4,4)
      
      real*8 Ce(ElemCountAct,4,4),ke(ElemCountAct,4,4)
      real*8 Cglobal(NodeCountAct,NodeCountAct), Kglobal(NodeCountAct,NodeCountAct)
      real*8  rRn(NodeCountAct),rTn(NodeCountAct)
      real*8  Rq(ElemCountAct,4) 
      real*8 thermalkd(ElemCountAct),DijSijd(ElemCountAct)!,rTn_1(NodeCountAct)
      
 
!!!!!!      real hycon     
      integer idt,i,idtPrintOutput,j,iNodes(4)
      real*8 rDiffCoeffElem,DtCondDiff
      integer node1,node2,node3,node4
      
      DtCondDiff=DtCurrent
      DtCondDiff=DtCurrent/idtimeStepsSol(1)
! ##################################################
      do idt=1,idtimeStepsSol(1) !---- this is the refinement loop, chops dt by n-sections
        Cglobal=0.0;        Kglobal=0.0
        DijSijd=0.0;        thermalkd=0.0
        rhocp=1.0
        Rq=0.0;
        rRn=0.0;rTn=0.0
        ke=0.0;  Ce=0.0
        !!!!!!!!!!!!!        pressnode=0.0
               
!!!!!!!!!! ##################################################               
        call CNGet_kd_RoCp_DijSij(thermalkd,DijSijd,rhocp)
        
!!!!!!!!!!!!        write(*,*) '#################################################'
!!!!!!!!!!!!        write(*,*) matp(1:100)
!!!!!!!!!!!!        write(*,*) int(matp(1:100))
!!!!!!!!!!!!        write(*,*) ElemCountAct
        CALL CNGet_press_kd(matp,thermalkd)
!!!  we can update this here b/c rTn is saved         
!!!!        rTn_1(1:NodeCountAct)=rTn(1:NodeCountAct)
      
        CALL CNCalShpFcnDeriv(ShpFcnDeriv)
! ##################################################
      
!c     starting the big element loop
        do i=1,ElemCountAct
!c      reading nodes associated with element i 
            iNodes(1:4)=iElemConnect(1:4,i)
            
            Coords(1:4,1)=y(iNodes(1:4))+u(id(1,iNodes(1:4)))          ! current configuration
            Coords(1:4,2)=z(iNodes(1:4))+u(id(2,iNodes(1:4)))          ! current configuration
            
!c     given shape function derivatives wrt isoparametric coords
!c     and coordinate matrix for the 4 Gauss points, find the B-matrix  
            do j=1,GaussIntPtsCount
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

            !c           thermal conductivity   / or diffusion   
            kappa(1,1)=thermalkd(i) 
            kappa(1,2)=0.
            kappa(2,1)=0.
            kappa(2,2)=thermalkd(i)

            ke1a=matmul(kappa,B1)
            ke2a=matmul(kappa,B2)
            ke3a=matmul(kappa,B3)
            ke4a=matmul(kappa,B4)

            ke1b=matmul(B1T,ke1a)
            ke2b=matmul(B2T,ke2a)        
            ke3b=matmul(B3T,ke3a)        
            ke4b=matmul(B4T,ke4a)  
            !c      Element k matrix for the 4 Gauss points            
            ke1(1:4,1:4)=ke1b(1:4,1:4)*Jacob(i,1)*GaussIntPtsWeight
            ke2(1:4,1:4)=ke2b(1:4,1:4)*Jacob(i,2)*GaussIntPtsWeight
            ke3(1:4,1:4)=ke3b(1:4,1:4)*Jacob(i,3)*GaussIntPtsWeight
            ke4(1:4,1:4)=ke4b(1:4,1:4)*Jacob(i,4)*GaussIntPtsWeight

            CALL CNGet_Rq_Ce(i,rhocp,Jacob,DijSijd,Rq,Ce)
            CALL CNGet_Ks(i,ke1,ke2,ke3,ke4,B1,B2,B3,B4,ks1,ks2,ks3,ks4)

            ke(i,1:4,1:4)=ke1(1:4,1:4)+ke2(1:4,1:4)+ke3(1:4,1:4)+ke4(1:4,1:4)!-ks1(1:4,1:4)-ks2(1:4,1:4) -ks3(1:4,1:4)-ks4(1:4,1:4)
!!!!!!            write(*,*)ke(i,1:4,1)
!!!!!!            write(*,*)ke(i,1:4,2)
!!!!!!            write(*,*)ke(i,1:4,3)
!!!!!!            write(*,*)ke(i,1:4,4)
!!!!!!c     Assembling global matrices
            Kglobal(iNodes(1:4),iNodes(1:4))=Kglobal(iNodes(1:4),iNodes(1:4))+ke(i,1:4,1:4)
            Cglobal(iNodes(1:4),iNodes(1:4))=Cglobal(iNodes(1:4),iNodes(1:4))+Ce(i,1:4,1:4)
            rRn(iNodes(1:4))=rRn(iNodes(1:4))+Rq(i,1:4)
            
        enddo
!!!!        do i=1,NodeCountAct
!!!!        write(2102,*)Kglobal(i,1:NodeCountAct)
!!!!        enddo

!!!!!!c     Assembling global matrices
!!!!!        do i=1,ElemCountAct
!!!!!!!!!            node1=iElemConnect(1,i)
!!!!!!!!!            node2=iElemConnect(2,i)
!!!!!!!!!            node3=iElemConnect(3,i)
!!!!!!!!!            node4=iElemConnect(4,i)
!!!!!
!!!!!            Kglobal(iElemConnect(1:4,i),iElemConnect(1:4,i))=Kglobal(iElemConnect(1:4,i),iElemConnect(1:4,i))+ke(i,1:4,1:4)
!!!!!            Cglobal(iElemConnect(1:4,i),iElemConnect(1:4,i))=Cglobal(iElemConnect(1:4,i),iElemConnect(1:4,i))+Ce(i,1:4,1:4)
!!!!!            rRn(iElemConnect(1:4,i))=rRn(iElemConnect(1:4,i))+Rq(i,1:4)
!!!!!        enddo
!!!!!     Crank-Nicolson 
        call CNmanagerGetSolution (DtCondDiff,Cglobal,Kglobal,rRn,rTn) 
!!!        write(*,*)rTn(13:15)
!!!  we can't update this here b/c rTn_1 is not saved in any block but rTn is saved   
!!!!!!!      rTn_1(1:NodeCount,1)=rTn(1:NodeCount,1)       ! this is the current step solution 
!!!!      write(*,*)'rTn(1:NodeCount,1)',rTn(1:NodeCount,1)
!!!!      write(*,*)'hycon(1:NodeCount,1)',hycon(1:ElemCount,1)
!!!!! ##################################################
        if(iPrintOutputFlag .eq. 1) then
            idtPrintOutput=mod(idt,idtimeStepsOutput(1))
            !!!        idtPrintOutput=mod(idt,1)
            if(idtPrintOutput == 0) then ! there is a reminder  ex. 2200/1000 = 200 but 3000/1000 =0
              ! if nstep is not divisible with iprint then do not print
                CALL CNWriteOutputFiles()
            endif

        endif      
!!!!!!                --print output for node 32 for testing
!!!!!      write(iFU_temperNodeID_out,26) TdSol(40,1) ! nodes ID
!!!!!      write(iFU_temperElemID_out,26) TdSol(26,1) !elem ID
        
!!!!!! ##################################################
        enddo !do idt=1,idtimeStepsSol(1) !---- this is the refinement loop, chops dt by n-sections
!!!!!! ##################################################
        if(iPrintOutputFlag .eq. 1) then
            CALL CNWriteOutputFilesDijSij(DijSijd,rhocp,thermalkd)   
            CALL CNmanager_WriteFiles21_2_3_4()
        endif
        
      
      
      
   26 format(e15.6)
      
      end
      
      
    
!!!            Kglobal(node1,iElemConnect(1:4,i))=Kglobal(node1,iElemConnect(1:4,i))+ke(i,1,1:4)
!!!            Kglobal(node2,iElemConnect(1:4,i))=Kglobal(node2,iElemConnect(1:4,i))+ke(i,2,1:4)
!!!            Kglobal(node3,iElemConnect(1:4,i))=Kglobal(node3,iElemConnect(1:4,i))+ke(i,3,1:4)
!!!            Kglobal(node4,iElemConnect(1:4,i))=Kglobal(node4,iElemConnect(1:4,i))+ke(i,4,1:4)
            
!!!!!            Kglobal(node1,node1)=Kglobal(node1,node1)+ke(i,1,1)
!!!!!            Kglobal(node1,node2)=Kglobal(node1,node2)+ke(i,1,2)
!!!!!            Kglobal(node1,node3)=Kglobal(node1,node3)+ke(i,1,3)
!!!!!            Kglobal(node1,node4)=Kglobal(node1,node4)+ke(i,1,4)
!!!!!            
!!!!!            Kglobal(node2,node1)=Kglobal(node2,node1)+ke(i,2,1)
!!!!!            Kglobal(node2,node2)=Kglobal(node2,node2)+ke(i,2,2)
!!!!!            Kglobal(node2,node3)=Kglobal(node2,node3)+ke(i,2,3)
!!!!!            Kglobal(node2,node4)=Kglobal(node2,node4)+ke(i,2,4)
!!!!!            
!!!!!            Kglobal(node3,node1)=Kglobal(node3,node1)+ke(i,3,1)
!!!!!            Kglobal(node3,node2)=Kglobal(node3,node2)+ke(i,3,2)
!!!!!            Kglobal(node3,node3)=Kglobal(node3,node3)+ke(i,3,3)
!!!!!            Kglobal(node3,node4)=Kglobal(node3,node4)+ke(i,3,4)
!!!!!            
!!!!!            Kglobal(node4,node1)=Kglobal(node4,node1)+ke(i,4,1)
!!!!!            Kglobal(node4,node2)=Kglobal(node4,node2)+ke(i,4,2)
!!!!!            Kglobal(node4,node3)=Kglobal(node4,node3)+ke(i,4,3)
!!!!!            Kglobal(node4,node4)=Kglobal(node4,node4)+ke(i,4,4)

      
      !!!!            Cglobal(node1,iElemConnect(1:4,i))=Cglobal(node1,iElemConnect(1:4,i))+Ce(i,1,1:4)
!!!!            Cglobal(node2,iElemConnect(1:4,i))=Cglobal(node2,iElemConnect(1:4,i))+Ce(i,2,1:4)
!!!!            Cglobal(node3,iElemConnect(1:4,i))=Cglobal(node3,iElemConnect(1:4,i))+Ce(i,3,1:4)
!!!!            Cglobal(node4,iElemConnect(1:4,i))=Cglobal(node4,iElemConnect(1:4,i))+Ce(i,4,1:4)
            
!!!!            Cglobal(node1,node1)=Cglobal(node1,node1)+Ce(i,1,1)
!!!!            Cglobal(node1,node2)=Cglobal(node1,node2)+Ce(i,1,2)
!!!!            Cglobal(node1,node3)=Cglobal(node1,node3)+Ce(i,1,3)
!!!!            Cglobal(node1,node4)=Cglobal(node1,node4)+Ce(i,1,4)
!!!!            
!!!!            Cglobal(node2,node1)=Cglobal(node2,node1)+Ce(i,2,1)
!!!!            Cglobal(node2,node2)=Cglobal(node2,node2)+Ce(i,2,2)
!!!!            Cglobal(node2,node3)=Cglobal(node2,node3)+Ce(i,2,3)
!!!!            Cglobal(node2,node4)=Cglobal(node2,node4)+Ce(i,2,4)
!!!!            
!!!!            Cglobal(node3,node1)=Cglobal(node3,node1)+Ce(i,3,1)
!!!!            Cglobal(node3,node2)=Cglobal(node3,node2)+Ce(i,3,2)
!!!!            Cglobal(node3,node3)=Cglobal(node3,node3)+Ce(i,3,3)
!!!!            Cglobal(node3,node4)=Cglobal(node3,node4)+Ce(i,3,4)
!!!!            
!!!!            Cglobal(node4,node1)=Cglobal(node4,node1)+Ce(i,4,1)
!!!!            Cglobal(node4,node2)=Cglobal(node4,node2)+Ce(i,4,2)
!!!!            Cglobal(node4,node3)=Cglobal(node4,node3)+Ce(i,4,3)
!!!!            Cglobal(node4,node4)=Cglobal(node4,node4)+Ce(i,4,4)

!!!            rRn(node1)=rRn(node1)+Rq(i,1)
!!!            rRn(node2)=rRn(node2)+Rq(i,2)
!!!            rRn(node3)=rRn(node3)+Rq(i,3)
!!!            rRn(node4)=rRn(node4)+Rq(i,4)

