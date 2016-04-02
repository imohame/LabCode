
! ###########################################################
module CN_Objects_manager
    use CN_Thermal
    use CN_Diffusion
    use CN_Consts
!    implicit none
    public 
    ! Declare a set geometric objects.
    type ( DiffusionClass ) :: CNObjDiffusion
    type ( ThermalClass )   :: CNObjThermal 
    real*8 , ALLOCATABLE :: rMatTH_k(:),rMatTH_h(:),rMatTH_Ro(:)
    real*8 , ALLOCATABLE :: rMatTH_cp(:),rMatTH_x(:),rMatDiff_D(:)
    real*8 , ALLOCATABLE :: gradpdata(:,:),sigalt(:,:),pressnode(:)
    integer , ALLOCATABLE :: iElemConnect(:,:)
    real*8 ::GaussIntPtsWeight=1.0d0
!!!        ---------------------------------------
    
    Interface CNSet_T_init;                       Module Procedure CNSet_T_init;                     End Interface   
    Interface CNCreateObjects;                    Module Procedure CNCreateObjects;                  End Interface     
    Interface CNCleanMem;                         Module Procedure CNCleanMem;                       End Interface     
    Interface CNSetBCdata;                        Module Procedure CNSetBCdata;                      End Interface    
    Interface CNTestPrint;                        Module Procedure CNTestPrint;                      End Interface    
    Interface CNSetMatProp;                       Module Procedure CNSetMatProp;                     End Interface   
    Interface CNmanagerGetSolution;               Module Procedure CNmanagerGetSolution;             End Interface   
    Interface CNmanagerGetElemTemp;               Module Procedure CNmanagerGetElemTemp;             End Interface   
    Interface CNmanagerSet_DijSij ;               Module Procedure CNmanagerSet_DijSij ;             End Interface   
    Interface CNmanagerCal_AdiabaticTemp ;        Module Procedure CNmanagerCal_AdiabaticTemp ;      End Interface   
    Interface CNmanager_CopyElemnt ;              Module Procedure CNmanager_CopyElemnt ;            End Interface   
    Interface CNmanager_CopyNode ;                Module Procedure CNmanager_CopyNode ;              End Interface   
    Interface CNmanager_Set_sigalt ;              Module Procedure CNmanager_Set_sigalt ;            End Interface   
    Interface CNmanager_Get_sigalt ;              Module Procedure CNmanager_Get_sigalt ;            End Interface   
    Interface CNmanager_WriteFiles21_2_3_4 ;      Module Procedure CNmanager_WriteFiles21_2_3_4 ;    End Interface   
    Interface CNmanager_Get_ElemConcentration ;   Module Procedure CNmanager_Get_ElemConcentration ; End Interface   
                
!##############################################################################        
!##############################################################################        
    contains 
        subroutine CNCreateObjects(mNodCountB,mElemCountB,mThermalFlag)
!            use CN_Consts
            use CN_Thermal
            use CN_Diffusion
            implicit none
            integer, intent(in) ::mNodCountB,mElemCountB,mThermalFlag
            integer ::ne
            type ( ThermalClass ) init_ThermalClass
            type ( DiffusionClass ) init_DiffusionClass
!!!!!           --calculate the constants
            call calconsts()
!!!!!        !      --- set the upper bound of the elem and nod count 30% extra
            NodCountB=1.5*mNodCountB
            ElemCountB=1.5*mElemCountB
            ElemCountAct=mElemCountB
            NodeCountAct=mNodCountB
!!-------------- allocate the mem for the material properties  
            ne=ElemCountB
            allocate(rMatTH_k(ne),rMatTH_h(ne),rMatTH_Ro(ne),rMatTH_cp(ne),rMatTH_x(ne),rMatDiff_D(ne))
            allocate(iElemConnect(4,ne))
            allocate(gradpdata(3,ne),sigalt(4,ne),pressnode(NodCountB))
            
            rMatTH_k=0.0
            rMatTH_h=0.0
            rMatTH_Ro=0.0
            rMatTH_cp=0.0
            rMatTH_x=0.0
            rMatDiff_D=0.0
            iElemConnect=0.0
            gradpdata=0.0
            sigalt=0.0
            pressnode=0.0
!!!!!           --create the objects according to the input options 
            select case (mThermalFlag)
                case (1) !-- thermal 
                        CNObjThermal=init_ThermalClass(NodCountB,ElemCountB)
                case (2) !-- diffusion
                        CNObjDiffusion=init_DiffusionClass(NodCountB,ElemCountB)
                        
                case (3) !-- both thermal and diffusion
                        CNObjThermal=init_ThermalClass(NodCountB,ElemCountB)
                        CNObjDiffusion=init_DiffusionClass(NodCountB,ElemCountB)
                case default !-- no thermal or diffusion
            end select
         end subroutine CNCreateObjects
!##############################################################################        
!##############################################################################        
        subroutine CNCleanMem()
            implicit none

            IF (ALLOCATED (rMatTH_k))              DEALLOCATE (rMatTH_k)
            IF (ALLOCATED (rMatTH_h))              DEALLOCATE (rMatTH_h)
            IF (ALLOCATED (rMatTH_Ro))             DEALLOCATE (rMatTH_Ro)
            IF (ALLOCATED (rMatTH_cp))             DEALLOCATE (rMatTH_cp)
            IF (ALLOCATED (rMatTH_x))              DEALLOCATE (rMatTH_x)
            IF (ALLOCATED (rMatDiff_D))            DEALLOCATE (rMatDiff_D)
            IF (ALLOCATED (iElemConnect))            DEALLOCATE (iElemConnect)
            IF (ALLOCATED (gradpdata))            DEALLOCATE (gradpdata)
            IF (ALLOCATED (sigalt))            DEALLOCATE (sigalt)
            IF (ALLOCATED (pressnode))            DEALLOCATE (pressnode)
            
            call CNObjThermal%CleanMem()
            call CNObjDiffusion%CleanMem()
            
         end subroutine CNCleanMem
!##############################################################################        
!##############################################################################        
        subroutine CNSetBCdata(NodeID,BC_thermal_t,BC_diffusion_c,BC_thermal_flag,BC_diffusion_flag)
!            use CN_Consts
            implicit none
            real,    intent(in) :: BC_thermal_t,BC_diffusion_c
            integer, intent(in) :: BC_thermal_flag,BC_diffusion_flag,NodeID

            call CNObjThermal%CNSetBCdataBase(NodeID,BC_thermal_t,BC_thermal_flag)
            call CNObjDiffusion%CNSetBCdataBase(NodeID,BC_diffusion_c,BC_diffusion_flag)

         end subroutine CNSetBCdata
!##############################################################################        
!##############################################################################        
        subroutine CNSet_T_init()
!!            use CN_Consts
!            implicit none
!!            integer,intent(in):: ElemCou-hntAct
!            real,intent(in):: rElemConnect(4,*)
!            real rElemConnect(4,*)
!            integer k02
!            common /main_block/ a
!            integer ix(4,*) 
!            real, dimension:: a(*)

            CALL CNObjThermal%CNBase_InitiateTn_1 (iElemconnect) 
            CALL CNObjDiffusion%CNBase_InitiateTn_1 (iElemconnect) 
            

         end subroutine CNSet_T_init
!##############################################################################        
!##############################################################################        
        subroutine CNTestPrint()
            implicit none
            integer i
            
            write(*,*) '------------------CNObjThermal%iSolutionActive',CNObjThermal%iSolutionActive
            if(CNObjThermal%iSolutionActive == 1)then
                do i=1,CNObjThermal%iBCcount
                    write(*,*) CNObjThermal%iBCNodeID(i)
                enddo
            endif
            write(*,*) '------------------CNObjDiffusion%iSolutionActive',CNObjDiffusion%iSolutionActive
            if(CNObjDiffusion%iSolutionActive == 1)then
                do i=1,CNObjDiffusion%iBCcount
                    write(*,*) CNObjDiffusion%iBCNodeID(i)
                enddo
            endif

         end subroutine CNTestPrint
!##############################################################################        
!##############################################################################        
        subroutine CNSetMatProp(nume,mat_type,MatCount,thermalki,thermalhi,thermalRoi,thermalcpi,thermalxi,thermalDi)
            implicit none
            integer ::nume,mat_type(nume)
            integer ::MatCount,MatID,i
            real ::thermalki(MatCount),thermalhi(MatCount),thermalRoi(MatCount)
            real ::thermalcpi(MatCount),thermalxi(MatCount),thermalDi(MatCount)

            do i = 1, ElemCountAct
                MatID = mat_type(i)         ! assign material properties to elements
                rMatTH_k(i)=thermalki(MatID)
                rMatTH_h(i)=thermalhi(MatID)
                rMatTH_Ro(i)=thermalRoi(MatID)
                rMatTH_cp(i)=thermalcpi(MatID)
                rMatTH_x(i)=thermalxi(MatID)
                rMatDiff_D(i)=thermalDi(MatID)
            enddo
            
         end subroutine CNSetMatProp
!##############################################################################        
!##############################################################################        
!!!        ---- This is used to calculate DijSije,thermalkd,DijSijd,rhocp 
!!!        ---------------------------------------
        subroutine CNGet_kd_RoCp_DijSij(thermalkd,DijSijd,rhocp)
            implicit none
!!!!!            integer , intent(in) ::nume
            real*8 , intent(out) ::thermalkd(ElemCountAct),DijSijd(ElemCountAct),rhocp(ElemCountAct)
            integer i
            
                if(CNObjThermal%iSolutionActive == 1) then !if thermal
                
                    thermalkd(1:ElemCountAct)=rMatTH_k(1:ElemCountAct)
                    !------------------------- this is the plastic work at each element as a heat source
                        DijSijd(1:ElemCountAct)=CNObjThermal%rThermal_DijSij(1:ElemCountAct)!    ! Q rate of heat generation, textbook by Cook
                    do i=1,ElemCountAct
                        rhocp(i)=(rMatTH_Ro(i)*rMatTH_cp(i)) 
                    enddo
                endif
            
         end subroutine CNGet_kd_RoCp_DijSij
!##############################################################################        
!##############################################################################        
!!!        ---- This is used to calculate DijSije,thermalkd,DijSijd,rhocp 
!!!        ---------------------------------------
        subroutine CNGet_press_kd(matp,thermalkd)
            use  mod_file_units
            use mod_DiffCoeffTable
            
!            implicit none
!!!            integer , intent(in) ::nume
!!!!!            real , intent(in) ::sigalt(4,nume)!,matp(*)
!!!!!            integer, intent(in) ::connect(ElemCountAct,4)
            dimension   ::matp(*)
            real*8 , intent(out) ::thermalkd(ElemCountAct)
            
            integer i,j
            integer nelenode(NodeCountAct)
            real*8 ::rElemHydroPress,rDiffCoeffElem
            
                if(CNObjDiffusion%iSolutionActive == 1) then !if diffusion                                
                    !------------------------- calculate the pressure at the nodes
                    pressnode=0.0
                    nelenode=0          ! number of elements connecting node i
!!!!!!!!!                    write(*,*)sigalt(1,1:ElemCountAct)
!!!!!!!!!                    write(*,*)sigalt(2,1:ElemCountAct)
!!!!!!!!!                    write(*,*)sigalt(3,1:ElemCountAct)
!!!!!!!!!                    write(*,*)iElemConnect(1:ElemCountAct,1)
!!!!!!!!!                    write(*,*)iElemConnect(1:ElemCountAct,2)
!!!!!!!!!                    write(*,*)iElemConnect(1:ElemCountAct,3)
!!!!!!!!!                    write(*,*)iElemConnect(1:ElemCountAct,4)
!!!!!!!!!                    write(*,*)pressnode(1:NodeCountAct)
!!!!!!!!!                    write(*,*)nelenode(1:NodeCountAct)
                    
                    do i=1,ElemCountAct
                        rElemHydroPress=(sigalt(1,i)+sigalt(2,i)+sigalt(3,i))/3.0
!!!!!                    write(*,*)sigalt(1,i),sigalt(2,i),sigalt(3,i)
                        do j=1,4
                            pressnode(iElemConnect(j,i))=pressnode(iElemConnect(j,i))+ rElemHydroPress
                             nelenode(iElemConnect(j,i))= nelenode(iElemConnect(j,i))+1
!!!!!                    write(*,*)iElemConnect(i,j),pressnode(iElemConnect(i,j))
!!!!!                    write(*,*)nelenode(iElemConnect(i,j))
                        end do
                    end do
                    
                    do i=1,ElemCountAct
                        if(nelenode(i)==0) then
                            write(*,*) 'warning a/0, check ThermalFEM.f'
                        end if
                        pressnode(i)=pressnode(i)/nelenode(i)   
                    end do
                    !------------------------- calculate the D which depends on temperature
                    thermalkd(1:ElemCountAct)=rMatDiff_D(1:ElemCountAct)
                    if(iDiffCoeffTableActive == 1)  then
                        do i=1,ElemCountAct
                            !------------- find the D from the table or equation
                            rDiffCoeffElem=0.0
    !!!!!!!!!!!!!!                        write(*,*) int(matp(1:ElemCountAct))
    !!!!!!!!!!!!!!                        write(*,*) '#######################################################################'
    !!!!!!!!!!!!!!                        write(*,*) CNObjDiffusion%rThermal_Diff_Ele
    !!!!!!!!!!!!!!                        write(*,*) i,ElemCountAct
                            call DiffCoeffTableGetD(i,rDiffCoeffElem,CNObjDiffusion%rThermal_Diff_Ele(i),matp)
                            thermalkd(i)=rDiffCoeffElem
                        enddo
                    endif
!!!!!                    write(iFU_DiffCoeff_out,26) thermalkd(1:ElemCountAct)
                endif
                
                26 format(e15.6)

         end subroutine CNGet_press_kd
!##############################################################################        
!##############################################################################        
        subroutine CNCalShpFcnDeriv(ShpFcnDeriv)
            implicit none
            real*8 , intent(out) :: ShpFcnDeriv(GaussIntPtsCount,2,4)
            integer i
            
                do i=1,GaussIntPtsCount

                    ShpFcnDeriv(i,1,1)=-(1-eta(i))/4.
                    ShpFcnDeriv(i,1,2)=(1-eta(i))/4.
                    ShpFcnDeriv(i,1,3)=(1+eta(i))/4.
                    ShpFcnDeriv(i,1,4)=-(1+eta(i))/4.
                    ShpFcnDeriv(i,2,1)=-(1-psi(i))/4.
                    ShpFcnDeriv(i,2,2)=-(1+psi(i))/4.
                    ShpFcnDeriv(i,2,3)=(1+psi(i))/4.
                    ShpFcnDeriv(i,2,4)=(1-psi(i))/4.
                enddo

         end subroutine CNCalShpFcnDeriv
!##############################################################################        
!##############################################################################        
        subroutine CNGet_Rq_Ce(ElemId,rhocp,Jacob,DijSijd,Rq,Ce)
            implicit none
            integer , intent(in) :: ElemId
            real*8 , intent(in) :: DijSijd(ElemCountAct),rhocp(ElemCountAct),Jacob(ElemCountAct,4)
            real*8 , intent(out) :: Rq(ElemCountAct,4,1) ,Ce(ElemCountAct,4,4)
            integer i,ii,jj
            real*8 Rq1(4,1),Rq2(4,1),Rq3(4,1),Rq4(4,1)
            real*8 ce1(4,4),ce2(4,4),ce3(4,4),ce4(4,4)
            !!--- no need to this check b/c this has to be done for both solvers
!!!            if(CNObjThermal%iSolutionActive == 1)then !if thermal
        !c     Element c matrix for the 4 Gauss points  
                i=ElemId
                do ii=1,4
                    do jj=1,4
                        ce1(ii,jj)=NtN1(ii,jj)*rhocp(i)*Jacob(i,1)*GaussIntPtsWeight
                        ce2(ii,jj)=NtN2(ii,jj)*rhocp(i)*Jacob(i,2)*GaussIntPtsWeight
                        ce3(ii,jj)=NtN3(ii,jj)*rhocp(i)*Jacob(i,3)*GaussIntPtsWeight
                        ce4(ii,jj)=NtN4(ii,jj)*rhocp(i)*Jacob(i,4)*GaussIntPtsWeight
                    enddo
                enddo
                Ce(i,1:4,1:4)=ce1(1:4,1:4)+ce2(1:4,1:4)+ce3(1:4,1:4)+ ce4(1:4,1:4)

        !c      Element Rq matrix given plastic work for each element 
                Rq1(1:4,1)=N1T(1:4,1)*DijSijd(i)*Jacob(i,1)*GaussIntPtsWeight
                Rq2(1:4,1)=N2T(1:4,1)*DijSijd(i)*Jacob(i,2)*GaussIntPtsWeight
                Rq3(1:4,1)=N3T(1:4,1)*DijSijd(i)*Jacob(i,3)*GaussIntPtsWeight
                Rq4(1:4,1)=N4T(1:4,1)*DijSijd(i)*Jacob(i,4)*GaussIntPtsWeight

                Rq(i,1:4,1)=Rq1(1:4,1)+Rq2(1:4,1)+Rq3(1:4,1)+Rq4(1:4,1)
!!!            endif

        end subroutine CNGet_Rq_Ce
!##############################################################################        
!##############################################################################        
        subroutine CNGet_Ks(ElemId,ke1,ke2,ke3,ke4,B1,B2,B3,B4,kes1,kes2,kes3,kes4)
            implicit none
            integer , intent(in) :: ElemId
!!!!!!!!            integer    , intent(in) :: connect(ElemCountAct,4)
            real*8 ,  intent(in) :: B1(2,4),B2(2,4),B3(2,4),B4(2,4)
            real*8 ,  intent(in) :: ke1(4,4),ke2(4,4),ke3(4,4),ke4(4,4)
            real*8 ,  intent(out):: kes1(4,4),kes2(4,4),kes3(4,4),kes4(4,4)
            integer i
            real*8 keP1(4,1),keP2(4,1),keP3(4,1),keP4(4,1)
            real*8 NodePress(4,1),gradp1a(2,1),gradp2a(2,1),gradp3a(2,1),gradp4a(2,1),gradp(2,1)
            real*8 gradpabs
            
            kes1=0.0;kes2=0.0;kes3=0.0;kes4=0.0
            
            if(CNObjDiffusion%iSolutionActive == 1)then !if diffusion
    !c      Ks matrix caused by stress assisted diffusion					
                i=ElemId
                NodePress(1,1)=pressnode(iElemConnect(1,i))
                NodePress(2,1)=pressnode(iElemConnect(2,i))
                NodePress(3,1)=pressnode(iElemConnect(3,i))
                NodePress(4,1)=pressnode(iElemConnect(4,i))

                keP1=matmul(ke1,NodePress)
                keP2=matmul(ke2,NodePress)
                keP3=matmul(ke3,NodePress)
                keP4=matmul(ke4,NodePress)

                kes1=matmul(keP1,N1s)*Vh_RT
                kes2=matmul(keP2,N2s)*Vh_RT
                kes3=matmul(keP3,N3s)*Vh_RT
                kes4=matmul(keP4,N4s)*Vh_RT

    !c       calculate pressure gradient
                gradp1a=matmul(B1, NodePress)
                gradp2a=matmul(B2, NodePress)
                gradp3a=matmul(B3, NodePress)
                gradp4a=matmul(B4, NodePress)
                gradp(1:2,1)=gradp1a(1:2,1)+gradp2a(1:2,1)+gradp3a(1:2,1)+gradp4a(1:2,1)
                gradp=gradp*0.25
                gradpabs=sqrt((gradp(1,1))**2+(gradp(2,1))**2)

                gradpdata(1,i)=gradp(1,1)
                gradpdata(2,i)=gradp(2,1)
                gradpdata(3,i)=gradpabs	
            endif
        end subroutine CNGet_Ks
!   ----------- this is a matlab code to check these cal
!!!!!!BtKB=[1 2 3 4; 5 6 7 8; 9 1 2 3]
!!!!!!j1=3
!!!!!!w=2
!!!!!!vhrt=5
!!!!!!p=[5;6;7;8]
!!!!!!N1=[ 6 5 4 3]
!!!!!!%-- with only BtKB
!!!!!!ks1a=BtKB*vhrt
!!!!!!ks1b=ks1a*p
!!!!!!ks1c=ks1b*N1
!!!!!!ks1_BtKB=ks1c*j1*w
!!!!!!%-- with Ke=BtKB * J1 * w
!!!!!!BtKB=BtKB*j1*w
!!!!!!ks1a=BtKB*vhrt
!!!!!!ks1b=ks1a*p
!!!!!!ks1c=ks1b*N1
!!!!!!ks1_Ke1=ks1c
!!!!!!
!!!!!!ks1_Ke1-ks1_BtKB
		
!##############################################################################        
!##############################################################################        
        subroutine CNmanagerGetSolution(DtCondDiff,Cglobal,Kglobal,rRn,rTn)
            implicit none
            real*8,     intent(in) :: DtCondDiff
            real*8,     intent(in) ::  Cglobal(NodeCountAct,NodeCountAct), Kglobal(NodeCountAct,NodeCountAct), rRn(NodeCountAct)
            real*8,     intent(out)::  rTn(NodeCountAct)
            
            
            if(CNObjDiffusion%iSolutionActive == 1)then !if diffusion
                CALL CNObjDiffusion%CN_Base_GetSolution (DtCondDiff,iElemConnect,Cglobal,Kglobal,rRn,rTn)
            endif
            if(CNObjThermal%iSolutionActive == 1)then !if thermal
                CALL CNObjThermal%CN_Base_GetSolution (DtCondDiff,iElemConnect,Cglobal,Kglobal,rRn,rTn)
            endif
            
         end subroutine CNmanagerGetSolution
!##############################################################################        
!##############################################################################        
         real function CNmanagerGetElemTemp(iElemID)
            implicit none
            integer,    intent(in) :: iElemID
            real retValue;
            !---- it works also when both are active
            if(CNObjDiffusion%iSolutionActive == 1)then !if diffusion
                retValue=0.0
            endif
            if(CNObjThermal%iSolutionActive == 1)then !if diffusion
                retValue = CNObjThermal%rThermal_Diff_Ele(iElemID)
            endif
            CNmanagerGetElemTemp=retValue
         end function CNmanagerGetElemTemp
!##############################################################################        
!##############################################################################        
         subroutine CNmanagerSet_DijSij(iElemID,mDijSij_val)
            implicit none
            integer,    intent(in) :: iElemID
            real,     intent(in) :: mDijSij_val
            
            if(CNObjThermal%iSolutionActive == 1)then !if thermal
                CNObjThermal%rThermal_DijSij(iElemID) = mDijSij_val*rMatTH_x(iElemID)
            endif
         end subroutine CNmanagerSet_DijSij
!##############################################################################        
!##############################################################################        
         subroutine CNmanagerCal_AdiabaticTemp(iElemID,mDijSij_val,mTemp)
            implicit none
            integer,    intent(in) :: iElemID
            real,     intent(in) :: mDijSij_val
            real,     intent(inout) :: mTemp
            real*8 dKa_RaoCp
            
            if(CNObjThermal%iSolutionActive == 1)then !if diffusion
!!!!!!!!!!!!!!!     dKa_RaoCp= thermalx(ink)/(thermalRo(ink)*thermalcp(ink))
                    dKa_RaoCp= rMatTH_x(iElemID)/(rMatTH_Ro(iElemID)*rMatTH_cp(iElemID))
                    
                 !This Line: Adiabatic Temp Update, Zikry: 1992-1993, pg 275, commented out for quasi-static
	             mTemp = mTemp + DtCurrent*mDijSij_val*dKa_RaoCp				
            
            endif
         end subroutine CNmanagerCal_AdiabaticTemp
!##############################################################################        
!##############################################################################        
         subroutine CNmanager_CopyNode(n,m)
            implicit none
            integer,    intent(in) :: n,m
                        
            
            if(CNObjThermal%iSolutionActive == 1)then !if thermal
                call CNObjThermal%CNBase_CopyNode (n,m) 
            endif
            if(CNObjDiffusion%iSolutionActive == 1)then !if diffusion
                call CNObjDiffusion%CNBase_CopyNode (n,m) 
            endif
         end subroutine CNmanager_CopyNode
!##############################################################################        
!##############################################################################        
         subroutine CNmanager_CopyElemnt(iElemId,numelt) 
            implicit none
            integer,    intent(in) :: iElemId,numelt
                        
             rMatTH_k(numelt+1)=rMatTH_k(iElemId)
             rMatTH_h(numelt+1)=rMatTH_h(iElemId)
             rMatTH_Ro(numelt+1)=rMatTH_Ro(iElemId)
             rMatTH_cp(numelt+1)=rMatTH_cp(iElemId)
             rMatTH_x(numelt+1)=rMatTH_x(iElemId)
             rMatDiff_D(numelt+1)=rMatDiff_D(iElemId)
             if(CNObjThermal%iSolutionActive == 1)then !if thermal
                CNObjThermal%rThermal_DijSij(numelt+1)=CNObjThermal%rThermal_DijSij(iElemId)
             endif
!!!!!!             thermalk(numelt+1)=thermalk(ele)
!!!!!!             thermalh(numelt+1)=thermalh(ele)
!!!!!!             thermalRo(numelt+1)=thermalRo(ele)
!!!!!!             thermalcp(numelt+1)=thermalcp(ele)
!!!!!!             thermalx(numelt+1)=thermalx(ele)             
!!!!!!             thermalD(numelt+1)=thermalD(ele)             
!!!!!!	           DijSije(numelt+1)=DijSije(ele)
            
         end subroutine CNmanager_CopyElemnt
!##############################################################################        
!##############################################################################        
         subroutine CNmanager_Get_ElemConcentration(iElemId,RetVal)
            implicit none
            integer,    intent(in) :: iElemId
            real,    intent(inout) :: RetVal  
            !!! column-major and pass by ref

            if(CNObjDiffusion%iSolutionActive == 1)then !if diffusion
                RetVal= CNObjDiffusion%rThermal_Diff_Ele(iElemId)
            endif
         end subroutine CNmanager_Get_ElemConcentration
!##############################################################################        
!##############################################################################        
         subroutine CNmanager_Get_sigalt(iElemId,RetVal)
            implicit none
            integer,    intent(in) :: iElemId
            real,    intent(inout) :: RetVal(4)  
            !!! column-major and pass by ref
            
            RetVal(1:4)= sigalt(1:4,iElemId)
!!!!    -------the following is nt needed because to get or set sigalt need not to be diffusion or thermal
!!!!!!!            if(CNObjDiffusion%iSolutionActive == 1)then !if diffusion
!!!!!!!                RetVal(1:4)= sigalt(1:4,iElemId)
!!!!!!!            endif
         end subroutine CNmanager_Get_sigalt
!##############################################################################        
!##############################################################################        
         subroutine CNmanager_Set_sigalt(iElemId,RetVal)
            implicit none
            integer,    intent(in) :: iElemId
            real,    intent(in) :: RetVal(4)                        
            
            sigalt(1:4,iElemId) = RetVal(1:4)
!!!!    -------the following is nt needed because to get or set sigalt need not to be diffusion or thermal
!!!!            if(CNObjDiffusion%iSolutionActive == 1)then !if diffusion
!!!!                sigalt(1:4,iElemId) = RetVal(1:4)
!!!!            endif
         end subroutine CNmanager_Set_sigalt
!##############################################################################        
!##############################################################################        
        subroutine CNWriteOutputFiles()
            use  mod_file_units
            implicit none
            !integer i
            
          ! if nstep is not divisible with iprint then do not print
            if(CNObjThermal%iSolutionActive == 1)then !if thermal
                !            write the temperature output for all nodes
                write(iFU_temperNodes_out,26) CNObjThermal%rTn_1(1:NodeCountAct)
                !            write the temperature output for all elements
                write(iFU_temperElem_out,26) CNObjThermal%rThermal_Diff_Ele(1:ElemCountAct)
            endif
            
            if(CNObjDiffusion%iSolutionActive == 1)then !if diffusion
                !            write the diff output for all nodes
                write(iFU_DiffNodes_out,26) CNObjDiffusion%rTn_1(1:NodeCountAct)
                !            write the diff output for all elements
                write(iFU_DiffElem_out,26) CNObjDiffusion%rThermal_Diff_Ele(1:ElemCountAct)
            endif

            26 format(e15.6)
            
         end subroutine CNWriteOutputFiles
!##############################################################################        
!##############################################################################        
        subroutine CNWriteOutputFilesDijSij(DijSijd,rhocp,thermalkd)
            use  mod_file_units
            implicit none
            real*8 ,  intent(inout) :: DijSijd(*),rhocp(*)
            real*8 , intent(out) ::thermalkd(ElemCountAct)
            
            integer i                        
            
!!!            if(CNObjDiffusion%iSolutionActive == 1)then !if diffusion
            if(CNObjThermal%iSolutionActive == 1)then !if diffusion
                do i=1,ElemCountAct
                    DijSijd(i)=DijSijd(i)/rhocp(i)
                enddo
!            write the DijSije output for all elements
                write(iFU_DijSijElem_out,26) DijSijd(1:ElemCountAct)
            endif
!!!            
            if(CNObjDiffusion%iSolutionActive == 1)then !if diffusion
                write(iFU_DiffCoeff_out,26) thermalkd(1:ElemCountAct)
            endif

            26 format(e15.6)
            
         end subroutine CNWriteOutputFilesDijSij
!##############################################################################        
!##############################################################################        
         subroutine CNmanager_WriteFiles21_2_3_4()
            implicit none
                        
            
            if(CNObjDiffusion%iSolutionActive == 1)then !if diffusion
!!!!!!!!            write(2101,*) hycon(ink,1)         !hydrogen_con.out
!!!!!!!!!!!c        print out pressure gradient
!!!!!!!!!!         write(2102,*) gradpdata(1,ink)
!!!!!!!!!!         write(2103,*) gradpdata(2,ink)
!!!!!!!!!!		 write(2104,*) gradpdata(3,ink)
!!!!                write(2101,26) CNObjDiffusion%rThermal_Diff_Ele(1:ElemCountAct)
                write(2102,26) gradpdata(1,1:ElemCountAct)
                write(2103,26) gradpdata(2,1:ElemCountAct)
                write(2104,26) gradpdata(3,1:ElemCountAct)
            endif
                26 format(e15.6)
         end subroutine CNmanager_WriteFiles21_2_3_4
!!!        ---------------------------------------
!!!        ---------------------------------------end module CN_Objects_manager
end module CN_Objects_manager
!###################################################################
!###################################################################
    subroutine CN_UpdateElemConnet(ix)
!!!!!!!!        implicit none
        use CN_Objects_manager
        dimension ix(4,*)
      iElemconnect(1:4,1:ElemCountAct)=ix(1:4,1:ElemCountAct)
!      do i=1,ElemCountAct
!          iElemconnect(1,i)=ix(1,i)
!          iElemconnect(2,i)=ix(2,i)
!          iElemconnect(3,i)=ix(3,i)
!          iElemconnect(4,i)=ix(4,i)
!      enddo
      
    end subroutine
!###################################################################
!###################################################################
