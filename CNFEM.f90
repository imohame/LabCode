
module CN_Consts
    use mod_parameters
    implicit none
    public
    integer:: ElemCountAct,NodeCountAct
    integer:: NodCountB,ElemCountB,iBCcount
    real*8 :: pid=3.14159265358979  !4.d0*atan(1.d0) ! a circle constant
    real*8 :: alpha=-.25d0
    real*8 :: beta=0.75d0
    real*8 :: C1,C2,C3,C4,C5
    real*8 :: sqr3_3=sqrt(3.d0)/3.d0
    real*8 :: Vh=2.0E-6 ! partial molar volume of hydrogen
    real*8 :: Temprd=293.0d0
    real*8 :: R=8.3142d0
    real*8 :: Vh_RT
    real*8 :: N1T(4,1),N2T(4,1),N3T(4,1),N4T(4,1)
    real*8 :: N1s(1,4),N2s(1,4),N3s(1,4),N4s(1,4)
    real*8 :: NtN1(4,4),NtN2(4,4),NtN3(4,4),NtN4(4,4)
    real*8 :: psi(4),eta(4)
    integer:: iPrintOutputFlag=0
    real*8 :: DtCurrent=1
    integer ::GaussIntPtsCount=4
    
    Interface CalConsts;                  Module Procedure CalConsts ;                   End Interface
    Interface CNCalcElemAverageFromNodes; Module Procedure CNCalcElemAverageFromNodes;   End Interface
!    ---------------------------
    contains
        subroutine CalConsts
            !--- to see only this variable from the entire module
!!!!!!!            use CN_Objects_manager ,only: GaussIntPtsCount
            implicit none
            real*8 ::N(GaussIntPtsCount,4)
            integer i
            
            C1=beta*(1+alpha)
            C2=beta*alpha
            C3=(1-beta)*(1+alpha)
            C4=(1-beta)*alpha
            C5=(1+alpha)*alpha        
            Vh_RT=Vh/(R*Tempr)        
            psi(1:4)=(/-sqr3_3, sqr3_3, sqr3_3, -sqr3_3/)
            eta(1:4)=(/-sqr3_3, -sqr3_3, sqr3_3, sqr3_3/)
            
            do i=1,GaussIntPtsCount
                N(i,1)=0.25d0*(1.-psi(i))*(1.-eta(i))
                N(i,2)=0.25d0*(1.+psi(i))*(1.-eta(i))
                N(i,3)=0.25d0*(1.+psi(i))*(1.+eta(i))
                N(i,4)=0.25d0*(1.-psi(i))*(1.+eta(i))
            enddo
            
            N1s(1,1:4)=N(1,1:4)
            N2s(1,1:4)=N(2,1:4)
            N3s(1,1:4)=N(3,1:4)
            N4s(1,1:4)=N(4,1:4)   

            call matrixtrans(N(1,1:4),N1T,1,4)
            call matrixtrans(N(2,1:4),N2T,1,4)
            call matrixtrans(N(3,1:4),N3T,1,4)
            call matrixtrans(N(4,1:4),N4T,1,4)  

            NtN1=matmul(N1T,N1s)
            NtN2=matmul(N2T,N2s)
            NtN3=matmul(N3T,N3s)
            NtN4=matmul(N4T,N4s)       
            
        end subroutine
!    ---------------------------
        subroutine printconsts
                write ( *,* )pi,alpha,beta,c1,c2,c3,c4,c5
                write ( *,* )psi
                write ( *,* )eta

        end subroutine
!!!        --------------------------------------- to cal. the average at the element knowing the nodes values
!!!        --------------------------------------- 
        subroutine CNCalcElemAverageFromNodes (ElemCount, NodeCount,ElemConnect,rNodesIn,rElemOut) 
            implicit none
            integer,   intent(in) :: ElemCount, NodeCount
            integer,    intent(in) :: ElemConnect(ElemCount,4)
            real*8,    intent(in) :: rNodesIn(NodeCount)
            real*8,    intent(out):: rElemOut(ElemCount)
            real*8 ::Valn1,Valn2,Valn3,Valn4
            integer i
            
            
            do i=1,ElemCount
                Valn1=rNodesIn(ElemConnect(i,1))
                Valn2=rNodesIn(ElemConnect(i,2))
                Valn3=rNodesIn(ElemConnect(i,3))
                Valn4=rNodesIn(ElemConnect(i,4))
                rElemOut(i)=0.25*(Valn1+Valn2+Valn3+Valn4)
            enddo
            
        end subroutine CNCalcElemAverageFromNodes
    
end module CN_Consts
! ###########################################################
! ###########################################################
module CN_BaseClass 
    implicit none
    TYPE CN_Base
!!!!        integer :: GaussIntPtsCount=4
        integer :: iSolutionActive=0
        integer :: iBCcount=0
!!!!!!        rRn=Rqg .... rRn_1=Rqo    .... rRn_2=Rqor ....
!!!!!!        rTn=Tinit .. rTn_1=Tinitd .... rTn_2=Tinitdold ....
        real*8 , ALLOCATABLE :: rRn_1(:),rRn_2(:)
        real*8 , ALLOCATABLE :: rTn_1(:),rTn_2(:)
        real*8 , ALLOCATABLE :: rThermal_Diff_Ele(:)
        real*8 , ALLOCATABLE :: rBCTempConc(:)
        integer, ALLOCATABLE :: iBCNodeID(:)
        
     contains ! Computation of area for rectangles.
        procedure ::AllocateMem => CN_Base_AllocateMem
        procedure ::CleanMem => CN_Base_CleanMem
        procedure ::CNSetBCdataBase => CN_Base_CNSetBCdata
        procedure ::CN_Base_GetSolution 
        procedure ::CN_Base_ForceBCs 
        procedure ::CNBase_InitiateTn_1 
        procedure ::CNBase_CopyNode 
        
    ENDTYPE
        Interface init_CN_Base;                    Module Procedure init_CN_Base;                 End Interface
!!        Interface AllocateMem;                Module Procedure CN_Base_AllocateMem;          End Interface
!!        Interface CleanMem;                   Module Procedure CN_Base_CleanMem;             End Interface
!!        Interface CNSetBCdataBase;            Module Procedure CN_Base_CNSetBCdata;          End Interface
        
        contains
            type ( CN_Base )function init_CN_Base ( mNodCountB,mElemCountB )
!            use CN_Consts
                integer, intent(in) :: mNodCountB,mElemCountB
                
!                NodCountB=mNodCountB
!                ElemCountB=mElemCountB
                init_CN_Base%iBCcount=0
                
            end function init_CN_Base
!!!        --------------------------------------- to allocate memory for the Rn's and CTn's
        subroutine CN_Base_AllocateMem (tCN_object) 
            use CN_Consts
            class ( CN_Base ), intent(inout) :: tCN_object
            integer ne,nn
            
            if(tCN_object%iSolutionActive == 0) return
            
            ne=ElemCountB
            nn=NodCountB
            allocate(tCN_object%rTn_1(nn),tCN_object%rTn_2(nn))
            allocate(tCN_object%rRn_1(nn),tCN_object%rRn_2(nn))
            allocate(tCN_object%rBCTempConc(nn),tCN_object%iBCNodeID(nn))
            allocate(tCN_object%rThermal_Diff_Ele(ne))
            
            tCN_object%rTn_1       = 0     
            tCN_object%rTn_2       = 0 
            tCN_object%rRn_1       = 0  
            tCN_object%rRn_2       = 0      
            tCN_object%rBCTempConc = 0
            tCN_object%iBCNodeID   = 0 
            tCN_object%rThermal_Diff_Ele      = 0
            
        end subroutine CN_Base_AllocateMem
!!!        --------------------------------------- to clean the memory for the Rn's and CTn's
        subroutine CN_Base_CleanMem (tCN_object) 
            class ( CN_Base ), intent(inout) :: tCN_object
            
            if(tCN_object%iSolutionActive == 0) return
            
            IF (ALLOCATED (tCN_object%rTn_1))           DEALLOCATE (tCN_object%rTn_1)
            IF (ALLOCATED (tCN_object%rTn_2))           DEALLOCATE (tCN_object%rTn_2)
            IF (ALLOCATED (tCN_object%rRn_1))           DEALLOCATE (tCN_object%rRn_1)
            IF (ALLOCATED (tCN_object%rRn_2))           DEALLOCATE (tCN_object%rRn_2)
            IF (ALLOCATED (tCN_object%rBCTempConc))     DEALLOCATE (tCN_object%rBCTempConc)
            IF (ALLOCATED (tCN_object%iBCNodeID))       DEALLOCATE (tCN_object%iBCNodeID)
            IF (ALLOCATED (tCN_object%rThermal_Diff_Ele))       DEALLOCATE (tCN_object%rThermal_Diff_Ele)
            
        end subroutine CN_Base_CleanMem
!!!        --------------------------------------- to clean the memory for the Rn's and CTn's
        subroutine CN_Base_CNSetBCdata (tCN_object,NodeID,BC_thermal_t,BC_thermal_flag) 
            implicit none
            class ( CN_Base ), intent(inout) :: tCN_object
            real,    intent(in) :: BC_thermal_t
            integer, intent(in) :: BC_thermal_flag,NodeID
            
            if(tCN_object%iSolutionActive == 0) return
            
            tCN_object%rBCTempConc(NodeID) = BC_thermal_t
            if(BC_thermal_flag == 1)then
                tCN_object%iBCcount = tCN_object%iBCcount+1
                tCN_object%iBCNodeID(tCN_object%iBCcount) = NodeID
            endif
            
        end subroutine CN_Base_CNSetBCdata
!!!        --------------------------------------- to clean the memory for the Rn's and CTn's
        subroutine CN_Base_ForceBCs (tCN_object,TdSol,rTn) 
            use CN_Consts,only:NodeCountAct
            implicit none
            class ( CN_Base ), intent(inout) :: tCN_object
            real*8,    intent(inout) :: TdSol(NodeCountAct)
            real*8,    intent(in) :: rTn(NodeCountAct)
            integer NodeID,i
            
            if(tCN_object%iSolutionActive == 0) return
            
            do i=1,tCN_object%iBCcount
                NodeID=tCN_object%iBCNodeID(i)
                TdSol(NodeID) =rTn(NodeID)
            enddo
            
        end subroutine CN_Base_ForceBCs
!!!        --------------------------------------- to clean the memory for the Rn's and CTn's
        subroutine CNBase_InitiateTn_1 (obj,ElemConnect) 
            use CN_Consts!,only:NodeCountAct,ElemCountAct
            implicit none
            class ( CN_Base ), intent(inout) :: obj
            integer,    intent(in) :: ElemConnect(ElemCountAct,4)
!!!!            real*8,    intent(out):: rThermal_Diff_Ele(ElemCountAct)
!!!!            integer i,n1,n2,n3,n4
            
            if(obj%iSolutionActive == 0) return ! if not active then return
            
            !---- initialize the Tn_1 from the input
            obj%rTn_1(1:NodeCountAct) = obj%rBCTempConc(1:NodeCountAct)
!!!!!            write(*,*)obj%rBCTempConc(1:NodeCountAct)
!!!!!            write(*,*)obj%rTn_1(1:NodeCountAct) 
!!!!!            write(*,*)ElemConnect
            
            !---- initialize the T or C at the element from the initial values at the nodes
            CALL CNCalcElemAverageFromNodes (ElemCountAct, NodeCountAct,ElemConnect,obj%rTn_1,obj%rThermal_Diff_Ele) 
              
        end subroutine CNBase_InitiateTn_1
!!!        --------------------------------------- to clean the memory for the Rn's and CTn's
        subroutine CNBase_CopyNode (obj,n,m) 
            implicit none
            class ( CN_Base ), intent(inout) :: obj
            integer,    intent(in) :: n,m
            
            obj%rTn_1(n)=obj%rTn_1(m)
            obj%rTn_2(n)=obj%rTn_2(m)
            obj%rRn_1(n)=obj%rRn_1(m)
            obj%rRn_2(n)=obj%rRn_2(m)
            
!!!!!!            Rqold(n)=Rqold(m)
!!!!!!            Rqolder(n)=Rqolder(m)
!!!!!!            Tinit(n)=Tinit(m)
!!!!!!            Tinitdold(n)=Tinitdold(m)
            
        end subroutine CNBase_CopyNode

!!!        --------------------------------------- to clean the memory for the Rn's and CTn's
        
        subroutine CN_Base_GetSolution (obj,DtCondDiff,ElemConnect,Cglobal,Kglobal,rRn,rTn)  
            use CN_Consts !!!!,only:NodeCountAct
            implicit none
            class ( CN_Base ), intent(inout) :: obj
            integer,    intent(in) :: ElemConnect(ElemCountAct,4)
            real*8, intent(in) :: DtCondDiff
            real*8, intent(in) ::  Cglobal(NodeCountAct,NodeCountAct), Kglobal(NodeCountAct,NodeCountAct), rRn(NodeCountAct)
            real*8, intent(out) ::  rTn(NodeCountAct)
            
            real*8 RHS(NodeCountAct)
            real*8 Tcoeff1(NodeCountAct,NodeCountAct),Rnew(NodeCountAct)
            real*8 Rold(NodeCountAct),Rolder(NodeCountAct)
            real*8 RHS2(NodeCountAct),LHS(NodeCountAct,NodeCountAct),TdSol(NodeCountAct)
            real*8 Tcoeff2(NodeCountAct,NodeCountAct)
            integer  IPIV(NodeCountAct),n
            integer ::INFO_flag=0

            if(obj%iSolutionActive == 0) return
            
            n=NodeCountAct
            LHS(1:n,1:n)=1./DtCondDiff*Cglobal(1:n,1:n) + (beta)*Kglobal(1:n,1:n)*(1.+alpha)
            Tcoeff1(1:n,1:n)=1./DtCondDiff*Cglobal(1:n,1:n) + (beta)*Kglobal(1:n,1:n)*(alpha)- (1.-beta)*(1.+alpha)*Kglobal(1:n,1:n)
            Tcoeff2(1:n,1:n)=(1.-beta)*Kglobal(1:n,1:n)*alpha
            
            Rold(1:n)=-(beta)*(alpha)*obj%rRn_1(1:n)+(1.-beta)*(1.+alpha)*obj%rRn_1(1:n)
            Rnew(1:n)=(beta)*rRn(1:n)*(1.+alpha)
            Rolder(1:n)=-(1.-beta)*(alpha)*obj%rRn_2(1:n)

!c     modifying C and K matrices for prescribed temperatures
            RHS2(1:n)=matmul(Tcoeff1(1:n,1:n),obj%rTn_1(1:n))+ matmul(Tcoeff2(1:n,1:n),obj%rTn_2(1:n))+ &
                        Rold(1:n)+Rnew(1:n)+Rolder(1:n)
!      to apply the effect of thermal loading
            call obj%CN_Base_ForceBCs (TdSol,rTn) 

            call ThermalSetLoading(n,RHS2,LHS,obj%rTn_1)
     
           call DGESV(n,1,LHS,n,IPIV,RHS2,n,INFO_flag)    
            if(INFO_flag.eq.0)then
                  TdSol(1:n)=RHS2(1:n)     
              write(*,*)'solve successful',DtCondDiff
            endif

            if(INFO_flag.ne.0)then    
              write(*,*)'solve unsuccessful',DtCondDiff
              call bye(2)
            endif
!------------------------------------------------	 
!      update only the nodes that has ICs not BCs
!   ----------------------------------testing

            call ThermalLdSetToSol(n,TdSol)

            rTn(1:n)=TdSol(1:n)

            obj%rRn_2(1:n)=obj%rRn_1(1:n)       ! save value at n step
            obj%rRn_1(1:n)=rRn(1:n)

            obj%rTn_2(1:n)=obj%rTn_1(1:n)
            obj%rTn_1(1:n)=rTn(1:n)
            
            !---- cal. the T or C at the element from the values at the nodes
            CALL CNCalcElemAverageFromNodes (ElemCountAct, NodeCountAct,ElemConnect,rTn,obj%rThermal_Diff_Ele) 
            
            
        end subroutine CN_Base_GetSolution      
        
        
end module CN_BaseClass
! ###########################################################
! ###########################################################
module CN_Thermal
    use CN_BaseClass
    implicit none
!!!        ---------------------------------------
    type,extends( CN_Base):: ThermalClass
        real*8 , ALLOCATABLE :: rThermal_DijSij(:)
        
        contains
!!!!!!!!!        procedure ::ThermalClass 
        procedure ::CN_Thermal_CleanMem 
!!!!!!!!!            Interface CN_Thermal_CleanMem;                Module Procedure CN_Thermal_CleanMem ;          End Interface
    end type ThermalClass
            Interface init_ThermalClass;                       Module Procedure init_ThermalClass ;            End Interface
!!!        ---------------------------------------          
!!!        ---------------------------------------
        contains
        type ( ThermalClass )function init_ThermalClass ( mNodCountB,mElemCountB )
            integer, intent(in) :: mNodCountB,mElemCountB
            
                init_ThermalClass%CN_Base=init_CN_Base(mNodCountB,mElemCountB)
                init_ThermalClass%iSolutionActive=1
                Call init_ThermalClass%AllocateMem()
                
                allocate(init_ThermalClass%rThermal_DijSij(mElemCountB))
                
        end function init_ThermalClass
! ###########################################################
! ###########################################################
        subroutine CN_Thermal_CleanMem (tCN_object)
            class ( ThermalClass ), intent(inout) :: tCN_object
            
            call tCN_object%CleanMem()
            IF (ALLOCATED (tCN_object%rThermal_DijSij))           DEALLOCATE (tCN_object%rThermal_DijSij)
            
        end subroutine CN_Thermal_CleanMem
end module CN_Thermal
! ###########################################################
! ###########################################################
module CN_Diffusion
    use CN_BaseClass
    implicit none
!!!        ---------------------------------------
    type,extends( CN_Base):: DiffusionClass
    end type DiffusionClass
!!!        ---------------------------------------
!!!        ---------------------------------------   
    Interface init_DiffusionClass;                Module Procedure init_DiffusionClass ;          End Interface
        
!!!        ---------------------------------------
        contains
        type ( DiffusionClass )function init_DiffusionClass ( mNodCountB,mElemCountB )
            integer, intent(in) :: mNodCountB,mElemCountB
            
                init_DiffusionClass%CN_Base=init_CN_Base(mNodCountB,mElemCountB)
                init_DiffusionClass%iSolutionActive=1
                Call init_DiffusionClass%AllocateMem()
            
        end function init_DiffusionClass
end module CN_Diffusion
! ###########################################################

      