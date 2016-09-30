
module EC_Consts
    use mod_parameters
    implicit none
    type typeEC_IntersectionPtsData
        INTEGER :: EC_NodeGlobalID=0    !--gobal id wrt all nodes
        real*8  :: EC_NodeLenRatio=0     !--ratio of the defined vector below
        INTEGER :: EC_NodeVectStart=0   !--start vector pt
        INTEGER :: EC_NodeVectEnd=0     !--end vector pt
    end type  typeEC_IntersectionPtsData


    public
    integer EC_eCrackNon,EC_eCrackInnerSide,EC_eCrackBothSides,EC_eCrackOutterSide
    parameter (EC_eCrackNon=0)          !--edge not cracked
    parameter (EC_eCrackInnerSide=1)    !--edge cracked on first side
    parameter (EC_eCrackBothSides=2)     !--edge cracked on both sides
    parameter (EC_eCrackOutterSide=3)   !--edge crack on second side

    integer EC_eElemTypeTri,EC_eElemTypeQuad,EC_eElemTypePenta
    parameter (EC_eElemTypeTri=3)   !--4-node elem triangle
    parameter (EC_eElemTypeQuad=4)   !--4-node elem quad
    parameter (EC_eElemTypePenta=5)   !--5-node elem pentagon

    integer EC_eElemOriginMain,EC_eElemOriginOverlap
    parameter (EC_eElemOriginMain=0)   !--original input elem
    parameter (EC_eElemOriginOverlap=1)   !--overlapping elem


    integer EC_bCracking !-- 1=allow fracture
    INTEGER EC_ZoneFactor !--- used to mark the crack front circle, do not crack within this area
    real*8  EC_ZoneRadius !--- used to mark the crack front circle, do not crack within this area
    INTEGER EC_DecayCount !--- used to mark the crack front circle, do not crack within this area
    real*8  EC_ElemEdgeMin !--- min elem edge
    real*8  EC_ElemEdgeMax !--- max elem edge
    real*8  EC_ElemAverageDia !--- average elem diagonal

    integer EC_PreCrackedElemCount !- number of pre-cracked elems
    integer EC_NodesCountAdded !-- nodes added as overlap nodes
    integer EC_ElemCountAdded !-- elems added as overlap elems
    integer EC_IntersectNodesCount !-- intersection nodes count

    integer EC_ElemCountInput !-- original elem count from the input before any phantom or overlap elements
    integer EC_ElemCountCurrent !-- current or updated or recent elem count that = original + overlapped

    integer EC_NodeCountInput !-- original node count from the input before any phantom or overlap elements
    integer EC_NodeCountCurrent !-- current or updated or recent node count that = original + overlapped

    integer EC_ElemEdgesConnect(2,4) !-- list the elem edges order 1-2 2-3 3-4 4-1
!!    integer EC_ElemNodesEdges(2,4) !-- list the elem nodes' edges order [1:2,1] [2:3,2] [3:4,3] [4:1,4]
    integer EC_ElemEdgesOpposite(4) !-- list the elem edges opposite edges 1-->3, 2-->4, 3-->1,4-->1
    integer EC_ElemEdgesBeforeAfter(2,4) !-- list the elem edges before and after edges  4-2,1-3,2-4,3-1


    Interface ECinitializeConsts;                  Module Procedure ECinitializeConsts;            End Interface
    Interface ECcalcAverageElemStatistics;         Module Procedure ECcalcAverageElemStatistics;   End Interface
    Interface ECprintTest;         Module Procedure ECprintTest;   End Interface
!    ---------------------------
    contains
    subroutine ECprintTest
        use mod_file_units
        write(iFU_crackinput_out,*)'EC_bCracking=',EC_bCracking
        write(iFU_crackinput_out,*)'EC_ZoneFactor=',EC_ZoneFactor
        write(iFU_crackinput_out,*)'EC_ZoneRadius=',EC_ZoneRadius
        write(iFU_crackinput_out,*)'EC_DecayCount=',EC_DecayCount
        write(iFU_crackinput_out,*)'EC_ElemEdgeMin=',EC_ElemEdgeMin
        write(iFU_crackinput_out,*)'EC_ElemEdgeMax=',EC_ElemEdgeMax
        write(iFU_crackinput_out,*)'EC_ElemAverageDia=',EC_ElemAverageDia
        write(iFU_crackinput_out,*)'EC_PreCrackedElemCount=',EC_PreCrackedElemCount
        write(iFU_crackinput_out,*)'EC_NodesCountAdded=',EC_NodesCountAdded
        write(iFU_crackinput_out,*)'EC_ElemCountAdded=',EC_ElemCountAdded
        write(iFU_crackinput_out,*)'EC_ElemCountInput=',EC_ElemCountInput
        write(iFU_crackinput_out,*)'EC_ElemCountCurrent=',EC_ElemCountCurrent
        write(iFU_crackinput_out,*)'EC_NodeCountInput=',EC_NodeCountInput
        write(iFU_crackinput_out,*)'EC_NodeCountCurrent=',EC_NodeCountCurrent
        write(iFU_crackinput_out,*)'EC_ElemEdgesConnect=',EC_ElemEdgesConnect
    end subroutine ECprintTest
!##############################################################################
!##############################################################################
    subroutine ECinitializeConsts(numnpo, numelto,ElemConnect,rNodesCoordx,rNodesCoordy,mEC_bCracking,mEC_ZoneFactor, &
                                mEC_DecayCount,mEC_PreCrackedElemCount)

        implicit none
        !--- set the input node and elem counts before cracking and overlap
        integer,    intent(in) :: ElemConnect(4,*)
        real,     intent(in) :: rNodesCoordx(*),rNodesCoordy(*)
        INTEGER numnpo, numelto
        INTEGER i
        integer mEC_ZoneFactor,mEC_ZoneRadius,mEC_DecayCount,mEC_PreCrackedElemCount,mEC_bCracking
        !--- if no fracture then no need to process this subroutine
        if (mEC_bCracking == 0) then
            return
        endif

        EC_bCracking=mEC_bCracking
        EC_ElemCountInput=numelto
        EC_NodeCountInput=numnpo
        do i=1,3
          EC_ElemEdgesConnect(1,i)=i
          EC_ElemEdgesConnect(2,i)=i+1
        enddo
        EC_ElemEdgesConnect(1,4)=4
        EC_ElemEdgesConnect(2,4)=1
!!        !-- list the elem nodes' edges order [1,2] [2,3] [3,4] [4,1]
!!        do i=1,3
!!            EC_ElemNodesEdges(1,i)=i
!!            EC_ElemNodesEdges(2,i)=i+1
!!        enddo
!!        EC_ElemNodesEdges(1,4)=4
!!        EC_ElemNodesEdges(2,4)=1

        !-- list the elem edges opposite edges 1-->3, 2-->4, 3-->1,4-->1
        do i=1,2
            EC_ElemEdgesOpposite(i)=i+2
        enddo
        do i=3,4
            EC_ElemEdgesOpposite(i)=i-2
        enddo
        !-- list the elem edges before and after edges  4-2,1-3,2-4,3-1
        EC_ElemEdgesBeforeAfter(1,1)=4
        EC_ElemEdgesBeforeAfter(2,1)=2
        do i=2,3
            EC_ElemEdgesBeforeAfter(1,i)=i-1
            EC_ElemEdgesBeforeAfter(2,i)=i+1
        enddo
        EC_ElemEdgesBeforeAfter(1,4)=3
        EC_ElemEdgesBeforeAfter(2,4)=1

        EC_ZoneFactor=mEC_ZoneFactor
        EC_DecayCount=mEC_DecayCount
        EC_PreCrackedElemCount=mEC_PreCrackedElemCount
        EC_NodesCountAdded=0
        EC_ElemCountAdded=0

        EC_ElemCountCurrent=EC_ElemCountInput
        EC_NodeCountCurrent=EC_NodeCountInput

        call ECcalcAverageElemStatistics (ElemConnect,rNodesCoordx,rNodesCoordy)
        EC_ZoneRadius=EC_ZoneFactor*EC_ElemEdgeMin

    end subroutine ECinitializeConsts
!##############################################################################
!##############################################################################
!!!        --------------------------------------- to cal. some statistics about elem geometry
    subroutine ECcalcAverageElemStatistics (ElemConnect,rNodesCoordx,rNodesCoordy)
        implicit none
        integer,    intent(in) :: ElemConnect(4,*)
        real,     intent(in) :: rNodesCoordx(*),rNodesCoordy(*)
        real*8 ::x1,y1,x2,y2, rEdgeMax,rEdgeMin,rElemDia,dist
        integer i,j,edgeNode1,edgeNode2

        rEdgeMax=-1e32
        rEdgeMin=1e32
        rElemDia=0

        do i=1,EC_ElemCountInput
            do j=1,4
                edgeNode1=ElemConnect(EC_ElemEdgesConnect(1,j),i)
                edgeNode2=ElemConnect(EC_ElemEdgesConnect(2,j),i)
                x1=rNodesCoordx(edgeNode1)
                y1=rNodesCoordy(edgeNode1)
                x2=rNodesCoordx(edgeNode2)
                y2=rNodesCoordy(edgeNode2)
                dist=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))
                !-- the max elem edge
                if(dist>rEdgeMax) then
                  rEdgeMax=dist
                endif
                !-- the min elem edge
                if((dist<rEdgeMin) .and. (dist>0)) then
                  rEdgeMin=dist
                endif
            enddo
            !-- calc the elem diagonal
            edgeNode1=ElemConnect(1,i)
            edgeNode2=ElemConnect(3,i)
            x1=rNodesCoordx(edgeNode1)
            y1=rNodesCoordy(edgeNode1)
            x2=rNodesCoordx(edgeNode2)
            y2=rNodesCoordy(edgeNode2)
            dist=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))
            rElemDia=rElemDia+dist
        enddo
        !-- calc the average diagonal
        rElemDia=rElemDia/EC_ElemCountInput

        !--- set the results to the CN_Consts
        EC_ElemEdgeMax=rEdgeMax
        EC_ElemEdgeMin=rEdgeMin
        EC_ElemAverageDia=rElemDia
    end subroutine ECcalcAverageElemStatistics

end module EC_Consts
! ###########################################################
! ###########################################################
module EC_ElemCrackingBaseClass
    use EC_Consts
    implicit none
    TYPE EC_ElemCrackingClass
        integer :: iElemType=EC_eElemTypeQuad !- default 4-nodes
        integer :: iElemStatus=0 !- 0=not cracked, 1:EC_DecayCount=decaying, >EC_DecayCount=cracked
        real*8  :: rCleavagePlane(3)
        real*8  :: rCoordRatioCracking(4) !- (-ve no cracking )for each edge, ratio of length for the cracking point
        integer :: EdgeStatus(4) !- for each edge status, see enum_EdgeStatus
        real*8  :: rAreaRatio,rArea !- for area ratio with respect of the original elem. also area value
        integer :: iElemOverlapping=EC_eElemOriginMain !- for the added new/overlapping/phantom elem;0=main, 1=overlap
        integer :: iElemEdgeNeighbors(6,4) !- each col lists the neighbors of this edge index
        integer :: iElemConnectivity(4) !- elem nodes


     contains
        procedure ::Initialize => EC_ElemCrackingBaseClass_Initialize
        procedure ::PrintTest => EC_ElemCrackingBaseClass_PrintTest
        procedure ::SetDataPre => EC_ElemCrackingBaseClass_SetDataPre
        procedure ::SetFailed => EC_ElemCrackingBaseClass_SetFailed
        procedure ::CalcElemMaxStress => EC_ElemCrackingBaseClass_CalcElemMaxStress
        procedure ::CheckDecaying => EC_ElemCrackingBaseClass_CheckDecaying
        procedure ::EC_CopyElem
        procedure ::EC_CalcAreaRatio
        procedure ::EC_GetElemAreaRatio
        procedure ::EC_GetElemUnloadingCount
        procedure ::EC_GetElemEdgeNeighbors
        procedure ::EC_ElemAddToNeighbors
        procedure ::EC_SetElemEdgeNeighbors
        procedure ::EC_CheckElemValid
        procedure ::EC_CalcCrackedElemAreaRatio

    ENDTYPE
    !-------------------
    contains
!!!        --------------------------------------- to allocate memory for the Rn's and CTn's
!##############################################################################
     subroutine EC_GetElemEdgeNeighbors (tEC_object,mEdgeId,mOutList)
        class ( EC_ElemCrackingClass ), intent(inout) :: tEC_object
        integer, intent(in) :: mEdgeId
        integer, intent(out) :: mOutList(6)
        mOutList(1:6)=tEC_object%iElemEdgeNeighbors(1:6,mEdgeId)
    end subroutine EC_GetElemEdgeNeighbors
!##############################################################################
     subroutine EC_SetElemEdgeNeighbors (tEC_object,mInList)
        class ( EC_ElemCrackingClass ), intent(inout) :: tEC_object
        integer, intent(in) :: mInList(8)
        integer i
        do i=1,4
            tEC_object%iElemEdgeNeighbors(1,i)=mInList(i)
            tEC_object%iElemEdgeNeighbors(2,i)=mInList(i+4)
        enddo
    end subroutine EC_SetElemEdgeNeighbors
!##############################################################################
     integer function EC_CheckElemValid (tEC_object)
        class ( EC_ElemCrackingClass ), intent(inout) :: tEC_object
        integer i,counter,valid(4)
        !- check that there are 2 consequtive nodes of the orignial elem
        valid=0
        counter=0
        EC_CheckElemValid=0
        do i=1,4
            if(tEC_object%iElemConnectivity(i) <= EC_NodeCountInput) then
                valid(i)=1
            endif
        enddo

        do i=1,3
            if(valid(i)+valid(i+1)>=2) then !-if there is 2 cosequtive nodes of
                EC_CheckElemValid=1
                return
            endif
        enddo
        if(valid(4)+valid(1)>=2) then !-if there is 2 cosequtive nodes of
            EC_CheckElemValid=1
            return
        endif

    end function EC_CheckElemValid
!##############################################################################
     subroutine EC_ElemAddToNeighbors (tEC_object,mMyEdgeId,mElemOther,mEdgeOther)
        class ( EC_ElemCrackingClass ), intent(inout) :: tEC_object
        integer, intent(in) :: mMyEdgeId,mElemOther,mEdgeOther
        integer i
        do i=1,3
            if(tEC_object%iElemEdgeNeighbors((i-1)*2+1,mMyEdgeId) == 0) then
                tEC_object%iElemEdgeNeighbors((i-1)*2+1,mMyEdgeId)=mElemOther
                tEC_object%iElemEdgeNeighbors((i-1)*2+2,mMyEdgeId)=mEdgeOther
                return
            endif
        enddo
    end subroutine EC_ElemAddToNeighbors
!##############################################################################
    real*8 function EC_GetElemAreaRatio (tEC_object)
        class ( EC_ElemCrackingClass ), intent(inout) :: tEC_object
        EC_GetElemAreaRatio=tEC_object%rAreaRatio
    end function EC_GetElemAreaRatio
!##############################################################################
    integer function EC_GetElemUnloadingCount (tEC_object)
        class ( EC_ElemCrackingClass ), intent(inout) :: tEC_object
        EC_GetElemUnloadingCount=tEC_object%iElemStatus
    end function EC_GetElemUnloadingCount
!##############################################################################
    subroutine EC_CopyElem (tEC_object,tEC_objectOUT)
        implicit none
        class ( EC_ElemCrackingClass ), intent(inout) :: tEC_object,tEC_objectOUT

        tEC_objectOUT%iElemStatus=tEC_object%iElemStatus
        tEC_objectOUT%rCleavagePlane=tEC_object%rCleavagePlane
        tEC_objectOUT%rCoordRatioCracking=tEC_object%rCoordRatioCracking
        tEC_objectOUT%EdgeStatus=tEC_object%EdgeStatus
        tEC_objectOUT%rAreaRatio=tEC_object%rAreaRatio
        tEC_objectOUT%rArea=tEC_object%rArea
        tEC_objectOUT%iElemType=tEC_object%iElemType
        tEC_objectOUT%iElemOverlapping=tEC_object%iElemOverlapping
        tEC_objectOUT%iElemEdgeNeighbors=tEC_object%iElemEdgeNeighbors
        tEC_objectOUT%iElemConnectivity=tEC_object%iElemConnectivity
    end subroutine EC_CopyElem
!##############################################################################

    subroutine EC_ElemCrackingBaseClass_Initialize (tEC_object)

        class ( EC_ElemCrackingClass ), intent(inout) :: tEC_object

        tEC_object%EdgeStatus=-1
        tEC_object%rCleavagePlane=[0,0,1]
        tEC_object%iElemStatus=0
        tEC_object%rCoordRatioCracking=-1
        tEC_object%rArea=0
        tEC_object%rAreaRatio=0
        tEC_object%iElemType=EC_eElemTypeQuad
        tEC_object%iElemOverlapping=EC_eElemOriginMain
        tEC_object%iElemEdgeNeighbors=-1
        tEC_object%iElemConnectivity=-1

    end subroutine EC_ElemCrackingBaseClass_Initialize
!##############################################################################
!##############################################################################
    subroutine EC_ElemCrackingBaseClass_PrintTest (tEC_object,mid)
        use mod_file_units
        class ( EC_ElemCrackingClass ), intent(inout) :: tEC_object
        integer mid

            write(iFU_crackinput_out,*)'================================ elem=',mid
            write(iFU_crackinput_out,*)'         ------------- iElemStatus='
            write(iFU_crackinput_out,*)tEC_object%iElemStatus
            write(iFU_crackinput_out,*)'         ------------- EdgeStatus='
            write(iFU_crackinput_out,*)tEC_object%EdgeStatus
            write(iFU_crackinput_out,*)'         ------------- rAreaRatio='
            write(iFU_crackinput_out,*)tEC_object%rAreaRatio
             write(iFU_crackinput_out,*)'         ------------- rArea='
            write(iFU_crackinput_out,*)tEC_object%rArea
            write(iFU_crackinput_out,*)'         ------------- iElemConnectivity='
            write(iFU_crackinput_out,*)tEC_object%iElemConnectivity
            write(iFU_crackinput_out,*)'         ------------- iElemEdgeNeighbors='
            write(iFU_crackinput_out,*)tEC_object%iElemEdgeNeighbors(1:2,1:4)
            write(iFU_crackinput_out,*)'         ------------- rCleavagePlane='
            write(iFU_crackinput_out,*)tEC_object%rCleavagePlane(1:3)

    end subroutine EC_ElemCrackingBaseClass_PrintTest
    !##############################################################################
    !##############################################################################
    subroutine EC_ElemCrackingBaseClass_SetDataPre (tEC_object,vn,ElemEdge1,rxc1,ryc1,r1,ElemEdge2,rxc2,ryc2,r2)

        implicit none
        class ( EC_ElemCrackingClass ), intent(inout) :: tEC_object
        real , intent(in) :: vn(3)
        integer , intent(in) ::ElemEdge1,ElemEdge2
        real , intent(in) ::rxc1,ryc1,r1,rxc2,ryc2,r2

        tEC_object%iElemStatus=EC_DecayCount+1
        tEC_object%rCleavagePlane=vn
        tEC_object%EdgeStatus=EC_eCrackInnerSide

    end subroutine EC_ElemCrackingBaseClass_SetDataPre
     !##############################################################################
     !##############################################################################
    subroutine EC_ElemCrackingBaseClass_SetFailed (tEC_object,miElemStatus,vn)

        implicit none
        class ( EC_ElemCrackingClass ), intent(inout) :: tEC_object
        real*8 , intent(in) :: vn(3)
        integer , intent(in) :: miElemStatus

        tEC_object%iElemStatus=miElemStatus
        tEC_object%rCleavagePlane=vn

    end subroutine EC_ElemCrackingBaseClass_SetFailed
     !##############################################################################
      !##############################################################################
    subroutine EC_ElemCrackingBaseClass_CalcElemMaxStress (tEC_object,ElemCleavagePlanes, &
                                                ElemStress,MaxStress100,iPlaneId)

        implicit none
        class ( EC_ElemCrackingClass ), intent(inout) :: tEC_object
        real*8 , intent(in) :: ElemCleavagePlanes(3,3),ElemStress(4)
        real*8 , intent(out) :: MaxStress100
        integer , intent(out) :: iPlaneId
        real*8 vy,vz,s11,s22,s44,StressComp
        integer j

           !-- obtain maximum normal component of the traction on cleavage planes and corresponding normal vector
        MaxStress100 = -100.0
        iPlaneId=-1
        do j = 1, 3
           vy=ElemCleavagePlanes(j,2)
           vz=ElemCleavagePlanes(j,3)
           s11=ElemStress(1)
           s22=ElemStress(2)
           s44=ElemStress(4)
           StressComp = s11*vy*vy+s22*vz*vz+s44*vy*vz*2.0
           !!!dum = sig(1,ele)*cleave(j,2)**2.0+sig(2,ele)*cleave(j,3)**2.0+sig(4,ele)*cleave(j,2)*cleave(j,3)*2.0
           if (abs(StressComp)>MaxStress100) then
               MaxStress100=abs(StressComp)
               iPlaneId=j
           end if
        end do
    end subroutine EC_ElemCrackingBaseClass_CalcElemMaxStress
    !##############################################################################
    !##############################################################################
    function EC_ElemCrackingBaseClass_CheckDecaying (tEC_object,ele,nstep)result(bDecaying)
        use mod_file_units
        implicit none
        class ( EC_ElemCrackingClass ), intent(inout) :: tEC_object
        integer, intent(in):: ele,nstep
        integer bDecaying

        !-- to check if the elem is failed and unloading
        !- if elem is unloading then increase the unloading steps
        if((tEC_object%iElemStatus <= EC_DecayCount) .and.  (tEC_object%iElemStatus>0)) then
          tEC_object%iElemStatus=tEC_object%iElemStatus+1
          write(iFU_crackprog_out,*) ele,nstep,tEC_object%iElemStatus,EC_DecayCount
          bDecaying= 1
        elseif( tEC_object%iElemStatus > EC_DecayCount) then !-for elem that unloaded but not split
          bDecaying= 1
        else
          bDecaying= 0
        endif
    end function EC_ElemCrackingBaseClass_CheckDecaying
    !##############################################################################
    subroutine EC_CalcAreaRatio (tEC_object,xs,ys,nPts)

        implicit none
        class ( EC_ElemCrackingClass ), intent(inout) :: tEC_object
        integer, intent(in):: nPts
        real*8 , intent(in)::xs(nPts),ys(nPts) !-input lines p-p2, q-q2
        real*8 CalcPolygonArea
        tEC_object%rArea=CalcPolygonArea(xs,ys,4)
        tEC_object%rAreaRatio=1

    end subroutine EC_CalcAreaRatio
    !##############################################################################
    subroutine EC_CalcCrackedElemAreaRatio (tEC_object,xs,ys,nPts)

        implicit none
        class ( EC_ElemCrackingClass ), intent(inout) :: tEC_object
        integer, intent(in):: nPts
        real*8 , intent(in)::xs(nPts),ys(nPts) !-input lines p-p2, q-q2
        real*8 CalcPolygonArea,DistRatio,AreaTemp
        real*8 rPtc(2,2),rPtc1(2),rPtc2(2),rP1(2),rP2(2)
        real*8 :: ElemNodesCoordx(6),ElemNodesCoordy(6),xNew(4),yNew(4)    
        integer :: i,counter
        integer :: nodesOrder(4),CrackedEdges(2)
        integer :: EdgeBefore,EdgeAfter,node1
        !----------------------------------- see sketch in P1
        nodesOrder=0
        CrackedEdges=0
        counter=1
        !-- find the cracked edges and the points
        do i=1,4
            if(tEC_object%rCoordRatioCracking(i)>0.0) then
                DistRatio=tEC_object%rCoordRatioCracking(i)
                rP1(1)=xs(EC_ElemEdgesConnect(1,i))
                rP1(2)=ys(EC_ElemEdgesConnect(1,i))
                rP2(1)=xs(EC_ElemEdgesConnect(2,i))
                rP2(2)=ys(EC_ElemEdgesConnect(2,i))
                CALL GetPtOnLine(rP1,rP2,DistRatio,rPtc1)
                rPtc(1:2,counter)=rPtc1(1:2)
                CrackedEdges(counter)=i
                counter=counter+1
            endif
            !- limit two cracked edges for now ....
            if ( counter>2 ) exit !-- break the loop
        enddo
        !- elem has nodes 1,2,3,4 and the edge crack points are 5, 6
        !- if they are opposite edges
        if (EC_ElemEdgesOpposite(CrackedEdges(1)) == CrackedEdges(2))then
            !- this is a 4 node elem
            tEC_object%iElemType=EC_eElemTypeQuad
            EdgeBefore=EC_ElemEdgesBeforeAfter(1,CrackedEdges(1))
            EdgeAfter=EC_ElemEdgesBeforeAfter(2,CrackedEdges(1))
            !- my first point
            Node1=EC_ElemEdgesConnect(1,CrackedEdges(1))
            if ( node1> EC_NodeCountInput ) then
                node1=EC_ElemEdgesConnect(2,CrackedEdges(1))
                !- my crack point
                nodesOrder(1)=5
                !- my second point
                nodesOrder(2)=node1
                !- the seond node of my after edge
                nodesOrder(3)=EC_ElemEdgesConnect(2,EdgeAfter)
                !- the other crack point
                nodesOrder(4)=6
            else
                !- my crack point
                nodesOrder(1)=5
                !- the other crack point
                nodesOrder(2)=6
                !- the first node of my before edge
                nodesOrder(3)=EC_ElemEdgesConnect(1,EdgeBefore)
                !- my first point
                nodesOrder(4)=node1
            end if
        end if
        ElemNodesCoordx(1:4)=xs(1:4)
        ElemNodesCoordy(1:4)=ys(1:4)
        ElemNodesCoordx(5:6)=rPtc(1,1:2)
        ElemNodesCoordy(5:6)=rPtc(2,1:2)
        do i = 1, 4
            xNew(i)=ElemNodesCoordx(nodesOrder(i))
            yNew(i)=ElemNodesCoordy(nodesOrder(i))
        end do

        AreaTemp=CalcPolygonArea(xNew,yNew,4)
!!!        write(*,*)tEC_object%rArea,AreaTemp
        tEC_object%rAreaRatio=AreaTemp/tEC_object%rArea

    end subroutine EC_CalcCrackedElemAreaRatio
    !##############################################################################

end module EC_ElemCrackingBaseClass
! ###########################################################
! ###########################################################
