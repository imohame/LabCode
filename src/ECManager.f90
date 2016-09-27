
! ###########################################################
module EC_Objects_manager
    use EC_ElemCrackingBaseClass
    use EC_Consts
!    implicit none
    public
    ! Declare a set geometric objects.
    type ( EC_ElemCrackingClass ) , ALLOCATABLE :: pEC_ElemData(:)
!!!!    integer, ALLOCATABLE :: EC_ElemNeighbors(:,:) !-(8,:) to hold the elem's edges neighboring elem/edges, update when overlapping
!-- EC_ElemNeighbors(e1,e2,e3,e4,ed1,ed2,ed3,ed4,:)

    type ( typeEC_IntersectionPtsData ) , ALLOCATABLE :: EC_NodesIntersectData(:)
    integer EC_NodesIntersectDataCount
!!!        ---------------------------------------

    Interface EC_AllocateMem;              Module Procedure EC_AllocateMem;                   End Interface
    Interface EC_CleanMem;                 Module Procedure EC_CleanMem;                      End Interface
    Interface EC_PrintTest;                Module Procedure EC_PrintTest;                     End Interface
    Interface EC_MarkElemForCrack;         Module Procedure EC_MarkElemForCrack;              End Interface
    Interface EC_SplitElem;                Module Procedure EC_SplitElem;                     End Interface
    Interface EC_CopyElemsData;            Module Procedure EC_CopyElemsData;                 End Interface
    Interface EC_CopyNodesData;            Module Procedure EC_CopyNodesData;                 End Interface


!##############################################################################
!##############################################################################
    contains
        subroutine EC_AllocateMem()
            use EC_Consts

            implicit none
            integer ne,i


            if(EC_bCracking == 0) then
              return
            endif
            ne=int(EC_ElemCountInput*1.25)
            Allocate(pEC_ElemData(ne))
!!!!!            Allocate(EC_ElemNeighbors(8,ne))
            ne=int(EC_NodeCountInput*.1) +100
            Allocate(EC_NodesIntersectData(ne))

            do i=1,ne
                pEC_ElemData(i)%EdgeStatus=-1
                pEC_ElemData(i)%rCleavagePlane=[0,0,1]
                pEC_ElemData(i)%iElemStatus=0
                pEC_ElemData(i)%rCoordRatioCracking=-1
                pEC_ElemData(i)%rArea=0
                pEC_ElemData(i)%rAreaRatio=0
                pEC_ElemData(i)%iElemType=EC_eElemTypeQuad
                pEC_ElemData(i)%iElemOverlapping=EC_eElemOriginMain
            enddo
!!!!            EC_ElemNeighbors=-1

            EC_NodesIntersectDataCount=0

         end subroutine EC_AllocateMem
!##############################################################################
!##############################################################################
        subroutine EC_CleanMem()
          use EC_Consts
          implicit none

          if(EC_bCracking == 0) then
            return
          endif
          IF (ALLOCATED (pEC_ElemData))              DEALLOCATE (pEC_ElemData)
!!!          IF (ALLOCATED (EC_ElemNeighbors))              DEALLOCATE (EC_ElemNeighbors)
          IF (ALLOCATED (EC_NodesIntersectData))              DEALLOCATE (EC_NodesIntersectData)
        end subroutine EC_CleanMem
!##############################################################################
!##############################################################################
        subroutine EC_PrintTest()
          use EC_Consts
          use EC_ElemCrackingBaseClass
          use mod_file_units
          implicit none
          integer i

          if(EC_bCracking == 0) then
            return
          endif

          CALL ECprintTest()

          do i=1,EC_ElemCountCurrent
              write(iFU_crackinput_out,*)'==========================='
              call pEC_ElemData(i)%PrintTest(i)
              write(iFU_crackinput_out,*)'==========================='
            !   write(*,*)pEC_ElemData(i)
            !   write(*,*)pEC_ElemData(i)%rCleavagePlane
            !   write(*,*)pEC_ElemData(i)%EdgeStatus
          enddo
        end subroutine EC_PrintTest
!##############################################################################
!##############################################################################
  subroutine EC_MarkElemForCrack(ElemId,rStartPt,ElemCoordx,ElemCoordy)
    use EC_Consts
    use EC_ElemCrackingBaseClass
    implicit none
    integer , intent(in)::ElemId
    real*8 , intent(in)::rStartPt(2) !- the crack staring point
    real*8 , intent(in)::ElemCoordx(4),ElemCoordy(4)
    integer i,ElemIdNeighbor,EdgeIdNeighbor,bIntersect,edgeNode1,edgeNode2
    real*8 Vx,Vy,p1(2),p2(2),q1(2),q2(2),rPtc(2),pRatio,qRatio
    real*8 ElemAreaTotal
    !-- calculate elem area
    call pEC_ElemData(ElemId)%CalcAreaRatio(ElemCoordx,ElemCoordy,4)
    !- this elem is cracked either on cleavage planes, or on the neighboring edges
    !-- set all edges to inner side crack
    pEC_ElemData(ElemId)%EdgeStatus(1:4)=EC_eCrackInnerSide
    !-- loop over each edge and check the neighbor for possible cracking
    do i=1,4
        ElemIdNeighbor=pEC_ElemData(ElemId)%iElemEdgeNeighbors(1,i)
        EdgeIdNeighbor=pEC_ElemData(ElemId)%iElemEdgeNeighbors(2,i)
        !-- if the edge is outter, has no neighbors then make it EC_eCrackBothSides
        if( EdgeIdNeighbor == 0) then
            pEC_ElemData(ElemId)%EdgeStatus(i)=EC_eCrackBothSides
            continue ! go to next edge
        endif
        !-if the neighboring edge is inner cracked then set this edge on both sides as cracked
        if(pEC_ElemData(ElemIdNeighbor)%EdgeStatus(EdgeIdNeighbor)==EC_eCrackInnerSide) then
            pEC_ElemData(ElemId)%EdgeStatus(i)=EC_eCrackBothSides
            pEC_ElemData(ElemIdNeighbor)%EdgeStatus(EdgeIdNeighbor)=EC_eCrackBothSides
        endif
    enddo
    !- this the crack direction
    !- since vx,vy are the components of the normal then the crack dir are vy,-vx
    !- the crack direction is normal to the cleavage dir.
    Vx=pEC_ElemData(ElemId)%rCleavagePlane(3)
    Vy=-pEC_ElemData(ElemId)%rCleavagePlane(2)

    !- make a crack line given a starting point and a direction
    CALL GetLineStartEndPts(rStartPt,Vx,Vy,EC_ElemAverageDia*2,p1,p2)
    !- intersect the crack line with the edges
    do i=1,4
        edgeNode1=EC_ElemEdgesConnect(1,i)
        edgeNode2=EC_ElemEdgesConnect(2,i)
        q1(1)=ElemCoordx(edgeNode1)
        q1(2)=ElemCoordy(edgeNode1)
        q2(1)=ElemCoordx(edgeNode2)
        q2(2)=ElemCoordy(edgeNode2)
        CALL LineSegIntersection(p1,p2,q1,q2,pRatio,qRatio,rPtc,bIntersect)
        if(bIntersect==1) then !- there is an intersection
            pEC_ElemData(ElemId)%rCoordRatioCracking(i)=qRatio
            !-- this for the neighboring edge properties
            if(pEC_ElemData(ElemId)%EdgeStatus(i)==EC_eCrackBothSides) then !- then set this value to the other edge
                ElemIdNeighbor=pEC_ElemData(ElemId)%iElemEdgeNeighbors(1,i)
                EdgeIdNeighbor=pEC_ElemData(ElemId)%iElemEdgeNeighbors(2,i)
                !-- if the edge is outter, has no neighbors then make it EC_eCrackBothSides
                if( EdgeIdNeighbor /= 0) then
                    !- 1-qRatio because the edge nodes order is reversed from elem to the neighboring elem
                    pEC_ElemData(ElemIdNeighbor)%rCoordRatioCracking(EdgeIdNeighbor)=1-qRatio
                endif
            endif
        endif
    enddo

  end subroutine EC_MarkElemForCrack

 !##############################################################################
 !##############################################################################
    subroutine EC_SplitElem(ElemId,NodesCoordx, NodesCoordy, ElemConnect, DofIds, NodesDispl, ElemMaterial, &
                            usi  , freep , ym,SolStepCount)
                !-a(k03)    ,a(k04)      ,a(k02)      ,a(k57)  ,a(k18)     , a(k08)      ,a(k20), a(k07), a(k09))
        use EC_Consts
        use EC_ElemCrackingBaseClass
        use mod_parameters
        implicit none

        real NodesCoordx(*), NodesCoordy(*),NodesDispl(*)
        real usi(*),freep(5,*), ym(4,*)
        INTEGER ElemConnect(4,*),ElemMaterial(*),DofIds(2,*)
        INTEGER SolStepCount

        common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
        integer nwebuf,ntime,numnp,neq,ibar,mthsol
        real dn1,dn2

        integer , intent(in)::ElemId
        integer i,j
        real*8 ElemAreaTotal
        integer EdgeStatus(4),mEdgeCount,mUniqueEdge
        INTEGER bDone

        EdgeStatus=-1
        mEdgeCount=1
        bDone=0
        do while(bDone==0)
            if(pEC_ElemData(ElemId)%EdgeStatus(mEdgeCount)==EC_eCrackBothSides) then
                CALL EC_SplitEdge(ElemId,mEdgeCount,NodesCoordx, NodesCoordy, ElemConnect, DofIds, NodesDispl,  &
                        ElemMaterial,usi  , freep , ym,SolStepCount)
                !-- split only one edge per elem per time step
                bDone=1
            endif
        enddo



    end subroutine EC_SplitElem
!##############################################################################
!##############################################################################
    subroutine EC_SplitEdge(ElemId,ElemEdgeId,NodesCoordx, NodesCoordy, ElemConnect, DofIds, NodesDispl, ElemMaterial, &
                            usi  , freep , ym,SolStepCount)
                !-a(k03)    ,a(k04)      ,a(k02)      ,a(k57)  ,a(k18)     , a(k08)      ,a(k20), a(k07), a(k09))
        use EC_Consts
        use EC_ElemCrackingBaseClass
        use mod_parameters
        implicit none

        real NodesCoordx(*), NodesCoordy(*),NodesDispl(*)
        real usi(*),freep(5,*), ym(4,*)
        INTEGER ElemConnect(4,*),ElemMaterial(*),DofIds(2,*)
        INTEGER SolStepCount

        common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
        integer nwebuf,ntime,numnp,neq,ibar,mthsol
        real dn1,dn2

        integer , intent(in)::ElemId,ElemEdgeId
        integer ::ElemId2,ElemEdgeId2,j
        INTEGER mElemEdgeNeighborsList(6)
        ! integer EdgeNode1,EdgeNode2,EdgeNode1N(4),EdgeNode2N(4),j
        type ( EC_ElemCrackingClass ) :: ElemCrackingClassMain(4),ElemCrackingClassOverlapping(4)
        INTEGER bElemValidMain(4),bElemValidOverlapping(4),ElemIdsNewMain(4),ElemIdsNewOverlapping(4)
        INTEGER edgeNodeOld1,edgeNodeOld2,edgeNodeNew1,edgeNodeNew2
        INTEGER bEdgeCracked !-- 0=not cracked, 1= cracked with elems added
        INTEGER nodeIdLocal,edge1,edge2,mMyEdgeId,mElemOther,mEdgeOther

        !- get new and old edge nodes
        edgeNodeOld1=ElemConnect(EC_ElemEdgesConnect(1,ElemEdgeId),ElemId)
        edgeNodeOld2=ElemConnect(EC_ElemEdgesConnect(2,ElemEdgeId),ElemId)
        edgeNodeNew1= EC_NodeCountCurrent+1
        edgeNodeNew2= EC_NodeCountCurrent+2

        bElemValidMain=0
        bElemValidOverlapping=0
        ElemIdsNewMain=0
        ElemIdsNewOverlapping=0
        !-- crack the edge of the current elem,create 2 nodes, add 1 elem
        call EC_CreateElemsTemp(ElemId,ElemEdgeId,ElemConnect,ElemCrackingClassMain(1),ElemCrackingClassOverlapping(1), &
                                edgeNodeNew1,edgeNodeNew2)

        !-- check that the elem are vaild wrt the connectivity
        bElemValidMain(1)= ElemCrackingClassMain(1)%EC_CheckElemValid ()
        !-- save the elem ids, the main and the overlapping
        if(bElemValidMain(1) ==1) ElemIdsNewMain(1)=ElemId
        bElemValidOverlapping(1)= ElemCrackingClassOverlapping(1)%EC_CheckElemValid ()
        if(bElemValidOverlapping(1) ==1) ElemIdsNewOverlapping(1)=EC_ElemCountCurrent+1

        !-- find the edge neighbors and crack their edges
        call pEC_ElemData(ElemId)%EC_GetElemEdgeNeighbors (ElemEdgeId,mElemEdgeNeighborsList)
        do j=1,3    !--each edge can be shared between more than one elem
            if(mElemEdgeNeighborsList((j-1)*2+1) /= 0) then
                ElemId2=mElemEdgeNeighborsList((j-1)*2+1)
                ElemEdgeId2=mElemEdgeNeighborsList((j-1)*2+2)
                call EC_CreateElemsTemp(ElemId2,ElemEdgeId2,ElemConnect,ElemCrackingClassMain(1+j), &
                                        ElemCrackingClassOverlapping(1+j),edgeNodeNew1,edgeNodeNew2)
                !-- check that the elem are vaild wrt the connectivity
                bElemValidMain(1+j)= ElemCrackingClassMain(1+j)%EC_CheckElemValid ()
                if(bElemValidMain(1+j) ==1) ElemIdsNewMain(1+j)=ElemId2

                bElemValidOverlapping(1+j)= ElemCrackingClassOverlapping(1+j)%EC_CheckElemValid ()
                if(bElemValidOverlapping(1+j) ==1) ElemIdsNewOverlapping(1+j)=EC_ElemCountCurrent+1+j
            endif
        enddo
        !- loop over all the temp elems and nodes and add them
        bEdgeCracked=0
        do j=1,4    !--each edge can be shared between more than one elem
            if(bElemValidMain(j) ==1) then
                bEdgeCracked=1
                ElemConnect(1:4,ElemIdsNewMain(j))=ElemCrackingClassMain(j)%iElemConnectivity(1:4)
                CALL EC_CopyElemsData(ElemIdsNewMain(j),ElemIdsNewMain(j),ElemCrackingClassMain(j),ElemCrackingClassMain(j), &
                                    ElemMaterial, freep,ym )
                pEC_ElemData(ElemIdsNewMain(j))%iElemOverlapping=ElemIdsNewOverlapping(j)
            !    !-- update edges neighbors
            !    nodeIdLocal=EC_ElemEdgesConnect(2,ElemEdgeId)
            !    edge1=EC_ElemNodesEdges(1,nodeIdLocal)
            !    edge2=EC_ElemNodesEdges(2,nodeIdLocal)
            !    ElemCrackingClassMain(j)%iElemEdgeNeighbors(:,edge1)=0
            !    ElemCrackingClassMain(j)%iElemEdgeNeighbors(:,edge2)=0
            !    !-- add the new overlapping elem as neighbor to this elem at the opposite edge
            !    !-- this is the opposite edge to the cracked edge.
            !    mMyEdgeId=EC_ElemEdgesOpposite(ElemEdgeId)
            !    mElemOther=EC_ElemCountCurrent+1
            !    mEdgeOther=mMyEdgeId
            !    call ElemCrackingClassMain(j)%EC_ElemAddToNeighbors (mMyEdgeId,mElemOther,mEdgeOther)
            endif
            if(bElemValidOverlapping(j) ==1) then
                ElemConnect(1:4,ElemIdsNewOverlapping(j))=ElemCrackingClassOverlapping(j)%iElemConnectivity(1:4)
                CALL EC_CopyElemsData(ElemIdsNewMain(j),ElemIdsNewOverlapping(j),ElemCrackingClassMain(j), &
                                        ElemCrackingClassOverlapping(j),ElemMaterial, freep,ym )
                pEC_ElemData(ElemIdsNewOverlapping(j))%iElemOverlapping=ElemIdsNewMain(j)
            !    !-- update edges neighbors
            !    nodeIdLocal=EC_ElemEdgesConnect(1,ElemEdgeId)
            !    edge1=EC_ElemNodesEdges(1,nodeIdLocal)
            !    edge2=EC_ElemNodesEdges(2,nodeIdLocal)
            !    ElemCrackingClassOverlapping(j)%iElemEdgeNeighbors(:,edge1)=0
            !    ElemCrackingClassOverlapping(j)%iElemEdgeNeighbors(:,edge2)=0
            !    !-- add the main elem as neighbor to new overlapping elem at the opposite edge
            !    !-- this is the opposite edge to the cracked edge.
            !    mElemOther=ElemId
            !    call ElemCrackingClassOverlapping(j)%EC_ElemAddToNeighbors (mMyEdgeId,mElemOther,mEdgeOther)
            !    !-- copy edge neighbors to the new edges
            !    EC_ElemCountCurrent=EC_ElemCountCurrent+1
            endif
        enddo
        if(bEdgeCracked ==1) then
            !- then add the new nodes
            CALL EC_AddOneNode(edgeNodeOld1,edgeNodeNew1,NodesCoordx, NodesCoordy, DofIds, NodesDispl,usi  , freep , ym)
            CALL EC_AddOneNode(edgeNodeOld2,edgeNodeNew2,NodesCoordx, NodesCoordy, DofIds, NodesDispl,usi  , freep , ym)
        endif

    end subroutine EC_SplitEdge
!##############################################################################
!##############################################################################
    subroutine EC_AddTwoElems(ElemId,ElemEdgeId,ElemCrackingClassMain,ElemCrackingClassOverlapping,EdgeNode1N,EdgeNode2N)

        use EC_ElemCrackingBaseClass
        implicit none
        type ( EC_ElemCrackingClass ), intent(out) :: ElemCrackingClassMain,ElemCrackingClassOverlapping
        integer , intent(in)::EdgeNode1N,EdgeNode2N
        integer , intent(in)::ElemId,ElemEdgeId

        integer nodeIdLocal,edge1,edge2
        integer mMyEdgeId,mElemOther,mEdgeOther

        !-------------------------- create new elems connectivity
        !-- this is the main elem
        !- e1: 1 - 2 - n+2 - 4
        ElemCrackingClassMain%iElemConnectivity(EC_ElemEdgesConnect(2,ElemEdgeId))=edgeNode2N
        ElemCrackingClassMain%iElemOverlapping=EC_eElemOriginMain
        !-- update edges neighbors
        nodeIdLocal=EC_ElemEdgesConnect(2,ElemEdgeId)
        edge1=EC_ElemNodesEdges(1,nodeIdLocal)
        edge2=EC_ElemNodesEdges(2,nodeIdLocal)
        ElemCrackingClassMain%iElemEdgeNeighbors(:,edge1)=0
        ElemCrackingClassMain%iElemEdgeNeighbors(:,edge2)=0

        !-- make this as the overlapping elem
        !- e2: 1 - n+1 - 3 - 4
        ElemCrackingClassOverlapping%iElemConnectivity(EC_ElemEdgesConnect(1,ElemEdgeId))=edgeNode1N
        ElemCrackingClassOverlapping%iElemOverlapping=EC_eElemOriginOverlap
        !-- update edges neighbors
        nodeIdLocal=EC_ElemEdgesConnect(1,ElemEdgeId)
        edge1=EC_ElemNodesEdges(1,nodeIdLocal)
        edge2=EC_ElemNodesEdges(2,nodeIdLocal)
        ElemCrackingClassOverlapping%iElemEdgeNeighbors(:,edge1)=0
        ElemCrackingClassOverlapping%iElemEdgeNeighbors(:,edge2)=0

    end subroutine EC_AddTwoElems
!##############################################################################
!##############################################################################
    subroutine EC_CreateElemsTemp(ElemId,ElemEdgeId,ElemConnect,ElemCrackingClassMain,ElemCrackingClassOverlapping, &
                                    EdgeNode1N,EdgeNode2N)
                !-a(k03)    ,a(k04)      ,a(k02)      ,a(k57)  ,a(k18)     , a(k08)      ,a(k20), a(k07), a(k09))
        use EC_Consts
        use EC_ElemCrackingBaseClass
        use mod_parameters

        implicit none
        type ( EC_ElemCrackingClass ), intent(inout) :: ElemCrackingClassMain,ElemCrackingClassOverlapping
        integer , intent(in)::EdgeNode1N,EdgeNode2N
        integer , intent(in)::ElemId,ElemEdgeId

        INTEGER ElemConnect(4,*)
        !-- add 2 temp elem objects
        !-- split into two quads -- copy the current connect to the new elem
        ElemCrackingClassMain%iElemConnectivity(1:4)=ElemConnect(1:4,ElemId)
        ElemCrackingClassOverlapping%iElemConnectivity(1:4)=ElemConnect(1:4,ElemId)
        CALL EC_AddTwoElems(ElemId,ElemEdgeId,ElemCrackingClassMain,ElemCrackingClassOverlapping,EdgeNode1N,EdgeNode2N)

    end subroutine EC_CreateElemsTemp
!##############################################################################
!##############################################################################
    subroutine EC_AddOneNode(NodeFrom,NodeTo,NodesCoordx, NodesCoordy, DofIds, NodesDispl,usi  , freep , ym)
                !-a(k03)    ,a(k04)      ,a(k02)      ,a(k57)  ,a(k18)     , a(k08)      ,a(k20), a(k07), a(k09))
        use EC_Consts
        use EC_ElemCrackingBaseClass
        use mod_parameters

        implicit none
        integer , intent(in)::NodeFrom,NodeTo

        real NodesCoordx(*), NodesCoordy(*),NodesDispl(*)
        real usi(*),freep(5,*), ym(4,*)
        INTEGER DofIds(2,*)
        !
        common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
        integer nwebuf,ntime,numnp,neq,ibar,mthsol
        real dn1,dn2

        !- copy nodes data to the new nodes
        CALL EC_CopyNodesData(NodeFrom,NodeTo,NodesCoordx, NodesCoordy, DofIds,NodesDispl, usi  , freep , ym)
        !---------------------------1- add Dof
        DofIds(1,NodeTo)=neq+1
        DofIds(2,NodeTo)=neq+2

        numnp=numnp+1
        neq=neq+1
        EC_NodeCountCurrent=EC_NodeCountCurrent+1

    end subroutine EC_AddOneNode
!##############################################################################
!##############################################################################
    subroutine EC_CopyNodesData(NodeFrom,NodeTo,NodesCoordx, NodesCoordy, DofIds,NodesDispl, usi  , freep , ym)
                             !-a(k03)    ,a(k04)      ,a(k02)      ,a(k57  ,a(k18)     , a(k08)      ,a(k20), a(k07), a(k09))
    use mod_parameters
    use CN_Objects_manager

    implicit none
    real NodesCoordx(*), NodesCoordy(*),NodesDispl(*)
    real usi(*),freep(5,*), ym(4,*)
    INTEGER DofIds(2,*)

    integer  , intent(in)::NodeFrom,NodeTo

    NodesCoordx(NodeTo)=NodesCoordx(NodeFrom)
    NodesCoordy(NodeTo)=NodesCoordy(NodeFrom)
    NodesDispl(DofIds(1,NodeTo))=NodesDispl(DofIds(1,NodeFrom))
    NodesDispl(DofIds(2,NodeTo))=NodesDispl(DofIds(2,NodeFrom))
    usi(DofIds(1,NodeTo))=usi(DofIds(1,NodeFrom))
    usi(DofIds(2,NodeTo))=usi(DofIds(2,NodeFrom))
    call CNmanager_CopyNode(NodeFrom,NodeTo)

    end subroutine EC_CopyNodesData
!##############################################################################
!##############################################################################
    subroutine EC_CopyElemsData(ElemFrom,ElemTo,ElemCrackingClassMain,ElemCrackingClassOverlapping,ElemMaterial, freep,ym )
                             !-a(k03)    ,a(k04)      ,a(k02)      ,a(k57  ,a(k18)     , a(k08)      ,a(k20), a(k07), a(k09))
    use mod_parameters
    use CN_Objects_manager
    use EC_ElemCrackingBaseClass

    implicit none
    type ( EC_ElemCrackingClass ), intent(inout) :: ElemCrackingClassMain,ElemCrackingClassOverlapping
    integer, intent(in) :: ElemFrom,ElemTo


    common/wblock8/  abc(573,nume,4), his(573,nume,4)
    real abc,his

    common /wblock12/ Y_modulus(nume),possion_ratio(nume),tau_y(nume)
    real Y_modulus,possion_ratio,tau_y

    common/hokao/    lst,nnn2(nume,4)
    integer lst,nnn2

    common/hourglass/fhg(nume,8),fhghis(nume,8),fhg1(nelemg), &
                   fhg2(nelemg),fhg3(nelemg),fhg4(nelemg),fhg5(nelemg), &
                   fhg6(nelemg),fhg7(nelemg),fhg8(nelemg)
    real fhg,fhghis
    real fhg1,fhg2,fhg3,fhg4,fhg5,fhg6,fhg7,fhg8

    common/hourglass2/hgsstore(nume),hgshis(nume)
    real hgsstore,hgshis

    common/totalenergy/totenerstore(nume),totenerhis(nume),inertener(nume)
    real totenerstore,totenerhis,inertener

    common/hgenergy/hgenerstore(nume),hgenerhis(nume)
    real hgenerstore,hgenerhis

    common/hgstress/hgstress1store(nume),hgstress2store(nume), &
                  hgstress1his(nume),hgstress2his(nume)
    real hgstress1store,hgstress2store,hgstress1his,hgstress2his

    common/hydroembrittle/critfrac(1000), sigfrac0(nume),sigfrac(nume),decfrac(nume)
    real sigfrac0,sigfrac,decfrac
    real critfrac

    real freep(5,*),ym(4,*)
    INTEGER ElemMaterial(*)
    integer j

    freep(1:5,ElemTo)=freep(1:5,ElemFrom)

    ElemMaterial(ElemTo)=ElemMaterial(ElemFrom)
    Y_modulus(ElemTo)=Y_modulus(ElemFrom)
    possion_ratio(ElemTo) = possion_ratio(ElemFrom)
    tau_y(ElemTo) = tau_y(ElemFrom)
    abc(1:573,ElemTo,1)=abc(1:573,ElemFrom,1)
    his(1:573,ElemTo,1)=his(1:573,ElemFrom,1)

    nnn2(ElemTo,1)=nnn2(ElemFrom,1)
    !!--- lumped mass matrix
    do j=1,4
        ym(j,ElemTo)=ym(j,ElemFrom)*ElemCrackingClassMain%EC_GetElemAreaRatio()
        ym(j,ElemFrom)=ym(j,ElemFrom)*ElemCrackingClassOverlapping%EC_GetElemAreaRatio()
    end do

    fhg(ElemTo,1:8)=fhg(ElemFrom,1:8)
    fhghis(ElemTo,1:8)=fhghis(ElemFrom,1:8)

    hgsstore(ElemTo)=hgsstore(ElemFrom)
    hgshis(ElemTo)=hgshis(ElemFrom)

    totenerstore(ElemTo)=totenerstore(ElemFrom)
    totenerhis(ElemTo)=totenerhis(ElemFrom)
    inertener(ElemTo)=inertener(ElemFrom)

    hgenerstore(ElemTo)=hgenerstore(ElemFrom)
    hgenerhis(ElemTo)=hgenerhis(ElemFrom)

    hgstress1store(ElemTo)=hgstress1store(ElemFrom)
    hgstress2store(ElemTo)=hgstress2store(ElemFrom)
    hgstress1his(ElemTo)=hgstress1his(ElemFrom)
    hgstress2his(ElemTo)=hgstress2his(ElemFrom)

    sigfrac0(ElemTo)=sigfrac0(ElemFrom)
    sigfrac(ElemTo)=sigfrac(ElemFrom)
    decfrac(ElemTo)=decfrac(ElemFrom)

    CALL CNmanager_CopyElemnt(ElemFrom,ElemTo)
    CALL pEC_ElemData(ElemFrom)%CopyElem (pEC_ElemData(ElemTo))


    end subroutine EC_CopyElemsData
!##############################################################################
!##############################################################################

end module EC_Objects_manager



! !##############################################################################
! !##############################################################################
!    subroutine EC_SplitElem2(ElemId,NodesCoordx, NodesCoordy, ElemConnect, DofIds, NodesDispl, ElemMaterial, &
!                            usi  , freep , ym,SolStepCount)
!                                 !-a(k03)    ,a(k04)      ,a(k02)      ,a(k57)  ,a(k18)     , a(k08)      ,a(k20), a(k07), a(k09))
!        use EC_Consts
!        use EC_ElemCrackingBaseClass
!        use mod_parameters
!        implicit none
!
!        real NodesCoordx(*), NodesCoordy(*),NodesDispl(*)
!        real usi(*),freep(5,*), ym(4,*)
!        INTEGER ElemConnect(4,*),ElemMaterial(*),DofIds(2,*)
!        INTEGER SolStepCount
!        common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
!        integer nwebuf,ntime,numnp,neq,ibar,mthsol
!        real dn1,dn2
!
!        integer , intent(in)::ElemId
!        integer i,j
!        integer EdgeStatus(4),nodesIds(8)
!        real*8 ElemAreaTotal
!
!        EdgeStatus=-1
!        do i=1,4
!            if(pEC_ElemData(ElemId)%EdgeStatus(i)==EC_eCrackBothSides) then
!                EdgeStatus(i)=1
!            endif
!        enddo
!        if(((EdgeStatus(1)==1) .and. (EdgeStatus(3)==1)) .or. &
!          ((EdgeStatus(2)==1) .and. (EdgeStatus(4)==1))) then
!            !- added nodes map to old nodes
!            !- 1 -- 2 -- 3 -- 4
!            !- 5 -- 6 -- 7 -- 8
!            nodesIds(1:4)=ElemConnect(1:4,ElemId)
!            do i=1,4
!                nodesIds(4+i)=nodesIds(i)+4+EC_NodeCountCurrent
!            enddo
!        endif
!
!        !- case 1: edges 1 & 3 or edges 2 & 4
!        if((EdgeStatus(1)==1) .and. (EdgeStatus(3)==1)) then
!            !-- split into two quads
!            !- e1: 1 - 6 - 7 - 4
!            ElemConnect(1,ElemId)=nodesIds(1)
!            ElemConnect(2,ElemId)=nodesIds(6)
!            ElemConnect(3,ElemId)=nodesIds(7)
!            ElemConnect(4,ElemId)=nodesIds(4)
!            !- e2: 5 - 2 - 3 - 8
!            ElemConnect(1,EC_ElemCountCurrent+1)=nodesIds(5)
!            ElemConnect(2,EC_ElemCountCurrent+1)=nodesIds(2)
!            ElemConnect(3,EC_ElemCountCurrent+1)=nodesIds(3)
!            ElemConnect(4,EC_ElemCountCurrent+1)=nodesIds(8)
!        endif
!        if((EdgeStatus(2)==1) .and. (EdgeStatus(4)==1)) then
!            !-- split into two quads
!            !- e1: 1 - 2 - 7 - 8
!            ElemConnect(1,ElemId)=nodesIds(1)
!            ElemConnect(2,ElemId)=nodesIds(2)
!            ElemConnect(3,ElemId)=nodesIds(7)
!            ElemConnect(4,ElemId)=nodesIds(8)
!            !- e2: 5 - 6 - 3 - 4
!            ElemConnect(1,EC_ElemCountCurrent+1)=nodesIds(5)
!            ElemConnect(2,EC_ElemCountCurrent+1)=nodesIds(6)
!            ElemConnect(3,EC_ElemCountCurrent+1)=nodesIds(3)
!            ElemConnect(4,EC_ElemCountCurrent+1)=nodesIds(4)
!        endif
!
!        !- copy nodes data to the new nodes
!        do i=1,4
!            call EC_CopyNodesData(nodesIds(i),nodesIds(i+4),NodesCoordx, NodesCoordy, DofIds, NodesDispl,usi,freep,ym)
!        enddo
!        !- copy elems data to the new nodes
!        ! call EC_CopyElemsData(ElemId,EC_ElemCountCurrent+1,NodesCoordx, NodesCoordy, ElemConnect, DofIds, NodesDispl, &
!        !                         ElemMaterial, usi  , freep , ym,SolStepCount)
!        ! CALL pEC_ElemData(ElemId)%CopyElem (pEC_ElemData(EC_ElemCountCurrent+1))
!
!        !---------------------------1- add Dof
!        do j=1,4
!            DofIds(1,EC_NodeCountCurrent+j)=neq+2*j-1
!            DofIds(2,EC_NodeCountCurrent+j)=neq+2*j
!        enddo
!        numnp=numnp+4
!        neq=neq+8
!        EC_ElemCountCurrent=EC_ElemCountCurrent+1
!        EC_NodeCountCurrent=EC_NodeCountCurrent+4
!
!    end subroutine EC_SplitElem2
!
