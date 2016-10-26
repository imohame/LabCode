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
    Interface EC_GetElemUnloadingCount;            Module Procedure EC_GetElemUnloadingCount;                 End Interface
    Interface EC_GetElemAreaRatio;            Module Procedure EC_GetElemAreaRatio;                 End Interface
    Interface EC_GetElemSplit;            Module Procedure EC_GetElemSplit;                 End Interface


!##############################################################################
!##############################################################################
    contains
        subroutine EC_AllocateMem(ElemConnect,rNodesCoordx,rNodesCoordy)
            use EC_Consts

            implicit none
            integer,    intent(in) :: ElemConnect(4,*)
            real,     intent(in) :: rNodesCoordx(*),rNodesCoordy(*)
            INTEGER numnpo, numelto
            integer ne,i,j
            real*8 x(4),y(4)

            if(EC_bCracking == 0) then
              return
            endif
            ne=int(EC_ElemCountInput*1.25)
            Allocate(pEC_ElemData(ne))
!!!!!            Allocate(EC_ElemNeighbors(8,ne))
!!!!!            EC_ElemNeighbors=-1
            !-- initialize the elems to default
            do i=1,ne
                CALL pEC_ElemData(i)%Initialize()
            enddo
            !-- calc the elems area
            do i=1,EC_ElemCountInput
                do j=1, 4
                    x(j)=rNodesCoordx(ElemConnect(j,i))
                    y(j)=rNodesCoordy(ElemConnect(j,i))
                    pEC_ElemData(i)%iElemConnectivity(j)=ElemConnect(j,i)
                enddo
                CALL pEC_ElemData(i)%EC_CalcAreaRatio(x,y,4)
            enddo

            ne=int(EC_NodeCountInput*.1) +100
            Allocate(EC_NodesIntersectData(ne))

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
          IF (ALLOCATED (pEC_ElemData))                         DEALLOCATE (pEC_ElemData)
          IF (ALLOCATED (EC_NodesIntersectData))                DEALLOCATE (EC_NodesIntersectData)
!!!          IF (ALLOCATED (EC_ElemNeighbors))              DEALLOCATE (EC_ElemNeighbors)
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
  subroutine EC_MarkElemForCrack(ElemId,rStartPt,ElemCoordx,ElemCoordy,mSolStepCount)
    use EC_Consts
    use EC_ElemCrackingBaseClass
    use mod_file_units

    implicit none
    integer , intent(in)::ElemId,mSolStepCount
    real*8 , intent(in)::rStartPt(2) !- the crack staring point
    real*8 , intent(in)::ElemCoordx(4),ElemCoordy(4)

    integer i,j,ElemIdNeighbor,EdgeIdNeighbor,bIntersect,edgeNode1,edgeNode2
    real*8 Vx,Vy,p1(2),p2(2),q1(2),q2(2),rPtc(2),pRatio,qRatio,rStartPtAct(2)
    real*8 edgeRatio
    integer EdgeCount,node1,node2
    real*8 rP2(2),rP1(2)

    !-- calculate elem area -- no need for it, it's done initialy
    ! call pEC_ElemData(ElemId)%EC_CalcAreaRatio(ElemCoordx,ElemCoordy,4)
    !- this elem is cracked either on cleavage planes, or on the neighboring edges
    !-- set all edges to inner side crack
    pEC_ElemData(ElemId)%EdgeStatus(1:4)=EC_eCrackInnerSide
    !--- set the starting point of the cracking line to the input of the subroutine
    rStartPtAct=rStartPt
    !-- loop over each edge and check the neighbor for possible cracking
    do i=1,4
        EdgeCount=pEC_ElemData(ElemId)%iElemEdgeNeighborsCount(i)
        if( EdgeCount == 0) then
            pEC_ElemData(ElemId)%EdgeStatus(i)=EC_eCrackBothSides
            rStartPtAct=rStartPt
            cycle ! go to next i, edge
        endif

        do j=1,EdgeCount    !--each edge can be shared between more than one elem
!!!---------------------------------debugging
!!write(*,*)ElemId,pEC_ElemData(ElemId)%EdgeStatus
!!write(*,*)pEC_ElemData(ElemId)%rCoordRatioCracking
!!!---------------------------------debugging
            ElemIdNeighbor=pEC_ElemData(ElemId)%iElemEdgeNeighbors((j-1)*2+1,i)
            EdgeIdNeighbor=pEC_ElemData(ElemId)%iElemEdgeNeighbors((j-1)*2+2,i)
!!!---------------------------------debugging
!!write(*,*)ElemIdNeighbor,pEC_ElemData(ElemIdNeighbor)%EdgeStatus
!!write(*,*)pEC_ElemData(ElemIdNeighbor)%rCoordRatioCracking
!!!---------------------------------debugging
            !-- if the edge is outter, has no neighbors then make it EC_eCrackBothSides
           !-if the neighboring edge is inner cracked then set this edge on both sides as cracked
            if(pEC_ElemData(ElemIdNeighbor)%EdgeStatus(EdgeIdNeighbor)==EC_eCrackInnerSide) then
                pEC_ElemData(ElemId)%EdgeStatus(i)=EC_eCrackBothSides
                pEC_ElemData(ElemIdNeighbor)%EdgeStatus(EdgeIdNeighbor)=EC_eCrackBothSides
                !-- check if there is edge with cracking pint then make it the starting point
                !-- when there is a semi cracked elem or elem has one cracked edge, one of the edges is shared among 3 elems
                edgeRatio=pEC_ElemData(ElemIdNeighbor)%rCoordRatioCracking(EdgeIdNeighbor)
                if(edgeRatio>0.0) then
                    !-- the ratio for the current edge is 1- neighbor edge ratio
                    edgeRatio=1-edgeRatio
                    node1=EC_ElemEdgesConnect(1,i)
                    node2=EC_ElemEdgesConnect(2,i)

                    rP1(1)=ElemCoordx(node1)
                    rP1(2)=ElemCoordy(node1)
                    rP2(1)=ElemCoordx(node2)
                    rP2(2)=ElemCoordy(node2)
                    CALL GetPtOnLine(rP1,rP2,edgeRatio,rStartPtAct)
                endif
            endif
        enddo
    enddo
    !- this the crack direction
    !- since vx,vy are the components of the normal then the crack dir are vy,-vx
    !- the crack direction is normal to the cleavage dir.
    Vx=pEC_ElemData(ElemId)%rCleavagePlane(3)
    Vy=-pEC_ElemData(ElemId)%rCleavagePlane(2)

    !- make a crack line given a starting point and a direction
    CALL GetLineStartEndPts(rStartPtAct,Vx,Vy,EC_ElemAverageDia*2,p1,p2)
    !- intersect the crack line with the edges
    write(iFU_frac_crackEdgeSplit,*) '------------------------------------------',ElemId
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
            !--- elem id , step no, edge no, edge node 1, edge node 2, intersection ratio, edge status
            write(iFU_frac_crackEdgeSplit,*) ElemId,mSolStepCount,i,edgeNode1,edgeNode2,qRatio,pEC_ElemData(ElemId)%EdgeStatus(i)
            !-- this for the neighboring edges
!!            if(pEC_ElemData(ElemId)%EdgeStatus(i)==EC_eCrackBothSides) then !- then set this value to the other edge
                EdgeCount=pEC_ElemData(ElemId)%iElemEdgeNeighborsCount(i)
                do j=1,EdgeCount    !--each edge can be shared between more than one elem
                    ElemIdNeighbor=pEC_ElemData(ElemId)%iElemEdgeNeighbors((j-1)*2+1,i)
                    EdgeIdNeighbor=pEC_ElemData(ElemId)%iElemEdgeNeighbors((j-1)*2+2,i)
                    !- 1-qRatio because the edge nodes order is reversed from elem to the neighboring elem
                    pEC_ElemData(ElemIdNeighbor)%rCoordRatioCracking(EdgeIdNeighbor)=1-qRatio
!!!---------------------------------debugging
!!write(*,*)ElemIdNeighbor,pEC_ElemData(ElemIdNeighbor)%EdgeStatus
!!write(*,*)pEC_ElemData(ElemIdNeighbor)%rCoordRatioCracking
!!!---------------------------------debugging
                    write(iFU_frac_crackEdgeSplit,*) ElemIdNeighbor,mSolStepCount,EdgeIdNeighbor, &
                                                1-qRatio,pEC_ElemData(ElemIdNeighbor)%EdgeStatus(EdgeIdNeighbor)
                enddo
!!            endif
        endif
    enddo

  end subroutine EC_MarkElemForCrack

 !##############################################################################
 !##############################################################################
    subroutine EC_SplitElem(ElemId,NodesCoordx, NodesCoordy, ElemConnect, DofIds, NodesDispl, ElemMaterial, &
                            usi  , freep , ym,mSolStepCount)
                !-a(k03)    ,a(k04)      ,a(k02)      ,a(k57)  ,a(k18)     , a(k08)      ,a(k20), a(k07), a(k09))
        use EC_Consts
        use EC_ElemCrackingBaseClass
        use mod_parameters
        use mod_file_units
        implicit none

        real NodesCoordx(*), NodesCoordy(*),NodesDispl(*)
        real usi(*),freep(5,*), ym(4,*)
        INTEGER ElemConnect(4,*),ElemMaterial(*),DofIds(2,*)
        INTEGER mSolStepCount

        common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
        INTEGER nprnt,mprint,itmpop,numelt,jprint,idump,locstr

        common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
        integer nwebuf,ntime,numnp,neq,ibar,mthsol
        real dn1,dn2

        integer , intent(in)::ElemId
        integer i,j,k
        real*8 ElemAreaTotal,xs(4),ys(4)
        integer EdgeStatus(4),mEdgeCount,mUniqueEdge
        INTEGER bDone,EC_ElemCountCurrent_old
        integer :: mElemIdsList1(16) !-- find the neighbors of mElemIdsList1 in mElemIdsList2
        integer :: mElemIdsList2(16) !-- assume any elem has 16 neghbors
        integer :: ElemIdsNewMain(4),ElemIdsNewOverlapping(4)
        integer :: mElemCount1,mElemCount2 !-- count the actual neighobrs
        integer :: e1,EdgeCount,ElemIdtemp

        mElemIdsList1=0
        mElemIdsList2=0
        mElemCount1=0
        mElemCount2=0


        EdgeStatus=-1
        mEdgeCount=0
        bDone=0
        EC_ElemCountCurrent_old=EC_ElemCountCurrent
        do while(bDone==0 .and. mEdgeCount<5)
            mEdgeCount=mEdgeCount+1
            if(pEC_ElemData(ElemId)%EdgeStatus(mEdgeCount)==EC_eCrackBothSides) then
!!                write(*,*)ElemId,pEC_ElemData(ElemId)%iElemStatus,pEC_ElemData(ElemId)%iElemSplit
                CALL EC_SplitEdge(ElemId,mEdgeCount,NodesCoordx, NodesCoordy, ElemConnect, DofIds, NodesDispl,  &
                        ElemMaterial,usi  , freep , ym,mSolStepCount,ElemIdsNewMain,ElemIdsNewOverlapping)
                !-- reset the elem
                pEC_ElemData(ElemId)%iElemSplit=1
!!                write(*,*)ElemId,pEC_ElemData(ElemId)%iElemStatus,pEC_ElemData(ElemId)%iElemSplit
                !-- split only one edge per elem per time step
                bDone=1
            endif
        enddo
        !--- this is to update the elem neighbors due to the creation of new elem
        EC_UpdateMesh=0 !-to update some parameters and redo the step again
        if(bDone==1) then
            EC_UpdateMesh=1 !-to update some parameters and redo the step again
            !-- find this elem neghbors
            !-- add the modified elems and the new elems
            do i=1,4
                ElemIdtemp=ElemIdsNewMain(i)
                if(ElemIdtemp>0) then
                    mElemCount1=mElemCount1+1
                    mElemIdsList1(mElemCount1)=ElemIdtemp
                    mElemCount2=mElemCount2+1
                    mElemIdsList2(mElemCount2)=ElemIdtemp
                endif
            enddo
            do i=1,4
                ElemIdtemp=ElemIdsNewOverlapping(i)
                if(ElemIdtemp>0) then
                    mElemCount1=mElemCount1+1
                    mElemIdsList1(mElemCount1)=ElemIdtemp
                    mElemCount2=mElemCount2+1
                    mElemIdsList2(mElemCount2)=ElemIdtemp
                endif
            enddo
!!!!            mElemCount1=mElemCount1+1
!!!!            mElemIdsList1(mElemCount1)=ElemId
!!!!            mElemCount2=mElemCount2+1
!!!!            mElemIdsList2(mElemCount2)=ElemId
            !-- get the niehgbors of the modified/new elems
            do k=1,mElemCount1    !--
                ElemIdtemp=mElemIdsList1(k)
                do i=1,4    !--4 edges
                    EdgeCount=pEC_ElemData(ElemIdtemp)%iElemEdgeNeighborsCount(i)
                    do j=1,EdgeCount    !--each edge can be shared between more than one elem
                        mElemCount2=mElemCount2+1
                        mElemIdsList2(mElemCount2)=pEC_ElemData(ElemIdtemp)%iElemEdgeNeighbors((j-1)*2+1,i)
                    enddo
                enddo
            enddo

!!!        !-- add the new overlapping elems to the list
!!!            do i=1,EC_ElemCountCurrent-EC_ElemCountCurrent_old
!!!                mElemCount1=mElemCount1+1
!!!                mElemIdsList1(mElemCount1)=EC_ElemCountCurrent_old+i
!!!                mElemCount2=mElemCount2+1
!!!                mElemIdsList2(mElemCount2)=EC_ElemCountCurrent_old+i
!!!            enddo
            !-- update edges neighbors
            CALL EC_CalcElemNeighborsInList(mElemIdsList1,mElemCount1,mElemIdsList2,mElemCount2)
            !-- add tracking info
            write(iFU_frac_crackOverlapInfo,*) '------------------------------------------step#',mSolStepCount
            !- elem no, step no, elem connect, elem neighbors
            do i=1,mElemCount1
                e1=mElemIdsList1(i)
                pEC_ElemData(e1)%iElemSplit=1
            enddo
            do i=1,mElemCount2
                e1=mElemIdsList2(i)
                write(iFU_frac_crackOverlapInfo,*) e1,pEC_ElemData(e1)%iElemConnectivity,pEC_ElemData(e1)%EC_GetMyAreaRatio() &
                                                    ,pEC_ElemData(e1)%iElemSplit
                do j=1,4
                    write(iFU_frac_crackOverlapInfo,*) pEC_ElemData(e1)%iElemEdgeNeighbors(:,j)
                enddo
                do j=1, 4
                    xs(j)=NodesCoordx(ElemConnect(j,e1))+NodesDispl(DofIds(1,ElemConnect(j,e1)))
                    ys(j)=NodesCoordy(ElemConnect(j,e1))+NodesDispl(DofIds(2,ElemConnect(j,e1)))
                enddo
                write(iFU_frac_crackOverlapInfo,*) xs
                write(iFU_frac_crackOverlapInfo,*) ys
                write(iFU_frac_crackOverlapInfo,*) '-------------------'
            enddo
        !!!================================================================
        !!!================================================================
        !!!========= update this global count b/c it's used in the rest of the code
        !!!================================================================
        !!!================================================================
            numelt=EC_ElemCountCurrent
        endif



    end subroutine EC_SplitElem
!##############################################################################
!##############################################################################
    subroutine EC_SplitEdge(ElemId,ElemEdgeId,NodesCoordx, NodesCoordy, ElemConnect, DofIds, NodesDispl, ElemMaterial, &
                            usi  , freep , ym,mSolStepCount,ElemIdsNewMain,ElemIdsNewOverlapping)
                !-a(k03)    ,a(k04)      ,a(k02)      ,a(k57)  ,a(k18)     , a(k08)      ,a(k20), a(k07), a(k09))
        use EC_Consts
        use EC_ElemCrackingBaseClass
        use mod_parameters
        use mod_file_units
        implicit none

        real NodesCoordx(*), NodesCoordy(*),NodesDispl(*)
        real usi(*),freep(5,*), ym(4,*)
        INTEGER ElemConnect(4,*),ElemMaterial(*),DofIds(2,*)
        INTEGER mSolStepCount

        common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
        integer nwebuf,ntime,numnp,neq,ibar,mthsol
        real dn1,dn2

        integer , intent(in)::ElemId,ElemEdgeId
        integer ::ElemId2,ElemEdgeId2,j,i,ElemIdtemp
        integer , intent(inout)::ElemIdsNewMain(4),ElemIdsNewOverlapping(4)
        ! integer EdgeNode1,EdgeNode2,EdgeNode1N(4),EdgeNode2N(4),j
        type ( EC_ElemCrackingClass ) :: ElemCrackingClassMain(4),ElemCrackingClassOverlapping(4)
        INTEGER bElemValidMain(4),bElemValidOverlapping(4)
        INTEGER edgeNodeOld1,edgeNodeOld2,edgeNodeNew1,edgeNodeNew2
        INTEGER bEdgeCracked !-- 0=not cracked, 1= cracked with elems added
        INTEGER nodeIdLocal,edge1,edge2,mMyEdgeId,mElemOther,mEdgeOther,e1,EdgeCount
        real*8 xs(4),ys(4)


        !- get new and old edge nodes
        edgeNodeOld1=ElemConnect(EC_ElemEdgesConnect(1,ElemEdgeId),ElemId)
        edgeNodeOld2=ElemConnect(EC_ElemEdgesConnect(2,ElemEdgeId),ElemId)
        edgeNodeNew1= EC_NodeCountCurrent+1
        edgeNodeNew2= EC_NodeCountCurrent+2

        bElemValidMain=0
        bElemValidOverlapping=0
        ElemIdsNewMain=0
        ElemIdsNewOverlapping=0
        !- initialize the temp. elems
        do j=1,4
            CALL ElemCrackingClassMain(j)%Initialize()
            call pEC_ElemData(ElemId)%EC_CopyElem (ElemCrackingClassMain(j))
            CALL ElemCrackingClassOverlapping(j)%Initialize()
            call pEC_ElemData(ElemId)%EC_CopyElem (ElemCrackingClassOverlapping(j))
        enddo
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
!!!        call pEC_ElemData(ElemId)%EC_GetElemEdgeNeighbors (ElemEdgeId,mElemEdgeNeighborsList)
        EdgeCount=pEC_ElemData(ElemId)%iElemEdgeNeighborsCount(ElemEdgeId)
        do j=1,EdgeCount    !--each edge can be shared between more than one elem
            ElemId2    =pEC_ElemData(ElemId)%iElemEdgeNeighbors((j-1)*2+1,ElemEdgeId)
            ElemEdgeId2=pEC_ElemData(ElemId)%iElemEdgeNeighbors((j-1)*2+2,ElemEdgeId)
            !-- note that the order of the edge nodes are the opposite for the neighbors
            call EC_CreateElemsTemp(ElemId2,ElemEdgeId2,ElemConnect,ElemCrackingClassMain(1+j), &
                                    ElemCrackingClassOverlapping(1+j),edgeNodeNew2,edgeNodeNew1)
            !-- check that the elem are valid wrt the connectivity
            bElemValidMain(1+j)= ElemCrackingClassMain(1+j)%EC_CheckElemValid ()
            if(bElemValidMain(1+j) ==1) ElemIdsNewMain(1+j)=ElemId2

            bElemValidOverlapping(1+j)= ElemCrackingClassOverlapping(1+j)%EC_CheckElemValid ()
            if(bElemValidOverlapping(1+j) ==1) ElemIdsNewOverlapping(1+j)=EC_ElemCountCurrent+1+j
            !-- in case the main is not valid but the overlapping is valid, swap them
            if((bElemValidOverlapping(1+j) ==1) .and. (bElemValidMain(1+j) ==0)) then
                call ElemCrackingClassOverlapping(1+j)%EC_CopyElem (ElemCrackingClassMain(1+j))
                ElemIdsNewMain(1+j)=ElemId2
                bElemValidMain(1+j)=1
                !----------
                ElemIdsNewOverlapping(1+j)=0
                bElemValidOverlapping(1+j)=0
            endif
        enddo
        !- loop over all the temp elems and nodes and add them
!!!        write(iFU_frac_crackOverlapInfo,*) '------------------------------------------step#',mSolStepCount
        !-- add the nodes first, b/c the connectivity uses it
        bEdgeCracked=0
        do j=1,4    !--each edge can be shared between more than one elem
            if(bElemValidMain(j) ==1) then
                bEdgeCracked=1
                exit
            endif
        enddo

        if(bEdgeCracked ==1) then
            !- then add the new nodes --- EC_NodeCountCurrent already incremented
            CALL EC_AddOneNode(edgeNodeOld1,edgeNodeNew1,NodesCoordx, NodesCoordy, DofIds, NodesDispl,usi  , freep , ym)
            !-- --- EC_NodeCountCurrent already incremented
            CALL EC_AddOneNode(edgeNodeOld2,edgeNodeNew2,NodesCoordx, NodesCoordy, DofIds, NodesDispl,usi  , freep , ym)
        endif

        do j=1,4    !--each edge can be shared between more than one elem
            if(bElemValidMain(j) ==1) then
                ElemIdtemp=ElemIdsNewMain(j)
                !- to get the elem nodes coordinates
                do i=1, 4
                    xs(i)=NodesCoordx(ElemConnect(i,ElemIdtemp))+NodesDispl(DofIds(1,ElemConnect(i,ElemIdtemp)))
                    ys(i)=NodesCoordy(ElemConnect(i,ElemIdtemp))+NodesDispl(DofIds(2,ElemConnect(i,ElemIdtemp)))
                enddo
                CALL ElemCrackingClassMain(j)%EC_CalcCrackedElemAreaRatio (xs,ys,4,ElemConnect(:,ElemIdtemp))

                !-- copy elem in EC_BaseClass copies also the connectivity, ... update that in the next two line
                CALL EC_CopyElemsData(ElemIdsNewMain(j),ElemIdsNewMain(j),ElemCrackingClassMain(j),ElemCrackingClassMain(j), &
                                    ElemMaterial, freep,ym )
                !-- ... update that in the next two line
                ElemConnect(1:4,ElemIdtemp)=ElemCrackingClassMain(j)%iElemConnectivity(1:4)
                pEC_ElemData(ElemIdsNewMain(j))%iElemConnectivity(1:4)=ElemCrackingClassMain(j)%iElemConnectivity(1:4)
                pEC_ElemData(ElemIdsNewMain(j))%rAreaRatio=ElemCrackingClassMain(j)%rAreaRatio

                !-- set the related elem to overlap
                pEC_ElemData(ElemIdsNewMain(j))%iElemOverlapping=ElemIdsNewOverlapping(j)
                !-- track the split elem in the file iFU_frac_crackOverlapInfo
!!!!                e1=ElemIdsNewMain(j)
!!!!                write(iFU_frac_crackOverlapInfo,*) e1,mSolStepCount,pEC_ElemData(e1)%iElemConnectivity, &
!!!!                                                    pEC_ElemData(e1)%iElemEdgeNeighbors

                !-- DO NOT increment EC_ElemCountCurrent ... it's the main elem
!!!!!!                EC_ElemCountCurrent=ElemIdsNewMain(j)
            endif
            if(bElemValidOverlapping(j) ==1) then
                ElemIdtemp=ElemIdsNewOverlapping(j)
                !-- copy elem in EC_BaseClass copies also the connectivity, ... update that in the next two line
                CALL EC_CopyElemsData(ElemIdsNewMain(j),ElemIdtemp,ElemCrackingClassMain(j), &
                                        ElemCrackingClassOverlapping(j),ElemMaterial, freep,ym )
                !-- ... update that in the next two line
                ElemConnect(1:4,ElemIdtemp)=ElemCrackingClassOverlapping(j)%iElemConnectivity(1:4)
                pEC_ElemData(ElemIdsNewOverlapping(j))%iElemConnectivity(1:4)=ElemCrackingClassOverlapping(j)%iElemConnectivity(1:4)
                !- has to be done here after the connectivity has been done
                !- to get the elem nodes coordinates
                do i=1, 4
                    xs(i)=NodesCoordx(ElemConnect(i,ElemIdtemp))+NodesDispl(DofIds(1,ElemConnect(i,ElemIdtemp)))
                    ys(i)=NodesCoordy(ElemConnect(i,ElemIdtemp))+NodesDispl(DofIds(2,ElemConnect(i,ElemIdtemp)))
                enddo
                !-- calc the new elem area ratio, arrange the new elem connectivity
                CALL ElemCrackingClassOverlapping(j)%EC_CalcCrackedElemAreaRatio (xs,ys,4,ElemConnect(:,ElemIdtemp))
                pEC_ElemData(ElemIdsNewOverlapping(j))%rAreaRatio=ElemCrackingClassOverlapping(j)%rAreaRatio
                !-- set the related elem to main
                pEC_ElemData(ElemIdsNewOverlapping(j))%iElemOverlapping=ElemIdsNewMain(j)
                !====================================================================
                !====================================================================
                !================ here where to update the EC_ElemCountCurrent
                !====================================================================
                !-- increment EC_ElemCountCurrent
                EC_ElemCountCurrent=ElemIdsNewOverlapping(j)
                !-- track the split elem in the file iFU_frac_crackOverlapInfo
!!!!                e1=ElemIdsNewOverlapping(j)
!!!!                write(iFU_frac_crackOverlapInfo,*) e1,mSolStepCount,pEC_ElemData(e1)%iElemConnectivity, &
!!!!                                                    pEC_ElemData(e1)%iElemEdgeNeighbors
            endif
        enddo

    end subroutine EC_SplitEdge
!##############################################################################
!##############################################################################
    subroutine EC_AddTwoElems(ElemId,ElemEdgeId,ElemCrackingClassMain,ElemCrackingClassOverlapping,EdgeNode1N,EdgeNode2N)

        use EC_ElemCrackingBaseClass
        implicit none
        type ( EC_ElemCrackingClass ), intent(inout) :: ElemCrackingClassMain,ElemCrackingClassOverlapping
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
        edge1=EC_ElemEdgesConnect(1,nodeIdLocal)
        edge2=EC_ElemEdgesConnect(2,nodeIdLocal)
!!        ElemCrackingClassMain%iElemEdgeNeighbors(:,edge1)=0
!!        ElemCrackingClassMain%iElemEdgeNeighbors(:,edge2)=0

        !-- make this as the overlapping elem
        !- e2: 1 - n+1 - 3 - 4
        ElemCrackingClassOverlapping%iElemConnectivity(EC_ElemEdgesConnect(1,ElemEdgeId))=edgeNode1N
        ElemCrackingClassOverlapping%iElemOverlapping=EC_eElemOriginOverlap
        !-- update edges neighbors
        nodeIdLocal=EC_ElemEdgesConnect(1,ElemEdgeId)
        edge1=EC_ElemEdgesConnect(1,nodeIdLocal)
        edge2=EC_ElemEdgesConnect(2,nodeIdLocal)
!!        ElemCrackingClassOverlapping%iElemEdgeNeighbors(:,edge1)=0
!!        ElemCrackingClassOverlapping%iElemEdgeNeighbors(:,edge2)=0

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
        common/fissn0/maxneq,mwspac,ntpe0,ntpe1,nfissl(3)
        integer maxneq,mwspac,ntpe0,ntpe1,nfissl

        common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
        integer nwebuf,ntime,numnp,neq,ibar,mthsol
        real dn1,dn2

        !- copy nodes data to the new nodes
        CALL EC_CopyNodesData(NodeFrom,NodeTo,NodesCoordx, NodesCoordy, DofIds,NodesDispl, usi  , freep , ym)
    !!!================================================================
    !!!================================================================
    !!!========= update this global count b/c it's used in the rest of the code
    !!!================================================================
    !!!================================================================
    !---------------------------1- add Dof
        DofIds(1,NodeTo)=neq+1
        DofIds(2,NodeTo)=neq+2

        numnp=numnp+1
        neq=neq+2
        maxneq=neq
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
    real*8 :: RetVal8(4)
    real*4 :: RetVal4(4)

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
        ym(j,ElemTo)=ym(j,ElemFrom)*ElemCrackingClassMain%rAreaRatio
        ym(j,ElemFrom)=ym(j,ElemFrom)*ElemCrackingClassOverlapping%rAreaRatio
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
    CALL pEC_ElemData(ElemFrom)%EC_CopyElem (pEC_ElemData(ElemTo))

    !--- copy stress
    RetVal4=0
    RetVal8=0
    CALL CNmanager_Get_sigalt(ElemFrom,RetVal8)
    RetVal4=RetVal8*ElemCrackingClassOverlapping%rAreaRatio
    CALL CNmanager_Set_sigalt(ElemTo,RetVal4)
    !--- try to set the stress to half its value --- for testing
    RetVal4=RetVal8*ElemCrackingClassMain%rAreaRatio
    CALL CNmanager_Set_sigalt(ElemFrom,RetVal4)

    end subroutine EC_CopyElemsData
!##############################################################################
    subroutine EC_CalcElemNeighbors(mmElemIdsList,mmElemCount,bAllElems)

    implicit none
    integer, intent(in) :: mmElemCount
    integer, intent(in) :: bAllElems !-1= all elem, 0=use the list
    integer, intent(in) :: mmElemIdsList(mmElemCount)

    integer i,j,mElemCount,ElemIdi,ElemIdj
    integer edgei,edgej,nodei1,nodei2,nodej1,nodej2
    integer, ALLOCATABLE::mElemIdsList(:)
    integer ierr,mlist(8)
    integer*4 setvbuf3f_local

    !-------------------------- allocate the new list if all elems
    if(bAllElems==1)then
        mElemCount=EC_ElemCountCurrent
        allocate(mElemIdsList(mElemCount))
        do i=1,mElemCount
            mElemIdsList(i)=i
        enddo
    endif
    !----------------- set all neighbors to zero
    do i=1,mElemCount
        if(bAllElems==1)then
            ElemIdi=mElemIdsList(i) !- use all elems, the new list
        else
            ElemIdi=mmElemIdsList(i)  !- use the list
        endif
        pEC_ElemData(ElemIdi)%iElemEdgeNeighbors=0
        pEC_ElemData(ElemIdi)%iElemEdgeNeighborsCount=0
    enddo

    !---------------------------- loop over the elems to find the neighbors
    do i=1,mElemCount
        if(bAllElems == 1)then
            ElemIdi=mElemIdsList(i) !- use all elems, the new list
        else
            ElemIdi=mmElemIdsList(i)  !- use the list
        endif
        !-- loop over the edges of the ElemIdi
        do edgei=1,4

            if(pEC_ElemData(ElemIdi)%iElemEdgeNeighborsCount(edgei) > 0) cycle !- goto next loop i

            nodei1=pEC_ElemData(ElemIdi)%iElemConnectivity(EC_ElemEdgesConnect(1,edgei))
            nodei2=pEC_ElemData(ElemIdi)%iElemConnectivity(EC_ElemEdgesConnect(2,edgei))
            !-- loop over the other elems
            do j=1,mElemCount

                if(bAllElems == 1)then
                    ElemIdj=mElemIdsList(j) !- use all elems, the new list
                else
                    ElemIdj=mmElemIdsList(j)  !- use the list
                endif

                if(ElemIdi == ElemIdj) cycle !- goto next loop j
                !-- loop over the edges of the ElemIdj
                do edgej=1,4
                    nodej1=pEC_ElemData(ElemIdj)%iElemConnectivity(EC_ElemEdgesConnect(1,edgej))
                    nodej2=pEC_ElemData(ElemIdj)%iElemConnectivity(EC_ElemEdgesConnect(2,edgej))
                    if((nodei1==nodej1 .and. nodei2==nodej2) .or.(nodei1==nodej2 .and. nodei2==nodej1))then
                        pEC_ElemData(ElemIdi)%iElemEdgeNeighbors(1,edgei)=ElemIdj
                        pEC_ElemData(ElemIdi)%iElemEdgeNeighbors(2,edgei)=edgej
                        pEC_ElemData(ElemIdi)%iElemEdgeNeighborsCount(edgei)= 1

                        pEC_ElemData(ElemIdj)%iElemEdgeNeighbors(1,edgej)=ElemIdi
                        pEC_ElemData(ElemIdj)%iElemEdgeNeighbors(2,edgej)=edgei
                        pEC_ElemData(ElemIdj)%iElemEdgeNeighborsCount(edgej)= 1
                    endif
                enddo !-do edgej=1,4
            enddo !-do j=1,mElemCount
        enddo !-do edgei=1,4
    enddo !-- do i=1,mElemCount
    !----------------- print output file ElemsNeighbors2.out
    if(bAllElems == 1)then
        open(7011, file = 'ElemsNeighbors3.out', status = 'unknown')
        ierr=setvbuf3f_local(7011,1,100)
        do i=1,mElemCount
            do j=1,4
                mlist(j)=pEC_ElemData(i)%iElemEdgeNeighbors(1,j)
                mlist(j+4)=pEC_ElemData(i)%iElemEdgeNeighbors(2,j)
            enddo
            write (7011, *) mlist(1:8)
        enddo
        close(7011)
    endif

    !------------------------------- clean the memory
    if(bAllElems==1)then
        deallocate(mElemIdsList)
    endif

    end subroutine EC_CalcElemNeighbors
!##############################################################################
    subroutine EC_CalcElemNeighborsInList(mElemIdsList1,mElemCount1,mElemIdsList2,mElemCount2)

    implicit none
    integer, intent(in) :: mElemCount1,mElemCount2
    integer, intent(in) :: mElemIdsList1(mElemCount1),mElemIdsList2(mElemCount2)

    integer i,j,mElemCount,ElemIdi,ElemIdj,e1,jj
    integer edgei,edgej,nodei1,nodei2,nodej1,nodej2



    do i=1,mElemCount1
        ElemIdi=mElemIdsList1(i) !- use all elems, the new list
        !----------------- set all neighbors to zero
        pEC_ElemData(ElemIdi)%iElemEdgeNeighbors=0
        pEC_ElemData(ElemIdi)%iElemEdgeNeighborsCount=0
        !----------------- remove this elem from all the neighbors reference
        do j=1,mElemCount2
            ElemIdj=mElemIdsList2(j) !- use all elems, the new list
            if(ElemIdi == ElemIdj) cycle !-- go to next loop item
            !- remove the elem i from the neighbors of elem j
            call pEC_ElemData(ElemIdj)%EC_ElemRemoveFromNeighbors(ElemIdi)
        enddo
    enddo
!!!!    write(*,*) '---------------------------------------------------------->start'
!!!!    write(*,*) '---------------------------------------------------------->'
!!!!    do i=1,mElemCount2
!!!!        e1=mElemIdsList2(i) !- use all elems, the new list
!!!!        write(*,*) e1,'--->',pEC_ElemData(e1)%iElemConnectivity
!!!!        do j=1,4
!!!!            write(*,*) pEC_ElemData(e1)%iElemEdgeNeighbors(:,j)
!!!!        enddo
!!!!    enddo
   !---------------------------- loop over the elems to find the neighbors
    do i=1,mElemCount1
        ElemIdi=mElemIdsList1(i) !- use all elems, the new list
        !-- loop over the edges of the ElemIdi
        do edgei=1,4
!            -- no need for this if b/c there could be multiple elems at this edge
!!            if(pEC_ElemData(ElemIdi)%iElemEdgeNeighbors(1,edgei) == 0) then
                nodei1=pEC_ElemData(ElemIdi)%iElemConnectivity(EC_ElemEdgesConnect(1,edgei))
                nodei2=pEC_ElemData(ElemIdi)%iElemConnectivity(EC_ElemEdgesConnect(2,edgei))
                !-- loop over the other elems
                do j=1,mElemCount2
                    ElemIdj=mElemIdsList2(j) !- use all elems, the new list
                    if(ElemIdi /= ElemIdj) then
!!!---------------------------------------------------------------------- debugging
!!write(*,*) ElemIdi,edgei,'--------------------------------------------->'
!!write(*,*) ElemIdi,'--->',pEC_ElemData(ElemIdi)%iElemConnectivity
!!write(*,*) ElemIdj,'--->',pEC_ElemData(ElemIdj)%iElemConnectivity
!!write(*,*) ElemIdi,'--->',pEC_ElemData(ElemIdi)%iElemEdgeNeighbors
!!write(*,*) ElemIdj,'--->',pEC_ElemData(ElemIdj)%iElemEdgeNeighbors
!!!---------------------------------------------------------------------- debugging
                        !-- loop over the edges of the ElemIdj
                        do edgej=1,4
                            nodej1=pEC_ElemData(ElemIdj)%iElemConnectivity(EC_ElemEdgesConnect(1,edgej))
                            nodej2=pEC_ElemData(ElemIdj)%iElemConnectivity(EC_ElemEdgesConnect(2,edgej))
                            if((nodei1==nodej1 .and. nodei2==nodej2) .or.(nodei1==nodej2 .and. nodei2==nodej1))then
                                call pEC_ElemData(ElemIdi)%EC_ElemAddToNeighbors(edgei,ElemIdj,edgej)
                                call pEC_ElemData(ElemIdj)%EC_ElemAddToNeighbors(edgej,ElemIdi,edgei)
!!!!!!---------------------------------------------------------------------- debugging
!!!                                write(*,*) '---------------------------------------------------------->mid'
!!!                                write(*,*) ElemIdi,'--->',pEC_ElemData(ElemIdi)%iElemConnectivity
!!!                                do jj=1,4
!!!                                    write(*,*) pEC_ElemData(ElemIdi)%iElemEdgeNeighbors(:,jj)
!!!                                enddo
!!!                                write(*,*) ElemIdj,'--->',pEC_ElemData(ElemIdj)%iElemConnectivity
!!!                                do jj=1,4
!!!                                    write(*,*) pEC_ElemData(ElemIdj)%iElemEdgeNeighbors(:,jj)
!!!                                enddo
!!!!!!---------------------------------------------------------------------- debugging
                                exit
                            endif
                        enddo !-do edgej=1,4
                    endif !-if(ElemIdi/=ElemIdj) then
                enddo !-do j=1,mElemCount
!            endif !- if(pEC_ElemData(ElemIdi)%iElemEdgeNeighbors(1,edgei)==0) then
        enddo !-do edgei=1,4
    enddo !-- do i=1,mElemCount
!!!---------------------------------------------------------------------- debugging
!!!---------------------------------------------------------------------- debugging
!!!!    write(*,*) '---------------------------------------------------------->'
!!!!    write(*,*) '---------------------------------------------------------->end'
!!!!    do i=1,mElemCount2
!!!!        e1=mElemIdsList2(i) !- use all elems, the new list
!!!!        write(*,*) e1,'--->',pEC_ElemData(e1)%iElemConnectivity
!!!!        do j=1,4
!!!!            write(*,*) pEC_ElemData(e1)%iElemEdgeNeighbors(:,j)
!!!!        enddo
!!!!    enddo
!!!---------------------------------------------------------------------- debugging
!!!---------------------------------------------------------------------- debugging

    end subroutine EC_CalcElemNeighborsInList
!##############################################################################
    integer function EC_GetElemUnloadingCount (ElemId)
        integer, intent(in) :: ElemId
        EC_GetElemUnloadingCount=pEC_ElemData(ElemId)%iElemStatus
    end function EC_GetElemUnloadingCount
!##############################################################################
    integer function EC_GetElemSplit (ElemId)
        integer, intent(in) :: ElemId
        EC_GetElemSplit=pEC_ElemData(ElemId)%iElemSplit
    end function EC_GetElemSplit
!##############################################################################
    real*8 function EC_GetElemAreaRatio (ElemId)
        integer, intent(in) :: ElemId
        EC_GetElemAreaRatio=pEC_ElemData(ElemId)%rAreaRatio
    end function EC_GetElemAreaRatio
!##############################################################################

end module EC_Objects_manager
