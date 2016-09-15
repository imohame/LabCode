
! ###########################################################
module EC_Objects_manager
    use EC_ElemCrackingBaseClass
    use EC_Consts
!    implicit none
    public
    ! Declare a set geometric objects.
    type ( EC_ElemCrackingClass ) , ALLOCATABLE :: pEC_ElemData(:)
    integer, ALLOCATABLE :: EC_ElemNeighbors(:,:) !-(8,:) to hold the elem's edges neighboring elem/edges, update when overlapping
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
            Allocate(EC_ElemNeighbors(8,ne))
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
                pEC_ElemData(i)%iElemOverlapping=0
            enddo
            EC_ElemNeighbors=-1

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
          IF (ALLOCATED (EC_ElemNeighbors))              DEALLOCATE (EC_ElemNeighbors)
          IF (ALLOCATED (EC_NodesIntersectData))              DEALLOCATE (EC_NodesIntersectData)
        end subroutine EC_CleanMem
!##############################################################################
!##############################################################################
        subroutine EC_PrintTest()
          use EC_Consts
          use EC_ElemCrackingBaseClass
          implicit none
          integer i

          if(EC_bCracking == 0) then
            return
          endif

          CALL ECprintTest()

          do i=1,5
              call pEC_ElemData(i)%PrintTest(i)
              write(*,*)'==========================='
              write(*,*)pEC_ElemData(i)
              write(*,*)pEC_ElemData(i)%rCleavagePlane
              write(*,*)pEC_ElemData(i)%EdgeStatus
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
        ElemIdNeighbor=EC_ElemNeighbors(i,ElemId)
        EdgeIdNeighbor=EC_ElemNeighbors(i+4,ElemId)
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
                ElemIdNeighbor=EC_ElemNeighbors(i,ElemId)
                EdgeIdNeighbor=EC_ElemNeighbors(i+4,ElemId)
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
    subroutine EC_SplitElem2(ElemId,NodesCoordx, NodesCoordy, ElemConnect, DofIds, NodesDispl, ElemMaterial, &
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
        integer EdgeStatus(4),nodesIds(8)
        real*8 ElemAreaTotal

        EdgeStatus=-1
        do i=1,4
            if(pEC_ElemData(ElemId)%EdgeStatus(i)==EC_eCrackBothSides) then
                EdgeStatus(i)=1
            endif
        enddo
        if(((EdgeStatus(1)==1) .and. (EdgeStatus(3)==1)) .or. &
          ((EdgeStatus(2)==1) .and. (EdgeStatus(4)==1))) then
            !- added nodes map to old nodes
            !- 1 -- 2 -- 3 -- 4
            !- 5 -- 6 -- 7 -- 8
            nodesIds(1:4)=ElemConnect(1:4,ElemId)
            do i=1,4
                nodesIds(4+i)=nodesIds(i)+4+EC_NodeCountCurrent
            enddo
        endif

        !- case 1: edges 1 & 3 or edges 2 & 4
        if((EdgeStatus(1)==1) .and. (EdgeStatus(3)==1)) then
            !-- split into two quads
            !- e1: 1 - 6 - 7 - 4
            ElemConnect(1,ElemId)=nodesIds(1)
            ElemConnect(2,ElemId)=nodesIds(6)
            ElemConnect(3,ElemId)=nodesIds(7)
            ElemConnect(4,ElemId)=nodesIds(4)
            !- e2: 5 - 2 - 3 - 8
            ElemConnect(1,EC_ElemCountCurrent+1)=nodesIds(5)
            ElemConnect(2,EC_ElemCountCurrent+1)=nodesIds(2)
            ElemConnect(3,EC_ElemCountCurrent+1)=nodesIds(3)
            ElemConnect(4,EC_ElemCountCurrent+1)=nodesIds(8)
        endif
        if((EdgeStatus(2)==1) .and. (EdgeStatus(4)==1)) then
            !-- split into two quads
            !- e1: 1 - 2 - 7 - 8
            ElemConnect(1,ElemId)=nodesIds(1)
            ElemConnect(2,ElemId)=nodesIds(2)
            ElemConnect(3,ElemId)=nodesIds(7)
            ElemConnect(4,ElemId)=nodesIds(8)
            !- e2: 5 - 6 - 3 - 4
            ElemConnect(1,EC_ElemCountCurrent+1)=nodesIds(5)
            ElemConnect(2,EC_ElemCountCurrent+1)=nodesIds(6)
            ElemConnect(3,EC_ElemCountCurrent+1)=nodesIds(3)
            ElemConnect(4,EC_ElemCountCurrent+1)=nodesIds(4)
        endif

        !- copy nodes data to the new nodes
        do i=1,4
            call EC_CopyNodesData(nodesIds(i),nodesIds(i+4),NodesCoordx, NodesCoordy, ElemConnect, DofIds, NodesDispl, &
                                ElemMaterial, usi  , freep , ym,SolStepCount)
        enddo
        !- copy elems data to the new nodes
        call EC_CopyElemsData(ElemId,EC_ElemCountCurrent+1,NodesCoordx, NodesCoordy, ElemConnect, DofIds, NodesDispl, &
                                ElemMaterial, usi  , freep , ym,SolStepCount)
        CALL pEC_ElemData(ElemId)%CopyElem (pEC_ElemData(EC_ElemCountCurrent+1))

        !---------------------------1- add Dof
        do j=1,4
            DofIds(1,EC_NodeCountCurrent+j)=neq+2*j-1
            DofIds(2,EC_NodeCountCurrent+j)=neq+2*j
        enddo
        numnp=numnp+4
        neq=neq+8
        EC_ElemCountCurrent=EC_ElemCountCurrent+1
        EC_NodeCountCurrent=EC_NodeCountCurrent+4

    end subroutine EC_SplitElem2
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
        INTEGER mElemEdgeNeighborsList(4)
        integer EdgeNode1,EdgeNode2,EdgeNode1N,EdgeNode2N,j


        edgeNode1=ElemConnect(EC_ElemEdgesConnect(1,ElemEdgeId),ElemId)
        edgeNode2=ElemConnect(EC_ElemEdgesConnect(2,ElemEdgeId),ElemId)
        edgeNode1N= EC_NodeCountCurrent+1
        edgeNode2N= edgeNode1N+1
        !- copy nodes data to the new nodes
        call EC_CopyNodesData(edgeNode1,edgeNode1N,NodesCoordx, NodesCoordy, ElemConnect, DofIds, NodesDispl, &
                        ElemMaterial, usi  , freep , ym,SolStepCount)
        call EC_CopyNodesData(edgeNode2,edgeNode2N,NodesCoordx, NodesCoordy, ElemConnect, DofIds, NodesDispl, &
                        ElemMaterial, usi  , freep , ym,SolStepCount)
        !---------------------------1- add Dof
        do j=1,2
            DofIds(1,EC_NodeCountCurrent+j)=neq+2*j-1
            DofIds(2,EC_NodeCountCurrent+j)=neq+2*j
        enddo
        numnp=numnp+2
        neq=neq+4
        EC_NodeCountCurrent=EC_NodeCountCurrent+2

        call pEC_ElemData(ElemId)%EC_GetElemEdgeNeighbors (ElemEdgeId,mElemEdgeNeighborsList)

        !-- split into two quads -- copy the current connect to the new elem
        ElemConnect(1:4,EC_ElemCountCurrent+1)=ElemConnect(1:4,ElemId)
        !- e1: 1 - 2 - n+2 - 4
        ElemConnect(EC_ElemEdgesConnect(2,ElemEdgeId),ElemId)=edgeNode2N
        !- e2: 1 - n+1 - 3 - 4
        ElemConnect(EC_ElemEdgesConnect(1,ElemEdgeId),EC_ElemCountCurrent+1)=edgeNode1N
        !- copy elems data to the new nodes
        call EC_CopyElemsData(ElemId,EC_ElemCountCurrent+1,NodesCoordx, NodesCoordy, ElemConnect, DofIds, NodesDispl, &
                                ElemMaterial, usi  , freep , ym,SolStepCount)
        EC_ElemCountCurrent=EC_ElemCountCurrent+1
        
        do j=1,2
            if(mElemEdgeNeighborsList(j) == 0) then

            elseif(mElemEdgeNeighborsList(j) /= 0) then
                !-- split into two quads -- copy the current connect to the new elem
                ElemConnect(1:4,EC_ElemCountCurrent+1)=ElemConnect(1:4,ElemId)
                !- e1: 1 - 2 - n+2 - 4
                ElemConnect(EC_ElemEdgesConnect(2,ElemEdgeId),ElemId)=edgeNode2N
                !- e2: 1 - n+1 - 3 - 4
                ElemConnect(EC_ElemEdgesConnect(1,ElemEdgeId),EC_ElemCountCurrent+1)=edgeNode1N
                !- copy elems data to the new nodes
                call EC_CopyElemsData(ElemId,EC_ElemCountCurrent+1,NodesCoordx, NodesCoordy, ElemConnect, DofIds, NodesDispl, &
                                        ElemMaterial, usi  , freep , ym,SolStepCount)
                EC_ElemCountCurrent=EC_ElemCountCurrent+1

            endif
        enddo

    end subroutine EC_SplitEdge
!##############################################################################
!##############################################################################
    subroutine EC_CopyNodesData(n1,n2,NodesCoordx, NodesCoordy, ElemConnect, DofIds, &
                            NodesDispl, ElemMaterial, usi  , freep , ym,SolStepCount)
                             !-a(k03)    ,a(k04)      ,a(k02)      ,a(k57  ,a(k18)     , a(k08)      ,a(k20), a(k07), a(k09))
    use mod_parameters
    use CN_Objects_manager

    implicit none
    real NodesCoordx(*), NodesCoordy(*),NodesDispl(*)
    real usi(*),freep(5,*), ym(4,*)
    INTEGER ElemConnect(4,*),ElemMaterial(*),DofIds(2,*)
    INTEGER SolStepCount
    common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
    integer nwebuf,ntime,numnp,neq,ibar,mthsol
    real dn1,dn2
    integer n1,n2

    NodesCoordx(n2)=NodesCoordx(n1)
    NodesCoordy(n2)=NodesCoordy(n1)
    NodesDispl(DofIds(1,n2))=NodesDispl(DofIds(1,n1))
    NodesDispl(DofIds(2,n2))=NodesDispl(DofIds(2,n1))
    usi(DofIds(1,n2))=usi(DofIds(1,n1))
    usi(DofIds(2,n2))=usi(DofIds(2,n1))
    call CNmanager_CopyNode(n1,n2)

    end subroutine EC_CopyNodesData
!##############################################################################
!##############################################################################
    subroutine EC_CopyElemsData(e1,e2,NodesCoordx, NodesCoordy, ElemConnect, DofIds, NodesDispl, &
                                ElemMaterial, usi  , freep , ym,SolStepCount)
                             !-a(k03)    ,a(k04)      ,a(k02)      ,a(k57  ,a(k18)     , a(k08)      ,a(k20), a(k07), a(k09))
    use mod_parameters
    use CN_Objects_manager

    implicit none
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

    real NodesCoordx(*), NodesCoordy(*),NodesDispl(*)
    real usi(*),freep(5,*), ym(4,*)
    INTEGER ElemConnect(4,*),ElemMaterial(*),DofIds(2,*)
    INTEGER SolStepCount
    integer e1,e2

    freep(1:5,e2)=freep(1:5,e1)

    Y_modulus(e2)=Y_modulus(e1)
    possion_ratio(e2) = possion_ratio(e1)
    tau_y(e2) = tau_y(e1)
    abc(1:573,e2,1)=abc(1:573,e1,1)
    his(1:573,e2,1)=his(1:573,e1,1)

    nnn2(e2,1)=nnn2(e1,1)

    fhg(e1,1:8)=fhg(e2,1:8)
    fhghis(e1,1:8)=fhghis(e2,1:8)

    hgsstore(e1)=hgsstore(e2)
    hgshis(e1)=hgshis(e2)

    totenerstore(e1)=totenerstore(e2)
    totenerhis(e1)=totenerhis(e2)
    inertener(e1)=inertener(e2)

    hgenerstore(e1)=hgenerstore(e2)
    hgenerhis(e1)=hgenerhis(e2)

    hgstress1store(e1)=hgstress1store(e2)
    hgstress2store(e1)=hgstress2store(e2)
    hgstress1his(e1)=hgstress1his(e2)
    hgstress2his(e1)=hgstress2his(e2)

    sigfrac0(e1)=sigfrac0(e2)
    sigfrac(e1)=sigfrac(e2)
    decfrac(e1)=decfrac(e2)

    CALL CNmanager_CopyElemnt(e1,e2)
    CALL pEC_ElemData(e1)%CopyElem (pEC_ElemData(e2))

    end subroutine EC_CopyElemsData
!##############################################################################
!##############################################################################

end module EC_Objects_manager
