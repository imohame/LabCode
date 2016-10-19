module EC_OutputData
    INTEGER, ALLOCATABLE:: bNodeFlag(:),iElemConnectivityNew(:,:)
    real*8, ALLOCATABLE:: rNodesCoords(:,:)
    INTEGER :: iNodeCountCurrentNew,iElemCountCurrentNew
    save
end module EC_OutputData

subroutine EC_WriteConnectivity(ElemConnect,NodesCoordx, NodesCoordy,NodesDispl,DofIds,ElemMaterial, &
                                udt,udtt,temp)
    use EC_Objects_manager   
    use EC_Consts
    use EC_OutputData

    implicit none
    INTEGER    :: DofIds(2,*),ElemConnect(4,*),ElemMaterial(*)
    real       :: NodesCoordx(*), NodesCoordy(*),NodesDispl(*)
    real       :: udt(*),udtt(*),temp(*)

    INTEGER :: i,j,k,node1,node2,mphantomNode,nodeOut
    INTEGER :: ElemId2,ElemEdgeId2,EdgeCount
    real*8 ::DistRatio,rPtc(2)


    iNodeCountCurrentNew=EC_NodeCountCurrent
    iElemCountCurrentNew=EC_ElemCountCurrent

    allocate(bNodeFlag(EC_NodeCountCurrent))
    allocate(rNodesCoords(2,int(EC_NodeCountCurrent*1.2)+100))
    allocate(iElemConnectivityNew(4,int(EC_ElemCountCurrent*1.2)+100))
    !- it's used to mark the node if it's already adjusted for the actual coordinates or not
    bNodeFlag=0 !- 0=not yet ,  1= adjusted / calcualted
!!!!!!!!    rNodesCoords(1,1:EC_NodeCountCurrent)=NodesCoordx(1:EC_NodeCountCurrent)
!!!!!!!!    rNodesCoords(2,1:EC_NodeCountCurrent)=NodesCoordy(1:EC_NodeCountCurrent)

    iElemConnectivityNew(1:4,1:EC_ElemCountCurrent)=ElemConnect(1:4,1:EC_ElemCountCurrent)
    call hspwr (NodesDispl,udt,udtt,DofIds,temp,NodesCoordx, NodesCoordy)
!!!!     hspwr (a(k18),a(k55),a(k56),a(k57),a(k81+1),a(k03),a(k04))
      !--          u  , udt  ,udtt  ,dof   ,temp     ,x,y

!-- loop over all the elems to change the nodes coordinates if the nodes are phantom
    do i=1, EC_ElemCountCurrent
!!        write(*,*)i,iElemConnectivityNew(:,i)
        do j=1,4
            !----- if it has a corresponding / related elem,  then it has phantom nodes
            DistRatio=pEC_ElemData(i)%rCoordRatioCracking(j)
            if(DistRatio>0.0) then
                node1=pEC_ElemData(i)%iElemConnectivity(EC_ElemEdgesConnect(1,j))
                node2=pEC_ElemData(i)%iElemConnectivity(EC_ElemEdgesConnect(2,j))
                CALL EC_GetPhantomNodeCoord(DistRatio,node1,node2,rPtc,NodesCoordx, NodesCoordy,NodesDispl,DofIds)
                mphantomNode=0
                if ((node1 > EC_NodeCountInput) .and. (bNodeFlag(node1)==0)) then !-- this is a phantom node
                    mphantomNode=node1
                endif
                if ((node2 > EC_NodeCountInput) .and. (bNodeFlag(node2)==0)) then !-- this is a phantom node
                    mphantomNode=node2
                endif
                !- case where there is one old node and one phantom node
                if (mphantomNode>0) then !-- this is a phantom node
                    !- set the new coordinates
                    rNodesCoords(1:2,mphantomNode)=rPtc(1:2)
                    !- mark this node as updated
                    bNodeFlag(mphantomNode)=1
                end if
                !-- elem has only one cracked edge,,!-- add a new node, reconnect the elem to it
                if ((node1 <= EC_NodeCountInput) .and. (node2 <= EC_NodeCountInput)) then
                    if ((bNodeFlag(node1)==1) .and. (bNodeFlag(node2)==1)) then
                        cycle
                    endif
                    iNodeCountCurrentNew=iNodeCountCurrentNew+1
                    rNodesCoords(1:2,iNodeCountCurrentNew)=rPtc(1:2)
                    bNodeFlag(iNodeCountCurrentNew)=1
                    bNodeFlag(node1)=1
                    bNodeFlag(node2)=1
                    !-- update the elems connectivity
                    !- find the edge that has a phantom node, then chnage the non-phantom node
                    CALL pEC_ElemData(i)%EC_FindNodeEdgeWithPhantomNode (nodeOut)
                    if ( nodeOut>0 ) then
                        iElemConnectivityNew(nodeOut,i)=iNodeCountCurrentNew
                    end if
!!                    write(*,*)i,iElemConnectivityNew(:,i)
                    !-- find the edge neighbors and update the elems connectivity
!!                    call pEC_ElemData(i)%EC_GetElemEdgeNeighbors (j,mElemEdgeNeighborsList)
                    EdgeCount=pEC_ElemData(i)%iElemEdgeNeighborsCount(j)
                    do k=1,EdgeCount    !--each edge can be shared between more than one elem
                        ElemId2    =pEC_ElemData(i)%iElemEdgeNeighbors((k-1)*2+1,j)
                        ElemEdgeId2=pEC_ElemData(i)%iElemEdgeNeighbors((k-1)*2+2,j)
                        CALL pEC_ElemData(ElemId2)%EC_FindNodeEdgeWithPhantomNode (nodeOut)
                        if ( nodeOut>0 ) then
!!!                                write(*,*)ElemId2,iElemConnectivityNew(:,ElemId2)
                            iElemConnectivityNew(nodeOut,ElemId2)=iNodeCountCurrentNew
!!!                                write(*,*)ElemId2,iElemConnectivityNew(:,ElemId2)
                        end if
                    enddo
                endif
            endif !- if(pEC_ElemData(i)%rCoordRatioCracking(j)>0.0) then
        enddo
    enddo !--do i=1, EC_ElemCountInput

    !---- write to file 29
    CALL EC_WriteFile29(ElemMaterial)
    !--- clean the memory
    DEALLOCATE(bNodeFlag,rNodesCoords,iElemConnectivityNew)
    iNodeCountCurrentNew=0
    iElemCountCurrentNew=0

END subroutine EC_WriteConnectivity
!##############################################################################
!##############################################################################
subroutine EC_GetPhantomNodeCoord(DistRatio,node1,node2,rPtc,NodesCoordx, NodesCoordy,NodesDispl,DofIds)
    implicit none
    integer,    intent(in) :: DofIds(2,*)
    real,     intent(in) :: NodesCoordx(*), NodesCoordy(*),NodesDispl(*)
    integer,    intent(in) :: node1,node2
    real*8,     intent(in) :: DistRatio
    real*8,     intent(out) :: rPtc(2)

    real*8 rP2(2),rP1(2)

    rP1(1)=NodesCoordx(node1)+NodesDispl(DofIds(1,node1))
    rP1(2)=NodesCoordy(node1)+NodesDispl(DofIds(2,node1))
    rP2(1)=NodesCoordx(node2)+NodesDispl(DofIds(1,node2))
    rP2(2)=NodesCoordy(node2)+NodesDispl(DofIds(2,node2))
    CALL GetPtOnLine(rP1,rP2,DistRatio,rPtc)

END subroutine EC_GetPhantomNodeCoord
!##############################################################################
!##############################################################################
subroutine EC_WriteFile29(ElemMaterial)
    use EC_Consts
    use EC_OutputData

    implicit none
    INTEGER,    intent(in) :: ElemMaterial(*)

    common/bk18/nummat,ityp2d,ako(31)
    integer nummat,ityp2d,ako

    common/bk08/kprint,nstep,ite,ilimit,newstf
    integer kprint,nstep,ite,ilimit,newstf

    integer ::numnptotal,numeltotal,i,ElemMatId

      write(29,40) nummat, EC_NodeCountCurrent, EC_ElemCountCurrent, iNodeCountCurrentNew,iElemCountCurrentNew
      if (nstep.eq.0) then
          numnptotal=0
          numeltotal=0
      endif

!----------------------------ismail20150715
      numnptotal=numnptotal+EC_NodeCountCurrent
      numeltotal=numeltotal+EC_ElemCountCurrent
      write(7017,*) numnptotal,numeltotal
!      -- dump the entire history matrix for post-processing
    !   write(7021)his(1:573,1:numelt,1)
!----------------------------------------------
    !-- write elem connectivity and material id
    do i=1,iElemCountCurrentNew
        ElemMatId=ElemMaterial(i)
        write(29,26) i, iElemConnectivityNew(1:4,i), ElemMatId
        write(7019,426) iElemConnectivityNew(1:4,i), ElemMatId
    end do
    !---- write nodes coordinates
    do i=1,iNodeCountCurrentNew
      write(29,100) i, rNodesCoords(1,i), rNodesCoords(2,i)
      write(7020,104) rNodesCoords(1,i), rNodesCoords(2,i)
    end do



   26 format(6i5)
  426 format(5i5)
   36 format(i5)
  100 format(i5,5x,e15.6,5x,e15.6)
  104 format(e15.6,5x,e15.6)
   40 format(5i5)
END subroutine EC_WriteFile29
