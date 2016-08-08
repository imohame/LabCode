subroutine FindElemEdgeNeighbor(ElemConnectivity, ElemId1, ElemEdge1, bflag, ElemId2, ElemEdge2)
!c     find the element which shares the same edge  neighbor element on this edge
!c	  parameter (nume=40000)
	  common/meshnum/ numnpo, numelto	  
	  dimension ElemConnectivity(4,*)
	  integer numnpo, numelto,ElemConnectivity,j
      
	  integer ElemId1, ElemEdge1, bflag, ElemId2, ElemEdge2
	  integer EdgeNode1, EdgeNode2
	  
	  bflag=0

	  
        if(ElemEdge1==1) then
            EdgeNode1=1
            EdgeNode2=2
        else if(ElemEdge1==2) then
            EdgeNode1=2
            EdgeNode2=3
        else if(ElemEdge1==3) then
            EdgeNode1=3
            EdgeNode2=4
        else if(ElemEdge1==4) then
            EdgeNode1=4
            EdgeNode2=1
        end if
	  
        do j=1, numelto
            if(ElemConnectivity(EdgeNode1,ElemId1)==ElemConnectivity(2,j) .and. &
               ElemConnectivity(EdgeNode2,ElemId1)==ElemConnectivity(1,j)) then
                bflag=1
                ElemId2=j
                ElemEdge2=1
            else if(ElemConnectivity(EdgeNode1,ElemId1)==ElemConnectivity(3,j) .and. &
                    ElemConnectivity(EdgeNode2,ElemId1)==ElemConnectivity(2,j)) then
                bflag=1
                ElemId2=j
                ElemEdge2=2
            else if(ElemConnectivity(EdgeNode1,ElemId1)==ElemConnectivity(4,j) .and.  &
                    ElemConnectivity(EdgeNode2,ElemId1)==ElemConnectivity(3,j)) then
                bflag=1
                ElemId2=j
                ElemEdge2=3
            else if(ElemConnectivity(EdgeNode1,ElemId1)==ElemConnectivity(1,j) .and. &
                    ElemConnectivity(EdgeNode2,ElemId1)==ElemConnectivity(4,j)) then
                bflag=1
                ElemId2=j
                ElemEdge2=4
            end if
        end do
	  
    end
	  