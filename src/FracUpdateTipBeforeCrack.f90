      subroutine FracUpdateTipBeforeCrack(ix, ele)
!!c     update crack tip variable for element behind crack tip	
!!c     element crack status from 31 to 33
  
	  parameter (nume=40000)
	  common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
	  common /crackline/ ncleave(3,nume), elecrack(4,nume), nodeflag(4,nume)
	  
	  dimension ix(4,*)
      integer nprnt,mprint,itmpop,jprint,idump,locstr,ix
      
      
	  integer elecrack, numelt, bflag, ele, ElemEdge1
	  integer ElemId2, ElemEdge2
	  
	  if(elecrack(2,ele)==3) then
	      ElemEdge1=elecrack(1,ele)
	      call FindElemEdgeNeighbor(ix, ele, ElemEdge1, bflag, ElemId2, ElemEdge2)
	      if (bflag==1) then
	        if(elecrack(1,ElemId2)==ElemEdge2 .and. elecrack(2,ElemId2)==1) then
		          elecrack(2,ElemId2)=3 
		    elseif(elecrack(3,ElemId2)==ElemEdge2.and.elecrack(4,ElemId2)==1)then
		          elecrack(4,ElemId2)=3  
		    end if
          end if
      end if
	  
      if(elecrack(4,ele)==3) then
	      ElemEdge1=elecrack(3,ele)
	      call FindElemEdgeNeighbor(ix, ele, ElemEdge1, bflag, ElemId2, ElemEdge2)
	      if (bflag==1) then
	        if(elecrack(1,ElemId2)==ElemEdge2 .and. elecrack(2,ElemId2)==1) then
		          elecrack(2,ElemId2)=3 
		    elseif(elecrack(3,ElemId2)==ElemEdge2.and.elecrack(4,ElemId2)==1)then
		          elecrack(4,ElemId2)=3  
		    end if
          end if
      end if
	  
	  end
	  