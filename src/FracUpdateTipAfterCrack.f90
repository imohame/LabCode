      subroutine FracUpdateTipAfterCrack(ix, ele, ElemEdge1, p)
!!!!!!!!c     update crack tip variable for element ahead of crack tip	
!!!!!!!!c     element crack status from 00 to 20, or from 20 to 22
!!!  
!!!	  parameter (nume=40000)
!!!	  common /crackline/ ncleave(3,nume), elecrack(4,nume),nodeflag(4,nume)
!!!	  common /overlapping/ intersec(4, nume), area_coeff(nume), update_flag
!!!	  common/meshnum/ numnpo, numelto
!!!	  
!!!	  dimension ix(4,*)
!!!      integer ix,numnpo, numelto,elecrack
!!!      real intersec
!!!	  
!!!	  integer bflag, ele, ElemEdge1, p
!!!	  
!!!      integer ElemId2, ElemEdge2
!!!	  
!!!	  call FindElemEdgeNeighbor(ix, ele,ElemEdge1, bflag, ElemId2, ElemEdge2)
!!!	  
!!!	  if (bflag==0) then !-- if not found then set it to opened/free
!!!	      elecrack(p+2,ele)=3
!!!      else if(bflag==1) then !-- if found
!!!	      if(elecrack(1,ElemId2)==0) then   !-- if this location empty
!!!	          elecrack(1,ElemId2)=ElemEdge2 !-- set the edge id
!!!			  elecrack(2,ElemId2)=2         !-- the code is cracked
!!!			  intersec(1,ElemId2)=1.0-intersec(p+1,ele)
!!!			  intersec(2,ElemId2)=1.0-intersec(p+2,ele)
!!!		  else if(elecrack(3,ElemId2)==0) then
!!!			  elecrack(3,ElemId2)=ElemEdge2
!!!			  elecrack(4,ElemId2)=2
!!!			  intersec(3,ElemId2)=1.0-intersec(p+1,ele)
!!!			  intersec(4,ElemId2)=1.0-intersec(p+2,ele)
!!!		  end if
!!!	  end if
!!!	  
	  end
	  