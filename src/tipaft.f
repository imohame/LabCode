      subroutine tipaft(ix, ele)
c     update crack tip variable for element behind crack tip	
c     element crack status from 31 to 33
  
	  parameter (nume=40000)
	  common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
	  common /crackline/ ncleave(3,nume), elecrack(4,nume), 
     >       nodeflag(4,nume)
	  
	  dimension ix(4,*)
	  
	  integer elecrack, numelt, bflag, ele, nced
	  integer nele, ne2
	  
	  if(elecrack(2,ele)==3) then
	      nced=elecrack(1,ele)
	      call FindElemEdgeNeighbor(ix, ele, nced, bflag, nele, ne2)
	      if (bflag==1) then
	        if(elecrack(1,nele)==ne2 .and. elecrack(2,nele)==1) then
		          elecrack(2,nele)=3 
		    elseif(elecrack(3,nele)==ne2.and.elecrack(4,nele)==1)then
		          elecrack(4,nele)=3  
		    end if
          end if
      end if
	  
      if(elecrack(4,ele)==3) then
	      nced=elecrack(3,ele)
	      call FindElemEdgeNeighbor(ix, ele, nced, bflag, nele, ne2)
	      if (bflag==1) then
	        if(elecrack(1,nele)==ne2 .and. elecrack(2,nele)==1) then
		          elecrack(2,nele)=3 
		    elseif(elecrack(3,nele)==ne2.and.elecrack(4,nele)==1)then
		          elecrack(4,nele)=3  
		    end if
          end if
      end if
	  
	  end
	  