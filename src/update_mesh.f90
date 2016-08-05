      subroutine update_mesh(ele, y, z, ix, id, u, usi)
	  
	  parameter (nume=40000)
	  common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
      common /crackline/ ncleave(3,nume), elecrack(4,nume), nodeflag(4,nume)	  
	  dimension y(*), z(*), ix(4,*), id(2,*), u(*), usi(*)	  
	  integer ele, numnp, neq, elecrack

!!     two crack tip for the element, update area_coeff, but do not add phantom nodes.		  
    if(elecrack(2,ele)==1 .and. elecrack(4,ele)==1) then
        call eleupdate11(ele, ix)		  
    else if((elecrack(2,ele)==3 .and. elecrack(4,ele)==1) .or.(elecrack(2,ele)==1 .and. elecrack(4,ele)==3)) then
!c     one line cracked, add two phantom nodes
        do j=1,2
            id(1,numnp+j)=neq+2*j-1
            id(2,numnp+j)=neq+2*j
        end do		  
        call eleupdate31(ele, y, z, ix, id, u, usi)
        numnp=numnp+2
        neq=neq+4		  
    else if(elecrack(2,ele)==3 .and. elecrack(4,ele)==3) then
!c     two line cracked, add four phantom nodes     
        do j=1,4
            id(1,numnp+j)=neq+2*j-1
            id(2,numnp+j)=neq+2*j
        end do
        call eleupdate33(ele, y, z, ix, id, u, usi)
        numnp=numnp+4
        neq=neq+8		  
    end if
	  
    call areacoeff(ele)	  		  
		    
end		  
		      