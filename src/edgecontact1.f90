      subroutine edgecontact1(ix, ele, ne1, bflag, nele, ne2)
c     find the element which shares the same edge 
  
c	  parameter (nume=40000)
	  common/meshnum/ numnpo, numelto
	  
	  dimension ix(4,*)
	  
	  integer ele, ne1, bflag, nele, ne2
	  integer nn1, nn2
	  
	  bflag=0

	  
        if(ne1==1) then
            nn1=1
            nn2=2
        else if(ne1==2) then
            nn1=2
            nn2=3
        else if(ne1==3) then
            nn1=3
            nn2=4
        else if(ne1==4) then
            nn1=4
            nn2=1
        end if
	  
        do j=1, numelto
            if(ix(nn1,ele)==ix(2,j) .and. ix(nn2,ele)==ix(1,j)) then
                bflag=1
                nele=j
                ne2=1
            else if(ix(nn1,ele)==ix(3,j) .and. ix(nn2,ele)==ix(2,j)) then
                bflag=1
                nele=j
                ne2=2
            else if(ix(nn1,ele)==ix(4,j) .and.  ix(nn2,ele)==ix(3,j)) then
                bflag=1
                nele=j
                ne2=3
            else if(ix(nn1,ele)==ix(1,j) .and. ix(nn2,ele)==ix(4,j)) then
                bflag=1
                nele=j
                ne2=4
            end if
        end do
	  
    end
	  