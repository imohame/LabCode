      subroutine nodrepla(m, n, ix)
	  
	  parameter (nume=40000)
	  common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
	  common /crackline/ ncleave(3,nume), elecrack(4,nume), 
     1       nodeflag(4,nume)
c	  common/meshnum/ numnpo, numelto
	  
	  dimension ix(4,*)
	  
	  integer elecrack, nodeflag, numelt, m, n, nd
	  
	  do i=1, numelt
	      do j=1, 4
		      if(ix(j,i)==m .and. nodeflag(j,i)==1) then
			      ix(j,i)=n     
				  nodeflag(j,i)=0  
			  end if
		  end do
	  end do
						
	  end  