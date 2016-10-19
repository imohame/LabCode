      subroutine eleupdate11(ele, ix)
!!!c     update nodeflag for two crack front tip  
!!!
!!!	  parameter (nume=40000)
!!!	  common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
!!!	  common /crackline/ ncleave(3,nume), elecrack(4,nume), 
!!!     1          nodeflag(4,nume)
!!!	  
!!!	  dimension ix(4,*)
!!!	  integer elecrack, nodeflag, nd, md, ele, numelt
!!!	  
!!!	  ix(1,numelt+1)=ix(1,ele)
!!!	  ix(2,numelt+1)=ix(2,ele)
!!!	  ix(3,numelt+1)=ix(3,ele)
!!!	  ix(4,numelt+1)=ix(4,ele)
!!!	  
!!!	  nd=elecrack(1,ele)
!!!	  md=elecrack(3,ele)
!!!	  if(abs(nd-md)==2) then
!!!c     two quadrilaterial
!!!          
!!!		  if(nd==1 .or. nd==3) then
!!!		      nodeflag(2,ele)=1
!!!			  nodeflag(3,ele)=1
!!!			  nodeflag(1,numelt+1)=1
!!!			  nodeflag(4,numelt+1)=1
!!!		  else 
!!!		      nodeflag(3,ele)=1
!!!			  nodeflag(4,ele)=1
!!!			  nodeflag(1,numelt+1)=1
!!!			  nodeflag(2,numelt+1)=1
!!!		  end if
!!!      
!!!	  else
!!!c     one pentagon and triangle
!!!          
!!!		  if((nd==1 .and. md==2) .or. (nd==2 .and. md==1)) then
!!!		      nodeflag(2,ele)=1
!!!			  nodeflag(1,numelt+1)=1
!!!			  nodeflag(3,numelt+1)=1
!!!			  nodeflag(4,numelt+1)=1
!!!		  else if((nd==2 .and. md==3) .or. (nd==3 .and. md==2)) then
!!!		      nodeflag(3,ele)=1
!!!			  nodeflag(1,numelt+1)=1
!!!			  nodeflag(2,numelt+1)=1
!!!			  nodeflag(4,numelt+1)=1
!!!		  else if((nd==3 .and. md==4) .or. (nd==4 .and. md==3)) then
!!!		      nodeflag(4,ele)=1
!!!			  nodeflag(1,numelt+1)=1
!!!			  nodeflag(2,numelt+1)=1
!!!			  nodeflag(3,numelt+1)=1
!!!		  else if((nd==4 .and. md==1) .or. (nd==1 .and. md==4)) then
!!!		      nodeflag(2,ele)=1
!!!		      nodeflag(3,ele)=1
!!!		      nodeflag(4,ele)=1
!!!			  nodeflag(1,numelt+1)=1
!!!		  end if
!!!		  
!!!	  end if

      end
	  
			  
			  
	  

      
		      
			  
	  
	  