      subroutine areacoeff(ele)
	  
	  parameter (nume=40000)	  
	  common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
	  common /overlapping/ intersec(4, nume), area_coeff(nume), update_flag
	  common /crackline/ ncleave(3,nume), elecrack(4,nume), 
     1       nodeflag(4,nume)
	  
	  real intersec, area_coeff
	  integer elecrack, ele, nd, md
	  
	  nd=elecrack(1,ele)
	  md=elecrack(3,ele)
	  if(abs(nd-md)==2) then
c     two quadrilaterial
          
		  if(nd==1) then
		      area_coeff(ele)=(intersec(1,ele)+1-intersec(3,ele))/2.0
		  else if(nd==3) then
		      area_coeff(ele)=(1-intersec(1,ele)+intersec(3,ele))/2.0
		  else if(nd==2) then
		      area_coeff(ele)=(intersec(2,ele)+1-intersec(4,ele))/2.0
		  else if(nd==4) then
		      area_coeff(ele)=(1-intersec(2,ele)+intersec(4,ele))/2.0
		  end if
      
	  else
c     one pentagon and triangle
          
		  if(nd==1 .and. md==2) then
		      area_coeff(ele)=1-(1-intersec(1,ele))*intersec(4,ele)/2.0
		  else if(nd==2 .and. md==1) then
		      area_coeff(ele)=1-(1-intersec(3,ele))*intersec(2,ele)/2.0
		  else if(nd==2 .and. md==3) then
		      area_coeff(ele)=1-(1-intersec(2,ele))*intersec(3,ele)/2.0
		  else if(nd==3 .and. md==2) then
		      area_coeff(ele)=1-(1-intersec(4,ele))*intersec(1,ele)/2.0
		  else if(nd==3 .and. md==4) then
		      area_coeff(ele)=1-(1-intersec(1,ele))*intersec(4,ele)/2.0
		  else if(nd==4 .and. md==3) then
		      area_coeff(ele)=1-(1-intersec(3,ele))*intersec(2,ele)/2.0
		  else if(nd==4 .and. md==1) then
		      area_coeff(ele)=(1-intersec(2,ele))*intersec(3,ele)/2.0
		  else if(nd==1 .and. md==4) then
		      area_coeff(ele)=(1-intersec(4,ele))*intersec(1,ele)/2.0
		  end if
		  
	  end if	  
	  
	  area_coeff(numelt+1)=1-area_coeff(ele)
	  
	  end
	  