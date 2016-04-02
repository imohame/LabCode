      subroutine tipbef(ix, ele, nced, p)
c     update crack tip variable for element ahead of crack tip	
c     element crack status from 00 to 20, or from 20 to 22
  
	  parameter (nume=40000)
	  common /crackline/ ncleave(3,nume), elecrack(4,nume), 
     >       nodeflag(4,nume)
	  common /overlapping/ intersec(4, nume), area_coeff(nume), update_flag
	  common/meshnum/ numnpo, numelto
	  
	  dimension ix(4,*)
	  
	  integer elecrack, numelto, bflag, ele, nced, p
	  real intersec
      integer nele, ne2
	  
	  call edgecontact1(ix, ele, nced, bflag, nele, ne2)
	  
	  if (bflag==0) then 
	      elecrack(p+2,ele)=3
	  else if(bflag==1) then
	      if(elecrack(1,nele)==0) then
	          elecrack(1,nele)=ne2
			  elecrack(2,nele)=2
			  intersec(1,nele)=1.0-intersec(p+1,ele)
			  intersec(2,nele)=1.0-intersec(p+2,ele)
		  else if(elecrack(3,nele)==0) then
			  elecrack(3,nele)=ne2
			  elecrack(4,nele)=2
			  intersec(3,nele)=1.0-intersec(p+1,ele)
			  intersec(4,nele)=1.0-intersec(p+2,ele)
		  end if
	  end if
	  
	  end
	  