      subroutine eleupdate33(ele, y, z, ix, id, u, usi)
c     update connectivity 

	  parameter (nume=40000)
	  common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
	  common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
	  common /crackline/ ncleave(3,nume), elecrack(4,nume), 
     1       nodeflag(4,nume)
	  common /overlapping/ intersec(4, nume), area_coeff(nume), update_flag
	  
	  dimension y(*), z(*), ix(4,*), id(2,*), u(*), usi(*)
	  
	  integer elecrack, nodeflag, nd, md, ele
	  real intersec
	 
	  nd=elecrack(1,ele)
	  md=elecrack(3,ele)
      
	  do i=1, 4
	      call cpynode(ix(i,ele), numnp+i, y, z, id, u, usi)
          call nodrepla(ix(i,ele), numnp+i, ix)
	  end do
	  
	  if(abs(nd-md)==2) then
c     two quadrilaterial
          
		  if(nd==1 .or. nd==3) then
			  
		      ix(1,numelt+1)=numnp+1
		      ix(2,numelt+1)=ix(2,ele)
		      ix(3,numelt+1)=ix(3,ele)
		      ix(4,numelt+1)=numnp+4
		       
		      ix(2,ele)=numnp+2
			  ix(3,ele)=numnp+3
				  
		  else if(nd==2 .or. nd==4) then
		      
			  ix(1,numelt+1)=numnp+1
			  ix(2,numelt+1)=numnp+2
			  ix(3,numelt+1)=ix(3,ele)
			  ix(4,numelt+1)=ix(4,ele)
				  
			  ix(3,ele)=numnp+3
			  ix(4,ele)=numnp+4

		  end if
      
	  else
c     one pentagon and triangle
          if((nd==1 .and. md==2) .or. (nd==2 .and. md==1)) then
		      ix(1,numelt+1)=numnp+1
			  ix(2,numelt+1)=ix(2,ele)
			  ix(3,numelt+1)=numnp+3
			  ix(4,numelt+1)=numnp+4
			  
			  ix(2,ele)=numnp+2
			  
		  else if((nd==2 .and. md==3) .or. (nd==3 .and. md==2)) then
		      ix(1,numelt+1)=numnp+1
			  ix(2,numelt+1)=numnp+2
			  ix(3,numelt+1)=ix(3,ele)
			  ix(4,numelt+1)=numnp+4
			  
			  ix(3,ele)=numnp+3
			  
		  else if((nd==3 .and. md==4) .or. (nd==4 .and. md==3)) then
		      ix(1,numelt+1)=numnp+1
			  ix(2,numelt+1)=numnp+2
			  ix(3,numelt+1)=numnp+3
			  ix(4,numelt+1)=ix(ele,4)
			  
			  ix(4,ele)=numnp+4
			  
		  else if((nd==4 .and. md==1) .or. (nd==1 .and. md==4)) then
		      ix(1,numelt+1)=numnp+1
			  ix(2,numelt+1)=ix(2,ele)
			  ix(3,numelt+1)=ix(3,ele)
			  ix(4,numelt+1)=ix(4,ele)
			  
			  ix(2,ele)=numnp+2
			  ix(3,ele)=numnp+3
			  ix(4,ele)=numnp+4
			  
		  end if
		  
	  end if
	  
	  end
	  
			  
			  
	  

      
		      
			  
	  
	  