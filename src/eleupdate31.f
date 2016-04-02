      subroutine eleupdate31(ele, y, z, ix, id, u, usi)
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
	  
	  if(elecrack(2,ele)==3) then
	      nd=elecrack(1,ele)
	      md=elecrack(3,ele)
	  else if(elecrack(4,ele)==3) then
	      nd=elecrack(3,ele)
	      md=elecrack(1,ele)
	  end if
	  
	  if(abs(nd-md)==2) then
c     two quadrilaterial
          
		  if(nd==1) then
              call cpynode(ix(nd,ele), numnp+1, y, z, id, u, usi)
			  call nodrepla(ix(nd,ele), numnp+1, ix)
				  
		      call cpynode(ix(nd+1,ele), numnp+2, y, z, id, u, usi)
		      call nodrepla(ix(nd+1,ele), numnp+2, ix)
		       
		      ix(1,numelt+1)=numnp+1
		      ix(2,numelt+1)=ix(2,ele)
		      ix(3,numelt+1)=ix(3,ele)
		      ix(4,numelt+1)=ix(4,ele)
		      nodeflag(4,numelt+1)=1
		       
		      ix(2,ele)=numnp+2
		      nodeflag(3,ele)=1
				  
		  else if(nd==3) then
		      call cpynode(ix(nd,ele), numnp+1, y, z, id, u, usi)
			  call nodrepla(ix(nd,ele), numnp+1, ix)
			  call cpynode(ix(nd+1,ele), numnp+2, y, z, id, u, usi)
			  call nodrepla(ix(nd+1,ele), numnp+2, ix)
				  
			  ix(1,numelt+1)=ix(1,ele)
			  ix(2,numelt+1)=ix(2,ele)
			  ix(3,numelt+1)=ix(3,ele)
			  ix(4,numelt+1)=numnp+2
			  nodeflag(1,numelt+1)=1
				  
			  ix(3,ele)=numnp+1
			  nodeflag(2,ele)=1
			  
		  else if(nd==2) then
		      call cpynode(ix(nd,ele), numnp+1, y, z, id, u, usi)
			  call nodrepla(ix(nd,ele), numnp+1, ix)
			  call cpynode(ix(nd+1,ele), numnp+2, y, z, id, u, usi)
			  call nodrepla(ix(nd+1,ele), numnp+2, ix)
			  
			  ix(1,numelt+1)=ix(1,ele)
			  ix(2,numelt+1)=numnp+1
			  ix(3,numelt+1)=ix(3,ele)
			  ix(4,numelt+1)=ix(4,ele)
			  nodeflag(1,numelt+1)=1
				  
			  ix(3,ele)=numnp+2
			  nodeflag(4,ele)=1
			  
		  else if(nd==4) then
		      call cpynode(ix(nd,ele), numnp+1, y, z, id, u, usi)
			  call nodrepla(ix(nd,ele), numnp+1, ix)
			  call cpynode(ix(1,ele), numnp+2, y, z, id, u, usi)
			  call nodrepla(ix(1,ele), numnp+2, ix)
			  
			  ix(1,numelt+1)=numnp+2
			  ix(2,numelt+1)=ix(2,ele)
			  ix(3,numelt+1)=ix(3,ele)
			  ix(4,numelt+1)=ix(4,ele)
			  nodeflag(2,numelt+1)=1
				  
			  ix(4,ele)=numnp+1
			  nodeflag(3,ele)=1
			  
		  end if
		  
	  else
c     one pentagon and triangle
          if(nd==1 .and. md==2) then
		      call cpynode(ix(4,ele), numnp+1, y, z, id, u, usi)
			  call nodrepla(ix(4,ele), numnp+1, ix)
			  call cpynode(ix(1,ele), numnp+2, y, z, id, u, usi)
			  call nodrepla(ix(1,ele), numnp+2, ix)
			  
			  ix(1,numelt+1)=numnp+2
			  ix(2,numelt+1)=ix(2,ele)
			  ix(3,numelt+1)=ix(3,ele)
			  ix(4,numelt+1)=numnp+1
			  nodeflag(3,numelt+1)=1
				  
			  nodeflag(2,ele)=1

		  else if(nd==2 .and. md==1) then
		      call cpynode(ix(3,ele), numnp+1, y, z, id, u, usi)
			  call nodrepla(ix(3,ele), numnp+1, ix)
			  call cpynode(ix(4,ele), numnp+2, y, z, id, u, usi)
			  call nodrepla(ix(4,ele), numnp+2, ix)
			  
			  ix(1,numelt+1)=ix(1,ele)
			  ix(2,numelt+1)=ix(2,ele)
			  ix(3,numelt+1)=numnp+1
			  ix(4,numelt+1)=numnp+2
			  nodeflag(1,numelt+1)=1
				  
			  nodeflag(2,ele)=1
		      
		  else if(nd==2 .and. md==3) then
		      call cpynode(ix(1,ele), numnp+1, y, z, id, u, usi)
			  call nodrepla(ix(1,ele), numnp+1, ix)
			  call cpynode(ix(2,ele), numnp+2, y, z, id, u, usi)
			  call nodrepla(ix(2,ele), numnp+2, ix)
			  
			  ix(1,numelt+1)=numnp+1
			  ix(2,numelt+1)=numnp+2
			  ix(3,numelt+1)=ix(3,ele)
			  ix(4,numelt+1)=ix(4,ele)
			  nodeflag(4,numelt+1)=1
				  
			  nodeflag(3,ele)=1
			  
		  else if(nd==3 .and. md==2) then
		      call cpynode(ix(4,ele), numnp+1, y, z, id, u, usi)
			  call nodrepla(ix(4,ele), numnp+1, ix)
			  call cpynode(ix(1,ele), numnp+2, y, z, id, u, usi)
			  call nodrepla(ix(1,ele), numnp+2, ix)
			  
			  ix(1,numelt+1)=numnp+2
			  ix(2,numelt+1)=ix(2,ele)
			  ix(3,numelt+1)=ix(3,ele)
			  ix(4,numelt+1)=numnp+1
			  nodeflag(2,numelt+1)=1
				  
			  nodeflag(3,ele)=1
			  
		  else if(nd==3 .and. md==4) then
		      call cpynode(ix(2,ele), numnp+1, y, z, id, u, usi)
			  call nodrepla(ix(2,ele), numnp+1, ix)
			  call cpynode(ix(3,ele), numnp+2, y, z, id, u, usi)
			  call nodrepla(ix(3,ele), numnp+2, ix)
			  
			  ix(1,numelt+1)=ix(1,ele)
			  ix(2,numelt+1)=numnp+1
			  ix(3,numelt+1)=numnp+2
			  ix(4,numelt+1)=ix(4,ele)
			  nodeflag(1,numelt+1)=1
				  
			  nodeflag(4,ele)=1
			  
		  else if(nd==4 .and. md==3) then
		      call cpynode(ix(1,ele), numnp+1, y, z, id, u, usi)
			  call nodrepla(ix(1,ele), numnp+1, ix)
			  call cpynode(ix(2,ele), numnp+2, y, z, id, u, usi)
			  call nodrepla(ix(2,ele), numnp+2, ix)
			  
			  ix(1,numelt+1)=numnp+1
			  ix(2,numelt+1)=numnp+2
			  ix(3,numelt+1)=ix(3,ele)
			  ix(4,numelt+1)=ix(4,ele)
			  nodeflag(3,numelt+1)=1
				  
			  nodeflag(4,ele)=1
			  
		  else if(nd==4 .and. md==1) then
		      call cpynode(ix(3,ele), numnp+1, y, z, id, u, usi)
			  call nodrepla(ix(3,ele), numnp+1, ix)
			  call cpynode(ix(4,ele), numnp+2, y, z, id, u, usi)
			  call nodrepla(ix(4,ele), numnp+2, ix)
			  
			  ix(1,numelt+1)=ix(1,ele)
			  ix(2,numelt+1)=ix(2,ele)
			  ix(3,numelt+1)=ix(3,ele)
			  ix(4,numelt+1)=ix(4,ele)
			  nodeflag(1,numelt+1)=1
				
              ix(3,ele)=numnp+1 
              ix(4,ele)=numnp+2			  
			  nodeflag(2,ele)=1
			  
		  else if(nd==1 .and. md==4) then
		      call cpynode(ix(2,ele), numnp+1, y, z, id, u, usi)
			  call nodrepla(ix(2,ele), numnp+1, ix)
			  call cpynode(ix(3,ele), numnp+2, y, z, id, u, usi)
			  call nodrepla(ix(3,ele), numnp+2, ix)
			  
			  ix(1,numelt+1)=ix(1,ele)
			  ix(2,numelt+1)=ix(2,ele)
			  ix(3,numelt+1)=ix(3,ele)
			  ix(4,numelt+1)=ix(4,ele)
			  nodeflag(1,numelt+1)=1
				
              ix(2,ele)=numnp+1 
              ix(3,ele)=numnp+2			  
			  nodeflag(4,ele)=1
			  
		  end if
		  
	  end if
	  
	  end
	  
			  
			  
	  

      
		      
			  
	  
	  