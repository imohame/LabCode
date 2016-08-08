      subroutine postele22(ix, ele)
	  
	  parameter (nume=40000)
	  parameter (nume2=20000)
	  common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
	  common /crackline/ ncleave(3,nume), elecrack(4,nume),nodeflag(4,nume)
	  common /pcracktip/ connect(4,nume2), node(2,nume2), &
     	     penta(nume2),ndflag(2,nume2), numnpt, numeltu, ndc
	  
	  dimension ix(4,*)
	  
	  integer numelt, elecrack
	  integer connect, numeltu, penta
	  integer ele, ndc, ndc1, ndc2
	  
!!c     divide as two quadrilaterals		  
	  if(abs(elecrack(1,ele)-elecrack(3,ele))==2) then

!!c         go through edge 1 and 3	  
	      if(elecrack(1,ele)==1 .or. elecrack(1,ele)==3) then
              call cracknode1(1, ele, ix)		  
			  ndc1=ndc
			  
			  call cracknode1(3, ele, ix)
			  ndc2=ndc
			  
			  numeltu=numeltu+1
		      connect(1,numeltu)=ndc1
		      connect(2,numeltu)=connect(2,ele)
		      connect(3,numeltu)=connect(3,ele)
		      connect(4,numeltu)=ndc2
		      penta(numeltu-numelt)=ele
			  
		      connect(2,ele)=ndc1
		      connect(3,ele)=ndc2
			  
!!!c        go throuth edge 2 and 4				  
	      else if(elecrack(1,ele)==2 .or. elecrack(1,ele)==4) then
		      call cracknode1(2, ele, ix)
			  ndc1=ndc
			  
			  call cracknode1(4, ele, ix)
			  ndc2=ndc
			  
			  numeltu=numeltu+1
		      connect(1,numeltu)=ndc2
		      connect(2,numeltu)=ndc1
		      connect(3,numeltu)=connect(3,ele)
		      connect(4,numeltu)=connect(4,ele)
		      penta(numeltu-numelt)=ele
			  
		      connect(3,ele)=ndc1
		      connect(4,ele)=ndc2
			  
		  end if
		  
	  end if
	  
!!!c     divide as one triangle and one pentagon
      if((elecrack(1,ele)==1 .and. elecrack(3,ele)==2) .or. &
         (elecrack(1,ele)==2 .and. elecrack(3,ele)==1)) then
	 
          call cracknode1(1, ele, ix)
		  ndc1=ndc
		
		  call cracknode1(2, ele, ix)
		  ndc2=ndc
		  
		  numeltu=numeltu+1
		  connect(1,numeltu)=ndc1
		  connect(2,numeltu)=connect(2,ele)
		  connect(3,numeltu)=ndc2
		  connect(4,numeltu)=ndc2
		  penta(numeltu-numelt)=ele
		  
		  numeltu=numeltu+1
		  connect(1,numeltu)=ndc1
		  connect(2,numeltu)=ndc2
		  connect(3,numeltu)=connect(3,ele)
		  connect(4,numeltu)=connect(3,ele)
		  penta(numeltu-numelt)=ele
			  
		  connect(2,ele)=ndc1
			  
	  else if((elecrack(1,ele)==2 .and. elecrack(3,ele)==3) .or. &
           	 (elecrack(1,ele)==3 .and. elecrack(3,ele)==2)) then
	 
	      call cracknode1(2, ele, ix)
		  ndc1=ndc
		
		  call cracknode1(3, ele, ix)
		  ndc2=ndc
		  
		  numeltu=numeltu+1
		  connect(1,numeltu)=ndc1
		  connect(2,numeltu)=ndc1
		  connect(3,numeltu)=connect(3,ele)
		  connect(4,numeltu)=ndc2
		  penta(numeltu-numelt)=ele
		  
		  numeltu=numeltu+1
		  connect(1,numeltu)=ndc1
		  connect(2,numeltu)=ndc1
		  connect(3,numeltu)=ndc2
		  connect(4,numeltu)=connect(4,ele)
		  penta(numeltu-numelt)=ele
			  
		  connect(3,ele)=ndc1
		  
	  else if((elecrack(1,ele)==3 .and. elecrack(3,ele)==4) .or. &
         	 (elecrack(1,ele)==4 .and. elecrack(3,ele)==3)) then
	 
	      call cracknode1(3, ele, ix)
		  ndc1=ndc
		
		  call cracknode1(4, ele, ix)
		  ndc2=ndc
		  
		  numeltu=numeltu+1
		  connect(1,numeltu)=ndc2
		  connect(2,numeltu)=ndc2
		  connect(3,numeltu)=ndc1
		  connect(4,numeltu)=connect(4,ele)
		  penta(numeltu-numelt)=ele
		  
		  numeltu=numeltu+1
		  connect(1,numeltu)=connect(1,ele)
		  connect(2,numeltu)=ndc1
		  connect(3,numeltu)=ndc1
		  connect(4,numeltu)=ndc2
		  penta(numeltu-numelt)=ele
			  
		  connect(4,ele)=ndc1
		  
	  else if((elecrack(1,ele)==4 .and. elecrack(3,ele)==1) .or. &
       	 (elecrack(1,ele)==1 .and. elecrack(3,ele)==4)) then 
	 
	      call cracknode1(4, ele, ix)
		  ndc1=ndc
		
		  call cracknode1(1, ele, ix)
		  ndc2=ndc
		  
		  numeltu=numeltu+1
		  connect(1,numeltu)=ndc2
		  connect(2,numeltu)=connect(2,ele)
		  connect(3,numeltu)=connect(2,ele)
		  connect(4,numeltu)=ndc1
		  penta(numeltu-numelt)=ele
		  
		  numeltu=numeltu+1
		  connect(1,numeltu)=ndc1
		  connect(2,numeltu)=connect(2,ele)
		  connect(3,numeltu)=connect(3,ele)
		  connect(4,numeltu)=connect(4,ele)
		  penta(numeltu-numelt)=ele
			  
		  connect(2,ele)=ndc2
		  connect(3,ele)=ndc2
		  connect(4,ele)=ndc1
	  
	  end if
	  
	  end
	  