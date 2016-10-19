      subroutine postele3x(ix, ele)
	  
!!!	  parameter (nume=40000)
!!!	  parameter (nume2=20000)
!!!	  common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
!!!	  common /crackline/ ncleave(3,nume), elecrack(4,nume), nodeflag(4,nume)
!!!	  common /pcracktip/ connect(4,nume2), node(2,nume2), &
!!!     	      penta(nume2),ndflag(2,nume2), numnpt, numeltu, ndc
!!!	  common/meshnum/ numnpo, numelto
!!!	  
!!!	  dimension ix(4,*)
!!!      integer ix
!!!	  
!!!	  integer numelt, numelto, elecrack
!!!	  integer connect, numeltu, penta
!!!	  integer ele, ndc
!!!	  
!!!!!c     divide as two quadrilaterals		  
!!!	  if(abs(elecrack(1,ele)-elecrack(3,ele))==2) then
!!!
!!!!!c         go through edge 1 and 3	  
!!!	      if(elecrack(1,ele)==1 .or. elecrack(1,ele)==3) then
!!!!!c             left part of element		  
!!!	          if(ele<=numelto) then
!!!			      call cracknode1(1,ele,ix)
!!!				  connect(2,ele)=ndc
!!!				  
!!!				  call cracknode2(3,ele,ix)
!!!				  connect(3,ele)=ndc
!!!				  
!!!!!c            right part of element					  
!!!			  else if(ele>numelto) then
!!!			      call cracknode2(1,ele,ix)
!!!				  connect(1,ele)=ndc
!!!				  
!!!				  call cracknode1(3,ele,ix)
!!!				  connect(4,ele)=ndc
!!!			  end if
!!!			  
!!!!!c        go throuth edge 2 and 4				  
!!!	      else if(elecrack(1,ele)==2 .or. elecrack(1,ele)==4) then
!!!!!c             lower part of element
!!!		      if(ele<=numelto) then
!!!                  call cracknode1(2,ele,ix)
!!!				  connect(3,ele)=ndc
!!!				  
!!!				  call cracknode2(4,ele,ix)
!!!				  connect(4,ele)=ndc
!!!				  
!!!!!c             upper part of element				  
!!!		      else if(ele>numelto) then
!!!			      call cracknode2(2,ele,ix)
!!!				  connect(2,ele)=ndc
!!!				  
!!!				  call cracknode1(4,ele,ix)
!!!                  connect(1,ele)=ndc				  
!!!			  end if
!!!			  
!!!		  end if
!!!		  
!!!	  end if
!!!	  
!!!!!c     divide as one triangle and one pentagon
!!!      if((elecrack(1,ele)==1 .and. elecrack(3,ele)==2) .or.  &
!!!         (elecrack(1,ele)==2 .and. elecrack(3,ele)==1)) then
!!!	      if(ele<=numelto) then
!!!		      call cracknode1(1,ele,ix)
!!!			  connect(2,ele)=ndc
!!!           
!!!              call cracknode2(2,ele,ix)
!!!			  
!!!			  numeltu=numeltu+1
!!!			  connect(1,numeltu)=connect(2,ele)
!!!			  connect(2,numeltu)=ndc
!!!			  connect(3,numeltu)=connect(3,ele)
!!!			  connect(4,numeltu)=connect(3,ele)
!!!			  penta(numeltu-numelt)=ele
!!!			  
!!!		  else if(ele>numelto) then
!!!		      call cracknode2(1,ele,ix)
!!!			  connect(1,ele)=ndc
!!!			  
!!!			  call cracknode1(2,ele,ix)
!!!			  
!!!			  connect(3,ele)=ndc
!!!			  connect(4,ele)=connect(3,ele)
!!!		  end if
!!!			  
!!!	  else if((elecrack(1,ele)==2 .and. elecrack(3,ele)==3) .or. &
!!!        	 (elecrack(1,ele)==3 .and. elecrack(3,ele)==2)) then
!!!	     if(ele<=numelto) then
!!!		      call cracknode1(2,ele,ix)
!!!			  connect(3,ele)=ndc
!!!           
!!!              call cracknode2(3,ele,ix)
!!!			  
!!!			  numeltu=numeltu+1
!!!			  connect(1,numeltu)=connect(3,ele)
!!!			  connect(2,numeltu)=ndc
!!!			  connect(3,numeltu)=connect(4,ele)
!!!			  connect(4,numeltu)=connect(4,ele)
!!!			  penta(numeltu-numelt)=ele
!!!			  
!!!		  else if(ele>numelto) then
!!!		      call cracknode2(2,ele,ix)
!!!			  connect(2,ele)=ndc
!!!			  connect(1,ele)=connect(2,ele)
!!!			  
!!!			  call cracknode1(3,ele,ix)
!!!			  connect(4,ele)=ndc
!!!		  end if
!!!	 
!!!	  else if((elecrack(1,ele)==3 .and. elecrack(3,ele)==4) .or. &
!!!       	 (elecrack(1,ele)==4 .and. elecrack(3,ele)==3)) then
!!!	      if(ele<=numelto) then
!!!		      call cracknode1(3,ele,ix)
!!!			  connect(4,ele)=ndc
!!!           
!!!              call cracknode2(4,ele,ix)
!!!			  
!!!			  numeltu=numeltu+1
!!!			  connect(1,numeltu)=connect(1,ele)
!!!			  connect(2,numeltu)=connect(4,ele)
!!!			  connect(3,numeltu)=connect(4,ele)
!!!			  connect(4,numeltu)=ndc
!!!			  penta(numeltu-numelt)=ele
!!!			  
!!!		  else if(ele>numelto) then
!!!		      call cracknode2(3,ele,ix)
!!!			  connect(3,ele)=ndc
!!!			  connect(2,ele)=connect(3,ele)
!!!			  
!!!			  call cracknode1(4,ele,ix)
!!!			  connect(1,ele)=ndc
!!!		  end if
!!!	 
!!!	  else if((elecrack(1,ele)==4 .and. elecrack(3,ele)==1) .or. &
!!!         	 (elecrack(1,ele)==1 .and. elecrack(3,ele)==4)) then
!!!	      if(ele<=numelto) then
!!!		      call cracknode1(1,ele,ix)
!!!			  connect(2,ele)=ndc
!!!			  connect(3,ele)=connect(2,ele)
!!!			  
!!!			  call cracknode2(4,ele,ix)
!!!			  connect(4,ele)=ndc
!!!			  
!!!		  else if(ele>numelto) then
!!!		      call cracknode2(1,ele,ix)
!!!			  connect(1,ele)=ndc  
!!!           
!!!              call cracknode1(4,ele,ix)
!!!			  
!!!			  numeltu=numeltu+1
!!!			  connect(1,numeltu)=connect(1,ele)
!!!			  connect(2,numeltu)=connect(4,ele)
!!!			  connect(3,numeltu)=connect(4,ele)
!!!			  connect(4,numeltu)=ndc
!!!			  penta(numeltu-numelt)=ele
!!!		      
!!!		  end if
!!!	  
!!!	  end if
	  
	  end
	  