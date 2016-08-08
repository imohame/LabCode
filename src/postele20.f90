      subroutine postele20(ix, ele, ntip_edge)
	  
	  parameter (nume=40000)
	  parameter (nume2=20000)
	  common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
	  common /pcracktip/ connect(4,nume2), node(2,nume2), &
     	       penta(nume2),ndflag(2,nume2), numnpt, numeltu, ndc
	  
	  dimension ix(4,*)
      integer ix
	  
	  integer numelt
	  integer connect, numeltu, penta
	  integer ele, ntip_edge, ndc
	  
	  if(ntip_edge==1) then
	      call cracknode1(1,ele,ix)
		  
		  numeltu=numeltu+1
		  connect(1,numeltu)=ndc
		  connect(2,numeltu)=connect(2,ele)
		  connect(3,numeltu)=connect(3,ele)
		  connect(4,numeltu)=connect(3,ele)
		  penta(numeltu-numelt)=ele
		  
		  connect(2,ele)=ndc
		  
	  else if(ntip_edge==2) then
	      call cracknode1(2,ele,ix)
		  
		  numeltu=numeltu+1
		  connect(1,numeltu)=ndc
		  connect(2,numeltu)=ndc
		  connect(3,numeltu)=connect(3,ele)
		  connect(4,numeltu)=connect(4,ele)
		  penta(numeltu-numelt)=ele
		  
		  connect(3,ele)=ndc
		  
	  else if(ntip_edge==3) then
	      call cracknode1(3,ele,ix)
		  
		  numeltu=numeltu+1
		  connect(1,numeltu)=connect(1,ele)
		  connect(2,numeltu)=ndc
		  connect(3,numeltu)=ndc
		  connect(4,numeltu)=connect(4,ele)
		  penta(numeltu-numelt)=ele
		  
		  connect(4,ele)=ndc
		  
	  else if(ntip_edge==4) then
	      call cracknode1(4,ele,ix)
		  
		  numeltu=numeltu+1
		  connect(1,numeltu)=connect(1,ele)
		  connect(2,numeltu)=connect(2,ele)
		  connect(3,numeltu)=ndc
		  connect(4,numeltu)=ndc
		  penta(numeltu-numelt)=ele
		  
		  connect(1,ele)=ndc
		  
	  end if
	  
      end
	  