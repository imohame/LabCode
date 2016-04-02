      subroutine outconnect(matp)
	  
	  parameter (nume=40000)
	  parameter (nume2=20000)
	  common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
	  common /pcracktip/ connect(4,nume2), node(2,nume2), 
     >	      penta(nume2),ndflag(2,nume2), numnpt, numeltu, ndc
	  
	  dimension matp(*)
	  
	  integer numelt, numnpt, numeltu, connect, penta, nmat
	  
	  real node
	  
	  do i=1,numeltu
	      if(i<=numelt) then
		      nmat=matp(i)
		  else if(i>numelt) then
		      nmat=matp(penta(i-numelt))
		  end if 
		  
	      write(29,26) i, connect(1:4,i), nmat
	      write(7019,426) connect(1:4,i), nmat
		  
	  end do
	  
	  do i=1, numeltu-numelt
	      write(29,36) penta(i)
	      write(7018,36) penta(i)
	  end do
	  
	  do i=1,numnpt
	      write(29,100) i, node(1,i), node(2,i)
	      write(7020,104) node(1,i), node(2,i)
      end do		  
	  
   26 format(6i5)
  426 format(5i5)
   36 format(i5)
  100 format(i5,5x,e15.6,5x,e15.6)
  104 format(e15.6,5x,e15.6)
	  end
	  