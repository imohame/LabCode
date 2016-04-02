      subroutine edge_cross(ndnum, crflag)
	  
	  parameter (nume2=20000)
	  common /pcracktip/ connect(4,nume2),node(2,nume2),penta(nume2),
     >	                 ndflag(2,nume2), numnpt, numeltu, ndc
	  
	  integer ndnum(4), crflag
	  real ncord(2,4), node
	  real k1, b1, k2, b2, yc, zc, tol
	  
	  tol=1.0e-6
	  yc=0.0
	  zc=0.0
	  
	  do i=1,4
	      ncord(1,i)=node(1,ndnum(i))
	      ncord(2,i)=node(2,ndnum(i))
	  end do
	  
	  if(sqrt((ncord(1,1)-ncord(1,3))**2+
     >    (ncord(2,1)-ncord(2,3))**2)<1.0e-8) then
	      go to 20
	  else if(sqrt((ncord(1,1)-ncord(1,4))**2+
     >    (ncord(2,1)-ncord(2,4))**2)<1.0e-8) then
	      go to 20
      else if(sqrt((ncord(1,2)-ncord(1,3))**2+
     >    (ncord(2,2)-ncord(2,3))**2)<1.0e-8) then
	      go to 20  
	  else if(sqrt((ncord(1,2)-ncord(1,4))**2+
     >    (ncord(2,2)-ncord(2,4))**2)<1.0e-8) then
	      go to 20  
	  end if
	  
c     calculate intersection points
      if(ncord(1,1)==ncord(1,2) .and. ncord(1,3)/=ncord(1,4)) then
	      k2=(ncord(2,4)-ncord(2,3))/(ncord(1,4)-ncord(1,3))
		  b2=ncord(2,3)-k2*ncord(1,3)
		  yc=ncord(1,1)
		  zc=k2*yc+b2
		  
	  elseif(ncord(1,3)==ncord(1,4) .and. ncord(1,1)/=ncord(1,2)) then
	      k1=(ncord(2,2)-ncord(2,1))/(ncord(1,2)-ncord(1,1))
		  b1=ncord(2,1)-k1*ncord(1,1)
		  yc=ncord(1,3)
		  zc=k1*yc+b1
		  
	  elseif(ncord(1,3)==ncord(1,4) .and. ncord(1,1)==ncord(1,2)) then
	      go to 20
		  
	  else
	      k1=(ncord(2,2)-ncord(2,1))/(ncord(1,2)-ncord(1,1))
		  b1=ncord(2,1)-k1*ncord(1,1)
		  k2=(ncord(2,4)-ncord(2,3))/(ncord(1,4)-ncord(1,3))
		  b2=ncord(2,3)-k2*ncord(1,3)
		  if(k1 .ne. k2) then
	          yc=(b2-b1)/(k1-k2)
		      zc=k1*yc+b1
	      end if
		  
	  end if
	  
	  if (yc-min(ncord(1,1),ncord(1,2))>tol .and. 
     >    max(ncord(1,1),ncord(1,2))-yc>tol .and. 
     >    zc-min(ncord(2,1),ncord(2,2))>tol .and. 
     >    max(ncord(2,1),ncord(2,2))-zc>tol .and.
     >    yc-min(ncord(1,3),ncord(1,4))>tol .and. 
     >    max(ncord(1,3),ncord(1,4))-yc>tol .and. 
     >    zc-min(ncord(2,3),ncord(2,4))>tol .and. 
     >    max(ncord(2,3),ncord(2,4))-zc>tol) then
	      crflag=1
		  write(*,*) 'crack surface overlapped'
		  write(*,*) 'crack line intersection position:'
		  write(*,*) 'yc:', yc, 'zc:', zc
      end if
      
   20 end 	  
	  
	 
	  
	  
	  
	      
	  
	  