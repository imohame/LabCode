      subroutine postcoord(n1, n2, yc, zc)
!!!!	  
!!!!	  parameter (nume2=20000)
!!!!	  common/meshnum/ numnpo, numelto
!!!!	  common /pcracktip/ connect(4,nume2), node(2,nume2), &
!!!!     	     penta(nume2),ndflag(2,nume2), numnpt, numeltu, ndc
!!!!	  
!!!!	  integer n1, n2, ndflag, ndc, numnpt
!!!!	  real yc, zc, node
!!!!	  integer npd, numnpo
!!!!	  real ync, znc, ds
!!!!      real distTol
!!!!      distTol= 1.0E-11
!!!!	  
!!!!	  ync=(node(1,n2)-node(1,n1))*yc+node(1,n1)
!!!!	  znc=(node(2,n2)-node(2,n1))*zc+node(2,n1)
!!!!	  
!!!!	  do i=1,2
!!!!	      if(ndflag(i,n2)/=0) then
!!!!		      npd=ndflag(i,n2)
!!!!		      ds=sqrt((ync-node(1,npd))**2+(znc-node(2,npd))**2)
!!!!	          if(ds<distTol) then
!!!!	              ndc=npd
!!!!				  go to 20
!!!!			  end if
!!!!		  end if
!!!!	  end do
!!!!	  
!!!!	  if(n2<=numnpo) then
!!!!	      do i=1,2
!!!!	          if(ndflag(i,n1)/=0) then
!!!!		          npd=ndflag(i,n1)
!!!!		          ds=sqrt((ync-node(1,npd))**2+(znc-node(2,npd))**2)
!!!!	              if(ds<distTol) then
!!!!	                  ndc=npd
!!!!				      go to 20
!!!!			      end if
!!!!		      end if
!!!!	      end do
!!!!	  end if
!!!!	  
!!!!	  numnpt=numnpt+1
!!!!	  node(1,numnpt)=ync
!!!!	  node(2,numnpt)=znc
!!!!	  if(ndflag(1,n2)==0) then
!!!!	      ndflag(1,n2)=numnpt
!!!!	  else if(ndflag(2,n2)==0) then
!!!!	      ndflag(2,n2)=numnpt
!!!!	  else
!!!!	      write(*,*) 'Warning, check postcoord.f'
!!!!	  end if
!!!!	  ndc=numnpt
!!!!	  
   20 end