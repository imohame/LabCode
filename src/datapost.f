      subroutine datapost(ix)
	  
	  parameter (nume=40000)
	  parameter (nume2=20000)
	  common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
	  common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
	  common /stressflag/ strflag1(nume),ctr_flag(nume)
	  common /crackline/ ncleave(3,nume), elecrack(4,nume), 
     1       nodeflag(4,nume)
	  common /pcracktip/ connect(4,nume2),node(2,nume2),penta(nume2),
     >	                 ndflag(2,nume2), numnpt, numeltu, ndc
	  
	  dimension ix(4,*)
	  
	  integer numelt, numnp, strflag1, elecrack
	  integer connect, penta, numnpt, numeltu
	  integer ele, ne1, ele2, ne2, ndflag
	  integer bflag, cbflag, est, ntip_edge
	  
	  numnpt=numnp
	  numeltu=numelt
	  
	  do i=1, nume2
	      do j=1,4
		      connect(j,i)=0
		  end do
		  penta(i)=0
		  ndflag(1,i)=0
		  ndflag(2,i)=0
	  end do
	  
	  do ele=1, numelt
	      do j=1, 4
		      connect(j,ele)=ix(j,ele)
		  end do
		  
c     cracked element		  
		  if(strflag1(ele)==2) then
		      call postele3x(ix, ele)
		  end if
		  
c     crack tip element
          if(strflag1(ele)==1) then
		      cbflag=0
		      do est=2,4,2
			      if(elecrack(est,ele)==3) then
				      ne1=elecrack(est-1,ele)
				      call edgecontact1(ix, ele, ne1, bflag, ele2, ne2)
					  if(bflag==1 .and. strflag1(ele2)==2) then
					      cbflag=cbflag+1
						  if(cbflag==1) ntip_edge=ne1
					  end if
				  end if
			  end do
			  
			  if(cbflag==1) then
			      call postele20(ix, ele, ntip_edge)
			  else if(cbflag==2) then
			      call postele22(ix, ele)
			  end if
			  
		  end if
		  
		  if(strflag1(ele)==0) then
              cbflag=0
		      do est=2,4,2
			      if(elecrack(est,ele)==2) then
				      ne1=elecrack(est-1,ele)
				      call edgecontact1(ix, ele, ne1, bflag, ele2, ne2)
					  if(bflag==1 .and. strflag1(ele2)==2) then
					      cbflag=cbflag+1
						  if(cbflag==1) ntip_edge=ne1
					  end if
				  end if
			  end do
			  
			  if(cbflag==1) then
			      call postele20(ix, ele,ntip_edge)
			  else if(cbflag==2) then
			      call postele22(ix, ele)
			  end if
			  
		  end if
		  
      end do
	  
	  end
			   