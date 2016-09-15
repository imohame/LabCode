      subroutine datapost(ix)
	  
	  parameter (nume=40000)
	  parameter (nume2=20000)
	  common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
	  common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
	  common /stressflag/ ElemFractCode(nume),ElemDecayCount(nume)
	  common /crackline/ ncleave(3,nume), elecrack(4,nume), 
     1       nodeflag(4,nume)
	  common /pcracktip/ connect(4,nume2),node(2,nume2),penta(nume2),
     >	                 ndflag(2,nume2), numnpt, numeltu, ndc
	  
	  dimension ix(4,*)
	  
	  integer numelt, numnp, ElemFractCode, elecrack
	  integer connect, penta, numnpt, numeltu
	  integer ele, ne1, ele2, ne2, ndflag
	  integer bflag, cbflag, est, ntip_edge
      integer jj,NonzeroNodeId
	  
	  numnpt=numnp
	  numeltu=numelt
	  
!!!	  do i=1, nume2
!!!	      do j=1,4
!!!		      connect(j,i)=0
!!!		  end do
!!!		  penta(i)=0
!!!		  ndflag(1,i)=0
!!!		  ndflag(2,i)=0
!!!	  end do
	  
	  do ele=1, numelt
!!	      do j=1, 4
!!		      connect(j,ele)=ix(j,ele)
!!		  end do
          connect(1:4,ele)=ix(1:4,ele)
		  
!c     cracked element		  
		  if(ElemFractCode(ele)==2) then
		      call postele3x(ix, ele)
		  end if
		  
!c     crack tip element
          if(ElemFractCode(ele)==1) then
		      cbflag=0
		      do est=2,4,2
			      if(elecrack(est,ele)==3) then
				      ne1=elecrack(est-1,ele)
				   call FindElemEdgeNeighbor(ix, ele, ne1, bflag, ele2, ne2)
					  if(bflag==1 .and. ElemFractCode(ele2)==2) then
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
		  
		  if(ElemFractCode(ele)==0) then
              cbflag=0
		      do est=2,4,2
			      if(elecrack(est,ele)==2) then
				    ne1=elecrack(est-1,ele)
				    call FindElemEdgeNeighbor(ix, ele, ne1, bflag, ele2, ne2)
					  if(bflag==1 .and. ElemFractCode(ele2)==2) then
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
		!---------------------- add this check to prevent connect of having zeros for tec-plot
          NonzeroNodeId=0
          do jj=1,4
            if (connect(jj,ele) >0) then
                NonzeroNodeId=connect(jj,ele)
                exit !--- break the loop
             endif
          enddo
        !---- set any zero node to that non-zero node
          do jj=1,4
            if (connect(jj,ele) >0) then
                NonzeroNodeId=connect(jj,ele)
            else
                connect(jj,ele)=NonzeroNodeId
            endif
          enddo
        !-----------------------
          
      end do
	  
	  end
			   