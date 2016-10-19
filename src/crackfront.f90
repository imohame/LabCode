      subroutine crackfront(Xcoord, ix)
!!!!!!c     find elements at crack front
!!!!      use CN_Consts
!!!!!!!      parameter (nume=40000)
!!!!	  common /crackline/ ncleave(3,nume), elecrack(4,nume), nodeflag(4,nume)
!!!!	  common/meshnum/ numnpo, numelto 
!!!!	  common/cracktip/ tipele(2,nume), ntp
!!!!!!!!	  common /tipvelocity/ ncrack,nelefail(1000),tipelenum(1000,nume)
!!!!	  
!!!!	  dimension Xcoord(*), ix(4,*)
!!!!	  
!!!!	  integer elecrack, numelto
!!!!	  integer ElemEdge1, bflag, ElemId2, ElemEdge2
!!!!	  integer tipele, ntp
!!!!	  integer i, j, ele
!!!!!!!!	  integer ncrack, nelefail, tipelenum 
!!!!	  integer nlc, llast, ElemId1, ElemIdtip, nrc, rlast
!!!!	  real ElemId1X, ElemIdtipX
!!!!	  
!!!!	  ntp=0
!!!!      tipele(1:2,1:nume)=0
!!!!!!!!!	  do i=1,nume
!!!!!!!!!	      do j=1,2
!!!!!!!!!		      tipele(j,i)=0
!!!!!!!!!		  end do
!!!!!!!!!	  end do
!!!!	  !!-- loop over all the elements to update the cracked elements neighbors
!!!!    do ele=1,numelto	  
!!!!        if(elecrack(2,ele)==1) then !! if element edge is cracked
!!!!            ElemEdge1=elecrack(1,ele)     !! cracked edge id
!!!!            call FindElemEdgeNeighbor(ix, ele, ElemEdge1, bflag, ElemId2, ElemEdge2) !! find the neighbor element on this cracked edge
!!!!            if(bflag==1) then 
!!!!                if((elecrack(1,ElemId2)==ElemEdge2 .and. elecrack(2,ElemId2)==2) .or. &
!!!!                   (elecrack(3,ElemId2)==ElemEdge2 .and. elecrack(4,ElemId2)==2))then
!!!!                    ntp=ntp+1
!!!!                    tipele(1,ntp)=ele
!!!!                    tipele(2,ntp)=ElemId2
!!!!                end if
!!!!            end if
!!!!        end if
!!!!			
!!!!        if(elecrack(4,ele)==1) then
!!!!            ElemEdge1=elecrack(3,ele)
!!!!            call FindElemEdgeNeighbor(ix, ele, ElemEdge1, bflag, ElemId2, ElemEdge2)
!!!!            if(bflag==1) then 
!!!!                if((elecrack(1,ElemId2)==ElemEdge2 .and. elecrack(2,ElemId2)==2) .or. &
!!!!                   (elecrack(3,ElemId2)==ElemEdge2 .and. elecrack(4,ElemId2)==2))then
!!!!                    ntp=ntp+1
!!!!                    tipele(1,ntp)=ele
!!!!                    tipele(2,ntp)=ElemId2
!!!!                end if
!!!!            end if
!!!!        end if			
!!!!    end do
!!!!!!!!!!=================================================	  
!!!!!!!!!!=============== this only for crack velocity ====	  
!!!!!!!!!!=================================================	  
!!!!!!!!!    do i=1, ncrack   ! save cracked elements information  ! to calculate crack velocity	      
!!!!!!!!!        nlc=2*i-1    ! left crack tip
!!!!!!!!!        llast=nelefail(nlc)
!!!!!!!!!        if(llast==1) then
!!!!!!!!!            do j=1,ntp
!!!!!!!!!                if(tipelenum(nlc, llast)==tipele(1,j)) then
!!!!!!!!!                    ElemId1=tipele(1,j)
!!!!!!!!!                    ElemIdtip=tipele(2,j)
!!!!!!!!!                    ElemId1X=(Xcoord(ix(1,ElemId1))+Xcoord(ix(2,ElemId1))+Xcoord(ix(3,ElemId1))+Xcoord(ix(4,ElemId1)))/4.0
!!!!!!!!!                    ElemIdtipX=(Xcoord(ix(1,ElemIdtip))+Xcoord(ix(2,ElemIdtip))+Xcoord(ix(3,ElemIdtip))+Xcoord(ix(4,ElemIdtip)))/4.0
!!!!!!!!!                    if(ElemIdtipX<ElemId1X) then	  ! left crack tip
!!!!!!!!!                        nelefail(nlc)=nelefail(nlc)+1
!!!!!!!!!                        tipelenum(nlc,llast+1)=tipele(2,j)
!!!!!!!!!                    end if
!!!!!!!!!                end if
!!!!!!!!!            end do
!!!!!!!!!        else
!!!!!!!!!            do j=1,ntp
!!!!!!!!!                if(tipelenum(nlc, llast)==tipele(1,j)) then
!!!!!!!!!                    nelefail(nlc)=nelefail(nlc)+1
!!!!!!!!!                    tipelenum(nlc,llast+1)=tipele(2,j)
!!!!!!!!!                end if 
!!!!!!!!!            end do
!!!!!!!!!        end if
!!!!!!!!!		  
!!!!!!!!!        nrc=2*i    ! right crack tip
!!!!!!!!!        rlast=nelefail(nrc)
!!!!!!!!!        if(rlast==1) then
!!!!!!!!!            do j=1,ntp
!!!!!!!!!                if(tipelenum(nrc, rlast)==tipele(1,j)) then
!!!!!!!!!                    ElemId1=tipele(1,j)
!!!!!!!!!                    ElemIdtip=tipele(2,j)
!!!!!!!!!                    ElemId1X=(Xcoord(ix(1,ElemId1))+Xcoord(ix(2,ElemId1))+Xcoord(ix(3,ElemId1))+Xcoord(ix(4,ElemId1)))/4.0
!!!!!!!!!                    ElemIdtipX=(Xcoord(ix(1,ElemIdtip))+Xcoord(ix(2,ElemIdtip))+Xcoord(ix(3,ElemIdtip))+Xcoord(ix(4,ElemIdtip)))/4.0
!!!!!!!!!                    if(ElemIdtipX>ElemId1X) then	  ! right crack tip
!!!!!!!!!                        nelefail(nrc)=nelefail(nrc)+1
!!!!!!!!!                        tipelenum(nrc,rlast+1)=tipele(2,j)
!!!!!!!!!                    end if
!!!!!!!!!                end if
!!!!!!!!!            end do
!!!!!!!!!        else
!!!!!!!!!            do j=1,ntp
!!!!!!!!!                if(tipelenum(nrc, rlast)==tipele(1,j)) then
!!!!!!!!!                    nelefail(nrc)=nelefail(nrc)+1
!!!!!!!!!                    tipelenum(nrc,rlast+1)=tipele(2,j)
!!!!!!!!!                end if 
!!!!!!!!!            end do
!!!!!!!!!        end if
!!!!!!!!!		  
!!!!!!!!!    end do	  
!!!!!!!!!	  

end