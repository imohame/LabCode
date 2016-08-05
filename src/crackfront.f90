      subroutine crackfront(y, ix)
!!c     find elements at crack front
      use CN_Consts
!!!      parameter (nume=40000)
	  common /crackline/ ncleave(3,nume), elecrack(4,nume), nodeflag(4,nume)
	  common/meshnum/ numnpo, numelto 
	  common/cracktip/ tipele(2,nume), ntp
	  common /tipvelocity/ ncrack,nelefail(1000),tipelenum(1000,nume)
	  
	  dimension y(*), ix(4,*)
	  
	  integer elecrack, numelto
	  integer ne1, bflag, ele2, ne2
	  integer tipele, ntp
	  integer i, j, ele
	  integer ncrack, nelefail, tipelenum 
	  integer nlc, llast, ele1, ele22, nrc, rlast
	  real yele1, yele2
	  
	  ntp=0
	  do i=1,nume
	      do j=1,2
		      tipele(j,i)=0
		  end do
	  end do
	  
    do ele=1,numelto	  
        if(elecrack(2,ele)==1) then
            ne1=elecrack(1,ele)
            call edgecontact1(ix, ele, ne1, bflag, ele2, ne2)
            if(bflag==1) then 
                if((elecrack(1,ele2)==ne2 .and. elecrack(2,ele2)==2).or.(elecrack(3,ele2)==ne2 .and. elecrack(4,ele2)==2))then
                    ntp=ntp+1
                    tipele(1,ntp)=ele
                    tipele(2,ntp)=ele2
                end if
            end if
        end if
			
        if(elecrack(4,ele)==1) then
            ne1=elecrack(3,ele)
            call edgecontact1(ix, ele, ne1, bflag, ele2, ne2)
            if(bflag==1) then 
                if((elecrack(1,ele2)==ne2 .and. elecrack(2,ele2)==2).or.(elecrack(3,ele2)==ne2 .and. elecrack(4,ele2)==2))then
                    ntp=ntp+1
                    tipele(1,ntp)=ele
                    tipele(2,ntp)=ele2
                end if
            end if
        end if			
    end do
	  
    do i=1, ncrack   ! save cracked elements information  ! to calculate crack velocity	      
        nlc=2*i-1    ! left crack tip
        llast=nelefail(nlc)
        if(llast==1) then
            do j=1,ntp
                if(tipelenum(nlc, llast)==tipele(1,j)) then
                    ele1=tipele(1,j)
                    ele22=tipele(2,j)
                    yele1=(y(ix(1,ele1))+y(ix(2,ele1))+y(ix(3,ele1))+y(ix(4,ele1)))/4.0
                    yele2=(y(ix(1,ele22))+y(ix(2,ele22))+y(ix(3,ele22))+y(ix(4,ele22)))/4.0
                    if(yele2<yele1) then	  ! left crack tip
                        nelefail(nlc)=nelefail(nlc)+1
                        tipelenum(nlc,llast+1)=tipele(2,j)
                    end if
                end if
            end do
        else
            do j=1,ntp
                if(tipelenum(nlc, llast)==tipele(1,j)) then
                    nelefail(nlc)=nelefail(nlc)+1
                    tipelenum(nlc,llast+1)=tipele(2,j)
                end if 
            end do
        end if
		  
        nrc=2*i    ! right crack tip
        rlast=nelefail(nrc)
        if(rlast==1) then
            do j=1,ntp
                if(tipelenum(nrc, rlast)==tipele(1,j)) then
                    ele1=tipele(1,j)
                    ele22=tipele(2,j)
                    yele1=(y(ix(1,ele1))+y(ix(2,ele1))+y(ix(3,ele1))+y(ix(4,ele1)))/4.0
                    yele2=(y(ix(1,ele22))+y(ix(2,ele22))+y(ix(3,ele22))+y(ix(4,ele22)))/4.0
                    if(yele2>yele1) then	  ! right crack tip
                        nelefail(nrc)=nelefail(nrc)+1
                        tipelenum(nrc,rlast+1)=tipele(2,j)
                    end if
                end if
            end do
        else
            do j=1,ntp
                if(tipelenum(nrc, rlast)==tipele(1,j)) then
                    nelefail(nrc)=nelefail(nrc)+1
                    tipelenum(nrc,rlast+1)=tipele(2,j)
                end if 
            end do
        end if
		  
    end do	  
	  

end