      subroutine compdis(y, z, ix, id, u, ele, nd, cflag)
!!!!!            compdis(y, z, ix, id, u, ele, ndcircle, cflag)
!!c     compare distance from current element to cracked elment 
!!c     with preset distance 

      parameter (nume=40000)
	  common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
	  common /stressflag/ ElemFractCode(nume),ElemDecayCount(nume)
	  common/meshnum/ numnpo, numelto 
	  common/cracktip/ tipele(2,nume), ntp
	  
      dimension y(*), z(*), ix(4,*), id(2,*), u(*)
	 
	  integer tipele, ntp, ele, nd, cflag
	  real y3, z3, y1, z1, ds0, ds
	  real ye, ze, yt, zt
	  integer ntel, ElemFractCode, numelt, numelto
	  
	  cflag=1
	  
	  y3=y(ix(3,1))
	  z3=z(ix(3,1))
	  y1=y(ix(1,1))
	  z1=z(ix(1,1))
	  ds0=sqrt((y3-y1)**2+(z3-z1)**2)*nd   ! nd is a coefficient
	  
	  ye=(y(ix(1,ele))+y(ix(3,ele)))/2.0
	  ze=(z(ix(1,ele))+z(ix(3,ele)))/2.0
	  
	  do i=1, ntp
	      if(ele==tipele(2,i)) then 
		      go to 30
		  end if
	  end do
	  
!c	  do i=1, ntp
!c	      if(ele/=tipele(2,i)) then   !element not at the crack tip
!c		      ntel=tipele(1,i)
!c		      yt=(y(ix(1,ntel))+u(id(1,ix(1,ntel)))
!c     >            +y(ix(3,ntel))+u(id(1,ix(3,ntel))))/2.0
!c	          zt=(z(ix(1,ntel))+u(id(2,ix(1,ntel)))
!c     >            +z(ix(3,ntel))+u(id(2,ix(3,ntel))))/2.0
!c	          ds=sqrt((ye-yt)**2+(ze-zt)**2)
!c			  if(ds<ds0) then
!c			      cflag=0
!c			  end if
!c		  end if
!c	  end do
	  
	  do i=1,numelt
	      if (ElemFractCode(i) > 0) then
		      ntel=i
			  yt=(y(ix(1,ntel))+y(ix(3,ntel)))/2.0
	          zt=(z(ix(1,ntel))+z(ix(3,ntel)))/2.0
	          ds=sqrt((ye-yt)**2+(ze-zt)**2)
			  if(ds<ds0) then   !- if this elem i 
			      cflag=0
				  go to 30
			  end if
		  end if
	  end do
	  
   30 end