      subroutine total_Qe(y, z, ix)
	  integer, parameter :: nume       = 40000
	  common/bk08/kprint,nstep,ite,ilimit,newstf
	  common/wblock8/  abc(573,nume,4), his(573,nume,4) 
	  common /stressflag/ strflag1(nume),ElemDecayCount(nume)
	  common/meshnum/ numnpo, numelto
	  common /t_Qe/ totqe, totqe_last, ele_area
	  
	  dimension y(*), z(*), ix(4,*)
	  
	  integer numelto, strflag1, nstep, ele
	  real totqe, totqe_last, ele_area
	  real y3, z3, y1, z1, ds0
	  
	  if(nstep==0) then            ! only for square element
	      y3=y(ix(3,1))
	      z3=z(ix(3,1))
	      y1=y(ix(1,1))
	      z1=z(ix(1,1))
	      ds0=sqrt((y3-y1)**2+(z3-z1)**2)
		  ele_area=(ds0/sqrt(2.0))**2
	  end if
	  
	  totqe_last=totqe
	  totqe=0.0
	  do ele=1,numelto
	      if(strflag1(ele)<2) then
		      totqe=totqe+abc(3,ele,1)*ele_area
		  end if
      end do
	  
	  write(959,10) nstep, totqe, totqe-totqe_last
	  
   10 format(i8,2(e20.10))
   
      end
	  
	  