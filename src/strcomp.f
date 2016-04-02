      subroutine strcomp(ele, nstep)

        use CN_Objects_manager
!!!!!!!!!	  parameter (nume=40000)
!!!!!!!!!!!!!!!!!	  common /propag/ sigalt(4,40000)
	  common/wblock8/  abc(573,nume,4), his(573,nume,4)

      real abc, sig(4), cleave(3,3), slip(24,3)  !!!sigalt, 
	  real str, strcl, ncl(3), stro, strsl, nsl(3)
	  integer ele, nstep, nsp
	  
	  str=0.0
	  strcl=0.0
	  
	  stro=0.0
	  strsl=0.0
      
	  do i=1,3
	      ncl(i)=0.0
		  nsl(i)=0.0
	  end do
	  
	  nsp=0
	  
!!!!!!!!!!!!	  do i=1,4
!!!!!!!!!!!!	      sig(i)=sigalt(i,ele)
!!!!!!!!!!!!	  end do
	  CALL CNmanager_Get_sigalt(ele,sig)
c     cleavage plane
      do i=1,3
          cleave(1,i)=abc(362+i,ele,1)
          cleave(2,i)=abc(365+i,ele,1)
          cleave(3,i)=abc(368+i,ele,1)
      end do

	  do i=1,3
	      str = sig(1)*cleave(i,2)**2.0+sig(2)*cleave(i,3)**2.0+
     >          sig(4)*cleave(i,2)*cleave(i,3)*2.0
	      if(str>strcl) then
		      strcl=str
			  do j=1,3
			      ncl(j)=cleave(i,j)
			  end do 
		  end if
	  end do 
	  
c     slip plane
	  do i=1,24
	      slip(i,1)=abc(105+i,ele,1)
		  slip(i,2)=abc(129+i,ele,1)
		  slip(i,3)=abc(153+i,ele,1)
	  end do
	  
	  do i=1,24
	      stro = sig(1)*slip(i,2)**2.0+sig(2)*slip(i,3)**2.0+
     >          sig(4)*slip(i,2)*slip(i,3)*2.0
	      if(stro>strsl) then
		      nsp=i
		      strsl=stro
			  do j=1,3
			      nsl(j)=slip(i,j)
			  end do 
		  end if
	  end do
	  
	  write(969,*) ele, nstep, strcl, strsl, ncl(1), ncl(2), ncl(3),
     >	  nsl(1), nsl(2), nsl(3), nsp
	  
	  end