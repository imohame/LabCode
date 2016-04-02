      subroutine unit_vector(m,n,slipn,slips)

      dimension slipn(n,m), slips(n,m)
      real norm_n, norm_s
	  integer nn
	  !write(*,*) 'unit_vector beginning'
	  !write(*,*) nn
	  !write(*,*) slipn(1,1)
	  !write(*,*) slipn(1,2)
	  !write(*,*) slipn(1,3)
	  !write(*,*) slips(1,1)
	  !write(*,*) slips(1,2)
	  !write(*,*) slips(1,3)
      sum_n = 0.0
      sum_s = 0.0
	  nn=n
      do i = 1, nn
         do k = 1, m 
            sum_n = sum_n + (abs(slipn(i,k)))**2
            sum_s = sum_s + (abs(slips(i,k)))**2
         end do
         norm_n = sqrt(sum_n)        
         norm_s = sqrt(sum_s)
		 if(abs(norm_n)<1.0E-20) then
		     norm_n=1.0e-20
		 end if
		 if(abs(norm_s)<1.0E-20) then
		     norm_s=1.0e-20
		 end if
         do j = 1, m
            slipn(i,j) = (1.0/norm_n)*slipn(i,j)
            slips(i,j) = (1.0/norm_s)*slips(i,j)
         end do
         sum_n = 0.0
         sum_s = 0.0
      end do
	  !write(*,*) 'unit_vector has been called'
	  !write(*,*) nn
	  !write(*,*) slipn(1,1)
	  !write(*,*) slipn(1,2)
	  !write(*,*) slipn(1,3)
	  !write(*,*) slips(1,1)
	  !write(*,*) slips(1,2)
	  !write(*,*) slips(1,3)

      return 
      end 
     
