	subroutine edge(y,z,matp)
!!!!!!!        edge(a(k03),a(k04),a(k08))
	
	integer, parameter :: nume   = 40000
	real, parameter    :: tempr  = 293.1
	integer, parameter :: no_mat     = 1000 
	real   , parameter :: xi         = 0.3 
	  common /gbblock/ gbvec(nume,3), gbflag(nume,3)
      common/wblock8/  abc(573,nume,4), his(573,nume,4)
      common /wblock10/ng,grain_mo(1000,3),bv(no_mat,87),nssmat(1000)!changed to accomodate multiple bv lengths WML
     1                                 ,nssimmat(1000)
	  common /aij/ aijhcp(24,87), aijfcc(12,18), aijbcc(24,43),
	1         aijhcpt(24,86)
      common/wblock12/ Y_modulus(nume),possion_ratio(nume),tau_y(nume)
	  common/wblock2/  g_source, g_immob, g_minter, g_recov, b_v,
     1                             b_vvec(87),nmo,nim
	  common /stressflag/ strflag1(nume),ElemDecayCount(nume)
	  
	integer gbflag, i, j, nssmat, ctr, matp
	real gbvec, normgb, abc, aijbcc, aijfcc, aijhcp, aijhcpt, g, bv,
	1         temp, tau_y
	real possion_ratio,thermal_factor,den_im(87),tau(24),taur(24)
	integer nintg, strflag1
	
	dimension y(*), z(*), matp(*)
	
	nintg=1
	
	do i = 1, nume
	  if(strflag1(i)==0) then     ! only for non-cracked elements
	    if (gbflag(i,1) .ne. 0) then
	        gbvec(i,1) = -(z(gbflag(i,2))-z(gbflag(i,3)))
	        gbvec(i,2) = (y(gbflag(i,2))-y(gbflag(i,3)))
	        normgb = sqrt(gbvec(i,1)**2.0+gbvec(i,2)**2.0)
	        gbvec(i,1) = gbvec(i,1)/normgb
	        gbvec(i,2) = gbvec(i,2)/normgb
	        gbvec(i,3) = normgb
	
	        g  = Y_modulus(i)/(2.0*(1.0 + possion_ratio(i)))
	        temp = abc(4,i,1)
	        temp = max(temp,tempr)
	        thermal_factor = (tempr/temp)**xi
	
	        do j = 1, nssmat(matp(i))
	            den_im(j) = abc(81+j,i,1)
	            tau(j) = abc(9+j,i,1)
	        end do
	
	        if (nssmat(matp(i)) == 12 .and. nim==18) then
	            ctr = 1
	            do j = 1, 18
	                if (j .ge. 13) then
		                den_im(j) = abc(397+ctr,i,nintg)
		                ctr = ctr + 1
	                end if
	            end do
	        end if
	  
	        if (nssmat(matp(i)) == 24 .and. nim==43) then
	            ctr = 1
	            do j = 1, 43
	                if (j .ge. 25) then
		                den_im(j) = abc(397+ctr,i,nintg)
		                ctr = ctr + 1
	                end if
	            end do
	        end if

	        if (nssmat(matp(i)) == 24 .and. nim==86) then
	            ctr = 1
	            do j = 1, 86
	                if (j .ge. 25) then
		                den_im(j) = abc(397+ctr,i,nintg)
		                ctr = ctr + 1
	                end if
	            end do
	        end if

	        if (nssmat(matp(i)) == 24 .and. nim==87) then
	            ctr = 1
	            do j = 1, 87
	                if (j .ge. 25) then
		                den_im(j) = abc(397+ctr,i,nintg)
		                ctr = ctr + 1
	                end if
	            end do
	        end if
	  
	        if (nssmat(matp(i)) == 12 .and. nim==18) then
	            do j = 1, 12
	                taur(j) = tau_y(i)
	                do k = 1, 18
		                taur(j) = taur(j) + 
     >						 aijfcc(j,k)*g*b_v*sqrt(den_im(k))
	                end do
	                taur(j) = taur(j)*thermal_factor
	            end do
	        end if
	   
	        if (nssmat(matp(i)) == 24 .and. nim==43) then
	            do j = 1, 24
	                taur(j) = tau_y(i)
	                do k = 1, 43
		                taur(j) = taur(j) + 
     >					     aijbcc(j,k)*g*b_v*sqrt(den_im(k))
	                end do
	                taur(j) = taur(j)*thermal_factor
	            end do
	        end if

	        if (nssmat(matp(i)) == 24 .and. nim==86) then
	            do j = 1, 24
	                taur(j) = tau_y(i)
	                do k = 1, 86
		                taur(j) = taur(j) + 
     >					     aijhcpt(j,k)*g*b_v*sqrt(den_im(k))
	                end do
	                taur(j) = taur(j)*thermal_factor
	            end do
	        end if	
	        
	        if (nssmat(matp(i)) == 24 .and. nim==87) then
	            do j = 1, 24
	                taur(j) = tau_y(i)
	                do k = 1, 87
		                taur(j) = taur(j) + 
     >					     aijhcp(j,k)*g*b_v*sqrt(den_im(k))
	                end do
	                taur(j) = taur(j)*thermal_factor
	            end do
	        end if	        
	           
	        do j = 1, nssmat(matp(i))
	            abc(337+j,i,1) = abs(tau(j)/taur(j))
	            his(337+j,i,1) = abs(tau(j)/taur(j))
	        end do
	
	    end if
	  end if
	  
	end do
	
	return 
	end
	
	