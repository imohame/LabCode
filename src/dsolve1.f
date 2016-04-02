      subroutine dsolve1 (dt,time,NVAR1)
	  ! This Subroutine: 
	  ! ===============================================================
	  !
	  ! Computes Rho_m(alpha) & Rho_im(alpha) 
	  !
	  ! ===============================================================
	  ! LOCAL PARAMETERS & VARIABLES
	  ! ++++++++++++++++++++++++++++
	   integer NVAR1,nmo,nim
       parameter (nss = 24)
	   dimension ystart1(NVAR1)
	  ! GLOBAL PARAMETERS & VARIABLES
	  ! +++++++++++++++++++++++++++++
      common/wblock2/g_source, g_immob, g_minter, g_recov, b_v,
     1                b_vvec(87),nmo,nim
      common/wblock4/ den_m(nss), den_im(nss), gdot(nss), nnns, nnne
	  common /nab/ gn_bcc(3,43,43), n_bcc(3,43,43), np(43)
	  common /nab1/ gn_fcc(3,18,18), n_fcc(3,18,18), np1(18)	  
	  common /nab2/ gn_hcpt(3,86,86), n_hcpt(3,86,86), np2(86)
	  common /nab3/gn_hcp(3,87,87),n_hcp(3,87,87),den_im2(87),np3(87)
	  common /aij/ aijhcp(24,87), aijfcc(12,18), aijbcc(24,43),
	1         aijhcpt(24,86)
	  
	  ! ===============================================================
      ! Initialize Variables for ODEINT1
	  !++++++++++++++++++++++++++++++++++
	  ! ** ODEINT1(DENSTY,NVAR1,T1,T2,EPS1,H11,HMIN1,NOK1,NBAD1) **
	  ! NVAR1      = 2													!This Line: No. Variables to solve for
	  ! DENSTY(1)  = den_im(nnns)										!This Line: 1st Variable Starting Value (Non Self-Starting Method Used)
	  ! DENSTY(2)  = den_m (nnns)										!This Line: 2nd Variable Starting Value (Non Self-Starting Method Used)
	  ! T1         = time												!This Line: Time Point to move from
	  ! T2         = time+dt											    !This Line: Time Point to move to
	  ! EPS1       = 1.00e-05											!This Line: Tolerance on Convergnece Accuracy
	  ! H11        = dt													!This Line: Time Step Size
	  ! HMIN1      = H11/1.0e+4											!This Line: Minumum dt Size
	  ! ========================
	 
	   
		 

		  do i = 1, nmo
			  ystart1(i) = den_m(i)
		  end do
		  do i = 1, nim
			  ystart1(nmo+i) = den_im2(i)
		  end do
		  

	 ! ===============================================================

      CALL ODEINT1(ystart1,NVAR1,time,time+dt,1.00e-4,dt,dt/1000,
	1                   NOK1,NBAD1)!This Line: Call Nested Sequence of Subroutines for Adaptive RK4 Integration

		  


		  do i = 1, nim
              den_im2(i) = ystart1(nmo+i)
              if (i .le. nmo) then
                  den_m(i) = ystart1(i)
                  den_im(i) = ystart1(nmo+i)
              endif
		  end do


      RETURN
      END
