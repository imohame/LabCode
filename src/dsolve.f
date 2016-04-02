        subroutine dsolve(dt,time,nssm)

*  ____________________________________________________________________
* |                                                                    |
* |...  Define the number of the intial value systems                  |
* |...  Set-up for explicit-implicit method based on stiffness algothm |
* |                                                                    |
* |           ystart : intial values for resolved shear-stress         |
* |           x1     : time step                                       |
* |           x2     : t+dt                                            |
* |                                                                    |
* |                                     Modified by .................  |
* |                                                 Waeil M. Ashmawi   |
* |                                                 NC State University|
* |                                                 April 07, 1998     |
* |____________________________________________________________________|
*
      parameter (nss = 24)  

      common/wblock6/ xmhigh(nss),ref_gamma_dot(nss),Pij_Dijdev(nss), 
     >                tau(nss), p(nss,4), taur(nss), twomu, g
       
      dimension ystart(nss)

*
* ======================== S T A R T   C O D E ========================
* 
	  !write(*,*) 'beginning of dsolve'
	  !write(*,*) nssm
      do i = 1, nssm
         ystart(i) = tau(i)
      end do
      x1   = time
      x2   = time + dt

* [Set accuracy requirements]
* ...........................
      eps  = 1.00E-04
      h1   = dt

* [Set min.  time-step]
* .....................
      hmin = h1/1000

* [Call solver, get tau(i) back]
* ..............................
      call odeint(ystart,nssm,x1,x2,eps,h1,hmin,nok,nbad)
      !write(*,*) 'in dsolve after odeint'
	  !write(*,*) nssm
      do i = 1, nssm
         tau(i) = ystart(i)
      end do
	  !write(*,*) 'end of dsolve'
	  !write(*,*) nssm
*
* ========================== E N D   C O D E ==========================
* 
        
      return
      end


