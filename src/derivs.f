      subroutine derivs (x,y,yprime,nvar)

*  ____________________________________________________________________
* |                                                                    |
* |...  Define the system of nonlinear differential equations for each |
* |...  active slip-system (k)                                         |
* |                                                                    |
* |           tau_dot(k)   = 2*mu*Pij(k)*[Dij_dev - Dij_p_dev]         |
* |           Dij_p_dev    = Pij(n)*gamma_dot(n)                       |
* |           gamma_dot(n) = ref_gamma_dot(n)*[tau(n)/tau_r(n)]**(1/m) |
* |                                                                    |
* |____________________________________________________________________|
*

      parameter (nss = 24)

      common/wblock6/ xmhigh(nss),ref_gamma_dot(nss),Pij_Dijdev(nss), 
     >                tau(nss), p(nss,4), taur(nss), twomu, g
 
      dimension y(nss), yprime(nss), r(nss), g(nss)
	  integer nn,nvar

      dummy = 0.0
      sat   = 1.0/1.0
	  nn=nvar
      !write(*,*) 'beginning of derivs, nn then nvar'
	  !write(*,*) nn
	  !write(*,*) nvar
      do i = 1, nn
         r(i) = ((abs(y(i)/taur(i))*sat)**(xmhigh(i)-1.0))*
     >                                 (y(i)/taur(i))*sat
         g(i) = 0.0
      end do

      do i = 1, nn
         do j = 1, nn
            dummy = r(j)*ref_gamma_dot(j)*(p(i,1)*p(j,1) +
     >                                     p(i,2)*p(j,2) +
     >                                   2*p(i,3)*p(j,3) +
     >                                     p(i,4)*p(j,4))
            g(i) = g(i) + dummy
         end do
      end do

      do i = 1, nn
         yprime(i) = twomu*(Pij_Dijdev(i) - g(i))
c		 write(*,*) 'r,g,dummy,nn,Pij_Dijdev(i),rgd',r(i),g(i),
c	1        dummy,nn,Pij_Dijdev(i),ref_gamma_dot(i)
      end do 
	  nvar=nn
	  !write(*,*) 'in derivs'
	  
	  !write(*,*) nvar


      return
      end


































