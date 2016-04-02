      subroutine FCN(x,f,n,par)

      parameter (nss = 24)
      dimension x(n), f(n), par(n+1), r(nss), g(nss)
      common/wblock6/ xmhigh(nss),ref_gamma_dot(nss),Pij_Dijdev(nss), 
     >                tau(nss), p(nss,4), taur(nss), twomu, g
	  integer m
	  
	  m=n+1
	  !write(*,*) 'in fcn n size par xmhigh1 ref_gamma_dot1'
	  !write(*,*) n
	  !write(*,*) size(par)
	  !write(*,*) xmhigh(1)
	  !write(*,*) ref_gamma_dot(1)
	  !write(*,*) x(1)
      do i=1,m
c      write(*,*) 'in fcn,par',par(i)
	  enddo
      dummy = 0.0
      sat   = 1.0/1.0
	  !write(*,*) taur
      do i = 1, n
         r(i) = ((abs(x(i)/taur(i))*sat)**(xmhigh(i)-1.0))*
     >                                 (x(i)/taur(i))*sat
         g(i) = 0.0
      end do

      do i = 1, n
         do j = 1, n
            dummy = r(j)*ref_gamma_dot(j)*(p(i,1)*p(j,1) +
     >                                     p(i,2)*p(j,2) +
     >                                  2* p(i,3)*p(j,3) +
     >                                     p(i,4)*p(j,4))
            g(i) = g(i) + dummy
         end do
      end do

      do i = 1, n
         f(i) = x(i)-par(i)-par(m)*twomu*(Pij_Dijdev(i)-g(i))
c		 write(*,*) 'in FCN, f, x, par',f(i),x(i),par(i)
      end do 


       return
       end
