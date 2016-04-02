      subroutine FCN1(x1,f1,n1a,par1)
	  integer nmo, nim
      parameter (nss = 24)
      common /wblock2/ g_source, g_immob, g_minter, g_recov, b_v,
	1      b_vvec(87),nmo,nim
      common /wblock4/ den_m(nss), den_im(nss), gdot(nss), nnns, nnne
      common /wblock5/ enthalpy_coef, thermal_coef, temp
	  common /nab/ gn_bcc(3,43,43), n_bcc(3,43,43), np(43)
	  common /nab1/ gn_fcc(3,18,18), n_fcc(3,18,18), np1(18)	  
	  common /nab2/ gn_hcpt(3,86,86), n_hcpt(3,86,86), np2(86)
	  common /nab3/gn_hcp(3,87,87),n_hcp(3,87,87),den_im2(87),np3(87)
	  common /aij/ aijhcp(24,87), aijfcc(12,18), aijbcc(24,43),
	1         aijhcpt(24,86)

	  
	  real ll, rmin(87), rimin(87), srm(24), sgdot(87), c4, r1, r2
	  real yim(24), dum, gdum(24), y1(67), ysat, fphi
	  integer n, n1, n2, nssim, nssm, n1a, ni, nj, nk
      dimension x1(n1a),f1(n1a),par1(n1a+1),yprime1(n1a)

      do i = 1, n1a
	  y1(i) = x1(i)
	  end do
	  
	  ll = 0.0											!Mean free path
	  c4 = 0.0											!Gsource
	  
	  if (nim == 18) then
		  nssim = 18
		  nssm = 12
		  ysat = 1e14
		  fphi = 0.1
	  else if (nim == 43) then
		  nssim = 43
		  nssm = 24
		  ysat = 1e16
		  fphi = 0.05
	  else if (nim == 86) then
		  nssim = 86
		  nssm = 24
		  ysat = 1e16
		  fphi = 0.05
	  else
		  nssim = 87
		  nssm = 24
		  ysat = 1e16
		  fphi = 0.05

	  endif
	  
	  do i = 1, nssim
		  ni = nssm + i
		  ll = ll + y1(ni)
		  if (i .le. nssm) then
			  gdum(i) = abs(gdot(i))/b_vvec(i)
			  srm(i) = 0.7746*(y1(i)+y1(ni))
			  sgdot(i) = 0.7746*gdum(i)  
			  rimin(i) = 0.7746*(y1(i) + y1(ni))*gdum(i) 
		  else
			  rimin(i) = 0.0
			  sgdot(i) = 0.0
		  endif
	  end do
	  
	  
	  
	  ll = 1.0/sqrt(ll)
	  fphi = fphi*ll
	  
	  if (nim == 43) then
	  do j = 1, nssm
	  nj = nssm + j
	  do k = j+1, nssim 
	  nk = nssm + k
	  if (k .le. nssm) then
		  r1 = aijbcc(j,k)*(y1(k)+y1(nk))
		  r2 = aijbcc(j,k)*(y1(j)+y1(nj))
		  srm(j) = srm(j) + r1
		  srm(k) = srm(k) + r2
		  sgdot(j) = sgdot(j) + aijbcc(j,k)*gdum(k)
		  sgdot(k) = sgdot(k) + aijbcc(j,k)*gdum(j)
		  r1 = r2*gdum(k) + r1*gdum(j)
		  rimin(n_bcc(1,j,k)) = rimin(n_bcc(1,j,k)) + gn_bcc(1,j,k)*r1
		  rimin(n_bcc(2,j,k)) = rimin(n_bcc(2,j,k)) + gn_bcc(2,j,k)*r1
		  rimin(n_bcc(3,j,k)) = rimin(n_bcc(3,j,k)) + gn_bcc(3,j,k)*r1
	  else
		  r1 = y1(nk)*aijbcc(j,k)
		  srm(j) = srm(j) + r1
		  sgdot(k) = sgdot(k) + aijbcc(j,k)*gdum(j)
		  r1 = r1*gdum(j)
		  rimin(n_bcc(1,j,k)) = rimin(n_bcc(1,j,k)) + gn_bcc(1,j,k)*r1
		  rimin(n_bcc(2,j,k)) = rimin(n_bcc(2,j,k)) + gn_bcc(2,j,k)*r1
		  rimin(n_bcc(3,j,k)) = rimin(n_bcc(3,j,k)) + gn_bcc(3,j,k)*r1
	  endif
	  end do
	  end do
	  else if (nim==18) then
	  do j = 1, nssm
	  nj = nssm + j
	  do k = j+1, nssim 
	  nk = nssm + k
	  if (k .le. nssm) then
		  r1 = aijfcc(j,k)*(y1(k)+y1(nk))
		  r2 = aijfcc(j,k)*(y1(j)+y1(nj))
		  srm(j) = srm(j) + r1
		  srm(k) = srm(k) + r2
		  sgdot(j) = sgdot(j) + aijfcc(j,k)*gdum(k)
		  sgdot(k) = sgdot(k) + aijfcc(j,k)*gdum(j)
		  r1 = r2*gdum(k) + r1*gdum(j)
		  rimin(n_fcc(1,j,k)) = rimin(n_fcc(1,j,k)) + gn_fcc(1,j,k)*r1
		  rimin(n_fcc(2,j,k)) = rimin(n_fcc(2,j,k)) + gn_fcc(2,j,k)*r1
		  rimin(n_fcc(3,j,k)) = rimin(n_fcc(3,j,k)) + gn_fcc(3,j,k)*r1
	  else
		  r1 = y1(nk)*aijfcc(j,k)
		  srm(j) = srm(j) + r1
		  sgdot(k) = sgdot(k) + aijfcc(j,k)*gdum(j)
		  r1 = r1*gdum(j)
		  rimin(n_fcc(1,j,k)) = rimin(n_fcc(1,j,k)) + gn_fcc(1,j,k)*r1
		  rimin(n_fcc(2,j,k)) = rimin(n_fcc(2,j,k)) + gn_fcc(2,j,k)*r1
		  rimin(n_fcc(3,j,k)) = rimin(n_fcc(3,j,k)) + gn_fcc(3,j,k)*r1
	  endif
	  end do
	  end do
	  else if (nim==86) then
	  do j = 1, nssm
	  nj = nssm + j
	  do k = j+1, nssim 
	  nk = nssm + k
	  if (k .le. nssm) then
		  r1 = aijhcpt(j,k)*(y1(k)+y1(nk))
		  r2 = aijhcpt(j,k)*(y1(j)+y1(nj))
		  srm(j) = srm(j) + r1
		  srm(k) = srm(k) + r2
		  sgdot(j) = sgdot(j) + aijhcpt(j,k)*gdum(k)
		  sgdot(k) = sgdot(k) + aijhcpt(j,k)*gdum(j)
		  r1 = r2*gdum(k) + r1*gdum(j)
		  rimin(n_hcpt(1,j,k)) = rimin(n_hcpt(1,j,k)) + gn_hcpt(1,j,k)*r1
		  rimin(n_hcpt(2,j,k)) = rimin(n_hcpt(2,j,k)) + gn_hcpt(2,j,k)*r1
		  rimin(n_hcpt(3,j,k)) = rimin(n_hcpt(3,j,k)) + gn_hcpt(3,j,k)*r1
	  else
		  r1 = y1(nk)*aijhcpt(j,k)
		  srm(j) = srm(j) + r1
		  sgdot(k) = sgdot(k) + aijhcpt(j,k)*gdum(j)
		  r1 = r1*gdum(j)
		  rimin(n_hcpt(1,j,k)) = rimin(n_hcpt(1,j,k)) + gn_hcpt(1,j,k)*r1
		  rimin(n_hcpt(2,j,k)) = rimin(n_hcpt(2,j,k)) + gn_hcpt(2,j,k)*r1
		  rimin(n_hcpt(3,j,k)) = rimin(n_hcpt(3,j,k)) + gn_hcpt(3,j,k)*r1
	  endif
	  end do
	  end do
	  else if (nim==87) then
	  do j = 1, nssm
	  nj = nssm + j
	  do k = j+1, nssim 
	  nk = nssm + k
	  if (k .le. nssm) then
		  r1 = aijhcp(j,k)*(y1(k)+y1(nk))
		  r2 = aijhcp(j,k)*(y1(j)+y1(nj))
		  srm(j) = srm(j) + r1
		  srm(k) = srm(k) + r2
		  sgdot(j) = sgdot(j) + aijhcp(j,k)*gdum(k)
		  sgdot(k) = sgdot(k) + aijhcp(j,k)*gdum(j)
		  r1 = r2*gdum(k) + r1*gdum(j)
		  rimin(n_hcp(1,j,k)) = rimin(n_hcp(1,j,k)) + gn_hcp(1,j,k)*r1
		  rimin(n_hcp(2,j,k)) = rimin(n_hcp(2,j,k)) + gn_hcp(2,j,k)*r1
		  rimin(n_hcp(3,j,k)) = rimin(n_hcp(3,j,k)) + gn_hcp(3,j,k)*r1
	  else
		  r1 = y1(nk)*aijhcp(j,k)
		  srm(j) = srm(j) + r1
		  sgdot(k) = sgdot(k) + aijhcp(j,k)*gdum(j)
		  r1 = r1*gdum(j)
		  rimin(n_hcp(1,j,k)) = rimin(n_hcp(1,j,k)) + gn_hcp(1,j,k)*r1
		  rimin(n_hcp(2,j,k)) = rimin(n_hcp(2,j,k)) + gn_hcp(2,j,k)*r1
		  rimin(n_hcp(3,j,k)) = rimin(n_hcp(3,j,k)) + gn_hcp(3,j,k)*r1
	  endif
	  end do
	  end do	  	  
	  endif
	  
	  c4 = 0.1/(ll*ll)
	  
*------------------------------------
	  
	  do i = 1, nssim
		  ni = nssm + i
		  if (i .le. nssm) then
		  yprime1(i) = (gdum(i)*(c4*y1(ni)/y1(i)-srm(i)) 
     >                 - y1(i)*sgdot(i))*fphi
		  endif
		  thermal_coef = 3376.9112279916*(1.0 - sqrt(y1(ni)/ysat))
		  c5 = -1 * thermal_coef/temp
		  yprime1(ni) = fphi*(rimin(i) 
     >                - g_recov*exp(c5)*y1(ni)*sgdot(i))
	  end do

      do i = 1, n1a								!why is this only looped to n1?
	  f1(i) = x1(i) - par1(i) - par1(n1a+1)*yprime1(i)
	  end do
      
      

      
      return
      end
