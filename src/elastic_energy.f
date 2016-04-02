      subroutine elastic_energy(ele, qe)
	  
       use mod_parameters
	  common/custr/   sign1(nelemg),sign2(nelemg),
     > sign3(nelemg),sign4(nelemg)
	  common/wblock12/ Y_modulus(nume),possion_ratio(nume),tau_y(nume)
	  common/range/   mft,mlt,lft,llt,nftm1
	  
	  integer nftm1, ink, ele
	  real E, v, Y_modulus, possion_ratio
	  real epx, epy, epxy, qe
	  
	  qe=0.0
	  ink=ele+nftm1
	  E=Y_modulus(ink)
	  v=possion_ratio(ink)
	  
	  epx=(1+v)/E*((1-v)*sign1(ele)-v*sign2(ele))
	  epy=(1+v)/E*(-v*sign1(ele)+(1-v)*sign2(ele))
	  epxy=(1+v)/E*2*sign4(ele)
	  
	  qe=(sign1(ele)*epx+sign2(ele)*epy+sign4(ele)*epxy)/2.0
	  
	  end
	  
	  
	  
	  
	  
	  
	  
	  