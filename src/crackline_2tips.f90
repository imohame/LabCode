      subroutine crackline_2tips(ele, ix)
	  
	  parameter (nume=40000)
	  common /crackline/ ncleave(3,nume), elecrack(4,nume), nodeflag(4,nume)
	  
	  dimension ix(4,*)
	  
	  integer elecrack, ele,ix
	  
	  elecrack(2,ele)=3
!!!c	  call FracUpdateTipBeforeCrack(ix, ele, elecrack(1,ele))
	  elecrack(4,ele)=3
!!!c	  call FracUpdateTipBeforeCrack(ix, ele, elecrack(3,ele))
	  
	  end
	  