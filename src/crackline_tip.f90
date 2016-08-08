      subroutine crackline_tip(ele, y, z, ix, id, u)
	  
	  parameter (nume=40000)
	  common /crackline/ ncleave(3,nume), elecrack(4,nume),nodeflag(4,nume)
	  common /overlapping/ intersec(4, nume), area_coeff(nume), update_flag
	  common/meshnum/ numnpo, numelto
	  
	  dimension y(*), z(*), ix(4,*), id(2,*), u(*)
	  
	  real ncleave, k1, b1, ny, nz, yi, zi
	  real k2, b2, y1, z1, y2, z2
	  real yc, zc, intersec
	  integer p, elecrack, ele, numelto, bflag
	  
	  p=2
	  
	  ny=ncleave(2,ele)
	  nz=ncleave(3,ele)

      yi=intersec(1,ele)
	  zi=intersec(2,ele)
	  if(elecrack(1,ele)==1) then
	      y1=y(ix(1,ele))+u(id(1,ix(1,ele)))
	      z1=z(ix(1,ele))+u(id(2,ix(1,ele)))
	      y2=y(ix(2,ele))+u(id(1,ix(2,ele)))
	      z2=z(ix(2,ele))+u(id(2,ix(2,ele)))
	  else if(elecrack(1,ele)==2) then
	      y1=y(ix(2,ele))+u(id(1,ix(2,ele)))
	      z1=z(ix(2,ele))+u(id(2,ix(2,ele)))
	      y2=y(ix(3,ele))+u(id(1,ix(3,ele)))
	      z2=z(ix(3,ele))+u(id(2,ix(3,ele)))
	  else if(elecrack(1,ele)==3) then
	      y1=y(ix(3,ele))+u(id(1,ix(3,ele)))
	      z1=z(ix(3,ele))+u(id(2,ix(3,ele)))
	      y2=y(ix(4,ele))+u(id(1,ix(4,ele)))
	      z2=z(ix(4,ele))+u(id(2,ix(4,ele)))
	  else
	      y1=y(ix(4,ele))+u(id(1,ix(4,ele)))
	      z1=z(ix(4,ele))+u(id(2,ix(4,ele)))
	      y2=y(ix(1,ele))+u(id(1,ix(1,ele)))
	      z2=z(ix(1,ele))+u(id(2,ix(1,ele)))
	  end if
	  
      if(nz==0) then
!!!c   crack go through edge 1 and 3 perpendicularly
!!!c	      if(elecrack(1,ele)==1) then
!!!c              elecrack(3,ele)=3
!!!c			  elecrack(4,ele)=1
!!!c			  elecrack(2,ele)=3   ! edge status from 2 to 3(cracked)
!!!c			  call tipaft(ix, ele, 1)
!!!c			  call FracUpdateTipAfterCrack(ix, ele, 3, 2)
!!!c			  intersec(3,ele)=1.0-intersec(1,ele)
!!!c			  intersec(4,ele)=1.0-intersec(2,ele)	  
!!!c          else if(elecrack(1,ele)==3) then
!!!c              elecrack(3,ele)=1
!!!c			  elecrack(4,ele)=1
!!!c			  elecrack(2,ele)=3   ! edge status from 2 to 3(cracked)
!!!c			  call tipaft(ix, ele, 3)
!!!c			  call FracUpdateTipAfterCrack(ix, ele, 1, 2)
!!!c			  intersec(3,ele)=1.0-intersec(1,ele)
!!!c			  intersec(4,ele)=1.0-intersec(2,ele)
!!!c          end if
          
	  
      else	  
	      yi=y1+(y2-y1)*yi
	      zi=z1+(z2-z1)*zi
	           
	      k1=-ny/nz
	      b1=zi-k1*yi
	      
!!c          call tipaft(ix, ele, elecrack(1,ele))		  
	      do i=1,4
		  
	          if(i/=elecrack(1,ele)) then
	              if (i==1) then
	                  y1=y(ix(1,ele))+u(id(1,ix(1,ele)))
	                  z1=z(ix(1,ele))+u(id(2,ix(1,ele)))
	                  y2=y(ix(2,ele))+u(id(1,ix(2,ele)))
	                  z2=z(ix(2,ele))+u(id(2,ix(2,ele)))
	              else if(i==2) then
	                  y1=y(ix(2,ele))+u(id(1,ix(2,ele)))
	                  z1=z(ix(2,ele))+u(id(2,ix(2,ele)))
	                  y2=y(ix(3,ele))+u(id(1,ix(3,ele)))
	                  z2=z(ix(3,ele))+u(id(2,ix(3,ele)))
	              else if(i==3) then
	                  y1=y(ix(3,ele))+u(id(1,ix(3,ele)))
	                  z1=z(ix(3,ele))+u(id(2,ix(3,ele)))
	                  y2=y(ix(4,ele))+u(id(1,ix(4,ele)))
	                  z2=z(ix(4,ele))+u(id(2,ix(4,ele)))
	              else
	                  y1=y(ix(4,ele))+u(id(1,ix(4,ele)))
	                  z1=z(ix(4,ele))+u(id(2,ix(4,ele)))
	                  y2=y(ix(1,ele))+u(id(1,ix(1,ele)))
	                  z2=z(ix(1,ele))+u(id(2,ix(1,ele)))
	              end if
          
                  if (y1==y2) then
                      yc=y1
 	                  zc=k1*yc+b1
                  else 
	                 k2=(z2-z1)/(y2-y1)
	                 b2=z1-k2*y1
	                 if(k1 .ne. k2) then
	                     yc=(b2-b1)/(k1-k2)
	                     zc=k1*yc+b1
	                 end if
                  end if

!!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!!c			  
!!!c            estimate how crack intersect with element edge
!!!c
!!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc			  
	             if (yc>=min(y1,y2) .and. yc<=max(y1,y2) .and. &
                     zc>=min(z1,z2) .and. zc<=max(z1,z2)) then
	                 if(y1==y2) then
			  	         intersec(p+1, ele)=0.0
			  	     else
	                     intersec(p+1, ele)=(yc-y1)/(y2-y1)
			  	     end if
			  	     if(z1==z2) then
			  	         intersec(p+2, ele)=0.0
			  	     else
		                 intersec(p+2, ele)=(zc-z1)/(z2-z1)
			  	     end if
	                 elecrack(p+1, ele)=i
		             elecrack(p+2, ele)=1
					 elecrack(2,ele)=3   ! edge status from 2 to 3(cracked)

!!!c   update element crack variable for element ahead of crack front				  
                     call FracUpdateTipAfterCrack(ix, ele, i, p)
					 
	             end if
			 
			  end if
		
          end do	
		  
      end if
	  
	  end
	  