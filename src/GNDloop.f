      subroutine GNDloop(y, z, ix, id, u, numelt, numnp)
 
      parameter(nume = 40000) 
	  parameter(maxint = 4) 
	  common/wblock2/  g_source, g_immob, g_minter, g_recov, b_v,
     1                             b_vvec(87),nmo,nim
	  common/wblock8/  abc(421,nume,4), his(421,nume,4)
	  common /GND_loop/ gradslip(2,24,40000), rho_gnd(2,24,40000)
	  real*8 eta(maxint),psi(maxint)
	  real*8 N(maxint,4), ShpFcnDeriv(maxint,2,4)
	  real*8 Jmatrix(2,2), Jacob(numelt,maxint), Gamma(2,2),
     > B(maxint,2,4)
	  real*8 B1(2,4),B2(2,4),B3(2,4),B4(2,4),Coords(maxint,2)
	  real*8 gradslip1a(2,1), gradslip2a(2,1), gradslip3a(2,1)
	  real*8 gradslip4a(2,1), gradslipa(2,1)
	  integer connect(40000,4), numelt, numnp
	  integer nelenode2(40000), nsi
	  real slipnode(40000,24), slipp(4,1), ly, lz
	  real gradslip, rho_gnd
	  
	  dimension y(*), z(*), ix(4,*), id(2,*), u(*)
	        
      do i=1,numelt
	      do j=1,4
		      connect(i,j)=ix(j,i)
		  end do
	  end do	  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c               distribute shear slip to node
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do i=1,numnp
	      nelenode2(i)=0
		  do nsi=1,24
		      slipnode(i,nsi)=0.0
		  end do
	  end do
	  
      do i=1,numelt
          do j=1,4
              nelenode2(connect(i,j))=nelenode2(connect(i,j))+1 
          end do
      end do
	  
      do nsi=1,24     
	      do i=1,numelt
		      do j=1,4
		          slipnode(connect(i,j),nsi)=slipnode(connect(i,j),nsi)
     >				                    +abc(549+nsi,i,1)
	          end do
		  end do
		  do i=1,numnp
		      if(nelenode2(i)==0) then
			      write(*,*) 'warning a/0, check GNDloop.f'
			  end if
			  slipnode(i,nsi)=slipnode(i,nsi)/nelenode2(i)
		  end do
	  end do	   

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c              calculate B matrix
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
      psi(1)=-.57735
      eta(1)=-.57735
      psi(2)=.57735
      eta(2)=-.57735
      psi(3)=.57735
      eta(3)=.57735
      psi(4)=-.57735
      eta(4)=.57735            
      call shapefunctions(psi,eta,N,ShpFcnDeriv,maxint)      
      
      do i=1,numelt
c      reading nodes associated with element i 
            node1=connect(i,1)
            node2=connect(i,2)
       	    node3=connect(i,3)
            node4=connect(i,4)
            Coords(1,1)=y(node1)+u(id(1,node1))          ! current configuration
            Coords(1,2)=z(node1)+u(id(2,node1))
            Coords(2,1)=y(node2)+u(id(1,node2))
            Coords(2,2)=z(node2)+u(id(2,node2))
            Coords(3,1)=y(node3)+u(id(1,node3))
            Coords(3,2)=z(node3)+u(id(2,node3))
            Coords(4,1)=y(node4)+u(id(1,node4))
            Coords(4,2)=z(node4)+u(id(2,node4))
c     given shape function derivatives wrt isoparametric coords
c     and coordinate matrix for the 4 Gauss points, find the B-matrix  
       do j=1,maxint
          Jmatrix=matmul(ShpFcnDeriv(j,1:2,1:4),Coords)
          Jacob(i,j)=Jmatrix(1,1)*Jmatrix(2,2)-Jmatrix(2,1)*Jmatrix(1,2)
          Gamma(1,1)=1/Jacob(i,j)*Jmatrix(2,2)
          Gamma(1,2)=-1/Jacob(i,j)*Jmatrix(1,2)
          Gamma(2,1)=-1/Jacob(i,j)*Jmatrix(2,1)
          Gamma(2,2)=1/Jacob(i,j)*Jmatrix(2,2)
          B(j,1:2,1:4)=matmul(Gamma,ShpFcnDeriv(j,1:2,1:4))
       enddo
        B1(1:2,1:4)=B(1,1:2,1:4)
        B2(1:2,1:4)=B(2,1:2,1:4)
        B3(1:2,1:4)=B(3,1:2,1:4)
        B4(1:2,1:4)=B(4,1:2,1:4)        
               
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         calculate shear slip gradient and GND
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do nsi=1,24
            slipp(1,1)=slipnode(connect(i,1),nsi)
			slipp(2,1)=slipnode(connect(i,2),nsi)
			slipp(3,1)=slipnode(connect(i,3),nsi)
			slipp(4,1)=slipnode(connect(i,4),nsi)
		    gradslip1a=matmul(B1, slipp)
			gradslip2a=matmul(B2, slipp)
			gradslip3a=matmul(B3, slipp)
			gradslip4a=matmul(B4, slipp)
		    gradslipa(1:2,1)=gradslip1a(1:2,1)+gradslip2a(1:2,1)
     >                   +gradslip3a(1:2,1)+gradslip4a(1:2,1)
	        gradslipa(1,1)=gradslipa(1,1)/4.0
            gradslipa(2,1)=gradslipa(2,1)/4.0
				
		    gradslip(1,nsi,i)=gradslipa(1,1)
		    gradslip(2,nsi,i)=gradslipa(2,1)
c           dislocation line vector
            ly=-(abc(177+nsi,i,1)*abc(153+nsi,i,1)
     >          -abc(105+nsi,i,1)*abc(225+nsi,i,1))
            lz=	abc(177+nsi,i,1)*abc(129+nsi,i,1)
     >          -abc(201+nsi,i,1)*abc(105+nsi,i,1)				

c               calculate GND	
        rho_gnd(1,nsi,i)=-(gradslip(1,nsi,i)*abc(201+nsi,i,1)    ! edge GND
     >          +gradslip(2,nsi,i)*abc(225+nsi,i,1))/b_v
            rho_gnd(2,nsi,i)=-(gradslip(1,nsi,i)*ly                  ! screw GND
     >                       +gradslip(2,nsi,i)*lz)/b_v			
	 
			rho_gnd(1,nsi,i)=abs(rho_gnd(1,nsi,i))
			rho_gnd(2,nsi,i)=abs(rho_gnd(2,nsi,i))
		end do			
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc			
		
      enddo

      return
      end
      
      
    
