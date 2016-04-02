!!!      subroutine thermal(ix)
!!!          use CN_Objects_manager
!!!          dimension ix(4,*)
!!!        do i=1,ElemCountAct
!!!            iElemconnect(i,1)=ix(1,i)
!!!            iElemconnect(i,2)=ix(2,i)
!!!            iElemconnect(i,3)=ix(3,i)
!!!            iElemconnect(i,4)=ix(4,i)
!!!        enddo
!!!          
!!!!!!!!!!!!!c     this subroutine is to solve the thermal FEM problem
!!!!!!      integer ele,nod,ii,connect2(99999,4)
!!!!!!      dimension u(2*nod), yold(nod), zold(nod), y(nod), z(nod),
!!!!!!     1   yzold(2*nod),yz(2*nod),ix(4,*)
!!!!!!
!!!!!!      integer id(2*nod),maxint!,ix(4,ele)
!!!!!!      common /WMLthermali/connect(99999,4)
!!!!!!      common /WMLthermal/thermalflag,thermalconstraint(99999)
!!!!!!      common/bk08/kprint,nstep,ite,ilimit,newstf
!!!!!!      integer thermalconstraint,nstep
!!!!!!      integer connect
!!!!!!
!!!!!!!!!!!!!!!c      ele=fourele/4
!!!!!!!!!!!!!!!c      write(*,*)'in thermal.f'
!!!!!!!!!!!!!!!c      write(*,*)'ele=',ele
!!!!!!!!!!!!!!!c      write(*,*)'nod=',nod
!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!c      do i=1,ele
!!!!!!!!!!!!!!!c      	 connect(i,1)=ix(1,i)
!!!!!!!!!!!!!!!c        connect(i,2)=ix(2,i)
!!!!!!!!!!!!!!!c	 connect(i,3)=ix(3,i)
!!!!!!!!!!!!!!!c	 connect(i,4)=ix(4,i)
!!!!!!!!!!!!!!!c      enddo
!!!!!!      
!!!!!!      do i=1,nod
!!!!!!      yzold(2*i-1)=yold(i)
!!!!!!      yzold(2*i)=zold(i)
!!!!!!      enddo
!!!!!!      
!!!!!!      ii=1
!!!!!!      do i=1,2*nod
!!!!!!      if (id(i).eq.0) then
!!!!!!      yz(i)=yzold(i)
!!!!!!      ii=ii+1
!!!!!!      else
!!!!!!      yz(i)=yzold(i)+u(id(ii))   !ii could be replaced with do loop variable i
!!!!!!      ii=ii+1
!!!!!!      endif
!!!!!!      enddo
!!!!!!      
!!!!!!      do i=1,nod
!!!!!!      y(i)=yz(2*i-1)
!!!!!!      z(i)=yz(2*i)
!!!!!!      enddo
!!!!!!
!!!!!!      
!!!!!!
!!!!!!!!!!!!!      do i=1,2*nod
!!!!!!!!!!!!!c      write(*,*) u(i),y(i),z(i),id(i)
!!!!!!!!!!!!!      enddo
!!!!!!!!!!!!!c      write(*,*) 'before thermalFEMcall'
!!!!!!      maxint=4
!!!!!!      call thermalFEM(y,z,ele,nod,maxint)
!!!!!!     
!!!
!!!      return
!!!      end
