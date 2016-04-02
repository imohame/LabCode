      subroutine thermalBC(ele,nod)
!!!!c     this is used to set the BC's for the problem
!!!      integer, parameter :: nume   = 40000
!!!      integer ele,nod
!!!
!!!
!!!      common /WMLthermali/connect(nume,4),Telei(nume)
!!!
!!!      common/WMLthermal/thermalflag,thermalconstraint(nume),
!!!     1     Tinit(nume)    !,Tfl(nume)
!!!      common /WMLthermalBC/const,constraint(nume)
!!!      integer const,constraint
!!!!      common /WMLthermalBC/hconstraint(40000),sizeh,ele12size,ele23size,
!!!!     1       ele34size,ele41size,const,constraint(40000),ele12(40000),
!!!!     1       ele23(40000),ele34(40000),ele41(40000)
!!!!      integer hconstraint,sizeh,ele12size,ele23size,ele34size,ele41size
!!!!     1      ,ele12,ele23,ele34,ele41,
!!!      integer aa,bb,cc,dd
!!!      integer thermalconstraint
!!!      integer connect
!!!!c      ele=fourele/4
!!!!c      need to finish dimension statements and conversion from matlab
!!!!c      write(*,*) ele,nod
!!!      do i=1,ele
!!!!c      write(*,*) connect(i,1),connect(i,2),connect(i,3),connect(i,4)
!!!       Telei(i)=0.25*Tinit(connect(i,1))+0.25*Tinit(connect(i,2))+
!!!     1     0.25*Tinit(connect(i,3))+0.25*Tinit(connect(i,4))
!!!!c	  write(*,*) Telei(i,1)
!!!      enddo
!!!      
!!!
!!!      const=0
!!!      sizeh=0
!!!      do i=1,nod
!!!!c     eq 1 constant temperature 
!!!!      if (thermalconstraint(i).eq.1.or.thermalconstraint(i).eq.3)then
!!!        if (thermalconstraint(i).eq.1)then
!!!            const=const+1
!!!            constraint(const)=i
!!!            
!!!        endif
!!!!!c     eq 2 convection
!!!!!      if (thermalconstraint(i).eq.2.or.thermalconstraint(i).eq.3) then
!!!!!      sizeh=sizeh+1
!!!!!      hconstraint(sizeh)=i
!!!!!      endif
!!!      enddo
!!!      
!!!!!!      write(*,*)'thermalconstraint(1:nod)'
!!!!!!      write(*,*)thermalconstraint(1:nod)
!!!!!!      
!!!!!!      write(*,*)'constraint(1:const)'
!!!!!!      write(*,*)const,constraint(1:const)
!!!!!!      
!!!!      do i=1,const
!!!!c      write(*,*) constraint(i)
!!!!      enddo
!!!      
!!!!!c     dealing with convection BC's
!!!!!         ele12size=0
!!!!!         ele23size=0
!!!!!         ele34size=0
!!!!!         ele41size=0
!!!!!      do i=1,ele
!!!!!
!!!!!        do j=1,sizeh
!!!!!          aa=connect(i,1)
!!!!!          bb=connect(i,2)
!!!!!          cc=connect(i,3)
!!!!!          dd=connect(i,4)
!!!!!          if (aa.eq.int(hconstraint(j))) then
!!!!!             do kk=1,sizeh
!!!!!                if (bb.eq.int(hconstraint(kk))) then
!!!!!                ele12size=ele12size+1
!!!!!                ele12(ele12size)=i
!!!!!                
!!!!!                elseif (dd.eq.int(hconstraint(kk))) then
!!!!!                ele41size=ele41size+1
!!!!!                ele41(ele41size)=i
!!!!!                endif
!!!!!             enddo
!!!!!          endif
!!!!!          if (cc.eq.int(hconstraint(j))) then 
!!!!!             do kk=1,sizeh         
!!!!!                if(dd.eq.int(hconstraint(kk))) then
!!!!!                ele34size=ele34size+1
!!!!!                ele34(ele34size)=i
!!!!!                
!!!!!                elseif (bb.eq.int(hconstraint(kk))) then 
!!!!!                ele23size=ele23size+1
!!!!!                ele23(ele23size)=i
!!!!!                endif
!!!!!             enddo
!!!!!          endif
!!!!!          
!!!!!       enddo    
!!!!!      enddo
!!!!!
!!!!!c      do i=1,sizeh
!!!!!c      write(*,*) hconstraint(i)
!!!!!c      enddo
!!!!!c      do i=1,ele12size
!!!!!c      write(*,*) '12',ele12(i)
!!!!!c      enddo
!!!!!c      do i=1,ele23size
!!!!!c      write(*,*) '23',ele23(i)
!!!!!c      enddo
!!!!!c      do i=1,ele34size
!!!!!c      write(*,*) '34', ele34(i)
!!!!!c      enddo
!!!!!c      do i=1,ele41size
!!!!!c      write(*,*) '41', ele41(i)
!!!!!c     enddo
!!!      
!!!         
!!!      return
      end