      subroutine cpynode(m, n, y, z, id, u, usi)
	  use CN_Objects_manager

	  integer m, n
	  dimension y(*), z(*), id(2,*), u(*), usi(*)
	  
!!!!!!!!!	  common/WMLthermal/thermalflag !!!!!!!!!,thermalconstraint(40000),Tinit(40000),Rqold(40000)
!!!!!!!!!      common/WMLthermalSolve/Rqolder(40000),Tinitdold(40000)!,dummy(40000,1)

	  y(n)=y(m)
	  z(n)=z(m)
	  u(id(1,n))=u(id(1,m))
	  u(id(2,n))=u(id(2,m))
	  usi(id(1,n))=usi(id(1,m))
	  usi(id(2,n))=usi(id(2,m))
      
!!!!!!!!      Tinit(n)=Tinit(m)
!!!!!!!!	  Rqold(n)=Rqold(m)
!!!!!!!!	  Rqolder(n)=Rqolder(m)
!!!!!!!!	  Tinitdold(n)=Tinitdold(m)
	  call CNmanager_CopyNode(n,m)
	  end
	  