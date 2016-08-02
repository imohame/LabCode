      subroutine initial_values(properties,numelt)

!!!      integer, parameter :: nss    = 24
!!!      integer, parameter :: nume   = 40000
!!!      integer, parameter :: no_mat = 1000
      use  mod_parameters
      use CN_Objects_manager

      integer :: numelt
      integer :: i, j, ig
      dimension :: properties(48,*)
      common /wblock8/ abc(573,nume,4),his(573,nume,4)
      common /wblock9/ slip_n0_t(no_mat,nss,3), slip_s0_t(no_mat,nss,3)
      common /wblock10/ng,grain_mo(1000,3),bv(no_mat,87),nssmat(1000)
	1       ,nssimmat(1000)
      common /wblock11/ pd_counter,rhoim0(1000,nss),rhomo0(1000,nss)
      common /wblock12/ Y_modulus(nume),possion_ratio(nume),tau_y(nume)
      common /wblock20/ mat_type(nume)

      common/wblock3/  density_ms, density_ims,thermalEnthalpy(1000)
!!!!!!      common/WMLthermal/thermalflag !!!!!!,thermalconstraint(nume),Tinit(nume),Rqold(nume)
      common/WMLthermal3/thermalki(1000),thermalhi(1000),
     >      thermalRoi(1000),thermalcpi(1000),thermalxi(1000)
     >      ,thermalDi(1000)
!!!!!      common /WMLthermal2/thermalk(nume),thermalh(nume),
!!!!!     >      thermalRo(nume),thermalcp(nume),thermalx(nume)
!!!!!     >      ,thermalD(nume)
!!!
!!!!      common/wblock3/  density_ms, density_ims,etain(1000),ecin(1000)
!!!!      common /WMLthermal/thermalflag,thermalconstraint(40000),
!!!!     1     Tinit(nume),      Tfl(nume)
!!!!      common /WMLthermal2/thermalk(nume),thermalh(nume),etae(nume)
!!!!      common /WMLthermal3/thermalki(1000),thermalhi(1000)

	  common /cleavage_plane/ cleave_n(1000,3,3)
	  common /cleavage_plane0/ cleave_n0(3,3)
	  common/hydroembrittle/critfrac(1000), sigfrac0(40000),
     >       sigfrac(40000),decfrac(40000)
	  common/hydroembrittle110/critfrac110(1000), sigfrac0110(40000)

!!!!!!      integer thermalflag
	  real cleave_n, cleave_n0, critfrac, sigfrac0, sigfrac, decfrac
	  real sigfrac0110

!      --- call this to allocate the needed arrays for CN solver
        CALL CNSetMatProp(nume,mat_type,1000,thermalki,thermalhi,
     >                    thermalRoi,thermalcpi,thermalxi,thermalDi)


      do i = 1, numelt
       ig = mat_type(i)         ! assign material properties to elements
       Y_modulus(i)= properties(1,ig)
       possion_ratio(i) = properties(2,ig)
       tau_y(i) = properties(3,ig)
	   sigfrac0(i)=critfrac(ig)    ! assign critical fracture stress to elements
	   sigfrac0110(i)=critfrac110(ig)

!!------------------------------------------ these prop are saved in CN manager
!!!!!!!           thermalk(i)=thermalki(ig)
!!!!!!!           thermalh(i)=thermalhi(ig)
!!!!!!!           thermalRo(i)=thermalRoi(ig)
!!!!!!!           thermalcp(i)=thermalcpi(ig)
!!!!!!!           thermalx(i)=thermalxi(ig)
!!!!!!!           thermalD(i)=thermalDi(ig)
!!------------------------------------------ these prop are saved in CN manager

	   do j = 1, 20
		   abc(397+j,i,1) = 0.0
	   end do

	   abc(362,i,1) = 0.0

	   do j = 1, 3
		   abc(362+j,i,1) = cleave_n(ig,1,j)
		   abc(365+j,i,1) = cleave_n(ig,2,j)
		   abc(368+j,i,1) = cleave_n(ig,3,j)
	   end do

       do j = 1, nssmat(ig)
        j1=57+j
         abc(j1,i,1)= rhomo0(ig,j)
        j2=81+j
         abc(j2,i,1)= rhoim0(ig,j)
        j3=105+j
         abc(j3,i,1)=slip_n0_t(ig,j,1)
        j4=129+j
         abc(j4,i,1)=slip_n0_t(ig,j,2)
        j5 = 153 + j
         abc(j5,i,1) = slip_n0_t(ig,j,3)
        j6=177 + j
         abc(j6,i,1)=slip_s0_t(ig,j,1)
        j7=201+j
         abc(j7,i,1)=slip_s0_t(ig,j,2)
        j8=225+j
         abc(j8,i,1)=slip_s0_t(ig,j,3)
		 j9 = 337+j
		 abc(j9,i,1) = 0.0
       end do
!c------------------------
      end do
! #####################################################
! #####################################################
!! -- call the GB functions to apply any GB cleave planes to the elements
!! -- this is has to be done after applying the cleave for all elements
      call GBApplyCleavagePlanes()
! #####################################################
! #####################################################
!      to read precracked elements and apply proper cleavage planes
       call FractReadApplyPreCrackCleavagePlanes()!a(k03),a(k04),a(k02))
! #####################################################
! #####################################################
        write(*, *) '-->>>>>>>>initial_values --done'

!!!!           write(*,*)Y_modulus(1:numelt)
!!!!           write(*,*)possion_ratio(1:numelt)
!!!!           write(*,*)tau_y(1:numelt)
!!!!           write(*,*)thermalk(1:numelt)
!!!!           write(*,*)thermalh(1:numelt)
!!!!           write(*,*)thermalRo(1:numelt)
!!!!           write(*,*)thermalcp(1:numelt)
!!!!           write(*,*)thermalKa(1:numelt)

	  !write(*,*) 'in initial_values' !WMLWRITE51210
	  !write(*,*) nssmat(1)           !WMLWRITE51210
      return
      end
