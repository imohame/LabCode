      subroutine initial_values(properties,numelt)

!!!      integer, parameter :: nss    = 24
!!!      integer, parameter :: nume   = 40000
!!!      integer, parameter :: no_mat = 1000
      use  mod_parameters
      use EC_Objects_manager
      use CN_Objects_manager
      use mod_file_units

      integer :: numelt
      integer :: i, j, ig
	  common/bk00/ &
            k01,k02,k03,k04,k05,k06,k07,k08,k09,k10,k11,k12, &
            k13,k14,k15,k16,k17,k18,k19,k20,k21,k22,k23,k24, &
            k25,k26,k27,k28,k29,k30,k31,k32,k33,k34,k35,k36, &
            k37,k38,k39,k40,k41,k42,k43,k44,k45,k46,k47,k48, &
            k49,k50,k51,k52,k53,k54,k55,k56,k57,k58,k59,k60, &
            k61,k62,k63,k64,k65,k66,k67,k68,k69,k70,k71,k72, &
            k73,k74,k75,k76,k77,k78,k79,k80,k81,k82,k83,k84, &
            k85,k86,k87,k88,k89,k90,k91,k92,k93,k94,k95,k96
	  common /main_block/ a(1)
      dimension :: properties(48,*)
      common /wblock8/ abc(573,nume,4),his(573,nume,4)
      common /wblock9/ slip_n0_t(no_mat,nss,3), slip_s0_t(no_mat,nss,3)
      common /wblock10/ng,grain_mo(1000,3),bv(no_mat,87),nssmat(1000),nssimmat(1000)
      common /wblock11/ pd_counter,rhoim0(1000,nss),rhomo0(1000,nss)
      common /wblock12/ Y_modulus(nume),possion_ratio(nume),tau_y(nume)
      common /wblock20/ mat_type(nume)

      common/wblock3/  density_ms, density_ims,thermalEnthalpy(1000)
!!!!!!      common/WMLthermal/thermalflag !!!!!!,thermalconstraint(nume),Tinit(nume),Rqold(nume)
      common/WMLthermal3/thermalki(1000),thermalhi(1000), &
           thermalRoi(1000),thermalcpi(1000),thermalxi(1000),thermalDi(1000)
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
	  common/hydroembrittle/critfrac(1000), sigfrac0(40000),sigfrac(40000),decfrac(40000)
	  common/hydroembrittle110/critfrac110(1000), sigfrac0110(40000)

!!!!!!      integer thermalflag
	  real cleave_n, cleave_n0, critfrac, sigfrac0, sigfrac, decfrac
	  real sigfrac0110

!      --- call this to allocate the needed arrays for CN solver
        CALL CNSetMatProp(nume,mat_type,1000,thermalki,thermalhi,thermalRoi,thermalcpi,thermalxi,thermalDi)


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
    end do
!!!---------  for debugging                  
!!    write(iFU_check_transform_out,*) '---- elements cleavage planes before GB and crack.in   '
!!    do i = 1, numelt
!!        write(iFU_check_transform_out,*) '---- ele id= ',i
!!        do j = 1, 3
!!            write(iFU_check_transform_out,*)abc(362+j,i,1),abc(365+j,i,1),abc(368+j,i,1)
!!        end do
!!    end do      
!!!---------  for debugging          

! #####################################################
! #####################################################
!! -- call the GB functions to apply any GB cleave planes to the elements
!! -- this is has to be done after applying the cleave for all elements
      call GBApplyCleavagePlanes()
! #####################################################
! #####################################################
!      this is has to be called first b/c the pre-crack needs it
       call FractReadElemNeighbors()
!      to read precracked elements and apply proper cleavage planes
       call FractReadApplyPreCrackCleavagePlanes(a(k03),a(k04),a(k02),a(k08),a(k57),a(k18),a(k20), a(k07), a(k09)) 
                                              !!(y     ,z     ,ix    , matp , id   , u    , usi  , freep , ym)
       
!!!       call EC_PrintTest()
       
! #####################################################
! #####################################################
!!!!---------  for debugging                  
!!!    write(iFU_check_transform_out,*) '---- elements cleavage planes after GB and crack.in   '
!!!    do i = 1, numelt
!!!        write(iFU_check_transform_out,*) '---- ele id= ',i
!!!        do j = 1, 3
!!!            write(iFU_check_transform_out,*)abc(362+j,i,1),abc(365+j,i,1),abc(368+j,i,1)
!!!        end do
!!!    end do      
!!!!---------  for debugging          

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
       
       
    close(iFU_check_transform_out)
    write(*, *) '-->>>>>>>>initial_values --done'
!!    stop
end
