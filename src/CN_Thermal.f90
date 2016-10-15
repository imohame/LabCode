
! ###########################################################
! ###########################################################
module CN_Thermal
    use CN_BaseClass 
    implicit none
!!!        ---------------------------------------
    type,extends( CN_Base):: ThermalClass
        real*8 , ALLOCATABLE :: rThermal_DijSij(:)

        contains
!!!!!!!!!        procedure ::ThermalClass
        procedure ::CN_Thermal_CleanMem
!!!!!!!!!            Interface CN_Thermal_CleanMem;                Module Procedure CN_Thermal_CleanMem ;          End Interface
    end type ThermalClass
!            Interface init_ThermalClass;                       Module Procedure init_ThermalClass ;            End Interface
!!!        ---------------------------------------
!!!        ---------------------------------------
        contains
! ###########################################################
! ###########################################################
        subroutine CN_Thermal_CleanMem (tCN_object)
            class ( ThermalClass ), intent(inout) :: tCN_object

            call tCN_object%CleanMem()
            IF (ALLOCATED (tCN_object%rThermal_DijSij))           DEALLOCATE (tCN_object%rThermal_DijSij)

        end subroutine CN_Thermal_CleanMem
end module CN_Thermal
! ###########################################################
! ###########################################################
type ( ThermalClass )function init_ThermalClass ( mNodCountB,mElemCountB )
    use CN_Thermal
    integer, intent(in) :: mNodCountB,mElemCountB

        init_ThermalClass%CN_Base=init_CN_Base(mNodCountB,mElemCountB)
        init_ThermalClass%iSolutionActive=1
        Call init_ThermalClass%AllocateMem()

        allocate(init_ThermalClass%rThermal_DijSij(mElemCountB))

end function init_ThermalClass
