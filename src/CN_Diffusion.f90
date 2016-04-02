
! ###########################################################
! ###########################################################
module CN_Diffusion
    use CN_BaseClass
    implicit none
!!!        ---------------------------------------
    type,extends( CN_Base):: DiffusionClass
    end type DiffusionClass
!!!        ---------------------------------------
!!!        ---------------------------------------   
!    Interface init_DiffusionClass;                Module Procedure init_DiffusionClass ;          End Interface
        
!!!        ---------------------------------------
        contains
end module CN_Diffusion
! ###########################################################

type ( DiffusionClass )function init_DiffusionClass ( mNodCountB,mElemCountB )
    use CN_Diffusion
    integer, intent(in) :: mNodCountB,mElemCountB

        init_DiffusionClass%CN_Base=init_CN_Base(mNodCountB,mElemCountB)
        init_DiffusionClass%iSolutionActive=1
        Call init_DiffusionClass%AllocateMem()

end function init_DiffusionClass
      