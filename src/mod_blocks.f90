
!###########################################################
!###########################################################
module mod_parameters
!      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
!      IMPLICIT INTEGER(I-N)
      integer, parameter :: nelemg     =128 !512
	  integer, parameter :: nume       = 40000
      integer, parameter :: no_mat     = 1000
      integer, parameter :: nss        = 24
      integer, parameter :: print_intv = 10000
      real, parameter  :: pi=3.14159265358979 !4.d0*atan(1.d0) ! a circle constant
      real, parameter  :: Tempr=293.0d0
!!!!!!!      real   , parameter :: xi         = 0.3 ! it was 0.3
!!!!      integer ElemCountB,NodCountB

      save
end module mod_parameters

!#####################################################
module mod_DiffCoeffTable
    TYPE type_DiffCoeffTable
        integer iPtscount,iMatcount
        real*8, ALLOCATABLE :: rPts(:,:)
        INTEGER, ALLOCATABLE :: iMatIds(:)
    ENDTYPE
    integer iDiffCoeffTableCount,iDiffCoeffTableActive
    TYPE (type_DiffCoeffTable),ALLOCATABLE :: tDiffCoeffTables(:)
!!!!    integer iDiffCoeffTableCount,iDiffCoeffTableTotalptsCount,iDiffCoeffTableMatIDsCount,iDiffCoeffTableActive
!!!!    real*8, ALLOCATABLE :: rDiffCoeffTablepts(:,:)
!!!!    INTEGER, ALLOCATABLE :: iDiffCoeffTableptsCount(:), iDiffCoeffTableMatCount(:),iDiffCoeffTableMatIDs(:)
    real*8 rDiffCoeffDefault
    save

end module mod_DiffCoeffTable
!#####################################################
module mod_dtimeSpecs
    integer idtimeSpecsCount,idtimeSpecsActive
    real*8, ALLOCATABLE :: dtimeValues(:)
    INTEGER, ALLOCATABLE :: idtimeStepsSol(:),idtimeStepsOutput(:)

    save

end module mod_dtimeSpecs
!#####################################################
!#####################################################
module mod_ThermalLd
    integer ThermalLdNodesIdsCount,ThermalLdPtsCount,ThermalLdActive
    real*8, ALLOCATABLE :: ThermalLdPts(:,:)
    INTEGER, ALLOCATABLE :: ThermalLdNodesIds(:)
    real*8 mTemperature
    save

end module mod_ThermalLd
!#####################################################
!#####################################################
module mod_GBdata
    integer GBAreaCount,GBElemCount
    real*8, ALLOCATABLE :: GBNormals(:,:)
    INTEGER, ALLOCATABLE :: GBElemAreaMap(:,:)

    save

end module mod_GBdata
!!!#####################################################
!!!#####################################################
!!module mod_Fract
!!    ENUM (KIND=1):: enum_EdgeStatus
!!          ENUMERATOR :: CrackNon=0         !--edge not cracked
!!          ENUMERATOR :: CrackFirstSide=1   !--edge cracked on first side
!!          ENUMERATOR :: CrackBothSide=2    !--edge cracked on both sides
!!          ENUMERATOR :: CrackSecondSide=3  !--edge crack on second side
!!    END ENUM
!!    
!!    integer FractElemCountPreCracked !- number of pre-cracked elems
!!    integer bFractFlag !-- 1=allow fracture
!!    integer FractNodesCountAdded !-- nodes added as overlap nodes
!!    integer FractElemCountAdded !-- elems added as overlap elems
!!    integer FractIntersectNodesCount !-- intersection nodes count
!!    integer ElemCountInput !-- original elem count from the input before any phantom or overlap elements
!!    integer ElemCountCurrent !-- current or updated or recent elem count that = original + overlapped
!!    integer NodeCountInput !-- original node count from the input before any phantom or overlap elements
!!    integer NodeCountCurrent !-- current or updated or recent node count that = original + overlapped
!!    
!!    INTEGER FractZoneFactor !--- used to mark the crack front circle, do not crack within this area
!!    INTEGER FractDecayCount !--- used to mark the crack front circle, do not crack within this area
!!!!!    INTEGER, ALLOCATABLE :: FractElemIds(:,:) !--pre-cracked
!!    INTEGER, ALLOCATABLE :: FractElemNeighbors(:,:) !(8,nume) each edge has a neighbor elem and edge; e1:e4,ed1:ed4
!!    
!!    TYPE (enum_EdgeStatus) ,ALLOCATABLE :: FractElemEdgeStatus(:,:) !(4,nume) for each edge status
!!    real*8, ALLOCATABLE ::  FractNodesCoordCracking(:,:) !-- (4,:) for each edge, ratio of length for the cracking point
!!    real*8, ALLOCATABLE ::  FractElemCleavagePlanes(:,:) !-- (3,:) x,y,z elem cleavage plane
!!    INTEGER, ALLOCATABLE :: FractElemStatus(:) !(nume) 0=not cracked, 1:FractDecayCount=decaying, >FractDecayCount=cracked
!!    
!!    save
!!
!!end module mod_Fract
!#####################################################
!#####################################################
   !------ module for the files in/out unit numbers
module mod_file_units
!    integer :: ifile_unit_mazout           =18
!    integer :: ifile_unit_mazout_mirror    =12
!    integer :: ifile_unit_result           =17  !-- lfnt(2)
!    integer :: ifile_unit_order            =777 !write(777
!    integer :: ifile_unit_order_flag       =1 !0= do not print; 1= print the order
!    integer :: ifile_unit_forout           =126 !126,file='forout.plt'
!    integer :: ifile_unit_ge               =29 !29,file='ge.f'
!    integer, parameter :: LCHARS = 80
    integer :: iFU_solsteps_out             =7015
    integer :: iFU_times_out                =7016
    integer :: iFU_meshcount_out            =7017
    integer :: iFU_meshextra_out            =7018
    integer :: iFU_meshElems_out            =7019
    integer :: iFU_meshNodes_out            =7020
    integer :: iFU_abcDumpBIN_out           =7021
    integer :: iFU_temperNodes_out          =7022
    integer :: iFU_temperElem_out           =7023
    integer :: iFU_DijSijElem_out           =7024
    integer :: iFU_DiffNodes_out            =7025
    integer :: iFU_DiffElem_out             =7026
    integer :: iFU_DiffCoeff_out            =7027
    integer :: iFU_temperElemID_out         =55555
    integer :: iFU_temperNodeID_out         =66666
    !--- fracture files
    integer :: iFU_crackprog_out            =989
    integer :: iFU_crack_out                =979
    integer :: iFU_cracktip_out             =5001
    !--- input checking files
    integer :: iFU_check_transform_out      =66
    save
end module mod_file_units
!!!!!!!#####################################################
!!!!!!!#####################################################
!!!!!!module mod_CrankNicolson
!!!!!!    IMPLICIT NONE
!!!!!!!    public :: type_CrankNicolson_data
!!!!!!    TYPE type_CrankNicolson_data
!!!!!!        real*8 C1,C2,C3,C4,C5
!!!!!!        real*8 Vh,Tempr,R,Vh_RT
!!!!!!        integer maxint
!!!!!!        real*8 psi(4),eta(4)
!!!!!!        integer iSolutionActive
!!!!!!!        rRn=Rqg .... rRn_1=Rqo    .... rRn_2=Rqor ....
!!!!!!!        rTn=Tinit .. rTn_1=Tinitd .... rTn_2=Tinitdold ....
!!!!!!        real*8 , ALLOCATABLE :: rRn(:),rRn_1(:),rRn_2(:)
!!!!!!        real*8 , ALLOCATABLE :: rTn(:),rTn_1(:),rTn_2(:)
!!!!!!        real*8 , ALLOCATABLE :: rConcentration(:)
!!!!!!
!!!!!!        contains
!!!!!!             procedure,NOPASS :: initialize
!!!!!!    ENDTYPE
!!!!!!
!!!!!!    TYPE (type_CrankNicolson_data) tCNdataDiff,tCNdataThermal
!!!!!!    real*8 , ALLOCATABLE :: thermalk(:),thermalh(:),thermalRo(:),thermalcp(:),thermalx(:),thermalD(:)
!!!!!!
!!!!!!    contains
!!!!!!    subroutine initialize(this)
!!!!!!        TYPE (type_CrankNicolson_data) this
!!!!!!        real*8 alpha,beta
!!!!!!        real , parameter :: sqroot=0.577350269189626
!!!!!!        !----      at psi,eta=sqrt(3)/3, w=1 and maxint=4 and going ccw
!!!!!!        alpha=-.25
!!!!!!        beta=0.75
!!!!!!
!!!!!!        this%c1=beta*(1+alpha)
!!!!!!        this%c2=beta*alpha
!!!!!!        this%c3=(1-beta)*(1+alpha)
!!!!!!        this%c4=(1-beta)*alpha
!!!!!!        this%c5=(1+alpha)*alpha
!!!!!!
!!!!!!        this%Vh=2.0E-6           ! partial molar volume of hydrogen
!!!!!!        this%Tempr=293.0
!!!!!!        this%R=8.3142
!!!!!!        this%Vh_RT=this%Vh/(this%R*this%Tempr)
!!!!!!
!!!!!!        this%maxint=4
!!!!!!
!!!!!!        this%psi(1:4)=(/-sqroot, sqroot, sqroot, -sqroot/)
!!!!!!        this%eta(1:4)=(/-sqroot, -sqroot, sqroot, sqroot/)
!!!!!!        this%rRn=0.0
!!!!!!        this%rRn_1=0.0
!!!!!!        this%rRn_2=0.0
!!!!!!        this%rTn=0.0
!!!!!!        this%rTn_1=0.0
!!!!!!        this%rTn_2=0.0
!!!!!!        this%rConcentration=0.0
!!!!!!
!!!!!!    end subroutine initialize
!!!!!!
!!!!!!
!!!!!!end module mod_CrankNicolson

!#####################################################
!!!###########################################################
!module mod_main_block
!    !error #5520: A common block or variable may not exceed 2,147,483,647 bytes
!    !common /main_block/ a(0:400000000) ..so max= 2,147,483,647/8 = 200,000,000
!    real a(400000000)
!    save
!end module mod_main_block
!!!###########################################################
!!###########################################################
!  module mod_blocks_bk00
!
!  common/bk00/ioff(96)
!     common/bk00/ &
!     k01,k02,k03,k04,k05,k06,k07,k08,k09,k10,k11,k12, &
!     k13,k14,k15,k16,k17,k18,k19,k20,k21,k22,k23,k24, &
!     k25,k26,k27,k28,k29,k30,k31,k32,k33,k34,k35,k36, &
!     k37,k38,k39,k40,k41,k42,k43,k44,k45,k46,k47,k48, &
!     k49,k50,k51,k52,k53,k54,k55,k56,k57,k58,k59,k60, &
!     k61,k62,k63,k64,k65,k66,k67,k68,k69,k70,k71,k72, &
!     k73,k74,k75,k76,k77,k78,k79,k80,k81,k82,k83,k84, &
!     k85,k86,k87,k88,k89,k90,k91,k92,k93,k94,k95,k96
!
!    save
!  end module mod_blocks_bk00
!
