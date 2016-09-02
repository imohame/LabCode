
! ###########################################################
module EC_Objects_manager
    use EC_ElemCrackingBaseClass
    use EC_Consts
!    implicit none
    public
    ! Declare a set geometric objects.
    type ( EC_ElemCrackingClass ) , ALLOCATABLE :: pEC_ElemData(:)
    integer, ALLOCATABLE :: EC_ElemNeighbors(:,:) !-(8,:) to hold the elem's edges neighboring elem/edges, update when overlapping
!!!        ---------------------------------------

    Interface EC_AllocateMem;              Module Procedure EC_AllocateMem;                   End Interface
    Interface EC_CleanMem;                 Module Procedure EC_CleanMem;                      End Interface
    Interface EC_PrintTest;                 Module Procedure EC_PrintTest;                      End Interface

!##############################################################################
!##############################################################################
    contains
        subroutine EC_AllocateMem()
            use EC_Consts

            implicit none
            integer ne,i

            
            if(EC_bCracking == 0) then
              return
            endif
            ne=int(ElemCountInput*1.25)
            Allocate(pEC_ElemData(ne))
            Allocate(EC_ElemNeighbors(8,ne))
            
            do i=1,ne
                pEC_ElemData(i)%EdgeStatus=-1
                pEC_ElemData(i)%rCleavagePlane=[0,0,1]
                pEC_ElemData(i)%iElemStatus=0
                pEC_ElemData(i)%rCoordRatioCracking=-1
            enddo

         end subroutine EC_AllocateMem
!##############################################################################
!##############################################################################
        subroutine EC_CleanMem()
          use EC_Consts
          implicit none

          if(EC_bCracking == 0) then
            return
          endif
          IF (ALLOCATED (pEC_ElemData))              DEALLOCATE (pEC_ElemData)
          IF (ALLOCATED (EC_ElemNeighbors))              DEALLOCATE (EC_ElemNeighbors)
        end subroutine EC_CleanMem
!##############################################################################
!##############################################################################
        subroutine EC_PrintTest()
          use EC_Consts
          use EC_ElemCrackingBaseClass
          implicit none
          integer i

          if(EC_bCracking == 0) then
            return
          endif

          CALL ECprintTest()

          do i=1,5
              call pEC_ElemData(i)%PrintTest(i)
              write(*,*)'==========================='
              write(*,*)pEC_ElemData(i)
              write(*,*)pEC_ElemData(i)%rCleavagePlane
              write(*,*)pEC_ElemData(i)%EdgeStatus
          enddo
        end subroutine EC_PrintTest
!##############################################################################
!##############################################################################
end module EC_Objects_manager
