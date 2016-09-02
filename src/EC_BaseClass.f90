
module EC_Consts
    use mod_parameters
    implicit none

    public

    INTEGER :: EC_eCrackNon=0         !--edge not cracked
    INTEGER :: EC_eCrackFirstSide=1   !--edge cracked on first side
    INTEGER :: EC_eCrackBothSide=2    !--edge cracked on both sides
    INTEGER :: EC_eCrackSecondSide=3  !--edge crack on second side


    integer EC_bCracking !-- 1=allow fracture
    INTEGER EC_ZoneFactor !--- used to mark the crack front circle, do not crack within this area
    real*8  EC_ZoneRadius !--- used to mark the crack front circle, do not crack within this area
    INTEGER EC_DecayCount !--- used to mark the crack front circle, do not crack within this area
    real*8  EC_ElemEdgeMin !--- min elem edge
    real*8  EC_ElemEdgeMax !--- max elem edge
    real*8  EC_ElemAverageDia !--- average elem diagonal

    integer EC_PreCrackedElemCount !- number of pre-cracked elems
    integer EC_NodesCountAdded !-- nodes added as overlap nodes
    integer EC_ElemCountAdded !-- elems added as overlap elems
    integer EC_IntersectNodesCount !-- intersection nodes count

    integer ElemCountInput !-- original elem count from the input before any phantom or overlap elements
    integer ElemCountCurrent !-- current or updated or recent elem count that = original + overlapped

    integer NodeCountInput !-- original node count from the input before any phantom or overlap elements
    integer NodeCountCurrent !-- current or updated or recent node count that = original + overlapped

    integer EC_ElemEdgesConnect(2,4) !-- list the elem edges order 1-2 2-3 3-4 4-1


    Interface ECinitializeConsts;                  Module Procedure ECinitializeConsts;            End Interface
    Interface ECcalcAverageElemStatistics;         Module Procedure ECcalcAverageElemStatistics;   End Interface
    Interface ECprintTest;         Module Procedure ECprintTest;   End Interface
!    ---------------------------
    contains
    subroutine ECprintTest
        write(*,*)'EC_bCracking=',EC_bCracking
        write(*,*)'EC_ZoneFactor=',EC_ZoneFactor
        write(*,*)'EC_ZoneRadius=',EC_ZoneRadius
        write(*,*)'EC_DecayCount=',EC_DecayCount
        write(*,*)'EC_DecayCount=',EC_DecayCount
        write(*,*)'EC_ElemEdgeMin=',EC_ElemEdgeMin
        write(*,*)'EC_ElemEdgeMax=',EC_ElemEdgeMax
        write(*,*)'EC_ElemAverageDia=',EC_ElemAverageDia
        write(*,*)'EC_PreCrackedElemCount=',EC_PreCrackedElemCount
        write(*,*)'EC_NodesCountAdded=',EC_NodesCountAdded
        write(*,*)'EC_ElemCountAdded=',EC_ElemCountAdded
        write(*,*)'ElemCountInput=',ElemCountInput
        write(*,*)'ElemCountCurrent=',ElemCountCurrent
        write(*,*)'NodeCountInput=',NodeCountInput
        write(*,*)'NodeCountCurrent=',NodeCountCurrent
        write(*,*)'EC_ElemEdgesConnect=',EC_ElemEdgesConnect
    end subroutine ECprintTest
!##############################################################################
!##############################################################################
    subroutine ECinitializeConsts(numnpo, numelto,ElemConnect,rNodesCoordx,rNodesCoordy,mEC_bCracking,mEC_ZoneFactor, &
                                mEC_DecayCount,mEC_PreCrackedElemCount)

        implicit none
        !--- set the input node and elem counts before cracking and overlap
        integer,    intent(in) :: ElemConnect(4,*)
        real,     intent(in) :: rNodesCoordx(*),rNodesCoordy(*)
        INTEGER numnpo, numelto
        INTEGER i
        integer mEC_ZoneFactor,mEC_ZoneRadius,mEC_DecayCount,mEC_PreCrackedElemCount,mEC_bCracking
        !--- if no fracture then no need to process this subroutine
        if (mEC_bCracking == 0) then
            return
        endif

        EC_bCracking=mEC_bCracking
        ElemCountInput=numelto
        NodeCountInput=numnpo
        do i=1,3
          EC_ElemEdgesConnect(1,i)=i
          EC_ElemEdgesConnect(2,i)=i+1
        enddo
        EC_ElemEdgesConnect(1,4)=4
        EC_ElemEdgesConnect(2,4)=1

        EC_ZoneFactor=mEC_ZoneFactor
        EC_DecayCount=mEC_DecayCount
        EC_PreCrackedElemCount=mEC_PreCrackedElemCount
        EC_NodesCountAdded=0
        EC_ElemCountAdded=0

        ElemCountCurrent=ElemCountInput
        NodeCountCurrent=NodeCountInput

        call ECcalcAverageElemStatistics (ElemConnect,rNodesCoordx,rNodesCoordy)
        EC_ZoneRadius=EC_ZoneFactor*EC_ElemEdgeMin

    end subroutine ECinitializeConsts
!##############################################################################
!##############################################################################
!!!        --------------------------------------- to cal. some statistics about elem geometry
    subroutine ECcalcAverageElemStatistics (ElemConnect,rNodesCoordx,rNodesCoordy)
!            USE EC_Consts
        implicit none
        integer,    intent(in) :: ElemConnect(4,*)
        real,     intent(in) :: rNodesCoordx(*),rNodesCoordy(*)
        real*8 ::x1,y1,x2,y2, rEdgeMax,rEdgeMin,rElemDia,dist
        integer i,j,edgeNode1,edgeNode2

        rEdgeMax=-1e32
        rEdgeMin=1e32
        rElemDia=0

        do i=1,ElemCountInput
          do j=1,4
            edgeNode1=ElemConnect(EC_ElemEdgesConnect(1,j),i)
            edgeNode2=ElemConnect(EC_ElemEdgesConnect(2,j),i)
            x1=rNodesCoordx(edgeNode1)
            y1=rNodesCoordy(edgeNode1)
            x2=rNodesCoordx(edgeNode2)
            y2=rNodesCoordy(edgeNode2)
            dist=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))
            !-- the max elem edge
            if(dist>rEdgeMax) then
              rEdgeMax=dist
            endif
            !-- the min elem edge
            if((dist<rEdgeMin) .and. (dist>0)) then
              rEdgeMin=dist
            endif
          enddo
          !-- calc the elem diagonal
          edgeNode1=ElemConnect(1,i)
          edgeNode2=ElemConnect(3,i)
          x1=rNodesCoordx(edgeNode1)
          y1=rNodesCoordy(edgeNode1)
          x2=rNodesCoordx(edgeNode2)
          y2=rNodesCoordy(edgeNode2)
          dist=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))
          rElemDia=rElemDia+dist
        enddo
        !-- calc the average diagonal
        rElemDia=rElemDia/ElemCountInput

        !--- set the results to the CN_Consts
        EC_ElemEdgeMax=rEdgeMax
        EC_ElemEdgeMin=rEdgeMin
        EC_ElemAverageDia=rElemDia
    end subroutine ECcalcAverageElemStatistics

end module EC_Consts
! ###########################################################
! ###########################################################
module EC_ElemCrackingBaseClass
    use EC_Consts
    implicit none
    TYPE EC_ElemCrackingClass
        integer :: iElemStatus=0 !- 0=not cracked, 1:EC_DecayCount=decaying, >EC_DecayCount=cracked
        real*8  :: rCleavagePlane(3)
        real*8  :: rCoordRatioCracking(4) !- (-ve no cracking )for each edge, ratio of length for the cracking point
        integer :: EdgeStatus(4) !- for each edge status, see enum_EdgeStatus

     contains
        procedure ::Initialize => EC_ElemCrackingBaseClass_Initialize
        procedure ::PrintTest => EC_ElemCrackingBaseClass_PrintTest
        procedure ::SetDataPre => EC_ElemCrackingBaseClass_SetDataPre
        procedure ::SetFailed => EC_ElemCrackingBaseClass_SetFailed
        procedure ::CalcElemMaxStress => EC_ElemCrackingBaseClass_CalcElemMaxStress
        procedure ::CheckDecaying => EC_ElemCrackingBaseClass_CheckDecaying

    ENDTYPE
    !-------------------
    contains
!!!        --------------------------------------- to allocate memory for the Rn's and CTn's
    subroutine EC_ElemCrackingBaseClass_Initialize (tEC_object)

        class ( EC_ElemCrackingClass ), intent(inout) :: tEC_object

        tEC_object%iElemStatus=0
        tEC_object%rCleavagePlane=[0,0,1]
!            tEC_object%rCleavagePlane(2)=0
!            tEC_object%rCleavagePlane(3)=1
        tEC_object%rCoordRatioCracking=-1
        tEC_object%EdgeStatus=0

    end subroutine EC_ElemCrackingBaseClass_Initialize
!##############################################################################
!##############################################################################
    subroutine EC_ElemCrackingBaseClass_PrintTest (tEC_object,mid)
        class ( EC_ElemCrackingClass ), intent(inout) :: tEC_object
        integer mid

            write(*,*)'================================ elem=',mid
            write(*,*)tEC_object%iElemStatus
            write(*,*)tEC_object%rCleavagePlane(1)
            write(*,*)tEC_object%rCleavagePlane(2)
            write(*,*)tEC_object%rCleavagePlane(3)
            write(*,*)tEC_object%EdgeStatus

    end subroutine EC_ElemCrackingBaseClass_PrintTest
    !##############################################################################
    !##############################################################################
     subroutine EC_ElemCrackingBaseClass_SetDataPre (tEC_object,vn,ElemEdge1,rxc1,ryc1,r1,ElemEdge2,rxc2,ryc2,r2)

       implicit none
         class ( EC_ElemCrackingClass ), intent(inout) :: tEC_object
         real , intent(in) :: vn(3)
         integer , intent(in) ::ElemEdge1,ElemEdge2
         real , intent(in) ::rxc1,ryc1,r1,rxc2,ryc2,r2

          tEC_object%iElemStatus=EC_DecayCount+1
          tEC_object%rCleavagePlane=vn
          tEC_object%EdgeStatus=EC_eCrackFirstSide

     end subroutine EC_ElemCrackingBaseClass_SetDataPre
     !##############################################################################
     !##############################################################################
      subroutine EC_ElemCrackingBaseClass_SetFailed (tEC_object,miElemStatus,vn)

        implicit none
          class ( EC_ElemCrackingClass ), intent(inout) :: tEC_object
          real*8 , intent(in) :: vn(3)
          integer , intent(in) :: miElemStatus

           tEC_object%iElemStatus=miElemStatus
           tEC_object%rCleavagePlane=vn

      end subroutine EC_ElemCrackingBaseClass_SetFailed
      !##############################################################################
      !##############################################################################
      subroutine EC_ElemCrackingBaseClass_CalcElemMaxStress (tEC_object,ElemCleavagePlanes, &
                                                    ElemStress,MaxStress100,iPlaneId)

        implicit none
        class ( EC_ElemCrackingClass ), intent(inout) :: tEC_object
        real*8 , intent(in) :: ElemCleavagePlanes(3,3),ElemStress(4)
        real*8 , intent(out) :: MaxStress100
        integer , intent(out) :: iPlaneId
        real*8 vy,vz,s11,s22,s44,StressComp
        integer j

           !-- obtain maximum normal component of the traction on cleavage planes and corresponding normal vector
        MaxStress100 = -100.0
        iPlaneId=-1
        do j = 1, 3
           vy=ElemCleavagePlanes(j,2)
           vz=ElemCleavagePlanes(j,3)
           s11=ElemStress(1)
           s22=ElemStress(2)
           s44=ElemStress(4)
           StressComp = s11*vy*vy+s22*vz*vz+s44*vy*vz*2.0
           !!!dum = sig(1,ele)*cleave(j,2)**2.0+sig(2,ele)*cleave(j,3)**2.0+sig(4,ele)*cleave(j,2)*cleave(j,3)*2.0
           if (abs(StressComp)>MaxStress100) then
               MaxStress100=abs(StressComp)
               iPlaneId=j
           end if
        end do
       end subroutine EC_ElemCrackingBaseClass_CalcElemMaxStress
      !##############################################################################
      !##############################################################################
      function EC_ElemCrackingBaseClass_CheckDecaying (tEC_object,ele,nstep)result(bDecaying)
        use mod_file_units
        implicit none
        class ( EC_ElemCrackingClass ), intent(inout) :: tEC_object
        integer, intent(in):: ele,nstep
        integer bDecaying
        !-- to check if the elem is failed and unloading
        !- if elem is unloading then increase the unloading steps
        if((tEC_object%iElemStatus <= EC_DecayCount) .and.  (tEC_object%iElemStatus>0)) then
          tEC_object%iElemStatus=tEC_object%iElemStatus+1
          write(iFU_crackprog_out,*) ele,nstep,tEC_object%iElemStatus,EC_DecayCount
          bDecaying= 1
        elseif((tEC_object%iElemStatus > EC_DecayCount) !-for elem that unloaded but not split
          bDecaying= 1
        else
          bDecaying= 0
        endif
      end function EC_ElemCrackingBaseClass_CheckDecaying
!##############################################################################
!##############################################################################


end module EC_ElemCrackingBaseClass
! ###########################################################
! ###########################################################
