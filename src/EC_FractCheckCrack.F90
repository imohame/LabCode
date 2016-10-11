!!subroutine fractCheckFailure(NodesCoordx, NodesCoordy, ElemConnect, DofIds, NodesDispl, ElemMaterial,SolStepCount)
!!
!!    use CN_Objects_manager
!!    use mod_file_units
!!
!!    use mod_parameters
!!!    use EC_ElemCrackingBaseClass
!!    use EC_Objects_manager
!!    implicit none
!!
!!    real NodesCoordx(*), NodesCoordy(*),NodesDispl(*)
!!    INTEGER ElemConnect(4,*),ElemMaterial(*),DofIds(2,*)
!!    INTEGER SolStepCount
!!
!!    common/wblock8/  abc(573,nume,4)
!!    real abc
!!    common/hydroembrittle/critfrac(1000), sigfrac0(nume), sigfrac(nume),decfrac(nume)
!!    real critfrac,sigfrac0, sigfrac,decfrac
!!
!!    INTEGER i,j,iPlaneId
!!    real*8 ElemStress(4),ElemCleavagePlanes(3,3)
!!    real*8 ElemCriticalStress100,ElemCriticalStress110,MaxStress100,StressComp
!!    real*8 vy,vz,s11,s22,s44
!!
!!
!!!--- if no fracture then no need to process this subroutine
!!    if (EC_bCracking == 0) then
!!        return
!!    endif
!!
!!    do i=1, EC_ElemCountInput
!!      !--- check if the elem is cracked and decaying, then go to next elem
!!      if( pEC_ElemData(i)%CheckDecaying(i,SolStepCount)==1) then
!!          cycle !- goto next i
!!      endif
!!      !---- remeber to add check for zone limit
!!      !----------------
!!      !----------------
!!
!!      !-- get the stress at the current elem
!!
!!      CALL CNmanager_Get_sigalt(i,ElemStress)
!!      !- get the critical stress at the current elem
!!      ElemCriticalStress100=sigfrac0(i)		! {100} planes
!!      ElemCriticalStress110=sigfrac(i)         ! {110} planes
!!      !-- get the cleavage planes direction from abc
!!      do j = 1, 3
!!        ElemCleavagePlanes(1,j)=abc(362+j,i,1)
!!        ElemCleavagePlanes(2,j)=abc(365+j,i,1)
!!        ElemCleavagePlanes(3,j)=abc(368+j,i,1)
!!      enddo
!!      !-- obtain maximum normal component of the traction on cleavage planes and corresponding normal vector
!!      MaxStress100 = -100.0
!!      iPlaneId=-1
!!      do j = 1, 3
!!          vy=ElemCleavagePlanes(j,2)
!!          vz=ElemCleavagePlanes(j,3)
!!          s11=ElemStress(1)
!!          s22=ElemStress(2)
!!          s44=ElemStress(4)
!!          StressComp = s11*vy*vy+s22*vz*vz+s44*vy*vz*2.0
!!          !!!dum = sig(1,ele)*cleave(j,2)**2.0+sig(2,ele)*cleave(j,3)**2.0+sig(4,ele)*cleave(j,2)*cleave(j,3)*2.0
!!          if (abs(StressComp)>MaxStress100) then
!!              MaxStress100=abs(StressComp)
!!              iPlaneId=j
!!          end if
!!      end do
!!      !-- estimate failure, change element status flag
!!      if((pEC_ElemData(i)%iElemStatus==0) .and. (MaxStress100 > ElemCriticalStress100)) then
!!          call pEC_ElemData(i)%SetFailed(1,ElemCleavagePlanes(iPlaneId,1:3))
!!          !-- write to crackprog.out
!!          write(iFU_frac_crackFailedInfo,*) i,SolStepCount,MaxStress100,ElemCriticalStress100,ElemCleavagePlanes(iPlaneId,1:3)
!!      end if
!!
!!    enddo !--do i=1, EC_ElemCountInput
!!
!!END
!!!##############################################################################
!!!##############################################################################

subroutine fractCheckCracking(NodesCoordx, NodesCoordy, ElemConnect, DofIds, NodesDispl, ElemMaterial, &
                                usi  , freep , ym,SolStepCount)
                             !-a(k03)    ,a(k04)      ,a(k02)      ,a(k57  ,a(k18)     , a(k08)      ,a(k20), a(k07), a(k09))
    use CN_Objects_manager
    use mod_file_units

    use mod_parameters
    use EC_Objects_manager
    
    implicit none
    common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
    INTEGER nprnt,mprint,itmpop,numelt,jprint,idump,locstr

    real NodesCoordx(*), NodesCoordy(*),NodesDispl(*)
    real usi(*),freep(5,*), ym(4,*)
    INTEGER ElemConnect(4,*),ElemMaterial(*),DofIds(2,*)
    INTEGER SolStepCount

    common/wblock8/  abc(573,nume,4)
    real abc
    common/hydroembrittle/critfrac(1000), sigfrac0(nume), sigfrac(nume),decfrac(nume)
    real critfrac,sigfrac0, sigfrac,decfrac

    INTEGER i,j,iPlaneId
    real*8 ElemStress(4),ElemCleavagePlanes(3,3)
    real*8 ElemCriticalStress100,MaxStress100,StressComp
    real*8 x(4),y(4),Ptc(2)

!--- if no fracture then no need to process this subroutine
    if (EC_bCracking == 0) then
        return
    endif

!--- this part checks for failed element and start unloading them
    do i=1, EC_ElemCountInput
      !--- check if the elem is failed and decaying, then go to next elem
      if( pEC_ElemData(i)%CheckDecaying(i,SolStepCount)==1) then
          cycle !- goto next i
      endif
      !---- remeber to add check for zone limit
      !----------------
      !----------------

      !-- get the stress at the current elem
      CALL CNmanager_Get_sigalt(i,ElemStress)
      !- get the critical stress at the current elem
      ElemCriticalStress100=sigfrac0(i)		! {100} planes
      !-- get the cleavage planes direction from abc
      do j = 1, 3
        ElemCleavagePlanes(1,j)=abc(362+j,i,1)
        ElemCleavagePlanes(2,j)=abc(365+j,i,1)
        ElemCleavagePlanes(3,j)=abc(368+j,i,1)
      enddo
      !-- calc. maximum normal component of the traction on cleavage planes and corresponding normal vector
      call pEC_ElemData(i)%CalcElemMaxStress(ElemCleavagePlanes,ElemStress,MaxStress100,iPlaneId)

      !-- estimate failure, change element status flag
      if((pEC_ElemData(i)%iElemStatus==0) .and. (MaxStress100 > ElemCriticalStress100)) then
          call pEC_ElemData(i)%SetFailed(1,ElemCleavagePlanes(iPlaneId,1:3))
          !-- write to crackprog.out
          write(iFU_frac_crackFailedInfo,*) i,SolStepCount,MaxStress100,ElemCriticalStress100,ElemCleavagePlanes(iPlaneId,1:3)
      end if
    enddo !--do i=1, EC_ElemCountInput

!-- this part checks for the completely unloaded elems and crack/split them
    do i=1, EC_ElemCountInput
        !-if the elem unloaded but the neighbors are not failed then just mark the edges
        !- as EC_eCrackInnerSide
        if(pEC_ElemData(i)%iElemStatus > EC_DecayCount) then
        if(pEC_ElemData(i)%iElemSplit == 0) then
            !- to calc the elem C.G. as the crack starting point
            do j=1, 4
                x(j)=NodesCoordx(ElemConnect(j,i))+NodesDispl(DofIds(1,ElemConnect(j,i)))
                y(j)=NodesCoordy(ElemConnect(j,i))+NodesDispl(DofIds(2,ElemConnect(j,i)))
            enddo
            Ptc(1)=0.25*(x(1)+x(2)+x(3)+x(4))
            Ptc(2)=0.25*(y(1)+y(2)+y(3)+y(4))
            !-- set all edges to EC_eCrackInnerSide and neighbors to EC_eCrackOutterSide
            !-- if neighbors are also EC_eCrackInnerSide, then set this edge to EC_eCrackBothSides
            CALL EC_MarkElemForCrack(i,Ptc,x,y,SolStepCount)
        endif
        endif
    enddo !--do i=1, EC_ElemCountInput
!-- this part creates the overlapping elems, if the elem is failed and unloaded and has any edge EC_eCrackBothSides
    do i=1, EC_ElemCountInput
        if(pEC_ElemData(i)%iElemStatus > EC_DecayCount) then
        if(pEC_ElemData(i)%iElemSplit == 0) then
            CALL EC_SplitElem(i,NodesCoordx, NodesCoordy, ElemConnect, DofIds, NodesDispl, ElemMaterial, &
                            usi  , freep , ym,SolStepCount)

        endif
        endif
    enddo !--do i=1, EC_ElemCountInput
!!!!!!!!================================================================
!!!!!!!!================================================================
!!!!!!!!========= update this global count b/c it's used in the rest of the code
!!!!!!!!================================================================
!!!!!!!!================================================================
!!!!!    numelt=EC_ElemCountCurrent

END subroutine fractCheckCracking
