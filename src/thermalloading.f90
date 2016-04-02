
!     
! File:   thermalloading.f90! Author: imohame!! Created on August 15, 2015, 2:04 AM
!
subroutine ReadThermalLoad()
    use mod_ThermalLd
    
    implicit none
    integer i,IERR
    
    write(*, *) 'beginning of ReadThermalLoad ... reading thermalload.in'

    open(60, file = 'thermalload.in', status = 'old',IOSTAT=IERR, ERR=90)
    

    read (60, *) ThermalLdPtsCount
    if(ThermalLdPtsCount .eq. 0) then
        ThermalLdActive=0
        return
    endif
    allocate(ThermalLdPts(ThermalLdPtsCount,2))
    do i = 1, ThermalLdPtsCount
        read(60, *) ThermalLdPts(i,1),ThermalLdPts(i,2)
    end do

    
    read (60, *) ThermalLdNodesIdsCount
    allocate(ThermalLdNodesIds(ThermalLdNodesIdsCount))
    do i = 1, ThermalLdNodesIdsCount
        read(60, *) ThermalLdNodesIds(i)
    end do
    
    close(60)
!    set the thermal load flag to be active
    ThermalLdActive=1
!    modify the Tini to include the first point in the curve
    call ThermalSetTinitial()
    return
!    -------------------------------------------
!    -------------------------------------------
!     this to handle file not found
!    set the thermal load flag to be active if the file does not exist
    ThermalLdActive=0
    90  IF (IERR .EQ. 29 ) THEN  !FOR$IOS_FILNOTFOU
        write(*, *) '------------------------ no input file .. reading thermalload.in'
      ELSE
        write(*, *) 'Unrecoverable error in .. reading thermalload.in, code =', IERR
        STOP
      END IF

    return
END

!###################################################################
!###################################################################
!###################################################################
! memory clean up
subroutine ThermalLoadCleanMemory()
    
    use mod_ThermalLd
    if(ThermalLdActive .eq.0)  return 
        
    IF (ALLOCATED (ThermalLdPts)) DEALLOCATE (ThermalLdPts)
    IF (ALLOCATED (ThermalLdNodesIds)) DEALLOCATE (ThermalLdNodesIds)
end subroutine ThermalLoadCleanMemory

!###################################################################
!###################################################################
!###################################################################
! thermal curve interpolation
 function ThermalCurveInterpolate(mtime)
    use mod_ThermalLd
    implicit none
    integer i
    real*8 ThermalCurveInterpolate, mtime, mslope
    
!    loop over all the lines to find which one has the time=mtime
    do i = 1, ThermalLdPtsCount-1
        if(mtime >= ThermalLdPts(i,1) .and. mtime <= ThermalLdPts(i+1,1))then
!            get the slope of the line
            mslope=(ThermalLdPts(i+1,2)-ThermalLdPts(i,2))/(ThermalLdPts(i+1,1)-ThermalLdPts(i,1))
!            calc the temp
            ThermalCurveInterpolate = ThermalLdPts(i,2)+(mtime-ThermalLdPts(i,1))*mslope
            return
        end if
    end do

    
 END function ThermalCurveInterpolate
!###################################################################
!###################################################################
!###################################################################
! Sets the initial for the loaded nodes
 subroutine ThermalSetTinitial()
    use mod_ThermalLd
    use mod_parameters
    common/WMLthermal/thermalflag,thermalconstraint(nume),Tinit(nume)
    
    integer i
    real*8  ThermalCurveInterpolate
!!!!    -------------- is not working now, need to be fixed
    return
    
    if(ThermalLdActive.eq.0) return 
!    find the temperature at time=0
    mTemperature=ThermalCurveInterpolate(0.0d0)
    
!    loop over all the nodes and get the temperature at time=0.0
    do i = 1, ThermalLdNodesIdsCount
        Tinit(ThermalLdNodesIds(i))=mTemperature
    end do
    
 END subroutine ThermalSetTinitial
!###################################################################
!###################################################################
!###################################################################
! modifies the RHS2 to account for the thermal load on nodes
 subroutine ThermalSetLoading(nod,RHS2,LHS,Tinitd)
    use mod_ThermalLd
    
    integer i,NodeID
    real*8  ThermalCurveInterpolate,mtime
    
    integer nod
    real*8 RHS2(nod,1),LHS(nod,nod),Tinitd(nod,1)
    common/bk26/dt,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
    common/bk08/kprint,nstep,ite,ilimit,newstf
!!!!    -------------- is not working now, need to be fixed
    return
    if(ThermalLdActive.eq.0)  return 
    
!    find the current temperature of the current time
    mtime=dt*nstep
    mTemperature=ThermalCurveInterpolate(mtime)
!    write(*,*)dt,nstep,mtime,mTemperature
    
!    loop over all the nodes and set the pre-defined temperature from the loading curve
    do i = 1, ThermalLdNodesIdsCount
        NodeID=ThermalLdNodesIds(i)
        LHS(NodeID,1:nod)=0.
        LHS(NodeID,NodeID)=1.0
        RHS2(NodeID,1)=mTemperature
!        Tinitd(NodeID,1)=mTemperature
    end do

    
 END subroutine ThermalSetLoading

!###################################################################
!###################################################################
!###################################################################
! modifies the RHS2 to account for the thermal load on nodes
 subroutine ThermalLdSetToSol(nod,TdSol)
    use mod_ThermalLd
    
    integer i,NodeID
    integer nod
    real*8 TdSol(nod,1)
    
!!!!    -------------- is not working now, need to be fixed
    return
    if(ThermalLdActive.eq.0)  return 
    
!    loop over all the nodes and set the pre-defined temperature from the loading curve
    do i = 1, ThermalLdNodesIdsCount
        NodeID=ThermalLdNodesIds(i)
        TdSol(NodeID,1)=mTemperature
    end do

    
 END subroutine ThermalLdSetToSol

