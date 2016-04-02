
!     
! File:   dtimesteps.f90! Author: imohame!! Created on Dec 10, 2015, 2:04 pM
!
subroutine ReaddtimestepsSpecs()
    use mod_dtimeSpecs
    
    implicit none
    integer i,IERR
    
    IERR=0
    i=0
    idtimeSpecsCount=0
    idtimeSpecsActive=0
    write(*, *) 'beginning of dtimesteps ... reading dtimesteps.in'

    open(60, file = 'dtimesteps.in', status = 'old',IOSTAT=IERR, ERR=90)
    
    idtimeSpecsActive=1
    write(*, *) '.......................... idtimeSpecsActive=1'
    read (60, *) idtimeSpecsCount
    if(idtimeSpecsCount == 0) then
        idtimeSpecsActive=0
        write(*, *) '.......................... idtimeSpecsActive=0'
        return
    endif
 
    allocate(dtimeValues(idtimeSpecsCount),idtimeStepsSol(idtimeSpecsCount),idtimeStepsOutput(idtimeSpecsCount))
    do i = 1, idtimeSpecsCount
        read(60, *) dtimeValues(i),idtimeStepsSol(i),idtimeStepsOutput(i)
    end do
    
        write(*, *) '------------------------ idtimeSpecsCount',idtimeSpecsCount
    do i = 1, idtimeSpecsCount
        write(*, *) '------------------------ dtimeValues(i),idtimeStepsSol(i),idtimeStepsOutput(i)', & 
        dtimeValues(i),idtimeStepsSol(i),idtimeStepsOutput(i)
    end do

    close(60)

    return
!    -------------------------------------------
!    -------------------------------------------
    idtimeSpecsActive=0
    90  IF (IERR .EQ. 29 ) THEN  
        write(*, *) '------------------------ no input file .. reading dtimesteps.in'
      ELSE
        write(*, *) 'Unrecoverable error in .. reading dtimesteps.in, code =', IERR
        STOP
      END IF

    return
END

!###################################################################
!###################################################################
!###################################################################
! memory clean up
subroutine dtimeStepsCleanMemory()
    
    use mod_dtimeSpecs
    if(idtimesSpecsActive == 0)  return 
        
    IF (ALLOCATED (dtimeValues))       DEALLOCATE (dtimeValues)
    IF (ALLOCATED (idtimeStepsSol))     DEALLOCATE (idtimeStepsSol)
    IF (ALLOCATED (idtimeStepsOutput))  DEALLOCATE (idtimeStepsOutput)
    
end subroutine dtimeStepsCleanMemory
