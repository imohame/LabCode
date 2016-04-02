
!     
! File:   DiffCoeffTable.f90! Author: imohame!! Created on feb 10, 2016
!
subroutine DiffCoeffTableRead()
    
    use mod_DiffCoeffTable
    use mod_parameters
    
    implicit none
    integer i,j,IERR
    
    write(*, *) 'beginning of ReadDiffCoeffTable ... reading DiffCoeffTable.in'

    open(60, file = 'DiffCoeffTable.in', status = 'old',IOSTAT=IERR, ERR=90)
    
!--- read the number of tables
    read (60, *) iDiffCoeffTableCount
    if(iDiffCoeffTableCount .eq. 0) then
        iDiffCoeffTableActive=0
        return
    endif
    
    allocate(tDiffCoeffTables(iDiffCoeffTableCount))

!---- for each table read the number of pts then the table
    do i = 1, iDiffCoeffTableCount
        !--read the number of pts in this table
        read (60, *) tDiffCoeffTables(i)%iPtscount
        allocate(tDiffCoeffTables(i)%rPts(tDiffCoeffTables(i)%iPtscount,2))
        !--read these points
        do j = 1, tDiffCoeffTables(i)%iPtscount
            read(60, *) tDiffCoeffTables(i)%rPts(j,1),tDiffCoeffTables(i)%rPts(j,2)
        end do
        !--read the number of mat ids in this table
        read (60, *) tDiffCoeffTables(i)%iMatcount
        allocate(tDiffCoeffTables(i)%iMatIds(tDiffCoeffTables(i)%iMatcount))
        !--read these ids
        do j = 1, tDiffCoeffTables(i)%iMatcount
            read(60, *) tDiffCoeffTables(i)%iMatIds(j)
        end do
        
    end do
    
    close(60)
!    set the thermal load flag to be active
    iDiffCoeffTableActive=1
    
    
    return
!    -------------------------------------------
!    -------------------------------------------
!     this to handle file not found
!    set the thermal load flag to be active if the file does not exist
    iDiffCoeffTableActive=0
    90  IF (IERR .EQ. 29 ) THEN  !FOR$IOS_FILNOTFOU
        write(*, *) '------------------------ no input file .. reading DiffCoeffTable.in'
      ELSE
        write(*, *) 'Unrecoverable error in .. reading DiffCoeffTable.in, code =', IERR
        STOP
      END IF

    return
END

!###################################################################
!###################################################################
!###################################################################
! memory clean up
subroutine DiffCoeffTableCleanMemory()
    
    use mod_DiffCoeffTable
    
    if(iDiffCoeffTableActive .eq.0)  return 
    
    do i = 1, iDiffCoeffTableCount
        IF (ALLOCATED (tDiffCoeffTables(i)%rPts))      DEALLOCATE (tDiffCoeffTables(i)%rPts)
        IF (ALLOCATED (tDiffCoeffTables(i)%iMatIds))      DEALLOCATE (tDiffCoeffTables(i)%iMatIds)
    end do        
    IF (ALLOCATED (tDiffCoeffTables))      DEALLOCATE (tDiffCoeffTables)
    
end subroutine DiffCoeffTableCleanMemory
!###################################################################
!###################################################################
!###################################################################
 subroutine DiffCoeffTableGetD(elemId,rDiffCoeff,rTemperature,matp)
    use mod_DiffCoeffTable
    use mod_parameters
    use CN_Objects_manager
    
    integer i,j,elemId,elemIdMat
    real*8  rDiffCoeff,rTemperature,DiffCoeffTableInterpolate
    dimension  matp(*)
    
    
    if(iDiffCoeffTableActive == 0)  then
!    -- there is no tables, use the mat input data 
        rDiffCoeff=rMatDiff_D(elemId)
!!!!!!        rDiffCoeff=thermalk(elemId)
!        write(*,*) 'iDiffCoeffTableActive, rDiffCoeff=thermalk(elemId),rTemperature',iDiffCoeffTableActive,rDiffCoeff,rTemperature
        return
    endif 
    
    elemIdMat=int(matp(elemId))
!    write(*,*) elemIdMat 
    do i = 1, iDiffCoeffTableCount
        !--read these ids
        do j = 1, tDiffCoeffTables(i)%iMatcount
!            -- check if the elem mat is in one of the tables
!            write(*,*) elemIdMat, tDiffCoeffTables(i)%iMatIds(j),rTemperature
            if (elemIdMat == tDiffCoeffTables(i)%iMatIds(j)) then
                rDiffCoeff=DiffCoeffTableInterpolate(i,rTemperature)
!                write(*,*) elemId,elemIdMat,rDiffCoeff
!                write(*,*) 'iDiffCoeffTableActive, rDiffCoeff=DiffCoeffTableInterpolate,rTemperature',iDiffCoeffTableActive,rDiffCoeff,rTemperature
                return
            end if
        end do
    end do

 END subroutine DiffCoeffTableGetD
!###################################################################
!###################################################################
!###################################################################
! thermal curve interpolation
 function DiffCoeffTableInterpolate(TableId,rTemperature)
    use mod_DiffCoeffTable
    implicit none
    integer i,TableId
    real*8 DiffCoeffTableInterpolate, rTemperature,mslope,t1,t2,d1,d2
    
        if (rTemperature < tDiffCoeffTables(TableId)%rPts(1,1)) then
            DiffCoeffTableInterpolate=tDiffCoeffTables(TableId)%rPts(1,2)
            return
        endif
        if (rTemperature > tDiffCoeffTables(TableId)%rPts(tDiffCoeffTables(TableId)%iPtscount,1)) then
            DiffCoeffTableInterpolate=tDiffCoeffTables(TableId)%rPts(tDiffCoeffTables(TableId)%iPtscount,2)
            return
        endif
    
!    loop over all the lines to find which one has the time=mtime
    do i = 1, tDiffCoeffTables(TableId)%iPtscount-1
        t1=tDiffCoeffTables(TableId)%rPts(i,1)
        t2=tDiffCoeffTables(TableId)%rPts(i+1,1)
        d1=tDiffCoeffTables(TableId)%rPts(i,2)
        d2=tDiffCoeffTables(TableId)%rPts(i+1,2)
        if(rTemperature >= t1 .and. rTemperature <= t2)then
!            get the slope of the line
            mslope=(d2-d1)/(t2-t1)
!            calc the temp
            DiffCoeffTableInterpolate = d1+(rTemperature-t1)*mslope
            return
        endif
    end do
    
    return
    
 END function DiffCoeffTableInterpolate