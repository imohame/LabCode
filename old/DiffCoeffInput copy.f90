
!     
! File:   DiffCoeffTable.f90! Author: imohame!! Created on feb 10, 2016
!
subroutine DiffCoeffTableRead()
    
    use mod_DiffCoeffTable
    
    implicit none
    integer i,j,IERR,counterpts,counterMat
    
    write(*, *) 'beginning of ReadDiffCoeffTable ... reading DiffCoeffTable.in'

    open(60, file = 'DiffCoeffTable.in', status = 'old',IOSTAT=IERR, ERR=90)
    
!--- read the number of tables
    read (60, *) iDiffCoeffTableCount
    if(iDiffCoeffTableCount .eq. 0) then
        iDiffCoeffTableActive=0
        return
    endif
    allocate(iDiffCoeffTableptsCount(iDiffCoeffTableCount))
    allocate(iDiffCoeffTableMatCount(iDiffCoeffTableCount))
!--- read the total number of points in all tables
    read (60, *) iDiffCoeffTableTotalptsCount
    allocate(rDiffCoeffTablepts(iDiffCoeffTableTotalptsCount,2))
!--- read the total number of mat ids in all tables
    read (60, *) iDiffCoeffTableMatIDsCount
    allocate(iDiffCoeffTableMatIDs(iDiffCoeffTableMatIDsCount))

!---- for each table read the number of pts then the table
    counterpts=0
    counterMat=0
    do i = 1, iDiffCoeffTableCount
        !--read the number of pts in this table
        read (60, *) iDiffCoeffTableptsCount(i)
!        --- accumulate the counts
        iDiffCoeffTableptsCount(i)=iDiffCoeffTableptsCount(i)+counterpts
        !--read these points
        do j = counterpts+1, iDiffCoeffTableptsCount(i)
            read(60, *) rDiffCoeffTablepts(j,1),rDiffCoeffTablepts(j,2)
        end do
        !--read the number of mat ids in this table
        read (60, *) iDiffCoeffTableMatCount(i)
!        --- accumulate the counts
        iDiffCoeffTableMatCount(i)=iDiffCoeffTableMatCount(i)+counterMat
        !--read these ids
        do j = counterMat+1, iDiffCoeffTableMatCount(i)
            read(60, *) iDiffCoeffTableMatIDs(j)
        end do
        
!        --- save the old ones for the next iteration
        counterpts=iDiffCoeffTableptsCount(i)
        counterMat=iDiffCoeffTableMatCount(i)
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
        
    IF (ALLOCATED (rDiffCoeffTablepts))      DEALLOCATE (rDiffCoeffTablepts)
    IF (ALLOCATED (iDiffCoeffTableptsCount)) DEALLOCATE (iDiffCoeffTableptsCount)
    IF (ALLOCATED (iDiffCoeffTableMatCount)) DEALLOCATE (iDiffCoeffTableMatCount)
    IF (ALLOCATED (iDiffCoeffTableMatIDs))   DEALLOCATE (iDiffCoeffTableMatIDs)
    
end subroutine DiffCoeffTableCleanMemory


!###################################################################
!###################################################################
!###################################################################
 subroutine DiffCoeffTableGetD(elemId,rDiffCoeff,rTemperature)
    use mod_DiffCoeffTable
    use  mod_parameters
    
    integer i,elemId,counterMat,iMatElemIdStart
    real*8  rDiffCoeff,rTemperature,DiffCoeffTableInterpolate
    
    common /WMLthermal2/thermalk(nume),thermalh(nume),thermalRo(nume) ,thermalcp(nume),thermalKa(nume)
    common /main_block/ a(1)
    common/bk00/ &
          k01,k02,k03,k04,k05,k06,k07,k08,k09,k10,k11,k12, &
          k13,k14,k15,k16,k17,k18,k19,k20,k21,k22,k23,k24, &
          k25,k26,k27,k28,k29,k30,k31,k32,k33,k34,k35,k36, &
          k37,k38,k39,k40,k41,k42,k43,k44,k45,k46,k47,k48, &
          k49,k50,k51,k52,k53,k54,k55,k56,k57,k58,k59,k60, &
          k61,k62,k63,k64,k65,k66,k67,k68,k69,k70,k71,k72, &
          k73,k74,k75,k76,k77,k78,k79,k80,k81,k82,k83,k84, &
          k85,k86,k87,k88,k89,k90,k91,k92,k93,k94,k95,k96

    
    
    if(iDiffCoeffTableActive.eq.0)  then
!    -- there is no tables, use the mat input data
        rDiffCoeff=thermalk(elemId)
        write(*,*) 'iDiffCoeffTableActive, rDiffCoeff=thermalk(elemId),rTemperature',iDiffCoeffTableActive,rDiffCoeff,rTemperature
        return
    endif 
    counterMat=0
    do i = 1, iDiffCoeffTableCount
        !--read these ids
        do j = counterMat+1, iDiffCoeffTableMatCount(i)
!            -- check if the elem mat is in one of the tables
            if (a(k08+elemId-1) == iDiffCoeffTableMatIDs(j)) then
                rDiffCoeff=DiffCoeffTableInterpolate(i,rTemperature)
            end if
        end do
    end do

         write(*,*) 'iDiffCoeffTableActive, rDiffCoeff=DiffCoeffTableInterpolate,rTemperature',iDiffCoeffTableActive,rDiffCoeff,rTemperature
   
 END subroutine DiffCoeffTableGetD
!###################################################################
!###################################################################
!###################################################################
! thermal curve interpolation
 function DiffCoeffTableInterpolate(TableId,rTemperature)
    use mod_DiffCoeffTable
    implicit none
    integer i,TableId,istart
    real*8 DiffCoeffTableInterpolate, rTemperature,mslope
    
    istart=iDiffCoeffTableptsCount(TableId)
        if (rTemperature < rDiffCoeffTablepts(istart,1)) then
            DiffCoeffTableInterpolate=rDiffCoeffTablepts(istart,2)
        endif
        if (rTemperature > rDiffCoeffTablepts(i+1,1)) then
            DiffCoeffTableInterpolate=rDiffCoeffTablepts(i+1,2)
        end if
    
!    loop over all the lines to find which one has the time=mtime
    do i = 1, iDiffCoeffTableptsCount(TableId)-1
        if(rTemperature >= rDiffCoeffTablepts(i,1) .and. rTemperature <= rDiffCoeffTablepts(i+1,1))then
!            get the slope of the line
            mslope=(rDiffCoeffTablepts(i+1,2)-rDiffCoeffTablepts(i,2))/(rDiffCoeffTablepts(i+1,1)-rDiffCoeffTablepts(i,1))
!            calc the temp
            DiffCoeffTableInterpolate = rDiffCoeffTablepts(i,2)+(rTemperature-rDiffCoeffTablepts(i,1))*mslope
        endif
    end do
    return
    
 END function DiffCoeffTableInterpolate