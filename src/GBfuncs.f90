
!
!!###################################################################
!###################################################################
!###################################################################
!! to read the GB data, normals, element map
subroutine GBReadInput()
  CALL GBReadNormals()
  CALL GBReadElemMap()
end subroutine GBReadInput

!!###################################################################
!###################################################################
!###################################################################
subroutine GBReadNormals()
    use mod_GBdata

    implicit none
    integer*4 setvbuf3f_local
    integer i,IERR

    IERR=0
    i=0
    GBAreaCount=0
    write(*, *) 'beginning of ReadGBNormals ... reading GBnormal.in'

    open(60, file = 'GBnormal.in', status = 'unknown',IOSTAT=IERR, ERR=90)
    open(6011, file = 'GBnormal2.in', status = 'unknown')
    ierr=setvbuf3f_local(6011,1,100)

    read (60, *) GBAreaCount
    write (6011, *) GBAreaCount
    if (GBAreaCount == 0 ) then
      goto 90
    endif

    allocate(GBNormals(GBAreaCount,2))
    do i = 1, GBAreaCount
      read(60, *) GBNormals(i,1),GBNormals(i,2)
      write(6011, *) GBNormals(i,1),GBNormals(i,2)
    end do

    close(60)
    close(6011)

    return
!    -------------------------------------------
!    -------------------------------------------
    90  IF (IERR .EQ. 29 ) THEN
        write(*, *) '------------------------ no input file .. reading GBnormal.in'
      ELSE
        write(*, *) 'Unrecoverable error in .. reading GBnormal.in, code =', IERR
        STOP
      END IF

    return
END

!###################################################################
!###################################################################
!###################################################################
subroutine GBReadElemMap()
    use mod_GBdata

    implicit none
    integer*4 setvbuf3f_local
    integer i,IERR

    IERR=0
    i=0
    GBElemCount=0
    write(*, *) 'beginning of ReadGBNormals ... reading GBelemGrainMap.in'

    open(60, file = 'GBelemGrainMap.in', status = 'old',IOSTAT=IERR, ERR=90)
    open(6011, file = 'GBelemGrainMap2.in', status = 'unknown')
    ierr=setvbuf3f_local(6011,1,100)

    read (60, *) GBElemCount
    write (6011, *) GBElemCount
    if (GBElemCount == 0 ) then
      goto 90
    endif

    allocate(GBElemAreaMap(GBElemCount,2))
    do i = 1, GBElemCount
      read(60, *) GBElemAreaMap(i,1),GBElemAreaMap(i,2)
      write(6011, *) GBElemAreaMap(i,1),GBElemAreaMap(i,2)
    end do

    close(60)
    close(6011)

    return
!    -------------------------------------------
!    -------------------------------------------
    90  IF (IERR .EQ. 29 ) THEN
        write(*, *) '------------------------ no input file .. reading GBelemGrainMap.in'
      ELSE
        write(*, *) 'Unrecoverable error in .. reading GBelemGrainMap.in, code =', IERR
        STOP
      END IF

    return
END

!###################################################################
!###################################################################
!###################################################################
! memory clean up
subroutine GBCleanMemory()

    use mod_GBdata
    implicit none

    IF (ALLOCATED (GBNormals))       DEALLOCATE (GBNormals)
    IF (ALLOCATED (GBElemAreaMap))     DEALLOCATE (GBElemAreaMap)

end subroutine GBCleanMemory
!###################################################################
!###################################################################
!###################################################################
subroutine GBApplyCleavagePlanes()
    use mod_GBdata
    use  mod_parameters
    common /wblock8/ abc(573,nume,4),his(573,nume,4)


    integer i,j,GBelemId,GBareaId

    i=0
    write(*, *) '----Apply the GB cleave  planes to the elements'
    do i = 1, GBElemCount
      GBelemId=GBElemAreaMap(i,1)
      GBareaId=GBElemAreaMap(i,2)
      do j = 1, 3
        abc(362+j,GBelemId,1) = 0.0
        abc(365+j,GBelemId,1) = GBNormals(GBareaId,1)
        abc(368+j,GBelemId,1) = GBNormals(GBareaId,2)
      end do
    end do

END
