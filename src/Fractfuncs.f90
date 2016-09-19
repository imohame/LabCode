

!!###################################################################
!###################################################################
!###################################################################
subroutine FractReadApplyPreCrackCleavagePlanes2(x,y,ix)
    use EC_Objects_manager
!    use EC_ElemCrackingBaseClass
    use  mod_parameters
    common /wblock8/ abc(573,nume,4),his(573,nume,4)
    common/meshnum/ numnpo, numelto

    ! common /precrack/ npc, elepc(100),npc_count
    ! common /stressflag/ ElemFractCode(nume),ElemDecayCount(nume)
    ! common /crackopen/ ElemDecayed(nume), overlapele(2,nume)
    ! common /crackline/ ncleave(3,nume), elecrack(4,nume),nodeflag(4,nume)
    ! common /overlapping/ intersec(4, nume), area_coeff(nume), update_flag

    dimension x(*), y(*), ix(4,*)
    real x,y
    integer ix

    ! integer elecrack,nodeflag
    ! real ncleave
    ! real intersec,area_coeff
    ! integer ElemFractCode, ElemDecayed,npc, elepc,npc_count
    integer*4 setvbuf3f_local
    integer i,j,elemId,IERR
    real Vcx,Vcy
    real abc,his
    real vn(3)
    integer ElemEdge1,ElemEdge2
    real rxc1,ryc1,r1,rxc2,ryc2,r2
    character strTitle(200)
    real*8 rxc2Ratio,ryc2Ratio,rxc1Ratio,ryc1Ratio

    integer mEC_bCracking,mEC_ZoneFactor,mEC_DecayCount,mEC_PreCrackedElemCount
    write(*, *) '-->>>>>>>>FractReadApplyPreCrackCleavagePlanes'
    IERR=0
    i=0
    mEC_PreCrackedElemCount=0
    mEC_bCracking=0
    write(*, *) 'beginning of CrackProp ... reading cracks.in'

    open(60, file = 'cracks.in', status = 'old',IOSTAT=IERR, ERR=90)
    read (60, *) mEC_bCracking !-- 1=allow fracture
    if(mEC_bCracking==0) then
      close(60)
      return
    endif
    open(6011, file = 'cracks2.out', status = 'unknown')
    ierr=setvbuf3f_local(7011,1,100)

    write(6011, *) mEC_bCracking !-- 1=allow fracture

    read (60, *)    mEC_ZoneFactor !-- 20 * min elem diagonal
    write (6011, *) mEC_ZoneFactor !-- 20 * min elem diagonal

    read (60, *)    mEC_DecayCount !-- unloading steps
    write (6011, *) mEC_DecayCount !-- unloading steps

    read (60, *)    mEC_PreCrackedElemCount  !-- number of pre-cracked elems
    write (6011, *) mEC_PreCrackedElemCount  !-- number of pre-cracked elems
    !-- set the input data to the module
    call ECinitializeConsts(numnpo, numelto,ix,x,y,mEC_bCracking,mEC_ZoneFactor, &
                                mEC_DecayCount,mEC_PreCrackedElemCount)
    call EC_AllocateMem()

!!! if the fracture is not active then exit
    if (mEC_PreCrackedElemCount == 0) then
        close(60)
        close(6011)
        return
    endif

    write (6011, *) mEC_PreCrackedElemCount
    read(60, *) Vcx,Vcy
    write(6011, *) Vcx,Vcy

    read(60, *) strTitle
    write(6011, *) strTitle !'!!!--------- id  ed1     rxc1    ryc1     r1  ed2     rxc2    ryc2     r2'
    do i = 1, mEC_PreCrackedElemCount
        !!!--------- id  ed1     rxc1    ryc1     r1  ed2     rxc2    ryc2     r2
        read(60, *) elemId,ElemEdge1,rxc1,ryc1,r1,ElemEdge2,rxc2,ryc2,r2,rxc2Ratio,ryc2Ratio,rxc1Ratio,ryc1Ratio
        ! elepc(i)=elemId
        write(6011, *) elemId,ElemEdge1,rxc1,ryc1,r1,ElemEdge2,rxc2,ryc2,r2,rxc2Ratio,ryc2Ratio,rxc1Ratio,ryc1Ratio
        write(*, *) 'Pre-exist crack, element', elemId, 'cracks.in'
        vn(1)=0.0
        vn(2)=Vcx
        vn(3)=Vcy
        ! ElemFractCode(elemId)=2
        ! ElemDecayed(elemId)=1

        do j = 1, 3
          abc(362+j,elemId,1) = vn(j)
          abc(365+j,elemId,1) = vn(j)
          abc(368+j,elemId,1) = vn(j)
        end do
        !!--- this important b/c if the ElemFractCode(ele)=2, then it does not execute propagate for this element
        !!--- but it will execute overlap for intersection, so it needs the correct direction
        ! ncleave(2,elemId)=Vcx
        ! ncleave(3,elemId)=Vcy
        ! !!---- this to fill all the prop for the cracked elements
        ! elecrack(1,elemId)=ElemEdge1
        ! elecrack(2,elemId)=3
        ! elecrack(3,elemId)=ElemEdge2
        ! elecrack(4,elemId)=3
        !
        ! intersec(1, elemId)=rxc1
        ! intersec(2, elemId)=ryc1
        ! intersec(3, elemId)=rxc2
        ! intersec(4, elemId)=ryc2
        call pEC_ElemData(elemId)%SetDataPre(vn,ElemEdge1,rxc1,ryc1,r1,ElemEdge2,rxc2,ryc2,r2)
        ! write(*,*)'-----------------------------------'
        call pEC_ElemData(elemId)%PrintTest(elemId)
        ! if (i < FractElemCountPreCracked) then
        ! !-- update the neighbor elements
        !     call FracUpdateTipBeforeCrack(ix, elemId)
        ! else
        !     call FracUpdateTipAfterCrack(ix, elemId, ElemEdge2, 2)
        ! endif

    end do

    close(60)
    close(6011)
    write(*, *) '-->>>>>>>>FractReadApplyPreCrackCleavagePlanes'
    return
!    -------------------------------------------
!    -------------------------------------------
    90  IF (IERR .EQ. 29 ) THEN
        write(*, *) '------------------------ no input file .. reading cracks.in'
      ELSE
        write(*, *) 'Unrecoverable error in .. reading cracks.in, code =', IERR
        STOP
      END IF

    return
END
!!###################################################################
!###################################################################
!###################################################################
subroutine FractReadApplyPreCrackCleavagePlanes(NodesCoordx,NodesCoordy,ElemConnect)

    use mod_parameters
    use EC_Objects_manager
    implicit none

    real NodesCoordx(*), NodesCoordy(*)
    INTEGER ElemConnect(4,*)


    common /wblock8/ abc(573,nume,4),his(573,nume,4)
    real abc,his
    common/meshnum/ numnpo, numelto
    integer numnpo, numelto

    integer*4 setvbuf3f_local
    integer i,j,elemId,IERR
    real Vcx,Vcy
    real vn(3)
    integer ElemEdge1,ElemEdge2
    real rxc1,ryc1,r1,rxc2,ryc2,r2
    character strTitle(200)
    real*8 rxc2Ratio,ryc2Ratio,rxc1Ratio,ryc1Ratio
    integer mEC_bCracking,mEC_ZoneFactor,mEC_DecayCount,mEC_PreCrackedElemCount
    real*8 x(4),y(4),Ptc(2)

    write(*, *) '-->>>>>>>>FractReadApplyPreCrackCleavagePlanes'
    IERR=0
    i=0
    mEC_PreCrackedElemCount=0
    mEC_bCracking=0
    write(*, *) 'beginning of CrackProp ... reading cracks.in'

    open(60, file = 'cracks.in', status = 'old',IOSTAT=IERR, ERR=90)
    read (60, *) mEC_bCracking !-- 1=allow fracture
    if(mEC_bCracking==0) then
      close(60)
      return
    endif
    open(6011, file = 'cracks2.out', status = 'unknown')
    ierr=setvbuf3f_local(7011,1,100)

    write(6011, *) mEC_bCracking !-- 1=allow fracture

    read (60, *)    mEC_ZoneFactor !-- 20 * min elem diagonal
    write (6011, *) mEC_ZoneFactor !-- 20 * min elem diagonal

    read (60, *)    mEC_DecayCount !-- unloading steps
    write (6011, *) mEC_DecayCount !-- unloading steps

    read (60, *)    mEC_PreCrackedElemCount  !-- number of pre-cracked elems
    write (6011, *) mEC_PreCrackedElemCount  !-- number of pre-cracked elems
    !-- set the input data to the module
    call ECinitializeConsts(numnpo, numelto,ElemConnect,NodesCoordx,NodesCoordy, &
                           mEC_bCracking,mEC_ZoneFactor,mEC_DecayCount,mEC_PreCrackedElemCount)
    call EC_AllocateMem()

!!! if the fracture is not active then exit
    if (mEC_PreCrackedElemCount == 0) then
        close(60)
        close(6011)
        return
    endif

    write (6011, *) mEC_PreCrackedElemCount
    read(60, *) Vcx,Vcy
    write(6011, *) Vcx,Vcy
    vn(1)=0.0
    vn(2)=Vcx
    vn(3)=Vcy

    read(60, *) strTitle
    write(6011, *) strTitle !'!!!--------- id  ed1     rxc1    ryc1     r1  ed2     rxc2    ryc2     r2'
    do i = 1, mEC_PreCrackedElemCount
        !!!--------- id  ed1     rxc1    ryc1     r1  ed2     rxc2    ryc2     r2
        read(60, *) elemId,ElemEdge1,rxc1,ryc1,r1,ElemEdge2,rxc2,ryc2,r2,rxc2Ratio,ryc2Ratio,rxc1Ratio,ryc1Ratio
        write(6011, *) elemId,ElemEdge1,rxc1,ryc1,r1,ElemEdge2,rxc2,ryc2,r2,rxc2Ratio,ryc2Ratio,rxc1Ratio,ryc1Ratio
        write(*, *) 'Pre-exist crack, element', elemId, 'cracks.in'

        do j = 1, 3
          abc(362+j,elemId,1) = vn(j)
          abc(365+j,elemId,1) = vn(j)
          abc(368+j,elemId,1) = vn(j)
        end do

        call pEC_ElemData(elemId)%SetDataPre(vn,ElemEdge1,rxc1,ryc1,r1,ElemEdge2,rxc2,ryc2,r2)
        Ptc(1)=rxc1
        Ptc(2)=ryc1
        !- to get the elem nodes
        do j=1, 4
            x(j)=NodesCoordx(ElemConnect(j,i))
            y(j)=NodesCoordy(ElemConnect(j,i))
        enddo
        call EC_MarkElemForCrack(elemId,Ptc,x,y)
        !--- for debugging
        ! call pEC_ElemData(elemId)%PrintTest(elemId)

    end do

    close(60)
    close(6011)
    write(*, *) '-->>>>>>>>FractReadApplyPreCrackCleavagePlanes'
    return
!    -------------------------------------------
!    -------------------------------------------
    90  IF (IERR .EQ. 29 ) THEN
        write(*, *) '------------------------ no input file .. reading cracks.in'
      ELSE
        write(*, *) 'Unrecoverable error in .. reading cracks.in, code =', IERR
        STOP
      END IF

    return
END

! !###################################################################
! !###################################################################
! !###################################################################
 subroutine FractReadElemNeighbors()
     use EC_Objects_manager
     implicit none

    integer i,IERR
    integer*4 setvbuf3f_local
    integer EC_ElemNeighbors(8)

    write(*, *) '-->>>>>>>>FractReadElemNeighbors -- start'
    IERR=0
    i=0
    write(*, *) 'beginning of CrackProp ... reading cracks.in'

    open(60, file = 'ElemsNeighbors.in', status = 'old',IOSTAT=IERR, ERR=90)
    open(6011, file = 'ElemsNeighbors2.out', status = 'unknown')
    ierr=setvbuf3f_local(6011,1,100)

    do i = 1, EC_ElemCountInput

      read (60, *)    EC_ElemNeighbors(1:8)
      write (6011, *)    EC_ElemNeighbors(1:8)
      call pEC_ElemData(i)%EC_SetElemEdgeNeighbors(EC_ElemNeighbors)
    end do

    close(60)
    close(6011)
    write(*, *) '-->>>>>>>>FractReadElemNeighbors -- done'
    return
!    -------------------------------------------
!    -------------------------------------------
    90  IF (IERR .EQ. 29 ) THEN
        write(*, *) '------------------------ no input file .. reading cracks.in'
      ELSE
        write(*, *) 'Unrecoverable error in .. reading cracks.in, code =', IERR
        STOP
      END IF

    return
 END
