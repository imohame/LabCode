

!!###################################################################
!###################################################################
!###################################################################
subroutine FractReadApplyPreCrackCleavagePlanes(y,z,ix)
    use mod_Fract
    use  mod_parameters
    common /wblock8/ abc(573,nume,4),his(573,nume,4)
    common /precrack/ npc, elepc(100),npc_count
    common /stressflag/ ElemFractCode(nume),ElemDecayCount(nume)
    common /crackopen/ ElemDecayed(nume), overlapele(2,nume)
    common /crackline/ ncleave(3,nume), elecrack(4,nume),nodeflag(4,nume)
    common /overlapping/ intersec(4, nume), area_coeff(nume), update_flag
    
    dimension y(*), z(*), ix(4,*)
    integer elecrack,nodeflag
    real ncleave
    real intersec,area_coeff
    real y,z
    integer ix
    integer ElemFractCode, ElemDecayed,npc, elepc,npc_count
    integer*4 setvbuf3f_local
    integer i,j,elemId,IERR
    real Vcx,Vcy
    real abc,his
    real vn(3)
    integer ElemEdge1,ElemEdge2
    real rxc1,ryc1,r1,rxc2,ryc2,r2

    write(*, *) '-->>>>>>>>FractReadApplyPreCrackCleavagePlanes'
    IERR=0
    i=0
    FractLinesCount=0
    bFractFlag=0
    write(*, *) 'beginning of CrackProp ... reading cracks.in'

    open(60, file = 'cracks.in', status = 'old',IOSTAT=IERR, ERR=90)
    read (60, *) FractLinesCount
    
!!! if the fracture is not active then exit
    if (FractLinesCount == 0) then
        close(60)
        return
    endif

    bFractFlag=1

    open(6011, file = 'cracks2.out', status = 'unknown')
    ierr=setvbuf3f_local(7011,1,100)

    write (6011, *) FractLinesCount
    read(60, *) Vcx,Vcy
    write(6011, *) Vcx,Vcy

    npc=FractLinesCount
    npc_count=npc

    do i = 1, FractLinesCount
        !!!--------- id  ed1     rxc1    ryc1     r1  ed2     rxc2    ryc2     r2
        read(60, *) elemId,ElemEdge1,rxc1,ryc1,r1,ElemEdge2,rxc2,ryc2,r2
        elepc(i)=elemId
        write(6011, *) elemId,ElemEdge1,rxc1,ryc1,r1,ElemEdge2,rxc2,ryc2,r2
        write(*, *) 'Pre-exist crack, element', elemId, 'cracks.in'
        vn(1)=0.0
        vn(2)=Vcx
        vn(3)=Vcy
        ElemFractCode(elemId)=2
        ElemDecayed(elemId)=1

        do j = 1, 3
          abc(362+j,elemId,1) = vn(j)
          abc(365+j,elemId,1) = vn(j)
          abc(368+j,elemId,1) = vn(j)
        end do
        !!--- this important b/c if the ElemFractCode(ele)=2, then it does not execute propagate for this element
        !!--- but it will execute overlap for intersection, so it needs the correct direction
        ncleave(2,elemId)=Vcx
        ncleave(3,elemId)=Vcy
        !!---- this to fill all the prop for the cracked elements
        elecrack(1,elemId)=ElemEdge1
        elecrack(2,elemId)=3
        elecrack(3,elemId)=ElemEdge2
        elecrack(4,elemId)=3
        
        intersec(1, elemId)=rxc1
        intersec(2, elemId)=ryc1
        intersec(3, elemId)=rxc2
        intersec(4, elemId)=ryc2
        if (i < FractLinesCount) then
        !-- update the neighbor elements
            call FracUpdateTipBeforeCrack(ix, elemId)
        else
            call FracUpdateTipAfterCrack(ix, elemId, ElemEdge2, 2)
        endif
        
    end do

    close(60)
    close(6011)
    !--- no need for the excrack to do anything, so set these to zero
    npc_count=0
    npc=0
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

!###################################################################
!###################################################################
!###################################################################
! memory clean up
subroutine FractCleanMemory()

    use mod_Fract
    implicit none

    IF (ALLOCATED (FractElemIds))       DEALLOCATE (FractElemIds)

end subroutine FractCleanMemory
! !###################################################################
! !###################################################################
! !###################################################################
! subroutine FractApplyCleavagePlanes(y,z,ix)
!     use mod_Fract
!     use  mod_parameters
!     common /wblock8/ abc(573,nume,4),his(573,nume,4)
!     dimension y(*), z(*), ix(4,*)
!
!     integer i,j,e1,e2,ej
!     real*8 ye1,ze1,ye2,ze2,CrackLen,Vcrack(2),Vncrack(2)
! !!! ------------ read the crack lines
!     CALL FractReadCrackLines()
!     if (bFractFlag == 0) return
!
! !
! !     i=0
! !     write(*, *) '----Apply the pre-cracked cleave  planes to the elements'
! !     do i = 1, FractLinesCount
! ! !! get element ids
! !       e1=FractElemIds(i,1)
! !       e2=GBElemAreaMap(i,2)
! ! !! --  cg of the first elem
! !       ye1=(y(ix(1,e1))+y(ix(2,e1)) +y(ix(3,e1))+y(ix(4,e1)))*0.25
! !       ze1=(z(ix(1,e1))+z(ix(2,e1)) +z(ix(3,e1))+z(ix(4,e1)))*0.25
! !       !! --  cg of the second elem
! !       ye2=(y(ix(1,e2))+y(ix(2,e2)) +y(ix(3,e2))+y(ix(4,e2)))*0.25
! !       ze2=(z(ix(1,e2))+z(ix(2,e2)) +z(ix(3,e2))+z(ix(4,e2)))*0.25
! !       !! --  crack length
! !       CrackLen=sqrt((ye2-ye1)*(ye2-ye1)+(ze2-ze1)*(ze2-ze1))
! ! !! -- if the line has zero length then go to next line
! !       if (CrackLen < 1e-8) continue
! ! !!! -- calculate the crack vector
! !       Vcrack(1)=(ye2-ye1)/CrackLen
! !       Vcrack(2)=(ze2-ze1)/CrackLen
! ! !!! -- calculate the crack normal
! !       Vncrack(1)=Vcrack(2)
! !       Vncrack(2)=-Vcrack(1)
! !
! !       do ej = 1, nume
! !         ye1=(y(ix(1,ej))+y(ix(2,ej)) +y(ix(3,ej))+y(ix(4,ej)))*0.25
! !         ze1=(z(ix(1,ej))+z(ix(2,ej)) +z(ix(3,ej))+z(ix(4,ej)))*0.25
! !
! !
! !       END do
! !       do j = 1, 3
! !         abc(362+j,GBelemId,1) = 0.0
! !         abc(365+j,GBelemId,1) = GBNormals(GBareaId,1)
! !         abc(368+j,GBelemId,1) = GBNormals(GBareaId,2)
! !       end do
! !     end do
!
!
! END
