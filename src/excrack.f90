      subroutine excrack(nstep, y, ix)
!!!!c     generate pre-exist crack
!!!
!!!	  parameter(nume = 40000)
!!!	  common /stressflag/ ElemFractCode(nume),ElemDecayCount(nume)
!!!	  common /crackopen/ ElemDecayed(nume), overlapele(2,nume)
!!!	  common /crackline/ ncleave(3,nume), elecrack(4,nume),nodeflag(4,nume)
!!!	  common /precrack/ npc, elepc(100),npc_count
!!!!!!!	  common /excflag/ excf
!!!!!!	  common /tipvelocity/ncrack,nelefail(1000),tipelenum(1000,nume)
!!!!!!	  common /precrack2/ rlflag
!!!      common /sigfrac/ sigmacrit0, sigmacrit1, sigmacrit2,sigmacrit3,DecayCount, f_decay, penalty,fractFlag
!!!!!!!      common /fractureplane/ planeflag(nume), planeflag0(nume)
!!!
!!!	  integer ElemFractCode, ElemDecayed, npc, elepc, nstep,npc_count
!!!!!!	  integer ncrack, nelefail, tipelenum, rlflag
!!!
!!!!!!!	  integer excf
!!!      integer  ele, i, j, dum
!!!
!!!	  real ncleave, ycoord(1000)
!!!
!!!	  dimension y(*), ix(4,*)
!!!      integer ix
!!!      real y
!!!
!!!!!	  if(nstep==0 .and. excf==1) then
!!!	  if(npc_count > 0) then
!!!        write(*, *) '-->>>>>>>>excrack --start'
!!!        write(*, *) 'Pre-exist crack, # element', npc
!!!!!!!!!      write(*, *) 'Pre-exist crack, # element', elepc
!!!!!!        do i=1, npc
!!!        i=npc-(npc_count-1)
!!!        ele=elepc(i)
!!!!!!          write(*, *) 'Pre-exist crack, element', elepc(i)
!!!        ElemFractCode(ele)=2
!!!!!            ElemDecayCount(ele) = DecayCount+1
!!!        ElemDecayed(ele)=1
!!!        npc_count=npc_count-1
!!!!!!!!            planeflag(ele)=planeflag0(ele)
!!!!!			  ncleave(2,ele)=0.0
!!!!!			  ncleave(3,ele)=1.0
!!!!!		      write(*, *) 'Pre-crack excrack.f90 with ncleave(2:3,ele)=[0 1] , element', ele, 'excrack'
!!!		      write(*, *) 'Pre-crack excrack.f90 without ncleave(2:3,ele)=[0 1] , element', ele, 'excrack'
!!!!!!        end do
!!!
!!!
!!!!!!!!!c         crack velocity part
!!!!!!!!        if(npc>0) then
!!!!!!!!            write(*, *) '-->>>>>>>> crack velocity part'
!!!!!!!!            ncrack=1
!!!!!!!!            if(rlflag==0) then   ! propagating to from left to right
!!!!!!!!                nelefail(1)=1
!!!!!!!!                tipelenum(1,1)=0
!!!!!!!!                nelefail(2)=npc
!!!!!!!!                do i=1,npc
!!!!!!!!                    tipelenum(2,i)=elepc(i)
!!!!!!!!                end do
!!!!!!!!            else if(rlflag==1) then  ! propagating to from right to left
!!!!!!!!!          write(*, *) '>>>>>>>>> nelefail(2)',nelefail(2)
!!!!!!!!!          write(*, *) '>>>>>>>>> tipelenum(2,1)',tipelenum(2,1)
!!!!!!!!                nelefail(2)=1
!!!!!!!!                tipelenum(2,1)=0
!!!!!!!!                nelefail(1)=npc
!!!!!!!!                do i=1,npc
!!!!!!!!                    tipelenum(1,i)=elepc(i)
!!!!!!!!                end do
!!!!!!!!            end if
!!!!!!!!        end if
!!!!!!        excf=0
!!!        write(*, *) '-->>>>>>>>excrack --end'
!!!    end if !! if(nstep==0 .and. excf==1) then

end
