      subroutine propagate(y, z, ix, id, u, matp, numelt,nstep)
!!!!!!!          propagate(a(k03),a(k04),a(k02),a(k57),a(k18), a(k08),numelt,nstep)

    use CN_Objects_manager 
    use mod_file_units
!!!!!!!!!!	  integer, parameter :: nume       = 40000
!!!!!!!!!!	  parameter (n_dim = 3, nss = 24, no_mat = 1000 )

	  common /stressflag/ ElemFractCode(nume),ElemDecayCount(nume)
	  common/wblock8/  abc(573,nume,4), his(573,nume,4)
	  common /wblock10/ng,grain_mo(1000,3),bv(no_mat,87),nssmat(1000),nssimmat(1000)

	  common/WMLthermal/thermalflag!!!!!,thermalconstraint(nume),Tinit(nume),Rqold(nume)

!!!	  common /gbblock/ gbvec(nume,3), gbflag(nume,3)
!!!!!!!!!!!!!	  common /propag/ sigalt(4,nume)
    common/main_block/  a(1)
    common /sigfrac/ sigmacrit0, sigmacrit1, sigmacrit2,sigmacrit3,DecayCount, f_decay, penalty,fractFlag
    common /crackopen/ ElemDecayed(nume), overlapele(2,nume)
    common /crackline/ ncleave(3,nume), elecrack(4,nume),nodeflag(4,nume)
    common/meshnum/ numnpo, numelto
!!!!!!!!!!!	  common/hydrodiffusion/ hycon(nume,1)
!!!!!!!!!	  common/couplinganalysis/ TDflag
    common/hydroembrittle/critfrac(1000), sigfrac0(nume), sigfrac(nume),decfrac(nume)
    common/hydroembrittle110/critfrac110(1000),sigfrac0110(nume)
    common/slipplane110/ nsp110(6)
    common /cracktipcircle/ ndcircle
!!    common /fractureplane/ planeflag(nume), planeflag0(nume)

    integer ElemFractCode,ElemDecayCount,ele,numelt,lst,nstep!!,gbflag
    real sig1, abc, stress, sig, sigmafrac1, sigmafrac110       !!!!!!!!,sigalt
    real dum, sigmacrit, ncleave, hycon
    integer DecayCount, ElemDecayed, crackele, oele, cflag,numelto,fractFlag
    real sigmacrit0, sigmacrit1, sigmacrit2, sigmacrit3
    integer nssmat, nssimmat,  thermalflag, ndcircle  !!TDflag,
    real sigfrac0, sigfrac, decfrac, sigfrac0110
    real sigmacrit100, sigmacrit110, cleave110(6,3)
    integer nsp110, iHE, jHE
!!!    integer planeflag, planeflag0
    real*8 MaxStress100,StressComp


    dimension y(*), z(*), ix(4,*), id(2,*), u(*), matp(*)
    dimension sig(4,numelt)
    real cleave(3,3),retval(4)

    data nsp110/1, 2, 3, 7, 8, 14/
    real rElemConcentration
    real diffExpo
    real vy,vz

    if(fractFlag==0) then !off
      return
    endif
!c     Hydrogen Embrittle Process
!c     decrease critical fracture stress on {110} planes,
!c     based on hydrogen concentration
!!!!!!      if(thermalflag==1 .and. TDflag==1) then  ---- ismailbug
!!!!      if(thermalflag > 0 .and. TDflag==1) then !! for diff only

    if(thermalflag >= 2) then !! for diff only and both;diff+thermal
        do i=1,numelto
            if(ElemFractCode(i)==0) then
                CALL CNmanager_Get_ElemConcentration(i,rElemConcentration)
                if(rElemConcentration<1.0) then
                    sigfrac(i)=sigfrac0110(i)
                    decfrac(i)=0.0
                else
                    diffExpo=0.0 !!-0.14
                    rElemConcentration=rElemConcentration**(diffExpo)
                    sigfrac(i)=sigfrac0110(i)*rElemConcentration
                    decfrac(i)=1.0-rElemConcentration
                !!!!!!!!			          sigfrac(i)=sigfrac0110(i)*(hycon(i,1)**(-0.14))
                !!!!!!!!				      decfrac(i)=1.0-hycon(i,1)**(-0.14)
                end if
            end if
        end do
    end if
!------------- to list all the decayed elems every step
    ElemDecayed(1:numelto)=0 

	call crackfront(y, ix)

    do 20 ele = 1, numelto
        !-- if the elem is decaying
        if(ElemFractCode(ele)==1 .and. ElemDecayCount(ele) .lt.DecayCount)then !-if the element is decaying
            ElemDecayCount(ele) = ElemDecayCount(ele) + 1 !-- increment the decay
            write(iFU_crackprog_out,*) ele,ElemFractCode(ele),ElemDecayCount(ele),nstep
            go to 20
        else if(ElemFractCode(ele)==1 .and. ElemDecayCount(ele)==DecayCount)then !-if the element decayed
            ElemFractCode(ele)=2   !-- decayed
            ElemDecayed(ele)=1
!!!            planeflag(ele)=planeflag0(ele)
            write(*,*) 'element ', ele, 'overlapping'
            go to 20
        else if(ElemFractCode(ele)==2) then !-- completely cracked
		        go to 20
        end if
        !!- be sure that the elem is not within circle distance of another crack
        call FracCheckElemNearestCrack(y, z, ix, id, u, ele, ndcircle, cflag)
        if(cflag==0) then
            go to 20
        end if

        retval=0.0
        CALL CNmanager_Get_sigalt(ele,retval)
        sig(1:4,ele)=retval(1:4)

!!c         assign critical fracture stress for current element
        sigmacrit100=sigfrac0(ele)		! {100} planes
        sigmacrit110=sigfrac(ele)         ! {110} planes

!!c     HE in lath martensite occurs on {110} planes ==  if(thermalflag > 0 .and. TDflag==1) then  !! for diff only
        if(thermalflag >= 2) then !! for diff only and both;diff+thermal
            do iHE=1,6
                jHE=nsp110(iHE)
                cleave110(iHE,1)=abc(105+jHE,ele,1)
                cleave110(iHE,2)=abc(129+jHE,ele,1)
                cleave110(iHE,3)=abc(153+jHE,ele,1)
            end do
            sigmafrac110=0.0
            do iHE=1,6
                vy=cleave110(iHE,2)
                vz=cleave110(iHE,3)
                s11=sig(1,ele)
                s22=sig(2,ele)
                s44=sig(4,ele)
                dum = s11*vy*vy+s22*vz*vz+s44*vy*vz*2.0
                !!!!!dum = sig(1,ele)*cleave110(iHE,2)**2.0+sig(2,ele)*cleave110(iHE,3)**2.0+sig(4,ele)*cleave110(iHE,2)*cleave110(iHE,3)*2.0
                if (abs(dum)>sigmafrac110) then
                    sigmafrac110=abs(dum)
                    ncleave(2,ele)=vy !cleave110(iHE,2)
                    ncleave(3,ele)=vz !cleave110(iHE,3)
                end if
            end do

            if ((ElemFractCode(ele)==0) .and. (sigmafrac110 .gt. sigmacrit110)) then !!.and. (gbflag(ele,1) == 0)) then
                ElemFractCode(ele) = 1 !-- just started the decay process
                ElemDecayCount(ele) = 1 !-- start the decay increment
!!                planeflag0(ele) = 2
                write(iFU_crackprog_out,*) ele,ElemFractCode(ele),ElemDecayCount(ele),nstep,'str', '110'
                go to 20
            end if

        end if

!!!c    only check failure criteria at the crack tip element
!!!c         cleavage fracture on {100} planes
        do j = 1, 3
            cleave(1,j) = abc(362+j,ele,1)
            cleave(2,j) = abc(365+j,ele,1)
            cleave(3,j) = abc(368+j,ele,1)
        end do

!!!c    obtain maximum normal component of the traction on cleavage planes and corresponding normal vector
        MaxStress100 = -100.0
        iPlaneId=-1
        do j = 1, 3
            vy=cleave(j,2)
            vz=cleave(j,3)
            s11=sig(1,ele)
            s22=sig(2,ele)
            s44=sig(4,ele)
            StressComp = s11*vy*vy+s22*vz*vz+s44*vy*vz*2.0
            !!!dum = sig(1,ele)*cleave(j,2)**2.0+sig(2,ele)*cleave(j,3)**2.0+sig(4,ele)*cleave(j,2)*cleave(j,3)*2.0
            if (abs(StressComp)>MaxStress100) then
                MaxStress100=abs(StressComp)
                iPlaneId=j
            end if
        end do

!!!c    compare stress on cleavage plane and slip plane, when elements fail
!!!c		  if(sigmafrac1>sigmacrit) then
!!!c		      call strcomp(ele, nstep)
!!!c          end if
!!!
!!!c     estimate failure, change element status flag
        if((ElemFractCode(ele)==0).and.(MaxStress100 .gt. sigmacrit100)) then !.and. (gbflag(ele,1) == 0)) then
            ncleave(2,ele)=cleave(iPlaneId,2)
            ncleave(3,ele)=cleave(iPlaneId,3)
        
            ElemFractCode(ele) = 1
            ElemDecayCount(ele) = 1
!!            planeflag0(ele) = -2
            !-- write to crackprog.out
            write(iFU_crackprog_out,*) ele,ElemFractCode(ele), ElemDecayCount(ele),nstep,MaxStress100,sigmacrit100,ncleave(2,ele),ncleave(3,ele)
        end if

    20 continue

	  call excrack(nstep, y, ix)

      return
	  end
