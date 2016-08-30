      subroutine overlap(y, z, ix, matp, id, u, usi, freep, ym)
	  
      use mod_parameters
      use CN_Objects_manager
      use mod_file_units
      
	  common/wblock8/  abc(573,nume,4), his(573,nume,4) 
	  common /wblock12/ Y_modulus(nume),possion_ratio(nume),tau_y(nume)
      common/bk02/ioofc,iphase,imass,lpar(9)
	  common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
	  common/bk08/kprint,nstep,ite,ilimit,newstf
	  common/bk10/npb,nodep(2,8)
	  common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
	  common/bk26/    dt,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
      common/bk32/nsref,nequit,time,timep,lprint,nprint
      common/fissn0/maxneq,mwspac,ntpe0,ntpe1,nfissl(3)
	  common/hokao/    lst,nnn2(nume,4)
      common/hourglass/fhg(40000,8),fhghis(40000,8),fhg1(nelemg), &
                       fhg2(nelemg),fhg3(nelemg),fhg4(nelemg),fhg5(nelemg), &
                       fhg6(nelemg),fhg7(nelemg),fhg8(nelemg)
      common/hourglass2/hgsstore(40000),hgshis(40000)
      common/totalenergy/totenerstore(40000),totenerhis(40000),inertener(40000)
      common/hgenergy/hgenerstore(40000),hgenerhis(40000)
      common/hgstress/hgstress1store(40000),hgstress2store(40000), &
                      hgstress1his(40000),hgstress2his(40000)
	 
	  common/meshnum/ numnpo, numelto  
	  common /stressflag/ ElemFractCode(nume),ElemDecayCount(nume)
	  common /overlapping/ intersec(4, nume), area_coeff(nume),update_flag
	  common /stroverlap/ sflg
	  common /crackopen/ ElemDecayed(nume), overlapele(2,nume)
      common /crackline/ ncleave(3,nume), elecrack(4,nume), nodeflag(4,nume)
	  common /ovsum/ ovs
!!!!!!!!!      common /WMLthermal2/thermalk(40000),h(40000),etae(40000)
!!!!!!      common /WMLthermal2/thermalk(nume),thermalh(nume),
!!!!!!     >      thermalRo(nume),thermalcp(nume),thermalx(nume)
!!!!!!     >      ,thermalD(nume)
      
!!!!!!!!!!	  common/WMLthermalmatpoly/DijSije(nume),Tele(nume)
!!!!!!      common/WMLthermal/thermalflag !!!!!!!,thermalconstraint(nume), Tinit(nume),Rqold(nume)
	  common/hydroembrittle/critfrac(1000), sigfrac0(nume),sigfrac(nume),decfrac(nume)
	  common /sigfrac/ sigmacrit0, sigmacrit1, sigmacrit2, &
                       sigmacrit3,DecayCount, f_decay, penalty,fractFlag
!!!!	  common /fractureplane/ planeflag(nume), planeflag0(nume)
!!!!	  common /tipvelocity/ ncrack, nelefail(1000),tipelenum(1000,nume)
	  
	  dimension y(*), z(*), ix(4,*), matp(*), id(2,*), u(*), usi(*),freep(5,*), ym(4,*)
	  
	  integer numelt, numnp, ElemFractCode, ele,update_flag,fractFlag
	  integer sflg, nnn2, ElemDecayed, overlapele
	  integer lprint, nprint, nstep, numelto
	  integer elecrack, ovs, ElemDecayCount
      integer planeflag
	  real intersec, ym, area_coeff, sigfrac, decfrac
!!!	  integer ncrack, nelefail, tipelenum 
	  integer ic, nlast, itn, nfile
	  real yt, zt
	  
      if(fractFlag==0) then !off
        return
      endif
      
	  update_flag=0
    do ele=1, numelto
        if(((ElemFractCode(ele)==1) .AND.(ElemDecayCount(ele)==1))  &
           .or. (nstep==0 .and. ElemFractCode(ele)==2)) then
            !!!---- set crackline based on cleavage plane
            if(elecrack(2,ele)==0 .and. elecrack(4,ele)==0) then
            !!----- based on integration point
                call crackline_int(ele, y, z, ix, id, u)
!!!!!!!!                !--- new crack nucleated, pre-existed crack updated in excrack.f
!!!!!!!!                if(nstep>0) then     
!!!!!!!!                    ncrack=ncrack+1           
!!!!!!!!                    tipelenum(2*ncrack-1,1)=ele
!!!!!!!!                    tipelenum(2*ncrack,1)=ele
!!!!!!!!                    nelefail(2*ncrack-1)=1
!!!!!!!!                    nelefail(2*ncrack)=1
!!!!!!!!                end if
            else if(elecrack(2,ele)==2 .and. elecrack(4,ele)==0) then
            !!---based on crack tip
                call crackline_tip(ele, y, z, ix, id, u) 
            else if(elecrack(2,ele)==2 .and. elecrack(4,ele)==2) then
            !!---based on two crack tips
                call crackline_2tips(ele, ix)
            end if  
        end if		  
        !-- if edge is decayed .. cracked and opened
        if ((ElemFractCode(ele)==2) .AND. (ElemDecayed(ele)==1)) then	
!           !!!------ update elem ridge neighbor edges from 1 to 3
            call FracUpdateTipBeforeCrack(ix, ele)
            !!!!!!!!!!!write(*,*)'numelt',numelt
            do i=1, 4
                elecrack(i, numelt+1)= elecrack(i, ele)
                intersec(i, numelt+1)= intersec(i, ele)
            end do

            call update_mesh(ele, y, z, ix, id, u, usi)
			!!--- copy material parameters and status variables to the new element
            matp(numelt+1)=matp(ele)

            do j=1,5
                freep(j,numelt+1)=freep(j,ele)
            end do

            !!--- lumped mass matrix			 
            do j=1,4
                ym(j,numelt+1)=ym(j,ele)*area_coeff(numelt+1)
                ym(j,ele)=ym(j,ele)*area_coeff(ele)
            end do

            do j = 1, 573
               his(j,numelt+1,1) = his(j,ele,1)
               abc(j,numelt+1,1) = abc(j,ele,1)
            end do

            Y_modulus(numelt+1)=Y_modulus(ele)
            possion_ratio(numelt+1) = possion_ratio(ele)
            tau_y(numelt+1) = tau_y(ele)
            ElemFractCode(numelt+1)=ElemFractCode(ele)
!!!!            planeflag(numelt+1)=planeflag(ele)

            nnn2(numelt+1,1)=nnn2(ele,1)
            maxneq=neq

            newstf=0
            update_flag=1
            sflg=0

            do j=1,8
              fhg(numelt+1,j)=fhg(ele,j)
            end do
            hgsstore(numelt+1)=hgsstore(ele)
            totenerstore(numelt+1)=totenerstore(ele)
            hgenerstore(numelt+1)=hgenerstore(ele)
            hgstress1store(numelt+1)=hgstress1store(ele)
            hgstress2store(numelt+1)=hgstress2store(ele)
             
            CALL CNmanager_CopyElemnt(ele,numelt)
			  
!!!!!!!!!!!!!!             rMatTH_k(numelt+1)=rMatTH_k(MatID)
!!!!!!!!!!!!!!             rMatTH_h(numelt+1)=rMatTH_h(MatID)
!!!!!!!!!!!!!!             rMatTH_Ro(numelt+1)=rMatTH_Ro(MatID)
!!!!!!!!!!!!!!             rMatTH_cp(numelt+1)=rMatTH_cp(MatID)
!!!!!!!!!!!!!!             rMatTH_x(numelt+1)=rMatTH_x(MatID)
!!!!!!!!!!!!!!             rMatDiff_D(numelt+1)=rMatDiff_D(MatID)
!!!!!!!!!!!!!!!!!!!!             thermalk(numelt+1)=thermalk(ele)
!!!!!!!!!!!!!!!!!!!!             thermalh(numelt+1)=thermalh(ele)
!!!!!!!!!!!!!!!!!!!!             thermalRo(numelt+1)=thermalRo(ele)
!!!!!!!!!!!!!!!!!!!!             thermalcp(numelt+1)=thermalcp(ele)
!!!!!!!!!!!!!!!!!!!!             thermalx(numelt+1)=thermalx(ele)             
!!!!!!!!!!!!!!!!!!!!             thermalD(numelt+1)=thermalD(ele)             
!!!!!!!!!!!!!!!!!!!!             
!!!!!!!!!!!!!!!!!!!             Tfl(numelt+1)=Tfl(ele)
!!!!!!!!!!!!!!	         DijSije(numelt+1)=DijSije(ele)
            sigfrac(numelt+1)=sigfrac(ele)
            decfrac(numelt+1)=decfrac(ele)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c               overlapping elements records
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            ovs=ovs+1
            overlapele(1,ovs)=ele
            overlapele(2,ovs)=numelt+1	
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
            numelt=numelt+1
            !-------  crack.out output  
            write(iFU_crack_out,*) 'ele:', ele, 'step:', nstep
            write(iFU_crack_out,*) 'ix: ', ix(1,ele), ix(2,ele),ix(3,ele),ix(4,ele)
            write(iFU_crack_out,*) 'edge:', elecrack(1,ele), elecrack(3,ele)
            write(iFU_crack_out,*) 'status:', elecrack(2,ele), elecrack(4,ele)
            write(iFU_crack_out,*) 'y_coeffi:', intersec(1,ele), intersec(3,ele)
            write(iFU_crack_out,*) 'z_coeffi:', intersec(2,ele), intersec(4,ele)
            write(iFU_crack_out,*) '    '
            flush(iFU_crack_out)
!!!!!!            !---- print information to calculate crack velocity
!!!!!!            do ic=1,2*ncrack
!!!!!!!!!                nlast=nelefail(ic)
!!!!!!                do itn=1,nlast
!!!!!!                    if(ele==tipelenum(ic,itn)) then
!!!!!!                        yt=(y(ix(1,ele))+u(id(1,ix(1,ele)))+y(ix(2,ele))+u(id(1,ix(2,ele))) &
!!!!!!                           +y(ix(3,ele))+u(id(1,ix(3,ele)))+y(ix(4,ele))+u(id(1,ix(4,ele))))/4.0
!!!!!!                        zt=(z(ix(1,ele))+u(id(2,ix(1,ele)))+z(ix(2,ele))+u(id(2,ix(2,ele))) &
!!!!!!                           +z(ix(3,ele))+u(id(2,ix(3,ele)))+z(ix(4,ele))+u(id(2,ix(4,ele))))/4.0
!!!!!!                        nfile=5000+ic
!!!!!!!                       !!!!---write(nfile,*) ele, nstep, nstep*dt, yt, zt
!!!!!!                        write(iFU_cracktip_out,*) 'ele, nstep, nstep*dt, yt, zt ',ele, nstep, nstep*dt, yt, zt
!!!!!!                    end if
!!!!!!                end do
!!!!!!            end do					  
        end if!!!if ((ElemFractCode(ele)==2) .AND. (ElemDecayed(ele)==1)) then
    end do !!! do ele=1, numelto
	  
    if (update_flag==1) then
        lpar(4)=numelt
        lpar(2)=lpar(4)
        nodep(2,1)=numnp
        lprint=lprint-1
        nprint=nprint-1
    end if
		    
end
   
	  
			  
			  
			  
			  
			  
			  
			  
			  
			  
			  
		      