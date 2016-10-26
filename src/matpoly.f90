       subroutine matpoly(matp,sig,prop,ln)
!=======================================================================
!C		MATERIAL SUBROUTINE FOR MUTLI_SLIP SINGLE CRYSTALS
!  This subrotine calculates the Cauchy stresses wrt the reference axes
!  given the velocity gradient information from stfnf.f.  
!
!  Loading Finite Element Information
!
!  Stresses are loaded from the previous stable configuration, and split 
!  into hydrostatic and deviatoric parts.  The deformation rate tensor
!  is split to obtain the deviatoric deformation rate tensor for the 
!  calculation of the plastic deformation rate tensor.  The deformation
!  rate tensor is an Eulerian (current configuration) tensor expressed 
!  with respect to the reference axes. 
!
!  Crystal Plasticity
!
!  The reference shear stress is calculated given the immobile 
!  dislocation densities and Burgers vector lengths of the material slip 
!  systems.  The symmetric and asymmetric parts of the Schmid tensor are 
!  calculated and used for subsequent calculations of the resolved shear 
!  stress and plastic deformation rate and spin.  
!
!  The subroutine dsolve.f calculates the resolved shear stress using an 
!  rk5 or back-Euler integration scheme.  The resolved shear stresses 
!  are used to compute the strain rates on each system.  
!
!  The plastic deformation rates are calculated using the symmetric part 
!  of the Schmid tensor and the strain rates.  
!  
!  Constitutive Relationship
!
!  The plastic deformation rate tensor is subtracted from the total 
!  deformation rate tensor, which gives the elastic deformation rate, 
!  or D*, which is used with the elastic modulus to form the 
!  constitutive relationship for the Jaumann stress rate corotational 
!  with the lattice.  The MATERIAL stress rate, sigmadot, is the one
!  that needs to be programmed, since this code is set up using the 
!  reference axes (updated Lagrangian).
!
!  Final Calculations
!  
!  End result calculations (not used in other updates), such as
!  temperature, plastic work, and shear slip are computed here.
!  Also, dsolve1.f solves for the mobile and immobile dislocation 
!  densities, of which only the immobile ones are used in further 
!  calculations.
!
!=======================================================================
!	LOCAL PARAMETERS AND VARIABLES
!+++++++++++++++++++++++++++++++++
    use mod_parameters
    use CN_Objects_manager
    use EC_Objects_manager

    integer			 :: elem_mat_no, ss_m,nssm
    integer			 :: pd_counter(nume), ipr
    real		 :: n_dot1, n_dot2, s_dot1, s_dot2
    real		 :: den_t, slip_rate_ssm, dgamma
    dimension				sig(ln,*),prop(48,*),matp(*)
    dimension				Dij_dev(4),w12(nss), w21(nss)
    dimension				slip_n(nss,3), slip_s(nss,3), cleave(3,3)
! ======================================================================
!	GLOBAL PARAMETERS AND VARIABLES
!+++++++++++++++++++++++++++++++++++
      common/wblock2/  g_source, g_immob, g_minter, g_recov, b_v,b_vvec(87),nmo,nim
!!!      common/wblock3/  density_ms, density_ims,etain(1000),ecin(1000)!changed to accomodate material dependent parameters eta and enthalpy_coef
      common/wblock3/  density_ms, density_ims,thermalEnthalpy(1000)
      common/wblock4/  den_m(nss), den_im(nss), gdot(nss), nnns, nnne
      common/wblock5/  enthalpy_coef, thermal_coef, temp 
      common/wblock6/  xmhigh(nss),ref_gamma_dot(nss),Pij_Dijdev(nss), &
                      tau(nss), p(nss,4), taur(nss), twomu, g
!    >                 tau(nss), p(nss,3), taur, twomu, g                   !OLD WAY
      common/wblock7/  slip_n0(1000,nss,3), slip_s0(1000,nss,3)!changed to accomodate multiple slip WML 91109
      common/wblock8/  abc(573,nume,4), his(573,nume,4)        !changed to accomodate 14 slip systems WML
      common/wblock9/  slip_n0_t(no_mat,nss,3),slip_s0_t(no_mat,nss,3)
      common /wblock10/ng,grain_mo(1000,3),bv(no_mat,87),nssmat(1000),nssimmat(1000)!changed to accomodate multiple bv lengths WML
                                    
      common/wblock11/ pd_counter
      common/wblock12/ Y_modulus(nume),possion_ratio(nume),tau_y(nume)
      common/result/   press, plastic_work
      common/intgrt/   nintg
      common/bk16/     maxint,hgc
      common/hokao/    lst,nnn2(nume,4)
      common/berg/     dels11(nelemg),dels22(nelemg),dels33(nelemg),dels12(nelemg)
      common/vect21/   sig11(nelemg),sig22(nelemg),sig33(nelemg),sig12(nelemg)      
      common/vect3/ &
           dgi1(nelemg,4),dgi2(nelemg,4),dgi3(nelemg,4),dgi4(nelemg,4),dgi5(nelemg,4), &
           dgt1(nelemg,4),dgt2(nelemg,4),dgt3(nelemg,4),dgt4(nelemg,4), &
           f11v(nelemg),f22v(nelemg),f12v(nelemg),f21v(nelemg),dsd5(nelemg), &
           sig11s(nelemg),sig22s(nelemg),sig33s(nelemg),sig12s(nelemg), &
           ddp1(nelemg,4),ddp2(nelemg,4),ddp3(nelemg,4),ddp4(nelemg,4),ddp5(nelemg,4) 
      common/bk08/    kprint,nstep,ite,ilimit,newstf
      common/bk32/    nsref,nequit,time,timep,lprint,nprint
      common/meng/ &
     		f11dd(nelemg),f22dd(nelemg),f12dd(nelemg),f21dd(nelemg),f33dd(nelemg), &     !F matrix
     		f11dv(nelemg),f22dv(nelemg),f12dv(nelemg),f21dv(nelemg),f33dv(nelemg),  &   !Fdot matrix
     	    vg11(nelemg),vg22(nelemg),vg12(nelemg),vg21(nelemg),vg33(nelemg),fdet(nelemg), &   !FdotFinverse=velocity Gradient
     		d1(nelemg),d2(nelemg),d3(nelemg),d4(nelemg),spin(nelemg)                    !Dij=sym(FdotFinverse)
      common/custr/   sign1(nelemg),sign2(nelemg),sign3(nelemg),sign4(nelemg)
      common/range/   mft,mlt,lft,llt,nftm1
      common/bk26/    dt,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
      common/vect8/   dsave(4,4,nelemg)
      common/variab/  aaa1(nume),aaa2(nume),aaa3(nume),aaa4(nume),aaa5(nume),aaa6(nume),aaa7(nume)
      common/cccc/    ntt
      common/poutt/   nst(print_intv)
      
      common/WMLthermal/thermalflag         !!!!!!!,thermalconstraint(nume), Tinit(nume) 
!!!!!      common /WMLthermal2/thermalk(nume),thermalh(nume),
!!!!!     >      thermalRo(nume),thermalcp(nume),thermalx(nume)
!!!!!     >      ,thermalD(nume)

      common/intener/intener(nelemg)
!!!!!!!!!!      common/WMLthermalmatpoly/DijSije(nume),Tele(nume)
!!!!!!!!!!      common /WMLthermali/Telei(nume)           !connect(nume,4),
	  integer thermalflag
	  integer nstep,nnn2,ctr            !!!!!connect,
	  common /nab/ gn_bcc(3,43,43), n_bcc(3,43,43), np(43)
	  common /nab1/ gn_fcc(3,18,18), n_fcc(3,18,18), np1(18)	  
	  common /nab2/ gn_hcpt(3,86,86), n_hcpt(3,86,86), np2(86)
	  common /nab3/gn_hcp(3,87,87), n_hcp(3,87,87), den_im2(87),np3(87)
	  common /aij/ aijhcp(24,87), aijfcc(12,18), aijbcc(24,43),aijhcpt(24,86)
	  common /intgr/ rgen(24), rrecov(87), rintr(87), rintr_n(87),rintr_p(87)
!!!	  common /stressflag/ ElemFractCode(nume),ElemDecayCount(nume)
!!!	  common /gbblock/ gbvec(nume,3), gbflag(nume,3)
!!!!!!!!!!!!!!!	  common /propag/ sigalt(4,40000)
      common/effmod/effmod(nelemg)
	  common /cleavage_plane0/ cleave_n0(3,3)
	  common /test11/ ink
!!!	  common /overlapping/ intersec(4, nume), area_coeff(nume),
!!!     > update_flag
!!	  common /stroverlap/ sflg
	  common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
!!!	  common /crackopen/ ElemDecayed(nume), overlapele(2,nume)
	  common /sigfrac/ sigmacrit0, sigmacrit1, sigmacrit2, &
      sigmacrit3, DecayCount, f_decay, penalty,fractFlag
!!	  common /ovsum/ ovs
	  common /strvect/ stra(nelemg,4)
!!!!	  common /elas_energ/ inie(nume), ecount(nume)
!!!	  common /energbox/ gbco
!!!!!!	  common/couplinganalysis/ TDflag
!!c	  common /gbtranr/ gbtr(14)
      real :: intener, ystart1(67), yprime1(67),  cleave        !!!!!!sigalt,
!!	  integer :: ElemFractCode, ElemDecayCount
      integer ElemDecayed,overlapele,ele,oele
	  integer ecount, DecayCount, fractFlag !!!!TDflag,ovs, 
	  real  tempgb,ngb,taurgb,taugb, gdotgb, den_im2gb !gbvec,
	  real den_imgb, den_mgb, slip_ngb, slip_sgb, thermal_factorgb
	  real bres, cb, sb, cnu, snu, psgb, ugb1, ugb, beff, magvec1
	  real magvec2,bresdot,den_gb(24),den_gbtot,den_ggb(24),dumabc
	  real stra, qe, qp, inie, elas_energy
!!!      real intersec, area_coeff, f_decay
	  real deltarm, deltarim
	  real gbtr(24), gbtr_tot, gbco, gsh, gsh2, shearslip(24)
	  real sigmacrit0, sigmacrit1, sigmacrit2, sigmacrit3
	  integer  inkgb, slipgb(24), update_flag, sflg, numelt!gbflag,
	  dimension ngb(2), taurgb(24), taugb(24),gdotgb(24),den_im2gb(87)
	  dimension den_imgb(24),den_mgb(24),slip_ngb(24,3),slip_sgb(24,3)
	  dimension bres(24,24),cb(24,24),sb(24,24),cnu(24,24),snu(24,24)
	  dimension psgb(24,24),ugb1(24,24),ugb(24),beff(24,24),RetVal(4)
      real DijSij,DijSij_e,temp
      integer ink
      dimension den_im2SQRT(87)
      integer mElemUnloadingCount,mElemSplitCode !1=split
!=======================================================================
! timep      : point in time
! dt		 : time step size
! nstep      : step number (1st = zero)
! nintg      : intergr. point number in an element
! ntt	     : iteration number
!--------------------------
! %% (elements are grouped by: material no., and limitiation on max no. elements/group = 128)
! ink		 : absolute element number
! mft	     : first element in this group (number begins where last element in previous group ends)
! mlt	     : last element in this group  (< = 128)
! nftm1      : first element number in this group, minus '1'
! lft	     : Not used here *used to get mft*
! llt		 : Not used here *used to get mlt*
! matp(ink)  : material number attached to element 'ink'
!=======================================================================
    if (nstep.eq.0) then
        nnn2(ink,nintg)=0
        do ipr = 1, print_intv
            nst(ipr)=ipr*prop(9,1)                             
        end do
    end if
!      return
      !-- no need for the following, it's done in copy elem
!!!!!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!!!!!c
!!!!!!c            update stress for new added element
!!!!!!c
!!!!!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
!!!!!      do oele=1,ovs
!!!!!	      ele=overlapele(1,oele)
!!!!!          if (ElemDecayed(ele)==1 .and. sflg==0) then
!!!!!		      do i=1,4
!!!!!	              sig(i,overlapele(2,oele))=sig(i,overlapele(1,oele))
!!!!!	          end do
!!!!!		  end if
!!!!!	  end do
!!!!!	  sflg=sflg+1
!!!!!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	 !-------------------------------------------------
	 !  Define The Intervals at which Output is Printed
	 !-------------------------------------------------

	  
	 !--------------------------------------------
	 !  Load Stresses from previous Stable Config
	 !--------------------------------------------
    do 5 i = mft ,mlt
        sign1(i)= sig(1,i)!11stress from previous
        sign2(i)= sig(2,i)!22stress from previous
        sign3(i)= sig(3,i)!33stress from previous
        sign4(i)= sig(4,i)!12stress from previous
    !c      write(*,*) 'in matpoly, input sigs',sig(1,i),sig(2,i),sig(3,i),
    !c     1            sig(4,i),i
    !c      write(*,*) 'in matpoly, ds',d1(i),d2(i),d3(i),d4(i)
    5 continue
    
	
!c----------------------------------------------------------------------  
!c----------------------------------------------------------------------  
!c--------------------MAJOR LOOP----------BEGIN-------------------------  
!c----------------------------------------------------------------------  
!c---------------------------------------------------------------------- 
     
    do 10 i = mft, mlt              !loooing over elements in group
       !------------------------------------------------------
       ! ** Loop to update Internal and Field variables, **
       ! **  for Elements in the set: (mft,mlt)          **
       !------------------------------------------------------

        ink = i + nftm1
        mElemSplitCode=EC_GetElemSplit(ink)
        mElemUnloadingCount=EC_GetElemUnloadingCount(ink)
!!!!---------------------------debugging      
!!!      if((ink == 197).and.(mElemSplitCode /= 0)) then
!!!          qwer=0
!!!      endif
!!!!---------------------------debugging      
!!--------------------------------------debugging
!      write(969, *) 'matpoly loop',mft, mlt 
!         matid= matp(ink)
!         write(969, *) ink,matid,nssmat(matid),nssimmat(matid)
!!--------------------------------------debugging
       !--------------------------------------------------
       !  Define deviatoric deformation rate 
       ! (work conjugate to Cauchy Stress, 
       !  i.e. Nominal Stresses rerefenced in Original Config)			!This Claim: http://www.engin.brown.edu/courses/en222/notes_frame.htm
       !--------------------------------------------------		 
        traced =d1(i)+d2(i)+ d3(i)                 !trace of Dij
        !c		 write(*,*) d1(i),d2(i),d3(i),d4(i)
        Dij_dev(1) = d1(i)-traced/3.00				!11 deviatoric of Dij
        Dij_dev(2) = d2(i)-traced/3.00				!22 deviatoric of Dij
        Dij_dev(3) = d3(i)-traced/3.00				!33 deviatoric of Dij
        Dij_dev(4) = d4(i)							!12 deviatoric of Dij

       !-----------------------------------------------
       !  define deviatoric Cauchy Stresses
       !-----------------------------------------------
        press= (sign1(i) + sign2(i) + sign3(i))/3.00 !pressure
        sdev1=sign1(i) - press						!11 deviatoric stress
        sdev2=sign2(i) - press						!22 deviatoric stress
        sdev3=sign3(i) - press						!33 deviatoric stress
        sdev4= sign4(i)							!12 deviatoric stress
        !	     write(*,*) Telei(ink,1)
        nmo=nssmat(matp(ink))
        nim=nssimmat(matp(ink)) 

        if (nstep.ne.nnn2(ink,nintg)) then			!If convergence has occurred, nstep will be higher than the previous step (nnn2)
            ! ------------------------------------------------- 
            ! Store History Variables from iterations  	
            ! =======================================
            deltarm=0.0
            deltarim=0.0
!			 do j = 1, 573
!	            abc(j,ink,nintg) = his(j,ink,nintg)
!	         end do
            abc(1:573,ink,nintg) = his(1:573,ink,nintg)

           ! -------------------------------------------------
           ! Print Output at stipulated intervals
           ! ============
            if (mod(nstep,int(prop(9,1))) == 0) then
    !!!!		         if (thermalflag==0 .or. 
    !!!!     >			     (thermalflag==1 .and. TDflag==1)) then     !if diffusion
                if (thermalflag==0 .or. thermalflag==2) then !if no thermal or diffusion only
                    temp = abc(4,ink,nintg)			!previous temperature
                    if (nstep.eq.0)then
                        temp=tempr
                    end if
    !!!		  	     else if((thermalflag==1 .and. TDflag==0)       !if thermal
    !!!     >				     .or.thermalflag==2) then
                else if(thermalflag==1 .or. thermalflag==3)then       !if there is thermal; thermal only or both
                    temp=CNmanagerGetElemTemp(ink)
    !!!!     		             temp=Tele(ink)
    !!!!		             if (nstep.eq.0) then
    !!!!                         temp=Telei(ink)
    !!!!                     end if
                end if

                call print_result(ink,nintg,slip_n,slip_s,tau_y(ink))
                call print_str(i,ink,matp)
            end if          
        end if

       !---------------------------------------------------------
       ! Load internal variables independent of slip systems
       !---------------------------------------------------------
        angle_Psi    = abc(1,ink,nintg)*pi/180.0	!lattice rotation
        gamma        = abc(2,ink,nintg)			!slip
        elas_energy  = abc(3,ink,nintg)			!reference shear stress, qwu: use it to save elastic energy
!!!!		         if (thermalflag==0 .or. 
!!!!     >			     (thermalflag==1 .and. TDflag==1)) then     !if diffusion
        if (thermalflag==0 .or. thermalflag==2) then !if no thermal or diffusion only
            temp = abc(4,ink,nintg)			!previous temerature
            if (nstep.eq.0)then
               temp=tempr
            end if
!!!		  	     else if((thermalflag==1 .and. TDflag==0)       !if thermal
!!!     >				     .or.thermalflag==2) then
        else if(thermalflag==1 .or. thermalflag==3)then       !if there is thermal; thermal only or both
            temp=CNmanagerGetElemTemp(ink)
!!!!!			 DijSije(ink)=0.0            
!!!!		     temp=Tele(ink)
!!!!		     if (nstep.eq.0) then
!!!!                 temp=Telei(ink)
!!!!             end if    
        end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         temp=0
!		 temp         = max(temp,tempr)				!max of previous or initial temperature
        enthalpy_coef= thermalEnthalpy(matp(ink))			!Enthalpy coefficent now material dependent WML
    !c		 thermal_coef = enthalpy_coef/tempr			!This Line: coefficient = (H/[kT]) for dPm/dt and dPim/dt calculartions
        porosity     = abc(5,ink,nintg)
        plastic_work = abc(6,ink,nintg)

        nssm=nssmat(matp(ink))
       !---------------------------------------------------------
       ! Load internal variables dependent on slip systems on first iteration
       !---------------------------------------------------------
!!!!!$OMP PARALLEL DO       
        do j = 1, nssm								!Loop over system alpha:
            j1= 9 + j
            tau(j)      = abc(j1,ink,nintg)		!This Line: shear stress     
            j2= 33 + j
            gdot(j)     = abc(j2,ink,nintg)		!This Line: shear strain-rate
            j3= 57 + j 
            den_m(j)    = abc(j3,ink,nintg)		!This Line: rho mobile
            j4= 81 + j
            den_im(j)   = abc(j4,ink,nintg)
            den_im2(j) = den_im(j)          		!This Line: rho immobile
            j5= 105 + j 
            slip_n(j,1) = abc(j5,ink,nintg)		!This Line: slip plane normal (x)
            j6= 129 + j 
            slip_n(j,2) = abc(j6,ink,nintg)		!This Line: slip plane normal (y)
            j7= 153 + j
            slip_n(j,3) = abc(j7,ink,nintg)		!This Line: slip plane normal (z)
            j8= 177 + j 
            slip_s(j,1) = abc(j8,ink,nintg)		!This Line: slip direction    (x)
            j9= 201 + j 
            slip_s(j,2) = abc(j9,ink,nintg)		!This Line: slip direction    (y) 
            j10 = 225 + j 
            slip_s(j,3) = abc(j10,ink,nintg)		!This Line: slip direction    (z)
            j11= 270 + j 
            rgen(j) = abc(j11,ink,nintg)
            j12 = 337 + j
            den_gb(j) = abc(j12,ink,nintg)    
            j13 = 549 + j
            shearslip(j) = abc(j13,ink,nintg)  
        end do

        do j = 1, 3
            cleave(1,j) = abc(362+j,ink,nintg)
            cleave(2,j) = abc(365+j,ink,nintg)
            cleave(3,j) = abc(368+j,ink,nintg)
        end do	     

        den_gbtot = abc(362,ink,nintg)

       !---------------------------------------------------------

        if (nssimmat(matp(ink))==18) then
            ctr = 1
            do j = 1, 18
                rrecov(j) = abc(460+j,ink,nintg)
                if (j .ge. 13) then
                    den_im2(j) = abc(397+ctr,ink,nintg)
                    ctr = ctr + 1
                end if
            end do
        end if

        if (nssimmat(matp(ink))==43) then
            ctr = 1
            do j = 1, 43
                rrecov(j) = abc(460+j,ink,nintg)
                if (j .ge. 25) then
                    den_im2(j) = abc(397+ctr,ink,nintg)
                    ctr = ctr + 1
                endif
            end do
        end if

        if (nssimmat(matp(ink))==86) then
            ctr = 1
            do j = 1, 86
                rrecov(j) = abc(460+j,ink,nintg)
                if (j .ge. 25) then
                    den_im2(j) = abc(397+ctr,ink,nintg)
                    ctr = ctr + 1
                endif
            end do
        end if

        if (nssimmat(matp(ink))==87) then
            rrecov(1:87) = abc(460+1:460+87,ink,nintg)
            den_im2(25:87) = abc(397+1:397+1+(87-25),ink,nintg)
    !!!!		     ctr = 1
    !!!!		     do j = 1, 87
    !!!!		         rrecov(j) = abc(460+j,ink,nintg)
    !!!!		         if (j .ge. 25) then
    !!!!		             den_im2(j) = abc(397+ctr,ink,nintg)
    !!!!		             ctr = ctr + 1
    !!!!		         endif
    !!!!		     end do
        end if

      !-----------------------------------------------
      !            Load Material Properties
      ! ==============================================
        elem_mat_no = matp(ink)		 
        do j = 1, nssm
           xmhigh(j)        = prop(4,elem_mat_no)	!load strain rate sensitivity parameter per system (alpha)
           ref_gamma_dot(j) = prop(7,elem_mat_no)	!load reference strain rate per system (alpha)
        end do

        gdotcr   =  prop(8,elem_mat_no)*1.E29		!load critical strain rate per system (alpha), DEACTIVATED, WML
    !         gdotcr   =  prop(8,elem_mat_no)	         !OLD WAY
        gsh      =  Y_modulus(ink)/(2.0*(1.0 +possion_ratio(ink)))!This Line: Modulus of rigidity
        twomu    =  2.0*gsh							!This Line: Lame parameter mu, multiplied by '2'

        proppp   =  Y_modulus(ink)*possion_ratio(ink)/((1.0 + &
                 possion_ratio(ink))*(1.-2.*possion_ratio(ink)))!This Line: Lame parameter Lambda CHANGED WML 9/19/08
        alamdt   =  proppp*dt						!This Line: Lame parameter Lambda, multiplied by 'dt'

      !----------------------------------------------- 
      ! These are the stiffness values used in stfnf and pnbtcb				
      ! ==============================================
        dsave(1,1,i) = prop(10,elem_mat_no)
        dsave(1,2,i) = prop(11,elem_mat_no)
        dsave(1,3,i) = prop(11,elem_mat_no)
        dsave(1,4,i) = 0.
        dsave(2,1,i) = prop(11,elem_mat_no)
        dsave(2,2,i) = prop(10,elem_mat_no)
        dsave(2,3,i) = prop(11,elem_mat_no)
        dsave(2,4,i) = 0.
        dsave(3,1,i) = prop(11,elem_mat_no)
        dsave(3,2,i) = prop(11,elem_mat_no)
        dsave(3,3,i) = prop(10,elem_mat_no)
        dsave(3,4,i) = 0.
        dsave(4,1,i) = 0.
        dsave(4,2,i) = 0.
        dsave(4,3,i) = 0.
        dsave(4,4,i) = prop(25,elem_mat_no)!?????

      !-----------------------------------------------

        do j = 1, nssm
         !-----------------------------------------------
         ! Slip: Drag controlled or Thermally Activated ? for the most part, commented out above
         ! ==============================================
            if (gdot(j).ge.gdotcr) then
                xmhigh(j) = 1.00
                ref_gamma_dot(j) = gdotcr
            endif
         !-----------------------------------------------
        end do

       !------------------------------------------------------------
       ! Compute Tau Ref 
       ! ======================
!c         taur = tau_y(ink)													 
!c		 do j = 1, nssm
!c			b_v=bv(matp(ink),j)						!added for variable Burger's vectors, WML 91009
!c            taur = taur + 0.5*gsh*b_v*sqrt(den_im(j))	!This Loop: (Self Hardening Model)
!c         end do
        rnu=prop(40,elem_mat_no)
!!!!!!!         write(*,*)prop(40,elem_mat_no),rnu
        thermal_factor = (tempr/temp)**rnu			!This Line: (Thermal Multiplier)
!c         taur           = taur * thermal_factor		!This Line: Update Tau Reference, commented out for quasi-static
       !------------------------------------------------------------		 


!!!!!$OMP PARALLEL DO       
        do j = 1, nssm
        !------------------------------------------------------------
        ! Pij and Omega_ij for:
        !          Tau, gdot, plastic Deformation, and plastic Spin  SYMMETRIC SCHMID TENSOR
        !                  **** Also for Pij_Dijdev **** 
        ! =========================================================
            p(j,1) = (slip_s(j,2)*slip_n(j,2))              !11 plastic symmetric		!These 3 Assignments: Pij,  Zikry & Nasser, 1990 pg 217
            p(j,2) = (slip_s(j,3)*slip_n(j,3))              !22 plastic symmetric
            p(j,3) = 0.5 * (slip_s(j,2)*slip_n(j,3) +slip_s(j,3)*slip_n(j,2))    !12 plastic symmetric
            p(j,4)=(slip_s(j,1)*slip_n(j,1))                !33 plastic symmetric
            w12(j) = 0.5 * (slip_s(j,2)*slip_n(j,3) - slip_s(j,3)*slip_n(j,2))   !These 2 Assignments: Wij,  Zikry & Nasser, 1990 pg 217
                                     !plastic asymmetric
            w21(j) = - w12(j)                               !plastic asymmetric
        !-------------------------------------------------------------
        end do	


!!!!!!!!$OMP PARALLEL DO       
       do j = 1, nssm
       !------------------------------------------------------------
       ! Calculated for the Computation of Tau_dot						
       ! ========		
        !This loop: Pij(alpha)*Dij',  Zikry & Nasser, 1990 pg 218
        !Calculates time derivative of resolved shear stress
        !UPDATED FOR (3,3) term, WML 9/19/08 FOR OLD WAY, JUST COMMENT OUT
            Pij_Dijdev(j) = p(j,1)*Dij_dev(1)+p(j,2)*Dij_dev(2)+ &
                            2*p(j,3)*Dij_dev(4)+p(j,4)*Dij_dev(3)	
      !-------------------------------------------------------------
        end do

      !-------------------------------------------------------------
      ! Slip System Dependent taur
      !-------------------------------------------------------------
        den_im2SQRT=0.0
        den_im2SQRT(1:87)=sqrt(den_im2(1:87))*gsh
        if (nssimmat(matp(ink)).eq.18) then
            do j = 1, nssm
                b_v=bv(matp(ink),j)						!added for variable Burger's vectors, WML 91009
                taur(j) = tau_y(ink)
                do k = 1, 18
                    taur(j) = taur(j) + aijfcc(j,k)*b_v*den_im2SQRT(k)
    !!!!     >                gsh*sqrt(den_im2(k))
                end do
                taur(j) = taur(j)*thermal_factor
            end do
        end if

        if (nssimmat(matp(ink)).eq.43) then
            do j = 1, nssm
                b_v=bv(matp(ink),j)	
                taur(j) = tau_y(ink)
                do k = 1, 43
                    taur(j) = taur(j) + aijbcc(j,k)*b_v*den_im2SQRT(k)
    !!!!     >                gsh*sqrt(den_im2(k))
                end do
                taur(j) = taur(j)*thermal_factor
            end do
        end if

        if (nssimmat(matp(ink)).eq.86) then
            do j = 1, nssm
                b_v=bv(matp(ink),j)	
                taur(j) = tau_y(ink)
                do k = 1, 86
                    taur(j) = taur(j) + aijhcpt(j,k)*b_v*den_im2SQRT(k)
    !!!!     >                gsh*sqrt(den_im2(k))
                end do
                taur(j) = taur(j)*thermal_factor
            end do
        end if

        if (nssimmat(matp(ink)).eq.87) then
            do j = 1, nssm
                b_v=bv(matp(ink),j)	
                taur(j) = tau_y(ink)
                do k = 1, 87
                    taur(j) = taur(j) + aijhcp(j,k)*b_v*den_im2SQRT(k)
    !!!!     >                gsh*sqrt(den_im2(k))
                end do
                taur(j) = taur(j)*thermal_factor
            end do
        end if

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                  dislocation-GB interaction
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc		 
!!!!!!!$OMP PARALLEL DO  
        ugb(1:24)=0.0
        slipgb(1:24)=0.0
        bres(1:24,1:24)=0.0


       gbtr_tot=0.0
    !c         call dsolve(dt,timexx,nssm)   ! line for updating using rk5, implicit algorithm
       if (mElemUnloadingCount==0) then
           call dsolve(dt,timexx,nssm)   ! line for updating using rk5, implicit algorithm
       elseif (mElemUnloadingCount > 0) then
           tau(1:nssm)=1.0
       end if
        do j = 1, nssm
            gbtr(j)=0.0
            b_v=bv(matp(ink),j)
            !------------------------------------------------------------
            ! Tau(alpha) = Pij(alpha)Sij
            ! ==========================
    !c             tau(j)  = p(j,1)*sdev1 + p(j,2)*sdev2 + 2*p(j,3)*sdev4		!This Line: Cauchy Tau(alpha), Zikry & Nasser, 1990 pg 217   commented out b/c using dsolve
            ! Strain Rate: Power Law
            ! ======================	 
            gdot(j) = ref_gamma_dot(j)*(abs(tau(j)/taur(j))**(xmhigh(j)-1.))*(tau(j)/taur(j)) !This Line: strain rate(alpha), Zikry & Nasser, 1990 pg 217
                                         !This Line: strain rate(alpha), Zikry & Nasser, 1990 pg 217
    !     >         *exp(-gbco*0.5*gsh**2*b_v**3*ugb(j)/(1.38E-23*temp))
            if (mElemUnloadingCount == 0) then			 
                shearslip(j)=shearslip(j)+dt*abs(gdot(j))
            end if	
        end do

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c		       Update Plastic Deformation and Spin Tensors
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        D11_p    = 0                 !11 plastic deformation 
        D22_p    = 0                 !22 plastic deformation
        D12_p    = 0                 !12 plastic deformation
        D33_p    = 0                 !33 plastic deformation
        spin_p12 = 0                 !12 plastic spin tensor
        do j=1,nssm
        ! Sum on (alpha)
           D11_p	 = D11_p	+ p(j,1)*gdot(j)						!These 3 Lines: Plastic Deformation Rate, Zikry & Nasser, 1990 pg 216
           D22_p	 = D22_p	+ p(j,2)*gdot(j)
           D12_p	 = D12_p	+ p(j,3)*gdot(j)
        !   D33_p    = D33_p        + p(j,4)*gdot(j)                   !Commented back out WML 7910
           spin_p12 = spin_p12 + w12(j)*gdot(j)						!This Line: Plastic Spin Rate, Zikry & Nasser, 1990 pg 217
        end do	

        spin_p21= -spin_p12
       !-------------------------------------------------------------
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c           Update Slip Normals and Directions from:
!c                                  - Elastic Spin Tensor
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	 
        spin_21= spin(i)                              !Total spin
        spin_e21= spin_21 + spin_p12                  !Elastic Spin
        spin_e12= -spin_e21       

        if (mElemUnloadingCount==0) then
    !!!!!!!!!!$OMP PARALLEL DO       
            do j = 1 , nssm
                n_dot1=spin_e21*slip_n(j,2)		           !rate of change of slip normals				!These 4 Lines: n_dot and s_dot, Zikry & Nasser, 1990 pg 217
                n_dot2=spin_e12*slip_n(j,3)
                s_dot1=spin_e21*slip_s(j,2)                 !rate of change of slip directions
                s_dot2=spin_e12*slip_s(j,3) 
                slip_n(j,3)=slip_n(j,3) + n_dot1*dt         !These 4 Lines: n and s updated
                slip_n(j,2)=slip_n(j,2) + n_dot2*dt
                slip_s(j,3)=slip_s(j,3) + s_dot1*dt
                slip_s(j,2)=slip_s(j,2) + s_dot2*dt         ! ******* I NEED TO TAP THESE OFF, TO OUTPUT EULER ANGLES
            end do
            do j = 1, 3
                n_dot1=spin_e21*cleave(j,2)
                n_dot2=spin_e12*cleave(j,3) 
                cleave(j,3)=cleave(j,3) + n_dot1*dt         !These 4 Lines: n and s updated
                cleave(j,2)=cleave(j,2) + n_dot2*dt
            end do
        end if

        call unit_vector(3,nss,slip_n,slip_s)
        call unit_vector(3,3,cleave,cleave_n0)							!This Line: Normalize slip vectors to unity:

      !-------------------------------------------------------------

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
!c           Update  Cauchy Stresses from Objective Rates		!Next 7 Lines: Zikry & Nasser, 1990, pg 217
!c           This follows the work done by WML and Khalil to show that we need sigma dot and not sigma hat to update the stresses
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        o12=       -dt*spin_e12*(sdev1 - sdev2)	           !Cauchy stress rate derived from Jaumann rate corotational with elastic lattice distortion	!OBJECTIVE COMPONENTS ADDED TO MAKE STRESS RATES OBJECTIVE
        o11= 2.00*dt*spin_e12*(sdev4)
        o22= -o11		                                       !o33=0, this is correct, no 13 terms	 
        ssdev4 = sdev4+twomu*dt*(Dij_dev(4)-D12_p) + o12      !12 deviatoric cauchy stress
        ssdev1 = sdev1+twomu*dt*(Dij_dev(1)-D11_p) + o11      !11 deviatoric cauchy stress
        ssdev2 = sdev2+twomu*dt*(Dij_dev(2)-D22_p) + o22      !22 deviatoric cauchy stress
        ssdev3 = sdev3+twomu*dt*(Dij_dev(3)-D33_p)            !33 deviatoric cauchy stress
        if (mElemUnloadingCount==0) then
            DijSij = ssdev1*D11_p + ssdev2*D22_p+2.00*ssdev4*D12_p+ssdev3*D33_p		       !This Line: StressPower, Zikry: 1992-1993, pg 275
                                                 !REMOVED ABS Values WML 9/19/08

            DijSij_e = ssdev1*(Dij_dev(1)-D11_p)+ ssdev2*(Dij_dev(2)-D22_p)+ &
                        2.00*ssdev4*(Dij_dev(4)-D12_p)+ssdev3*(Dij_dev(3)-D33_p)

            effmod(i)=1./3.*((twomu/3.0*dt+alamdt)/dt*3.+2.*((Dij_dev(1)* &
                    twomu*(Dij_dev(1)-D11_p)+Dij_dev(2)*twomu*(Dij_dev(2)-D22_p)+ &
                    Dij_dev(3)*twomu*(Dij_dev(3)-D33_p)+2.*Dij_dev(4)* &
                    twomu*(Dij_dev(4)-D12_p))/(Dij_dev(1)*Dij_dev(1)+Dij_dev(2)* &
                    Dij_dev(2)+Dij_dev(3)*Dij_dev(3)+2.*Dij_dev(4)*Dij_dev(4))))


            press    = press  + (twomu/3.0*dt+alamdt)*traced         
            if ((Dij_dev(1)*Dij_dev(1)+Dij_dev(2)*Dij_dev(2)+Dij_dev(3)*Dij_dev(3)+ &
               2.*Dij_dev(4)*Dij_dev(4)).eq.0.) then
               effmod(i)=prop(10,matp(ink))
            end if
            if (effmod(i).lt.0.) effmod(i)=0.
                  !This Line: Update Hydrastatic Stresses
            sign1(i) = ssdev1 + press                         !These 8 Lines: Update Total Stresses	
            sign2(i) = ssdev2 + press 
            sign3(i) = ssdev3 + press
            sign4(i) = ssdev4
!!!!!!!!!!!!!!!!!!!!!		 else if(ElemFractCode(ink)==1) then
        else if((mElemUnloadingCount>0).and.(mElemUnloadingCount <= EC_DecayCount)) then
!!!           doptimizedConst=(inie(ink)**(mElemUnloadingCount**f_decay))
            doptimizedConst=1- mElemUnloadingCount / EC_DecayCount
            sign1(i) = sig(1,i)*doptimizedConst
            sign2(i) = sig(2,i)*doptimizedConst
            sign3(i) = sig(3,i)*doptimizedConst
            sign4(i) = sig(4,i)*doptimizedConst
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         else if(ElemFractCode(ink)==2) then
!!!!!!!!!!!!!!!         else if(EC_GetElemSplit(ink)>0) then
        else if(mElemUnloadingCount > EC_DecayCount) then
!!!!!!!		 else if(mElemSplitCode>0) then !-- already decayed
            sign1(i) = 0.0
            sign2(i) = 0.0
            sign3(i) = 0.0
            sign4(i) = 0.0
        end if
        sig(1,i) = sign1(i)              !11 stress updated
        sig(2,i) = sign2(i)              !22 stress updated
        sig(3,i) = sign3(i)              !33 stress updated
        sig(4,i) = sign4(i)              !12 stress updated
!!!!!!		 sigalt(1,ink) = sig(1,i)
!!!!!!	     sigalt(2,ink) = sig(2,i)
!!!!!!	     sigalt(3,ink) = sig(3,i)
!!!!!!	     sigalt(4,ink) = sig(4,i)
        RetVal(1:4)=sig(1:4,i)
        CALL CNmanager_Set_sigalt(ink,RetVal)

        intener(i)=sig(1,i)*d1(i)+sig(2,i)*d2(i)+sig(3,i)*d3(i)+2.*sig(4,i)*d4(i)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                Update Temperature and Plastic Work
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if (mElemUnloadingCount==0) then
!!!!		         if (thermalflag==0 .or. 
!!!!     >			     (thermalflag==1 .and. TDflag==1)) then     !if diffusion
            if (thermalflag==0 .or. thermalflag==2) then !if no thermal or diffusion only
!!!!!                    dKa_RaoCp= rMatTH_x(ink)/
!!!!!     >                   (rMatTH_Ro(ink)*rMatTH_cp(ink))
!!!!!!!!!                    dKa_RaoCp= thermalx(ink)/
!!!!!!!!!     >                   (thermalRo(ink)*thermalcp(ink))
!!!!!	             temp         = temp         + dt*DijSij*dKa_RaoCp				!This Line: Adiabatic Temp Update, Zikry: 1992-1993, pg 275, commented out for quasi-static
                CALL CNmanagerCal_AdiabaticTemp(ink,DijSij,temp)
!!!		  	     else if((thermalflag==1 .and. TDflag==0)       !if thermal
!!!     >				     .or.thermalflag==2) then
            else if(thermalflag==1 .or. thermalflag==3) then      !if there is thermal; thermal only or both
!!!!!		         DijSije(ink)=DijSij*thermalx(ink)
!!!!!		         DijSije(ink)=DijSij*rMatTH_x(ink)
                CALL CNmanagerSet_DijSij(ink,DijSij)
            end if
            plastic_work = plastic_work + dt*DijSij						!This Line: Update plastic work from stress power
            elas_energy = elas_energy + dt*DijSij_e
        end if
       !-------------------------------------------------------------

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                Compute Rho_m(alpha) and Rho_im(alpha) 				 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if ((nstep /= 0).and.(ntt == 1)) then
           timexx = timep - dt
        else
           timexx = timep
        end if
        do j = 1, nim
           b_vvec(j)=bv(matp(ink),j)											!added for variable Burger's vector WML 91009
        end do
        nnne = ink													!This Line: Pass apposite element number


        if (mElemUnloadingCount==0) then
           call dsolve1 (dt,timexx,nmo+nim)									!This Line: Call Rho_mobile & Rho_immobile Updating subroutine

    !!!!!!$OMP PARALLEL DO       
           do j = 1, nmo
               ystart1(j) = den_m(j)
           end do
    !!!!!!$OMP PARALLEL DO       
           do j = 1, nim
               ystart1(nmo+j) = den_im2(j)	
           end do

           call derivs1post(ystart1,yprime1,nmo+nim)

    !!!!!!$OMP PARALLEL DO       
           do j = 1, nmo
               rgen(j) = rgen(j) + dt*yprime1(j)
           end do
    !!!!!$OMP PARALLEL DO       
           do j = 1, nim
               rrecov(j) = rrecov(j) + dt*yprime1(nmo+j)
           end do
        end if

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c              Update History /Internal Variables
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        dgamma = 0.00													!This Line: Change in shear slip. Initialized at "zero"
!!!!$OMP PARALLEL DO       
        do j = 1, nssm 
        !	  dgamma = dt*abs(gdot(j)) + dgamma							!This Line: Sum Total CHANGE in Shear Slip COMMENTED OUT WML 7910 
            j1=9+j   
            his(j1,ink,nintg)  = tau(j)								!This Line: Shear Stress per system (alpha) 
            j2=33+ j
            his(j2,ink,nintg)  = gdot(j)								!This Line: Shear rate per system (alpha) 
            j3=57+j 
            his(j3,ink,nintg)  = den_m(j)								!This Line: Rho_m per system (alpha)
            j4=81+j
            his(j4,ink,nintg)  = den_im(j)							!This Line: Rho_im per system (alpha)
            j5=105+j
            his(j5,ink,nintg)  = slip_n(j,1)							!These 3 Assignments: New Slip Plane Normal for system (alpha)
            j6=129+j
            his(j6,ink,nintg)  = slip_n(j,2)
            j7=153+j
            his(j7,ink,nintg)  = slip_n(j,3)
            j8=177+j
            his(j8,ink,nintg)  = slip_s(j,1)							!These 3 Assignments: New Slip Direction for system (alpha)
            j9=201+j
            his(j9,ink,nintg)  = slip_s(j,2)
            j10=225+j
            his(j10,ink,nintg) = slip_s(j,3)
            j11 = 270 + j
            his(j11,ink,nintg) = rgen(j)
            j12 = 337 + j
            his(j12,ink,nintg) = den_gb(j)
            j13 =371+j
            his(j13,ink,nintg) = gbtr(j)                       !qwu save gbtr
            j14 =549+j                                       
            his(j14,ink,nintg) = shearslip(j)                  !qwu save shearslip for each slip system
        end do

        do j = 1, 3
            his(362+j,ink,nintg) = cleave(1,j)
            his(365+j,ink,nintg) = cleave(2,j)
            his(368+j,ink,nintg) = cleave(3,j)
        end do	  

        his(362,ink,nintg) = den_gbtot
        his(396,ink,nintg) = gbtr_tot

        if (mElemUnloadingCount==0) then
            angle_Psi              = angle_Psi + spin_e21*dt										
            his(1,ink,nintg) = angle_Psi*(180.0/pi)								!This Line: Angle through which Slip Planes & Directions have rotated,
            dgamma=0.666666667*dt*sqrt(D11_p**2+D22_p**2+2*D12_p**2)
        end if
        his(2,ink,nintg) = gamma + dgamma								!This Line: New Shear Slip,
        his(3,ink,nintg) = elas_energy										!This Line: New Reference Shear, qwu: save elastic energy
        his(4,ink,nintg) = temp										!This Line: New Temperature

        his(5,ink,nintg) = porosity									!This Line: New Porosity
        his(6,ink,nintg) = plastic_work								!This Line: Updated Plastic Work
        his(7,ink,nintg) = spin_e12
        his(8,ink,nintg) = spin_p12
        if (nssimmat(matp(ink)) == 18) then
            ctr = 1
            do j = 1, 18
                his(460+j,ink,nintg) = rrecov(j)
                if (j .ge. 13) then
                    his(397+ctr,ink,nintg) = den_im2(j)
                    ctr = ctr + 1
                end if
            end do
        end if

        if (nssimmat(matp(ink)) == 43) then
            ctr = 1
            do j = 1, 43
                his(460+j,ink,nintg) = rrecov(j)
                if (j .ge. 25) then
                    his(397+ctr,ink,nintg) = den_im2(j)
                    ctr = ctr + 1
                end if
            end do
        end if

                if (nssimmat(matp(ink)) == 86) then
            ctr = 1
            do j = 1, 86
                his(460+j,ink,nintg) = rrecov(j)
                if (j .ge. 25) then
                    his(397+ctr,ink,nintg) = den_im2(j)
                    ctr = ctr + 1
                end if
            end do
        end if



        if (nssimmat(matp(ink)) == 87) then
            his(460+1:460+87,ink,nintg) = rrecov(1:87)
            his(397+1:397+1+(87-25),ink,nintg) = den_im2(25:87)
    !!!!		     ctr = 1
    !!!!		     do j = 1, 87
    !!!!		   	     his(460+j,ink,nintg) = rrecov(j)
    !!!!		         if (j .ge. 25) then
    !!!!		   	         his(397+ctr,ink,nintg) = den_im2(j)
    !!!!		   	         ctr = ctr + 1
    !!!!		         endif
    !!!!		     end do
        endif

        nnn2(ink,nintg)  = nstep										!This Line: Store current Step number for element
!!!!!!		 if (thermalflag.eq.2) then										!isothermal
!!!!!!		     Tele(ink,1)=temp
!!!!!!		 end if
       !-- for debugging
!!!		 if(ink==nelec) then
!!!		     qepower=DijSij_e
!!!			 qppower=DijSij
!!!	     end if
10  continue
!c----------------------------------------------------------------------  
!c----------------------------------------------------------------------  
!c--------------------MAJOR LOOP-------END------------------------------  
!c----------------------------------------------------------------------  
!c----------------------------------------------------------------------      
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c      			Update Cauchy Stresses 				    	
!c     Good to go for Cartesian Equilibrium Equations    CAUCHY STRESSES STORED BY MATPOLY
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!$OMP PARALLEL DO       
    do 440 i = mft, mlt
        dels11(i) = sign1(i) - sig(1,i)
        dels22(i) = sign2(i) - sig(2,i)
        dels33(i) = sign3(i) - sig(3,i)
        dels12(i) = sign4(i) - sig(4,i)
        sig11(i)  = sig(1,i)
        sig22(i)  = sig(2,i)
        sig33(i)  = sig(3,i)
        sig12(i)  = sig(4,i)
        sig11s(i) = sig(1,i)
        sig22s(i) = sig(2,i)
        sig33s(i) = sig(3,i)
        sig12s(i) = sig(4,i)
440  continue
!c      if(update_flag==1) then
!c          write(*,*) (effmod(i), i=1,2)
!c	   end if
    return
    end





