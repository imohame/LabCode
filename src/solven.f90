      subroutine solven
!c     iterative matrix solvers
!c     method 1 = BFGS
!c     method 2 = Broyden
!c     method 3 = DFP
!c     method 4 = Davidson Symmetric
!c     method 5 = modified Newton
!c     implicit double precision (a-h,o-z)                                    dp
!c
      use CN_Objects_manager
      use mod_dtimeSpecs
      use mod_file_units
      real*8 hed                                                        !vax750

      common/hourglass/fhg(40000,8),fhghis(40000,8)
      common/hourglass2/hgsstore(40000),hgshis(40000)
      common/totalenergy/totenerstore(40000),totenerhis(40000),inertener(40000)
      common/hgenergy/hgenerstore(40000),hgenerhis(40000)
      common/irdmp1/lendr,lenhr,irt,trt,ityprs
      common/bk00/ &
            k01,k02,k03,k04,k05,k06,k07,k08,k09,k10,k11,k12, &
            k13,k14,k15,k16,k17,k18,k19,k20,k21,k22,k23,k24, &
            k25,k26,k27,k28,k29,k30,k31,k32,k33,k34,k35,k36, &
            k37,k38,k39,k40,k41,k42,k43,k44,k45,k46,k47,k48, &
            k49,k50,k51,k52,k53,k54,k55,k56,k57,k58,k59,k60, &
            k61,k62,k63,k64,k65,k66,k67,k68,k69,k70,k71,k72, &
            k73,k74,k75,k76,k77,k78,k79,k80,k81,k82,k83,k84, &
            k85,k86,k87,k88,k89,k90,k91,k92,k93,k94,k95,k96
      common/bk02/ioofc,iphase,imass,lpar(9)
      common/bk05/ifil,iadd,maxsiz,head(12)
      common/bk06/nprnt,mprint,itmpop,numelt,jprint,idump,locstr
      common/bk07/mbfc,nelpg,hed(12)
      common/bk08/kprint,nstep,ite,ilimit,newstf
      common/bk09/maxref,rhsn,rhsvn,cvtl,iteref,ectl,tolls
      common/bk11/cnwmk(2),iequit,iprint,isref
      integer iequit,iprint,isref
      common/bk12/ntlen
      common/bk14/lfna(15),lfnt(6)
      common/bk15/cpuio(36),cpuip(36)
      common/bk17/dn1,dn2,nwebuf,ntime,numnp,neq,ibar,mthsol
      common/bk23/itemp,itherm,irtin

      common/bk26/dt,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
      common/bk27/nlcur,nptst,nthpy,nthpz,nthps,xmy,xmz,xms,nload,nptm
      common/bk32/nsref,nequit,time,timep,lprint,nprint
      common/bk33/irfreq,krfreq,iress
      common/bk34/aa(6001)
      common/slar3/nsl,nsntl,nmntl,nslnmx,sltol,slhrd
      common/fissn0/maxneq,mwspac,ntpe0,ntpe1,nfissl(3)
      common/fissn1/melemt,nnns,ntpe2,n2g,llls
      common/fissn2/nwpblk,numblk,mwsusd,mxnepb,maxch,matpr, mench,ifa(2)
      common/riksw2/rlnew,alfa0,dsx,iteopt,idctrl,riksf,numspu,mthunl
      common/automt/dtmin,dtmax,mxback,termtm
!!!      character*8 lkname
      logical rezone
      common/rezone/rezone,nrzn,nctr,irzcnt,nrezon
      common/total/itrlas,irflas,irhlas,itrtot,irftot,irhtot
!c     common/psa2/inkopt,itpopt,icstep,ncstep,ctime,cdltim                   pl
      common/double/iprec,ncpw,unit
      common/cccc/ntt
!!!      common/kkk/km,kkjj
      common/pvari/iaa,ibb
      common/meichu/modeltype

      common /main_block/ a(1)

      common /wblock7/ slip_n0(1000,nss,3), slip_s0(1000,nss,3)
      common /wblock8/ abc(573,nume,4),his(573,nume,4)
      common /wblock9/ slip_n0_t(1000,nss,3), slip_s0_t(1000,nss,3)
      common /wblock20/ mat_type(nume)
      common /wblock12/ Y_modulus(nume),possion_ratio(nume),tau_y(nume)
!!!!!      common/WMLthermal/thermalflag  !!!!,thermalconstraint(nume),Tinit(nume),Rqold(nume)
      common/hgstress/hgstress1store(40000),hgstress2store(40000),hgstress1his(40000),hgstress2his(40000)
	  common/mbsize/numelt2, numnp2, neq2
	  common /overlapping/ intersec(4, nume), area_coeff(nume),update_flag
	  common /ovsum/ ovs
!!	  common /excflag/ excf
!!!!	  common /tipvelocity/ ncrack, nelefail(1000),tipelenum(1000,nume)
!!!!!      character*8 namef
!!!!!      logical term
      logical dynlnk
      data isw3on/1/
      data isw7on/0/
      integer  update_flag, ovs
!!      integer excf       !!thermalflag,
!!!!	  integer ncrack, nelefail, tipelenum 
	  real intersec, area_coeff
      real timesteps,time0,time1,time2,elpt1,timebase,timetotal
      real*8 tim(2), tim2(2)
      real*4 elapsed, wdiff  
      external wdiff       
	  integer ierr
      INTEGER :: count_rate, count_max,t1,t2
        INTEGER   nb_ticks_initial ! initial value of the clock tick counter
        INTEGER   nb_ticks_final   ! final value of the clock tick counter
        INTEGER   nb_ticks_max     ! maximum value of the clock counter
        INTEGER   nb_ticks_sec     ! number of clock ticks per second
        INTEGER   nb_ticks         ! number of clock ticks of the code
        REAL :: elapsed_time       ! real time in seconds
       integer connect,idt
       real*8 DtCondDiff
 
!!!       call CNTestPrint()
!!!       return
      if (modeltype.ge.4.and.modeltype.le.8) then
        call initial_values(a(k12),numelt)
      endif

!!!!c.....TECPLOT write
!!!      write(126,4)
    4 format(5x,'ZONE')
      write(5701,41)
   41 format(5x,'ZONE')

      ovs=0
!!!!	  ncrack=0
!!!!	  excf=1
	  update_flag=0
      
      !-- this for writing the forplot disp/stress every 100 step
      iaa = 0
      ibb = 100
!c----------------------
      dynlnk=.false.
      nn    = 0
      a(k08+numelt2:k08+numelt2+numelt2)=1.0
      
! ---- nstep starts as -1
      if (nstep >= 0) go to 40                                               

      nprint=jprint+nstep
      if (iprint.lt.1000000) then
        lprint=iprint+nstep
      endif

      call CPU_TIME(timebase)   
      CALL SYSTEM_CLOCK(COUNT_RATE=nb_ticks_sec,COUNT_MAX=nb_ticks_max)
      CALL SYSTEM_CLOCK(COUNT=nb_ticks_initial)
      
   10 nstep=nstep+1                                                          
      nbcku=0
      ntt=0

      iPrintOutputFlag=mod(nstep,iprint)
      if(iPrintOutputFlag > 0) then ! there is a reminder  ex. 2200/1000 = 200 but 3000/1000 =0
        iPrintOutputFlag=0 ! if nstep is not divisible with iprint then do not print
      else
        iPrintOutputFlag=1 !  print output
      endif
      call CPU_TIME(time0)     

      krfreq=krfreq+1
      if (irfreq.le.krfreq) krfreq=0

   30 continue

   nsref=isref-1
      call chgint

   40 timep=time

   time =time+dt

   nstp1=nstep+1
      write(iFU_times_out,*)'- new step,time,iPrintOutputFlag',nstep,time,iPrintOutputFlag

      irflas=0
      nsref =nsref+1
      nequit=nequit+1
      newstf=isref-nsref
      iter  =iequit-nequit
      if (newstf.eq.0) nsref =0
      if (nstep.le.0)  newstf=0
      if (iter.eq.0)   nequit=0
      if (newstf.eq.0) then
        irftot=irftot+1
        irflas=1
      endif
!c      ITERATION SOVER
!c     compute element stiffness matrices, stress divergence, and rhs
!c
   80 call CPU_TIME(time1)     
      call elstf (a(k08),a(k08+numelt2))
      call CPU_TIME(time2)     
      elpt1=time2-time1     
      write(iFU_times_out,*)'---- time of elstf code=',elpt1
!!!!!!!!    ----------------------------testing  
!!!!!!        write(*,*)a(k17:k17+neq)
!!!!!!!        stop
      
      if(update_flag==0) then
	      call total_Qe(a(k03),a(k04),a(k02))   ! calculate total elastic energy
	  end if
      
         if (nstep.eq.ntime) dynlnk=.true.
         if (time.gt.termtm) dynlnk=.true.
      
      if (dynlnk) then
         go to 50
      endif
      

      call blkcpy (a(k19),a(k21),neq)
      call bsolvr (a(k19),a(ntlen),newstf+2,6)
      call blkcpy (a(k19),a(k20),neq)

      if (iter > 0) go to 50

      ierr=0
        call CPU_TIME(time1)  
!!!!!!        write(*,*)k19,k20,k20,k21,neq2,time
!!!!!!        write(*,*)a(k19:k19+10)
        call quasin (a(k19),a(k20),a(k21),a(k22+neq2),a(k18),a(k22),time,ierr)
        call CPU_TIME(time2)     
        elpt1=time2-time1     
        write(iFU_times_out,*)'---- time of quasin code=',elpt1

        if (ierr.eq.1) then ! if the quasin did not converge then terminate
            call bye(2)
        endif

!c
!c     write high speed printer file, model geometry, stress
!c
   50 lprint=lprint+1
      kprint=iprint-lprint
      if (kprint.gt.0) go to 60 
      
      call prtout
!c
!c     write plot data into binary file
!c
   60 nprint=nprint+1
      mprint=jprint-nprint
      
      
      write(iFU_times_out,*)'output-> iprint,jprint,kprint,nstep',iprint,jprint,kprint,nstep
      write(iFU_times_out,*)'output-> lprint,mprint,nprint',lprint,mprint,nprint
      if (mprint.gt.0) go to 70
      
      call prtout
!c
!c     update displacements, velocities, and accelerations
!c

   70 continue

      call CPU_TIME(time1)
      call updat (a(k18),a(k55),a(k56),a(k19))
!!!!!      updat (u     ,udt   ,udtt  ,ui)
!!!	  call edge(a(k03),a(k04),a(k08))
!!!!!!!!      edge(y     ,z     ,matp)
      
      call propagate(a(k03),a(k04),a(k02),a(k57),a(k18), a(k08),numelt,nstep)
!!!!!      propagate(y     , z    , ix   , id   , u    , matp  ,numelt,nstep)
      call overlap(a(k03),a(k04),a(k02),a(k08),a(k57),a(k18),a(k20), a(k07), a(k09))
!!!!!!!!!e overlap(y     , z    , ix   , matp , id   , u    , usi  , freep , ym)
     
     
	  if (update_flag==1) then
		  go to 80
	  end if
	  
	  call GNDloop(a(k03),a(k04),a(k02),a(k57),a(k18),numelt,numnp)
      call CPU_TIME(time2)     
      elpt1=time2-time1     
      write(iFU_times_out,*)'---- time of fracture code=',elpt1
      
!!!!!!!!!c     this is my thermal FEM code WML 82010
!!!!    if (thermalflag.eq.1) then
        call CPU_TIME(time1)
!!!!!!----- update the total number of elements and nodes -- just in case there is cracks
        ElemCountAct=numelt
        NodeCountAct=numnp
        call CN_UpdateElemConnet(a(k02))
         
!!!!!!!!!!!!!        call thermal(a(k02),a(k03),a(k04),a(k18),numelt,numnp,a(k57))
!!!!!!!!!!!!!!            thermal(ix,yold,zold,u,ele,nod,id)       
!!!!!        write(*,*)numelt,numnp
        
!!!!!!!!        write(*,*) a(k08:k08+100)
!!!!!!!!        write(*,*) int(a(k08:k08+100))
        itempActiveFlagThermal=CNObjThermal%iSolutionActive
        itempActiveFlagDiffusion=CNObjDiffusion%iSolutionActive
                
        if(CNObjDiffusion%iSolutionActive == 1)then
            CNObjThermal%iSolutionActive=0
            call CN_FEM_Solver(a(k03),a(k04),a(k57),a(k18),a(k08))
            CNObjThermal%iSolutionActive=itempActiveFlagThermal
        endif   
        
        if(CNObjThermal%iSolutionActive == 1)then
            CNObjDiffusion%iSolutionActive=0
            call CN_FEM_Solver(a(k03),a(k04),a(k57),a(k18),a(k08))
!!!!!!!!!!!!     Diff(y     , z    , id   , u    ,matp)
            CNObjDiffusion%iSolutionActive=itempActiveFlagDiffusion
        endif
        iPrintOutputFlag=0
        
        call CPU_TIME(time2)     
        elpt1=time2-time1     
        write(iFU_times_out,*)'---- time of conduction/diffusion code=',elpt1
!    endif
    
    do i = 1, numelt
        do j = 1, 8
            fhg(i, j) = fhghis(i, j)
        enddo
        hgsstore(i) = hgshis(i)
        hgenerstore(i) = hgenerhis(i)
        totenerstore(i) = totenerhis(i)
        hgstress1store(i) = hgstress1his(i)
        hgstress2store(i) = hgstress2his(i)
    enddo

    nn=nn+1
   
!---------------------------------------------------
!*..... This is to call the subroutine for calculating the 
!*..... nominal stress-strain curve for both :
!*.....     - Single crystal, and
!*.....     - Polycrystal
         call m2bdforce (a(k02),a(k72),a(k73),a(k18),nn)
!---------------------------------------------------

!--------------------------------------------------- timing info
      call CPU_TIME(timesteps)     
      elpt1=timesteps-time0    
      write(iFU_times_out,*)'---- time of a step =',elpt1
      elpt1=timesteps-timebase
      write(iFU_times_out,*)'---- time of a step sum=',elpt1
!--------------------------------------------------- timing info
      
!!!!!! ---this is always the case mthsol=1=bfgs
      write(iFU_times_out,*)'---- nstep,ntime,time,termtm',nstep,ntime,time,termtm
      rewind(iFU_solsteps_out)
      write(iFU_solsteps_out,*)nstep,iprint
      flush(iFU_solsteps_out)
      write(*,*)' step#= ',nstep
      if (nstep < ntime-1)then ! .or. time < termtm) then
          goto 10
      endif
!!!!        if (nstep.lt.ntime) go to 10
!!!!        if (time.le.termtm) go to 10
                                                           
!--------------------------------------------------- timing info
       call CPU_TIME(timetotal)     
       elpt1=timetotal-timebase  
       
	   CALL SYSTEM_CLOCK(COUNT=nb_ticks_final)
	   nb_ticks = nb_ticks_final - nb_ticks_initial
	
	   IF (nb_ticks_final < nb_ticks_initial) nb_ticks = nb_ticks + nb_ticks_max
          
	   elapsed_time   = REAL(nb_ticks) / nb_ticks_sec
  
       write(iFU_times_out,*)'--------------------------------------'
       write(iFU_times_out,*)'--------------------------------------'
       write(iFU_times_out,*)'--------------------------------------'
       write(iFU_times_out,*)'---- time of CPU_TIME',elpt1
	   write(iFU_times_out,*)'---- time of SYSTEM_CLOCK=',elapsed_time
!--------------------------------------------------- timing info
       
             
        call bye (1)

   90 format(' solution phase  time=',e12.5,'  step number=',i6)
  100 format(' enter desired time step size (e10.0) ')
!!!!!!!!!  110 format(e10.0)
  120 format('  total iterations (last step)     = ',i6,2x,'(',i5,')',/, &
            '  total stiffness reformations     = ',i6,2x,'(',i5,')',/, &
            '  total rhs evaluations            = ',i6,2x,'(',i5,')',/, &
            '  current step size                = ',e12.5,//)
  130 format (//' restrt step with step size reset to:',e12.4,//)
 1000 format(//1x,12a6//' rezoner entered on cycle number ',i5,', time=',1pe10.3)
 1001 format (' mesh adjustment finished remap phase beginning')
 1002 format (/'remap phase completed-return to execution phase'/)

      end



