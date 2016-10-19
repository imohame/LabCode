      subroutine input_data
*
      parameter (n_dim = 3, nss = 24, no_mat = 1000 ) !added WML 91109
      integer, parameter :: nume   = 40000
	  common /wblock1/ iplotDirxy, x_area, yield_stress !iplotDirxy=1 for x, 2=y
      common/wblock2/  g_source, g_immob, g_minter, g_recov, b_v,
     1                             b_vvec(87),nmo,nim
      common/wblock3/  density_ms, density_ims,thermalEnthalpy(1000)
      common /wblock5/ enthalpy_coef, thermal_coef, temp
      common /wblock7/ slip_n0(1000,nss,3), slip_s0(1000,nss,3)!Changed to accomodate multiple slip WML 91109
      common /wblock9/ slip_n0_t(1000,nss,3),slip_s0_t(1000,nss,3)
      common /wblock10/ng,grain_mo(1000,3),bv(no_mat,87),nssmat(1000)  !added WML 91109
     >        ,nssimmat(1000)
      common /wblock11/ pd_counter,rhoim0(1000,nss),rhomo0(1000,nss)
!!!!!!!      common/WMLthermal/thermalflag
!!!!!!!      integer thermalflag
	  common /nab/ gn_bcc(3,43,43), n_bcc(3,43,43), np(43)
	  common /nab1/ gn_fcc(3,18,18), n_fcc(3,18,18), np1(18)
	  common /nab2/ gn_hcpt(3,86,86), n_hcpt(3,86,86), np2(86)
	  common /nab3/gn_hcp(3,87,87), n_hcp(3,87,87), den_im2(87),
     > np3(87)
	  common /aij/ aijhcp(24,87), aijfcc(12,18), aijbcc(24,43),
     >          aijhcpt(24,86)
!!!	  common /stressflag/ ElemFractCode(nume),ElemDecayCount(nume)
!!!!	  common /gbblock/ gbvec(nume,3), gbflag(nume,3)
	  common/bk00/
     1 k01,k02,k03,k04,k05,k06,k07,k08,k09,k10,k11,k12,
     2 k13,k14,k15,k16,k17,k18,k19,k20,k21,k22,k23,k24,
     3 k25,k26,k27,k28,k29,k30,k31,k32,k33,k34,k35,k36,
     4 k37,k38,k39,k40,k41,k42,k43,k44,k45,k46,k47,k48,
     5 k49,k50,k51,k52,k53,k54,k55,k56,k57,k58,k59,k60,
     6 k61,k62,k63,k64,k65,k66,k67,k68,k69,k70,k71,k72,
     7 k73,k74,k75,k76,k77,k78,k79,k80,k81,k82,k83,k84,
     8 k85,k86,k87,k88,k89,k90,k91,k92,k93,k94,k95,k96
	  common /main_block/ a(1)
	  common /cleavage_plane/ cleave_n(1000,3,3)
	  common /cleavage_plane0/ cleave_n0(3,3)
	  common /sigfrac/ sigmacrit0, sigmacrit1,sigmacrit2, sigmacrit3,
     >      DecayCount, f_decay, penalty,fractFlag
	  common /precrack/ npc, elepc(100)
      common /elas_energ/ inie(nume), ecount(nume)
	  common /energbox/ gbco
	  common /cracktipcircle/ ndcircle
!!!!	  common /fractureplane/ planeflag(nume), planeflag0(nume)
!!	  common /precrack2/ rlflag

      character *80 txt
	  real dum, b_v, inie,  gbco!!, gbvec
!!!	  integer ElemFractCode, ElemDecayCount !!, gbflag
	  integer var_no(1000), ngbctr, DecayCount ,fractFlag
	  real slip_s00(4,24,3), slip_n00(4,24,3)
	  real cleave_n, cleave_n0, sigmacrit0, sigmacrit1, sigmacrit2,
     > 	  sigmacrit3, f_decay, penalty
	  integer npc, elepc, ecount, ndcircle
!!!      integer  planeflag,planeflag0
	  integer  rlflag


      data (((slip_n00(i,j,k), k = 1, 3), j = 1, 24), i = 1, 4)
     >      / 1 , 1 , 1,  1, 1, 1,  1, 1, 1,
     >       -1 ,-1 , 1, -1,-1, 1, -1,-1, 1,
     >       -1 , 1 , 1, -1, 1, 1, -1, 1, 1,
     >        1 ,-1 , 1,  1,-1, 1,  1,-1, 1,
     >        0 , 0 , 0,  0, 0, 0,  0, 0, 0,
     >        0 , 0 , 0,  0, 0, 0,  0, 0, 0,
     >        0 , 0 , 0,  0, 0, 0,  0, 0, 0,
     >        0 , 0 , 0,  0, 0, 0,  0, 0, 0,
     >       -1 , 1 , 0, -1, 0, 1,  0,-1, 1,
     >        1 , 1 ,-2,  1,-2, 1, -2, 1, 1,
     >        1 , 1 , 0,  1, 0, 1,  0,-1, 1,
     >        2 , 1 , 1,  1,-1, 2,  1, 2,-1,
     >        1 , 1 , 0,  0, 1, 1, -1, 0, 1,
     >        1 , 2 , 1, -1, 1, 2,  2, 1,-1,
     >        1 , 0 , 1,  0, 1, 1, -1, 1, 0,
     >        1 , 1 , 2,  2,-1, 1, -1, 2, 1,
     >        0,0 ,-1.594,0,0 ,-1.594,0,0 ,-1.594,
     >         1,1.232,1.594, -2,0.134,-1.594,1,-1.366,0,
     >        1,1.232,0,-2,0.134,-3.188,1,-1.366,-1.594,
     >        -1,-1.232,-3.188,2,-0.134,0,-1,1.366,-1.594,
     >         1,1.232,0,1,1.232,0,-2,0.134,-3.188,
     >         -2,0.134,-3.188,1,-1.366,-1.594,1,-1.366,-1.594,
     >        -1,-1.232,-3.188,-1,-1.232,-3.188,2,-0.134,0,
     >        2,-0.134,0, -1,1.366,-1.594,-1,1.366,-1.594,
     >        0,0,0.3137,0,0,0.3137,0,0,0.3137,
     >        1,1.1547,0, -0.5,-1.4434,0,-0.5,0.2887,0,
     >      1,1.1547,0.3137,-0.5,-1.4434,0.3137,-0.5,0.2887,0.3137,
     >      -1,-1.1547,0.3137,0.5,1.4434,0.3137,0.5,-0.2887,0.3137,
     >       1,1.1547,0.3137, 1,1.1547,0.3137,-0.5,-1.4434,0.3137,
     >   -0.5,-1.4434,0.3137,-0.5,0.2887,0.3137,-0.5,0.2887,0.3137,
     >       -1,-1.1547,0.3137,-1,-1.1547,0.3137,0.5,1.4434,0.3137,
     >     0.5,1.4434,0.3137,0.5,-0.2887,0.3137,0.5,-0.2887,0.3137/


      data (((slip_s00(i,j,k), k = 1, 3), j = 1, 24), i = 1, 4)
     >      /-1 , 0 , 1, -1, 1, 0,  0,-1, 1,
     >		  0 , 1 , 1, -1, 1, 0,  1, 0, 1,
     >		  1 , 0 , 1,  1, 1, 0,  0,-1, 1,
     >		  0 , 1 , 1,  1, 1, 0, -1, 0, 1,
     >        0 , 0 , 0,  0, 0, 0,  0, 0, 0,
     >        0 , 0 , 0,  0, 0, 0,  0, 0, 0,
     >        0 , 0 , 0,  0, 0, 0,  0, 0, 0,
     >        0 , 0 , 0,  0, 0, 0,  0, 0, 0,
     >        1 , 1 , 1,  1, 1, 1,  1, 1, 1,
     >        1 , 1 , 1,  1, 1, 1,  1, 1, 1,
     >       -1 , 1 , 1, -1, 1, 1, -1, 1, 1,
     >       -1 , 1 , 1, -1, 1, 1, -1, 1, 1,
     >        1 ,-1 , 1,  1,-1, 1,  1,-1, 1,
     >        1 ,-1 , 1,  1,-1, 1,  1,-1, 1,
     >        1 , 1 ,-1,  1, 1,-1,  1, 1,-1,
     >        1 , 1 ,-1,  1, 1,-1,  1, 1,-1,
     >         0,0.866,0,1,-0.5,0,-1,-0.366,0,
     >          1,-0.5,0,0,0.866,0,-1,-0.366,0,
     >         1,-0.5,0,0,0.866,0,-1,-0.366,0,
     >         1,-0.5,0,0,0.866,0,1,0.366,0,
     >          0,0,-0.531,1,-0.5,-1.594,2,0.732,-1.594,
     >         2,-0.134,-0.531,1,0.4106,-0.5313, 0 ,0.866,-1.594,
     >         1,0.4106,-0.5313,2,0.732,-1.594,0,0,-0.5313,
     >          0,0.866,-1.594,1,-0.5,-1.594,2,-0.0446,-0.5313,
     >		   1,0,0,-0.5,0.866,0,-0.5,-0.866,0,
     >		  -0.5,0.866,0,1,0,0,-0.5,-0.866,0,
     >        -0.5,0.866,0,1,0,0,-0.5,-0.866,0,
     >		  -0.5,0.866,0,1,0,0,0.5,0.866,0,
     >        -0.5,-0.866,1.594,-1,0,-1.594,0.5,0.866,1.594,
     >        -0.5,0.866,1.594,1,0,1.594,0.5,-0.866,1.594,
     >         1,0 ,1.594,0.5,0.866,1.594,-0.5,-0.866,1.594,
     >         0.5,-0.866,1.594,-1, 0,-1.594,-0.5,0.866,1.594 /

      write(*,*) 'beginning of input data'
      open(60,file = 'data.in', status = 'unknown')

	  open(61,file = 'a.txt', status = 'unknown')
	  open(62,file = 'b.txt', status = 'unknown')
	  open(63,file = 'c.txt', status = 'unknown')
	  open(85,file = 'd.txt', status = 'unknown')

	  open(64,file = 'a1.txt', status = 'unknown')
	  open(65,file = 'b1.txt', status = 'unknown')
	  open(66,file = 'c1.txt', status = 'unknown')
	  open(86,file = 'd1.txt', status = 'unknown')

	  open(67,file = 'a3.txt', status = 'unknown')
	  open(68,file = 'b3.txt', status = 'unknown')
	  open(69,file = 'c3.txt', status = 'unknown')
	  open(87,file = 'd3.txt', status = 'unknown')

	  open(82,file = 'a4.txt', status = 'unknown')
	  open(83,file = 'b4.txt', status = 'unknown')
	  open(84,file = 'c4.txt', status = 'unknown')
	  open(88,file = 'd4.txt', status = 'unknown')

!!!!	  open(89,file = 'gb.f', status = 'unknown')

	  cleave_n0(1,1) = 1.
	  cleave_n0(1,2) = 0.
	  cleave_n0(1,3) = 0.
	  cleave_n0(2,1) = 0.
	  cleave_n0(2,2) = 1.
	  cleave_n0(2,3) = 0.
	  cleave_n0(3,1) = 0.
	  cleave_n0(3,2) = 0.
	  cleave_n0(3,3) = 1.


	  do i = 1, nume
!!!!!!	  planeflag(i)=0
!!!!!!	  planeflag0(i)=0
!!	  ElemFractCode(i) = 0
!!	  ElemDecayCount(i) = 1
!!!!!!	  gbvec(i,1) = 0.0
!!!!!!	  gbvec(i,2) = 0.0
!!!!!!	  gbvec(i,3) = 0.0
!!!!!!	  gbflag(i,1) = 0
!!!!!!	  gbflag(i,2) = 0
!!!!!!	  gbflag(i,2) = 0
	  inie(i)=0.0
	  ecount(i)=0
	  end do

!!!	  read(89,*) ngbctr
!!!	  write(999,*) ngbctr
!!!	  do i = 1,ngbctr
!!!	  read(89,*) j1,j2,j3,j4
!!!	  write(999,*) j1,j2,j3,j4
!!!	  gbflag(j1,1)=j2
!!!	  gbflag(j1,2)=j3
!!!	  gbflag(j1,3)=j4
!!!	  gbflag(j2,1)=j1
!!!	  gbflag(j2,2)=j3
!!!	  gbflag(j2,3)=j4
!!!	  end do

!!!!!	  call edge(a(k03),a(k04),a(k08))

	  do i = 1, 43
	  do j = 1, 43
	  do k = 1, 3
		  read(61,*) gn_bcc(k,j,i)
		  read(62,*) n_bcc(k,j,i)
	  end do
	  end do
	  read(63,*) np(i)
	  end do

	  do i = 1, 18
	  do j = 1, 18
	  do k = 1, 3
		  read(64,*) gn_fcc(k,j,i)
		  read(65,*) n_fcc(k,j,i)
	  end do
	  end do
	  read(66,*) np1(i)
	  end do

	  do i = 1, 86
	  do j = 1, 86
	  do k = 1, 3
		  read(67,*) gn_hcpt(k,j,i)
		  read(68,*) n_hcpt(k,j,i)
	  end do
	  end do
	  read(69,*) np2(i)
	  end do

	  do i = 1, 87
	  do j = 1, 87
	  do k = 1, 3
		  read(82,*) gn_hcp(k,j,i)
		  read(83,*) n_hcp(k,j,i)
	  end do
	  end do
	  read(84,*) np3(i)
	  end do

	  do i = 1, 43
	  do j = 1, 24
		  read(85,*) dum
		  aijbcc(j,i) = sqrt(dum)
	  end do
	  end do

	  do i = 1, 18
	  do j = 1, 12
		  read(86,*) dum
		  aijfcc(j,i) = sqrt(dum)
	  end do
	  end do

	  do i = 1, 86
	  do j = 1, 24
		  read(87,*) dum
		  aijhcpt(j,i) = sqrt(dum)
	  end do
	  end do

	  do i = 1, 87
	  do j = 1, 24
		  read(88,*) dum
		  aijhcp(j,i) = sqrt(dum)
	  end do
	  end do



!!!!      read (60,*) sigmacrit0
!!!!	  read (60,*) sigmacrit1
!!!!      read (60,*) sigmacrit2
!!!!	  read (60,*) sigmacrit3
	  read (60,*) fractFlag
	  read (60,*) DecayCount
	  read (60,*) f_decay
	  read (60,*) ndcircle
	  read (60,*) iplotDirxy
      read (60,*) x_area
      read (60,*) yield_stress
      read (60,*) b_v
      read (60,*) g_source
      read (60,*) g_immob
      read (60,*) g_minter
      read (60,*) g_recov
      read (60,*) density_ms
      read (60,*) density_ims
      !read (60,*) enthalpy_coef

!!!c     Elements for pre-existing crack
	  read(60,*) npc
	  do i=1, npc
	      read(60,*) elepc(i)
	  end do
      !-- the following two numbers are not used any more
	  read(60,*) rlflag     ! pre-existing crack at left or right edge, 0 left, 1 right
	  read(60,*) gbco



      read (60,*) ng
      do i = 1, ng
	  read (60,*) var_no(i),grain_mo(i,1), grain_mo(i,2),
     >  grain_mo(i,3)
	  if (var_no(i) .lt. 1) then   !fcc
	  nssmat(i) = 12
	  nssimmat(i) = 18
	  else if ((var_no(i) .ge. 1) .and. (var_no(i) .le. 24)) then   !bcc
	  nssmat(i) = 24
	  nssimmat(i) = 43
	  else if ((var_no(i) .ge. 25) .and. (var_no(i) .le. 60)) then
	  nssmat(i) = 12
	  nssimmat(i) = 18
	  else if (var_no(i) == 61) then    !hcpt
	  nssmat(i) = 24
	  nssimmat(i) = 86
	  else if (var_no(i) == 62) then   !hcp
	  nssmat(i) = 24
	  nssimmat(i) = 87

	  end if
	  write(*,*) 'input_data', nssimmat(i)
!c	  if (thermalflag.eq.0) then
!c	  read (60,*) ecin(i)!enthalpy coefficent
!c	  read (60,*) etain(i)!nondim fractionofworktoheat/rhocp
!c	  endif
      end do
!c	  write(*,*) 'in input_data' !WMLWRITE51210
!c	  write(*,*) nssmat(1)



	  do i = 1, ng
	  do j = 1, nssimmat(i)

	  bv(i,j) = b_v
	  end do
	  end do
      write(*,*) bv(1,1),'bv'
      do i = 1, ng
	  if (var_no(i) .lt. 1) then
         do j = 1, nssmat(i)
		 do k = 1, 3
               slip_n0(i,j,k) = slip_n00(1,j,k)
			   slip_s0(i,j,k) = slip_s00(1,j,k)
         end do
		 end do
	 else if ((var_no(i) .ge. 1) .and. (var_no(i) .le. 24)) then
		 do j = 1, nssmat(i)
		 do k = 1, 3
               slip_n0(i,j,k) = slip_n00(2,j,k)
			   slip_s0(i,j,k) = slip_s00(2,j,k)
         end do
		 end do
	  else if ((var_no(i) .ge. 25) .and. (var_no(i) .le. 60)) then
         do j = 1, nssmat(i)
		 do k = 1, 3
               slip_n0(i,j,k) = slip_n00(1,j,k)
			   slip_s0(i,j,k) = slip_s00(1,j,k)
         end do
		 end do

		 else if (var_no(i) == 61) then
		 do j = 1, nssmat(i)
		 do k = 1, 3
               slip_n0(i,j,k) = slip_n00(3,j,k)
			   slip_s0(i,j,k) = slip_s00(3,j,k)
         end do
		 end do
		 else if (var_no(i) == 62) then
		 do j = 1, nssmat(i)
		 do k = 1, 3
               slip_n0(i,j,k) = slip_n00(4,j,k)
			   slip_s0(i,j,k) = slip_s00(4,j,k)
         end do
		 end do


	  end if
      end do

	  do i=1,ng
		    do j=1,nssmat(i)
			      rhoim0(i,j) = density_ims
			      rhomo0(i,j) = density_ms
			  end do
	  end do

      call transform(nss,slip_n0,slip_s0,slip_n0_t,slip_s0_t, var_no)
      ! #####################################################
      !      read the multiple time step specs ....dtimesteps.in
            call ReaddtimestepsSpecs()

! #####################################################
! #####################################################
!      to read the GB data, normals, element map
            call GBReadInput()
! #####################################################
! #####################################################
!C      to read precracked elements and apply proper cleavage planes
!       call FractReadApplyPreCrackCleavagePlanes()!a(k03),a(k04),a(k02))
!!!!      call FractReadApplyPreCrackCleavagePlanes(y,z,ix)
! #####################################################

      return
      end
