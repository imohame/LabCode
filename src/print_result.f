      subroutine print_result(ink,nintg,slip_n,slip_s,tauy)
      
      use CN_Objects_manager
      
!!      parameter(nume = 40000, nss = 24)
      common/hgenergy/hgenerstore(40000),hgenerhis(40000),
     1inertener(40000)
      common/totalenergy/totenerstore(40000),totenerhis(40000)
	  common /wblock2/ g_source, g_immob, g_minter, g_recov, b_v,
     1         b_vvec(87),nmo,nim
      common/wblock8/ abc(573,nume,4), his(573,nume,4)
	  common /wblock5/ enthalpy_coef, thermal_coef, temp
      common/wblock12/ Y_modulus(nume),possion_ratio(nume),tau_y(nume)
      common/result/  press, plastic_work
	  common /nab/ gn_bcc(3,43,43), n_bcc(3,43,43), np(43)
	  common /nab1/ gn_fcc(3,18,18), n_fcc(3,18,18), np1(18)	  
	  common /nab2/ gn_hcpt(3,86,86), n_hcpt(3,86,86), np2(86)
	  common /nab3/gn_hcp(3,87,87), n_hcp(3,87,87), den_im2(87),
     > np3(87)
	  common /aij/ aijhcp(24,87), aijfcc(12,18), aijbcc(24,43),
     >          aijhcpt(24,86)
	  common /intgr/ rgen(24), rrecov(87), rintr(87), rintr_n(87),
     >               rintr_p(87)
	  common /stressflag/ ElemFractCode(nume),ElemDecayCount(nume)
      common /energbox/ gbco
	  common /gbtranr/ gbtr(14)
!!!!!!	  common/hydrodiffusion/ hycon(nume,1)
!!!!!!!!!	  common /grad_pressure/ gradpdata(3,40000)
      
!!!!!!!!	  common /WMLthermal2/thermalk(nume),thermalh(nume),etae(nume)
!!!!!      common /WMLthermal2/thermalk(nume),thermalh(nume),
!!!!!     >      thermalRo(nume),thermalcp(nume),thermalKa(nume)
      
	  common/hydroembrittle/critfrac(1000), sigfrac0(40000), 
     >       sigfrac(40000),decfrac(40000)
	  common /fractureplane/ planeflag(nume), planeflag0(nume)
	  common /GND_loop/ gradslip(2,24,40000), rho_gnd(2,24,40000)
	  
      dimension slip_n(nss,3), slip_s(nss,3)


	  real ll, srm, sgdot, gd(24), dm(24), dim1(24), dim2(87)
	  real yintr(87), yrecov(87), ygen(24),fphi,rhomotot
	  integer n, n1, n2, nssim, nssm, ctr, ElemFractCode, planeflag
      integer    ::  s_s_a(nss),nass, nelec
	  real rhointn, rhointp, rhomi, rhoimi, gbtr
	  real rintr_ntot, rintr_ptot, etae, rho_gnd
	  g  = Y_modulus(ink)/(2.0*(1.0 + possion_ratio(ink)))
	  
	  rhomi=1.0E+7
	  rhoimi=1.0E+10
	  
	  ctr = 1
	  if (nim == 18) then         ! output for FCC
		  do j = 1, 12
			  rgen(j) = abc(270+j,ink,nintg)
		  end do
		  do j = 1, 18
			  rrecov(j) = abc(460+j,ink,nintg)
			  if (j .le. 12) then
			  rintr_n(j) = -1.0*(abc(57+j,ink,nintg)-rhomi - rgen(j))
				  rintr_p(j)=abc(81+j,ink,nintg)-rhoimi + rrecov(j)
				  rintr(j)=rintr_p(j)-rintr_n(j)
			  else
				  rintr_p(j) = abc(397+ctr,ink,nintg) + rrecov(j)
				  rintr(j)=rintr_p(j)
				  ctr = ctr + 1
			  end if
		  end do
	  end if
	  
	  ctr = 1
	  if (nim == 43) then          ! output for BCC
		  do j = 1, 24
			  rgen(j) = abc(270+j,ink,nintg)
		  end do
		  do j = 1, 43
			  rrecov(j) = abc(460+j,ink,nintg)
			  if (j .le. 24) then
			  rintr_n(j) = -1.0*(abc(57+j,ink,nintg)-rhomi - rgen(j))
				  rintr_p(j)=abc(81+j,ink,nintg)-rhoimi + rrecov(j)
				  rintr(j)=rintr_p(j)-rintr_n(j)
			  else
				  rintr_p(j) = abc(397+ctr,ink,nintg) + rrecov(j)
				  rintr(j)=rintr_p(j)
				  ctr = ctr + 1
			  end if
		  end do
	  end if
	  
	  	  if (nim == 86) then          ! output for HCPT
		  do j = 1, 24
			  rgen(j) = abc(270+j,ink,nintg)
		  end do
		  do j = 1, 86
			  rrecov(j) = abc(460+j,ink,nintg)
			  if (j .le. 24) then
			  rintr_n(j) = -1.0*(abc(57+j,ink,nintg)-rhomi - rgen(j))
				  rintr_p(j)=abc(81+j,ink,nintg)-rhoimi + rrecov(j)
				  rintr(j)=rintr_p(j)-rintr_n(j)
			  else
				  rintr_p(j) = abc(397+ctr,ink,nintg) + rrecov(j)
				  rintr(j)=rintr_p(j)
				  ctr = ctr + 1
			  end if
		  end do
	  end if
	  
	  	  if (nim == 87) then          ! output for HCP
		  do j = 1, 24
			  rgen(j) = abc(270+j,ink,nintg)
		  end do
		  do j = 1, 87
			  rrecov(j) = abc(460+j,ink,nintg)
			  if (j .le. 24) then
			   rintr_n(j) = -1.0*(abc(57+j,ink,nintg)-rhomi - rgen(j))
				  rintr_p(j)=abc(81+j,ink,nintg)-rhoimi + rrecov(j)
				  rintr(j)=rintr_p(j)-rintr_n(j)
			  else
				  rintr_p(j) = abc(397+ctr,ink,nintg) + rrecov(j)
				  rintr(j)=rintr_p(j)
				  ctr = ctr + 1
			  end if
		  end do
	  end if
	  
!!!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!!!c            write for one position
!!!ccccccccccccccccccccccccccccccccccccccccccccccccccc

		 
		 write(993,*) ElemFractCode(ink)
		 write(938,*) planeflag(ink)
!!!!!!!!!!!!!!!!!!		 write(2101,*) hycon(ink,1)         !hydrogen_con.out
!!!!!!!!!!!c        print out pressure gradient
!!!!!!!!!!         write(2102,*) gradpdata(1,ink)
!!!!!!!!!!         write(2103,*) gradpdata(2,ink)
!!!!!!!!!!		 write(2104,*) gradpdata(3,ink)
		 		 
		 write(20,1001) abc(1,ink,nintg)
         write(21,1001) abc(2,ink,nintg)
         write(22,1001) abc(3,ink,nintg)
         write(23,1001) abc(4,ink,nintg)    !temperature.out
         write(24,1001) press               !pressure.out
         write(25,1001) abc(6,ink,nintg)    ! plastic work
         
!!!!!!!!!!!!!!         elem_eta=1/(thermalRo(ink)*thermalcp(ink))
         elem_eta=1/(rMatTH_Ro(ink)*rMatTH_cp(ink))
		 write(2105,1001) abc(6,ink,nintg)*elem_eta+293.0    ! adiabatic temperature
!!!!!!!!!!!!!!		 write(2105,1001) abc(6,ink,nintg)*etae(ink)+293.0    ! adiabatic temperature
         
!c		 temperature change by conduction
		 write(2106,1001) abc(6,ink,nintg)*elem_eta+293.0
     >                    -abc(4,ink,nintg)  !temp_cd.out
!		 write(2106,1001) abc(6,ink,nintg)*etae(ink)+293.0
!     >                    -abc(4,ink,nintg) 
!     
	     write(2107,1001) sigfrac(ink)  !crit_frac.out
		 write(2108,1001) decfrac(ink)  !dec_frac.out
         write(270,1001) abc(5,ink,nintg)       !porosity.out
         if(abs(abc(6,ink,nintg))<1.0e-6) then            !qeqp.out Qe/Qp
		     write(937,1001) abs(abc(3,ink,nintg)/1.0e-6)
		 else
 		     write(937,1001) abs(abc(3,ink,nintg)/abc(6,ink,nintg))
		 end if
!*----------------------------------------------------------------------
!*.... resolved shear stress  res_stress1 --> res_stress24
!*----------------------------------------------------------------------
!         write(30,1001) abc(10,ink,nintg)
!         write(31,1001) abc(11,ink,nintg) 
!         write(32,1001) abc(12,ink,nintg)
!         write(33,1001) abc(13,ink,nintg) 
!         write(34,1001) abc(14,ink,nintg)
!         write(35,1001) abc(15,ink,nintg)
!         write(36,1001) abc(16,ink,nintg)
!         write(37,1001) abc(17,ink,nintg)
!         write(38,1001) abc(18,ink,nintg)
!         write(39,1001) abc(19,ink,nintg)
!         write(40,1001) abc(20,ink,nintg)
!         write(41,1001) abc(21,ink,nintg)
!         write(130,1001) abc(22,ink,nintg)
!         write(131,1001) abc(23,ink,nintg)
!         write(132,1001) abc(24,ink,nintg)
!         write(133,1001) abc(25,ink,nintg)
!         write(134,1001) abc(26,ink,nintg)
!         write(135,1001) abc(27,ink,nintg)
!         write(136,1001) abc(28,ink,nintg)
!         write(137,1001) abc(29,ink,nintg)
!         write(138,1001) abc(30,ink,nintg)
!         write(139,1001) abc(31,ink,nintg)
!         write(140,1001) abc(32,ink,nintg)
!         write(141,1001) abc(33,ink,nintg)


!****** strain-rate for each slip system
!--- slip_rate1 --> slip_rate24
         write(50,1001) abc(34,ink,nintg)
         write(51,1001) abc(35,ink,nintg) 
         write(52,1001) abc(36,ink,nintg) 
         write(53,1001) abc(37,ink,nintg) 
         write(54,1001) abc(38,ink,nintg) 
         write(55,1001) abc(39,ink,nintg) 
         write(56,1001) abc(40,ink,nintg) 
         write(57,1001) abc(41,ink,nintg) 
         write(58,1001) abc(42,ink,nintg) 
         write(59,1001) abc(43,ink,nintg) 
         write(160,1001) abc(44,ink,nintg) 
         write(161,1001) abc(45,ink,nintg) 
         write(150,1001) abc(46,ink,nintg) 
         write(151,1001) abc(47,ink,nintg) 
         write(152,1001) abc(48,ink,nintg) 
         write(153,1001) abc(49,ink,nintg) 
         write(154,1001) abc(40,ink,nintg) 
         write(155,1001) abc(51,ink,nintg) 
         write(156,1001) abc(52,ink,nintg) 
         write(157,1001) abc(53,ink,nintg) 
         write(158,1001) abc(54,ink,nintg) 
         write(159,1001) abc(55,ink,nintg) 
         write(260,1001) abc(56,ink,nintg) 
         write(261,1001) abc(57,ink,nintg)
!!******
!---densityim_1 --> densityim_43
         write(70,1001) abc(82,ink,nintg)
         write(71,1001) abc(83,ink,nintg)
         write(72,1001) abc(84,ink,nintg)
         write(73,1001) abc(85,ink,nintg)
         write(74,1001) abc(86,ink,nintg)
         write(75,1001) abc(87,ink,nintg)
         write(76,1001) abc(88,ink,nintg)
         write(77,1001) abc(89,ink,nintg)
         write(78,1001) abc(90,ink,nintg)
         write(79,1001) abc(91,ink,nintg)
         write(80,1001) abc(92,ink,nintg)
         write(81,1001) abc(93,ink,nintg)
         write(170,1001) abc(94,ink,nintg)
         write(171,1001) abc(95,ink,nintg)
         write(172,1001) abc(96,ink,nintg)
         write(173,1001) abc(97,ink,nintg)
         write(174,1001) abc(98,ink,nintg)
         write(175,1001) abc(99,ink,nintg)
         write(176,1001) abc(100,ink,nintg)
         write(177,1001) abc(101,ink,nintg)
         write(178,1001) abc(102,ink,nintg)
         write(179,1001) abc(103,ink,nintg)
         write(180,1001) abc(104,ink,nintg)
         write(181,1001) abc(105,ink,nintg)
		 write(550,999) abc(398,ink,nintg)
		 write(551,999) abc(399,ink,nintg)
		 write(552,999) abc(400,ink,nintg)
		 write(553,999) abc(401,ink,nintg)
		 write(554,999) abc(402,ink,nintg)
		 write(555,999) abc(403,ink,nintg)
		 write(556,999) abc(404,ink,nintg)
		 write(557,999) abc(405,ink,nintg)
		 write(558,999) abc(406,ink,nintg)
		 write(559,999) abc(407,ink,nintg)
		 write(560,999) abc(408,ink,nintg)
		 write(561,999) abc(409,ink,nintg)
		 write(562,999) abc(410,ink,nintg)
		 write(563,999) abc(411,ink,nintg)
		 write(564,999) abc(412,ink,nintg)
		 write(565,999) abc(413,ink,nintg)
		 write(566,999) abc(414,ink,nintg)
		 write(567,999) abc(415,ink,nintg)
		 write(568,999) abc(416,ink,nintg)
!		 write(569,999) abc(417,ink,nintg)
!		 write(570,999) abc(418,ink,nintg)
!		 write(571,999) abc(419,ink,nintg)
!		 write(572,999) abc(420,ink,nintg)
!		 write(573,999) abc(421,ink,nintg)
!		 write(574,999) abc(422,ink,nintg)
!		 write(575,999) abc(423,ink,nintg)
!		 write(576,999) abc(424,ink,nintg)
!		 write(577,999) abc(425,ink,nintg)
!		 write(578,999) abc(426,ink,nintg)
!		 write(579,999) abc(427,ink,nintg)
!		 write(580,999) abc(428,ink,nintg)
!		 write(581,999) abc(429,ink,nintg)
!		 write(582,999) abc(430,ink,nintg)
!		 write(583,999) abc(431,ink,nintg)
!		 write(584,999) abc(432,ink,nintg)
!		 write(585,999) abc(433,ink,nintg)
!		 write(586,999) abc(434,ink,nintg)
!		 write(587,999) abc(435,ink,nintg)		 
!		 write(588,999) abc(436,ink,nintg)
!		 write(589,999) abc(437,ink,nintg)
!		 write(590,999) abc(438,ink,nintg)
!		 write(591,999) abc(439,ink,nintg)
!		 write(592,999) abc(440,ink,nintg)
!		 write(593,999) abc(441,ink,nintg)
!		 write(594,999) abc(442,ink,nintg)
!		 write(595,999) abc(443,ink,nintg)
!		 write(596,999) abc(444,ink,nintg)
!		 write(597,999) abc(445,ink,nintg)
!		 write(598,999) abc(446,ink,nintg)
!		 write(599,999) abc(447,ink,nintg)
!		 write(600,999) abc(448,ink,nintg)
!		 write(601,999) abc(449,ink,nintg)
!		 write(602,999) abc(450,ink,nintg)
!		 write(603,999) abc(451,ink,nintg)
!		 write(604,999) abc(452,ink,nintg)
!		 write(605,999) abc(453,ink,nintg)
!		 write(606,999) abc(454,ink,nintg)	
!		 write(607,999) abc(455,ink,nintg)
!		 write(608,999) abc(456,ink,nintg)
!		 write(609,999) abc(457,ink,nintg)
!		 write(610,999) abc(458,ink,nintg)
!		 write(611,999) abc(459,ink,nintg)
!		 write(612,999) abc(460,ink,nintg)
         
!	 	 --- ygen1 --> ygen24
!		 write(613,999) rgen(1)
!		 write(614,999) rgen(2)
!		 write(615,999) rgen(3)
!		 write(616,999) rgen(4)
!		 write(617,999) rgen(5)
!		 write(618,999) rgen(6)
!		 write(619,999) rgen(7)
!		 write(620,999) rgen(8)
!		 write(621,999) rgen(9)
!		 write(622,999) rgen(10)
!		 write(623,999) rgen(11)
!		 write(624,999) rgen(12)
!		 write(625,999) rgen(13)
!		 write(626,999) rgen(14)
!		 write(627,999) rgen(15)
!		 write(628,999) rgen(16)
!		 write(629,999) rgen(17)
!		 write(630,999) rgen(18)
!		 write(631,999) rgen(19)
!		 write(632,999) rgen(20)
!		 write(633,999) rgen(21)
!		 write(634,999) rgen(22)
!		 write(635,999) rgen(23)
!		 write(636,999) rgen(24)
         
!		 -- yintr1 --> yintr87
!		 write(701,999) rintr(1)
!		 write(702,999) rintr(2)
!		 write(703,999) rintr(3)
!		 write(704,999) rintr(4)
!		 write(705,999) rintr(5)
!		 write(706,999) rintr(6)
!		 write(707,999) rintr(7)
!		 write(708,999) rintr(8)
!		 write(709,999) rintr(9)
!		 write(710,999) rintr(10)
!		 write(711,999) rintr(11)
!		 write(712,999) rintr(12)
!		 write(713,999) rintr(13)
!		 write(714,999) rintr(14)
!		 write(715,999) rintr(15)
!		 write(716,999) rintr(16)
!		 write(717,999) rintr(17)
!		 write(718,999) rintr(18)
!		 write(719,999) rintr(19)
!		 write(720,999) rintr(20)
!		 write(721,999) rintr(21)
!		 write(722,999) rintr(22)
!		 write(723,999) rintr(23)
!		 write(724,999) rintr(24)
!		 write(725,999) rintr(25)
!		 write(726,999) rintr(26)
!		 write(727,999) rintr(27)
!		 write(728,999) rintr(28)
!		 write(729,999) rintr(29)
!		 write(730,999) rintr(30)
!		 write(731,999) rintr(31)
!		 write(732,999) rintr(32)
!		 write(733,999) rintr(33)
!		 write(734,999) rintr(34)
!		 write(735,999) rintr(35)
!		 write(736,999) rintr(36)
!		 write(737,999) rintr(37)
!		 write(738,999) rintr(38)
!		 write(739,999) rintr(39)
!		 write(740,999) rintr(40)
!		 write(741,999) rintr(41)
!		 write(742,999) rintr(42)
!		 write(743,999) rintr(43)
!		 write(744,999) rintr(44)
!		 write(745,999) rintr(45)
!		 write(746,999) rintr(46)
!		 write(747,999) rintr(47)
!		 write(748,999) rintr(48)
!		 write(749,999) rintr(49)
!		 write(750,999) rintr(50)
!		 write(751,999) rintr(51)
!		 write(752,999) rintr(52)
!		 write(753,999) rintr(53)
!		 write(754,999) rintr(54)
!		 write(755,999) rintr(55)
!		 write(756,999) rintr(56)
!		 write(757,999) rintr(57)
!		 write(758,999) rintr(58)
!		 write(759,999) rintr(59)
!		 write(760,999) rintr(60)
!		 write(761,999) rintr(61)
!		 write(762,999) rintr(62)
!		 write(763,999) rintr(63)
!		 write(764,999) rintr(64)
!		 write(765,999) rintr(65)
!		 write(766,999) rintr(66)
!		 write(767,999) rintr(67)
!		 write(768,999) rintr(68)
!		 write(769,999) rintr(69)
!		 write(770,999) rintr(70)
!		 write(771,999) rintr(71)
!		 write(772,999) rintr(72)
!		 write(773,999) rintr(73)
!		 write(774,999) rintr(74)
!		 write(775,999) rintr(75)
!		 write(776,999) rintr(76)
!		 write(777,999) rintr(77)
!		 write(778,999) rintr(78)
!		 write(779,999) rintr(79)
!		 write(780,999) rintr(80)
!		 write(781,999) rintr(81)
!		 write(782,999) rintr(82)
!		 write(783,999) rintr(83)
!		 write(784,999) rintr(84)
!		 write(785,999) rintr(85)
!		 write(786,999) rintr(86)
!		 write(787,999) rintr(87)
         
!		 -- yrecov1 --> yrecov87
!		 write(801,999) rrecov(1)
!		 write(802,999) rrecov(2)
!		 write(803,999) rrecov(3)
!		 write(804,999) rrecov(4)
!		 write(805,999) rrecov(5)
!		 write(806,999) rrecov(6)
!		 write(807,999) rrecov(7)
!		 write(808,999) rrecov(8)
!		 write(809,999) rrecov(9)
!		 write(810,999) rrecov(10)
!		 write(811,999) rrecov(11)
!		 write(812,999) rrecov(12)
!		 write(813,999) rrecov(13)
!		 write(814,999) rrecov(14)
!		 write(815,999) rrecov(15)
!		 write(816,999) rrecov(16)
!		 write(817,999) rrecov(17)
!		 write(818,999) rrecov(18)
!		 write(819,999) rrecov(19)
!		 write(820,999) rrecov(20)
!		 write(821,999) rrecov(21)
!		 write(822,999) rrecov(22)
!		 write(823,999) rrecov(23)
!		 write(824,999) rrecov(24)
!		 write(825,999) rrecov(25)
!		 write(826,999) rrecov(26)
!		 write(827,999) rrecov(27)
!		 write(828,999) rrecov(28)
!		 write(829,999) rrecov(29)
!		 write(830,999) rrecov(30)
!		 write(831,999) rrecov(31)
!		 write(832,999) rrecov(32)
!		 write(833,999) rrecov(33)
!		 write(834,999) rrecov(34)
!		 write(835,999) rrecov(35)
!		 write(836,999) rrecov(36)
!		 write(837,999) rrecov(37)
!		 write(838,999) rrecov(38)
!		 write(839,999) rrecov(39)
!		 write(840,999) rrecov(40)
!		 write(841,999) rrecov(41)
!		 write(842,999) rrecov(42)
!		 write(843,999) rrecov(43)
!		 write(844,999) rrecov(44)
!		 write(845,999) rrecov(45)
!		 write(846,999) rrecov(46)
!		 write(847,999) rrecov(47)
!		 write(848,999) rrecov(48)
!		 write(849,999) rrecov(49)
!		 write(850,999) rrecov(50)
!		 write(851,999) rrecov(51)
!		 write(852,999) rrecov(52)
!		 write(853,999) rrecov(53)
!		 write(854,999) rrecov(54)
!		 write(855,999) rrecov(55)
!		 write(856,999) rrecov(56)
!		 write(857,999) rrecov(57)
!		 write(858,999) rrecov(58)
!		 write(859,999) rrecov(59)
!		 write(860,999) rrecov(60)
!		 write(861,999) rrecov(61)
!		 write(862,999) rrecov(62)
!		 write(863,999) rrecov(63)
!		 write(864,999) rrecov(64)
!		 write(865,999) rrecov(65)
!		 write(866,999) rrecov(66)
!		 write(867,999) rrecov(67)
!		 write(868,999) rrecov(68)
!		 write(869,999) rrecov(69)
!		 write(870,999) rrecov(70)
!		 write(871,999) rrecov(71)
!		 write(872,999) rrecov(72)
!		 write(873,999) rrecov(73)
!		 write(874,999) rrecov(74)
!		 write(875,999) rrecov(75)
!		 write(876,999) rrecov(76)
!		 write(877,999) rrecov(77)
!		 write(878,999) rrecov(78)
!		 write(879,999) rrecov(79)
!		 write(880,999) rrecov(80)
!		 write(881,999) rrecov(81)
!		 write(882,999) rrecov(82)
!		 write(883,999) rrecov(83)
!		 write(884,999) rrecov(84)
!		 write(885,999) rrecov(85)
!		 write(886,999) rrecov(86)
!		 write(887,999) rrecov(87)		 
		 
!	  yrecovtot  ygentot  yintrtot  yintr_ntot  yintr_ptot
!		 write(888,999) rrecovtot
!		 write(889,999) rgentot
!		 write(890,999) rintrtot
!		 write(891,999) rintr_ntot
!		 write(892,999) rintr_ptot

! ---- densitymo_1  --> densitymo_24
         write(90,1001) abc(58,ink,nintg)
         write(91,1001) abc(59,ink,nintg)
         write(92,1001) abc(60,ink,nintg)
         write(93,1001) abc(61,ink,nintg)
         write(94,1001) abc(62,ink,nintg)
         write(95,1001) abc(63,ink,nintg)
         write(96,1001) abc(64,ink,nintg)
         write(97,1001) abc(65,ink,nintg)
         write(98,1001) abc(66,ink,nintg)
         write(99,1001) abc(67,ink,nintg)
         write(100,1001) abc(68,ink,nintg)
         write(101,1001) abc(69,ink,nintg)
         write(190,1001) abc(70,ink,nintg)
         write(191,1001) abc(71,ink,nintg)
         write(192,1001) abc(72,ink,nintg)
         write(193,1001) abc(73,ink,nintg)
         write(194,1001) abc(74,ink,nintg)
         write(195,1001) abc(75,ink,nintg)
         write(196,1001) abc(76,ink,nintg)
         write(197,1001) abc(77,ink,nintg)
         write(198,1001) abc(78,ink,nintg)
         write(199,1001) abc(79,ink,nintg)
         write(200,1001) abc(80,ink,nintg)
         write(201,1001) abc(81,ink,nintg)
!!!**********

      nn = nintg
!      -slip_plane1 --> slip_plane24
!      write(111,999) abc(154,ink,nn)*abc(226,ink,nn)
!      write(112,999) abc(155,ink,nn)*abc(227,ink,nn)
!      write(113,999) abc(156,ink,nn)*abc(228,ink,nn)
!      write(114,999) abc(157,ink,nn)*abc(229,ink,nn)
!      write(115,999) abc(158,ink,nn)*abc(230,ink,nn)
!      write(116,999) abc(159,ink,nn)*abc(231,ink,nn)
!      write(117,999) abc(160,ink,nn)*abc(232,ink,nn)
!      write(118,999) abc(161,ink,nn)*abc(233,ink,nn)
!      write(119,999) abc(162,ink,nn)*abc(234,ink,nn)
!      write(120,999) abc(163,ink,nn)*abc(235,ink,nn)
!      write(121,999) abc(164,ink,nn)*abc(236,ink,nn)
!      write(122,999) abc(165,ink,nn)*abc(237,ink,nn)
!      write(211,999) abc(166,ink,nn)*abc(238,ink,nn)
!      write(212,999) abc(167,ink,nn)*abc(239,ink,nn)
!      write(213,999) abc(168,ink,nn)*abc(240,ink,nn)
!      write(214,999) abc(169,ink,nn)*abc(241,ink,nn)
!      write(215,999) abc(170,ink,nn)*abc(242,ink,nn)
!      write(216,999) abc(171,ink,nn)*abc(243,ink,nn)
!      write(217,999) abc(172,ink,nn)*abc(244,ink,nn)
!      write(218,999) abc(173,ink,nn)*abc(245,ink,nn)
!      write(219,999) abc(174,ink,nn)*abc(246,ink,nn)
!      write(220,999) abc(175,ink,nn)*abc(247,ink,nn)
!      write(221,999) abc(176,ink,nn)*abc(248,ink,nn)
!      write(222,999) abc(177,ink,nn)*abc(249,ink,nn)
!	  -- den_gb1 --> den_gb24
!	  write(3111,1001) abc(338,ink,nintg)
!      write(3112,1001) abc(339,ink,nintg)
!      write(3113,1001) abc(340,ink,nintg)
!      write(3114,1001) abc(341,ink,nintg)
!      write(3115,1001) abc(342,ink,nintg)
!      write(3116,1001) abc(343,ink,nintg)
!      write(3117,1001) abc(344,ink,nintg)
!      write(3118,1001) abc(345,ink,nintg)
!      write(3119,1001) abc(346,ink,nintg)
!      write(3120,1001) abc(347,ink,nintg)
!      write(3121,1001) abc(348,ink,nintg)
!      write(3122,1001) abc(349,ink,nintg)
!      write(3211,1001) abc(350,ink,nintg)
!      write(3212,1001) abc(351,ink,nintg)
!      write(3213,1001) abc(352,ink,nintg)
!      write(3214,1001) abc(353,ink,nintg)
!      write(3215,1001) abc(354,ink,nintg)
!      write(3216,1001) abc(355,ink,nintg)
!      write(3217,1001) abc(356,ink,nintg)
!      write(3218,1001) abc(357,ink,nintg)
!      write(3219,1001) abc(358,ink,nintg)
!      write(3220,1001) abc(359,ink,nintg)
!      write(3221,1001) abc(360,ink,nintg)
!      write(3222,1001) abc(362,ink,nintg)

!       --- gbtr1 --> gbtr24
!c     print results for grain boundary transmission factor	  
!	  write(3123,1001) abc(372,ink,nintg)
!      write(3124,1001) abc(373,ink,nintg)
!      write(3125,1001) abc(374,ink,nintg)
!      write(3126,1001) abc(375,ink,nintg)
!      write(3127,1001) abc(376,ink,nintg)
!      write(3128,1001) abc(377,ink,nintg)
!      write(3129,1001) abc(378,ink,nintg)
!      write(3130,1001) abc(379,ink,nintg)
!      write(3131,1001) abc(380,ink,nintg)
!      write(3132,1001) abc(381,ink,nintg)
!      write(3133,1001) abc(382,ink,nintg)
!      write(3134,1001) abc(383,ink,nintg)
!      write(3235,1001) abc(384,ink,nintg)
!      write(3236,1001) abc(385,ink,nintg)
!      write(3237,1001) abc(386,ink,nintg)
!      write(3238,1001) abc(387,ink,nintg)
!      write(3239,1001) abc(388,ink,nintg)
!      write(3240,1001) abc(389,ink,nintg)
!      write(3241,1001) abc(390,ink,nintg)
!      write(3242,1001) abc(391,ink,nintg)
!      write(3243,1001) abc(392,ink,nintg)
!      write(3244,1001) abc(393,ink,nintg)
!      write(3245,1001) abc(394,ink,nintg)
!      write(3246,1001) abc(395,ink,nintg)
!	   write(3247,1001) abc(396,ink,nintg)
!!	  -- n1 --> n4
!	  write(996,999) abc(106,ink,nn),abc(130,ink,nn),abc(154,ink,nn)
!	  write(997,999) abc(109,ink,nn),abc(133,ink,nn),abc(157,ink,nn)
!	  write(998,999) abc(112,ink,nn),abc(136,ink,nn),abc(160,ink,nn)
!	  write(999,999) abc(115,ink,nn),abc(139,ink,nn),abc(163,ink,nn)

!c     print out GND	  
!	  do j=1,24
!      --- gnd_edge1 --> gnd_edge24
!	      j1=3247+j
!          write(j1, 1001) rho_gnd(1,j,ink)
!      --gnd_screw1  --> gnd_screw24
!          j2=3271+j
!		  write(j2, 1001) rho_gnd(2,j,ink)
!	  end do

!!*----------------------------------------------------------------------
!!* .... FORMAT STATEMENTS ....
!!*----------------------------------------------------------------------

999   format (5x,e20.10)
1001  format(5x,e20.10)
1002  format(2(e20.10,5x))
1003  format(12(i2,2x))

      return
      end




















