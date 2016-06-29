module lo_match
use global_var
implicit none
real(8) :: Cri(1:6), Crf(1:6), Cri_out(1:6)
real(8) :: solR_temp
contains
	!Andersson 2011 Matching---------------------------------------------------
	subroutine lo_bc_Solve(H0F, KF, cY, k2_, Coef)
	! solving with Andersson 2011
	use lin_sol_gen_int
	use puls_eqt_set_opt03_ex1
	use lo_eqt
use Root_Finding
	!	2 components
	!	Quantity Q = {QQ(i)} [Y] {c(i)} ------ {vec} [2nd Tensor] {vec} inner product
	implicit none
	real(8) :: H0F, KF, cY(1:6,1:5), k2_, Coef(1:5)
	real(8) :: eNu, eLamb, eNu_hf, eLamb_hf, dLamb, dNu, dP, gamma_1, P_, rho_, m_, mu_,  nu_, r_
	real(8) :: pi4, l_1, l_2_hf, ddNu, cs2
	real(8) :: XX(1:6), HH(1:6), KK(1:6), H2_C(1:6), PP(1:6), TT(1:6), PPF
	integer :: mode, i
	real(8) :: YI(1:6, 1:5),  YS(1:6, 1:5), M(1:3, 1:3), Vec(1:3,1), Sol(1:3,1), YS_Comb(1:6), HS, dHS, yR_
real(8) :: YI_Comb(1:6)
		!¡¸ Crust Core interface condition
		P_ = Cri(1)
		rho_ = Cri(2)
		m_ = Cri(3)
		!mu_ = mu_Cr(ii) Andersson's mu seems to be off by 1/2, st it cannot be reduced to RCA case
		mu_ = Cri(4)
		nu_ = Cri(5)
		r_ =  Cri(6)
gamma_1 = lo_gamma(P_)

		pi4 = 4.d0*pi
		l_1 = l_0 * (l_0 + 1.d0)
		l_2_hf = (l_0 + 2.d0)*(l_0 - 1.d0)/2.d0

		eLamb = pes03_eLamb(rho_, m_, r_) ** 2
		eNu = dexp(2.d0*nu_)
		eLamb_hf = eLamb**(0.5d0)
		eNu_hf = eNu**(0.5d0)
		dLamb = pes03_dLamb(rho_, m_, r_) * 2.d0
		dNu = pes03_dnu(P_, rho_, m_, r_) * 2.d0
		ddNu = pes03_ddnu(P_, rho_, m_, r_) * 2.d0
		dP = -dNu/2.d0 * (P_ + rho_)
		cs2 = (P_/(rho_+P_))*gamma_1


		!¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸
		! Solving for WE as given in Andersson (2011) Eqt(53)-(55)
		call lo_solid_Algebraic_Relation(P_, rho_, m_, nu_, mu_, gamma_1, r_, PP, TT)
		PPF = r_**2/2.d0 * (P_ + rho_) * (H0F)		

		Coef(1) = H0F
		Coef(3) = KF

		YI = 0.d0
		YI(1, 1) = 1.d0
		YI(2, 2) = 1.d0
		YI(3, 3) = 1.d0
		YI(4, 4) = 1.d0
		YI(5, 5) = 1.d0

		M(1, 1:3) = 0.d0
		Vec(1,1) = PPF							! Fluid, Zero Freq: P_Eu = -0.5 * (rho+P) * H0_LD
		do i = 1, 6
			M(1, 1) = M(1, 1) + YI(i, 2) * (r_**2.d0*PP(i) + TT(i))
			M(1, 2) = M(1, 2) + YI(i, 4) * (r_**2.d0*PP(i) + TT(i))
			M(1, 3) = M(1, 3) + YI(i, 5) * (r_**2.d0*PP(i) + TT(i))
			Vec(1,1) = Vec(1,1) - (Coef(1)*YI(i, 1) + Coef(3)*YI(i, 3) )* (r_**2.d0*PP(i) + TT(i))		! Continuity of (r^2 P_Eu - T1) 
		enddo

		!¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸
		! Solving for Crust Surface BC
		
		YS = 0.d0
		YS(1:6, 1) = cY(1:6, 1)
		YS(1:6, 2) = cY(1:6, 2)
		YS(1:6, 3) = cY(1:6, 3)
		YS(1:6, 4) = cY(1:6, 4)
		YS(1:6, 5) = cY(1:6, 5)

		!¡¸ T1 + P continuous
		!¡¸ Crust Surface condition
		P_ = Crf(1)
		rho_ = Crf(2)
		m_ = Crf(3)
		!mu_ = mu_Cr(ii) Andersson's mu seems to be off by 1/2, st it cannot be reduced to RCA case
		mu_ = Crf(4)
		nu_ = Crf(5)
		r_ = Crf(6)
gamma_1 = lo_gamma(P_)

		eLamb = pes03_eLamb(rho_, m_, r_) ** 2
		eNu = dexp(2.d0*nu_)
		eLamb_hf = eLamb**(0.5d0)
		eNu_hf = eNu**(0.5d0)
		dLamb = pes03_dLamb(rho_, m_, r_) * 2.d0
		dNu = pes03_dnu(P_, rho_, m_, r_) * 2.d0
		ddNu = pes03_ddnu(P_, rho_, m_, r_) * 2.d0
		dP = -dNu/2.d0 * (P_ + rho_)
		cs2 = (P_/(rho_+P_))*gamma_1

		call lo_solid_Algebraic_Relation(P_, rho_, m_, nu_, mu_, gamma_1, r_, PP, TT)
		XX = 0.d0
		XX(4) = -0.5d0*(P_+rho_)*dNu/r_
		XX(1:6) = XX(1:6) + PP(1:6)

		Vec(2,1) = 0.d0
		M(2, 1:3) = 0.d0
		do i = 1, 6
			M(2, 1) = M(2, 1) + YS(i, 2) * (r_**2.d0*XX(i) + TT(i))
			M(2, 2) = M(2, 2) + YS(i, 4) * (r_**2.d0*XX(i) + TT(i))
			M(2, 3) = M(2, 3) + YS(i, 5) * (r_**2.d0*XX(i) + TT(i))
			Vec(2,1) = Vec(2,1) - (Coef(1) * YS(i, 1) + Coef(3) * YS(i, 3))* (r_**2.d0*XX(i) + TT(i))
		enddo

		!¡¸ T2 = 0
		M(3, 1) = YS(6, 2)
		M(3, 2) = YS(6, 4)
		M(3, 3) = YS(6, 5)
		Vec(3,1) = 0.d0 - (  Coef(1) * YS(6, 1) + Coef(3) * YS(6, 3)  )

		call LIN_SOL_GEN(M, Vec, Sol)					!IMSL linear equations solver
		coef(2) = Sol(1,1)
		coef(4) = Sol(2,1)
		coef(5) = Sol(3,1)
		
		YS_Comb = 0.d0
		do i = 1,5
			YS_Comb(1:6) = YS_Comb(1:6) + coef(i)*YS(1:6, i)
		enddo
		
		dHS = YS_Comb(2)
		HS = YS_Comb(1)

		yR_ = r_*dHS/HS
		call lo_k2(yR_, YS_Comb(5)/HS, P_, rho_, m_, mu_, r_, k2_)

solR_temp = yR_- 4.d0*pi*rho_*r_**3/m_ + 32.d0*pi*	(dNu* r_)* (mu_*YS_Comb(5)/HS)

! Confirm that the Linear Equation Solver is accurate
!write(*,*) M(1, 1)*coef(2) + M(1, 2)*coef(4) + M(1, 3)*coef(5)
!write(*,*) Vec(1,1)
!write(*,*) M(2, 1)*coef(2) + M(2, 2)*coef(4) + M(2, 3)*coef(5)
!write(*,*) Vec(2,1)
!write(*,*) M(3, 1)*coef(2) + M(3, 2)*coef(4) + M(3, 3)*coef(5)
!write(*,*) Vec(3,1)
!pause

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*** Test
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*) "surface W, H0", r_*YS_Comb(4)*dP, 0.5d0*(rho_+P_)*r_**2.d0*YS_Comb(1)	

P_ = Cri(1)
rho_ = Cri(2)
m_ = Cri(3)
mu_ = Cri(4)
nu_ = Cri(5)
r_ =  Cri(6)
gamma_1 = lo_gamma(P_)

pi4 = 4.d0*pi
l_1 = l_0 * (l_0 + 1.d0)
l_2_hf = (l_0 + 2.d0)*(l_0 - 1.d0)/2.d0

eLamb = pes03_eLamb(rho_, m_, r_) ** 2
eNu = dexp(2.d0*nu_)
eLamb_hf = eLamb**(0.5d0)
eNu_hf = eNu**(0.5d0)
dLamb = pes03_dLamb(rho_, m_, r_) * 2.d0
dNu = pes03_dnu(P_, rho_, m_, r_) * 2.d0
ddNu = pes03_ddnu(P_, rho_, m_, r_) * 2.d0
dP = -dNu/2.d0 * (P_ + rho_)
cs2 = (P_/(rho_+P_))*gamma_1

write(*,*) "FE interface W, H0", r_*Coef(4)*dP, 0.5d0*(rho_+P_)*r_**2.d0*Coef(1)
call lo_solid_Algebraic_Relation(P_, rho_, m_, nu_, mu_, gamma_1, r_, PP, TT)
YI_Comb(1:5) = Coef(1:5)
YI_Comb(6) = 0.d0
write(*,*) "FE interface r^2 dP + T1, H0 terms", r_**2 * vec_dot(PP(1:6), YI_Comb(1:6)) + vec_dot(TT(1:6), YI_Comb(1:6)), - r_**2 /2.d0*(rho_ + P_) * Coef(1)
write(*,*) "FE interface T2", YI_Comb(6)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*** Test
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	endsubroutine lo_bc_Solve

!Andersson 2011 Matching---------------------------------------------------
	subroutine lo_bc_Solve_2Side(H0F, KF, cY, k2_, Coef)
	! solving with Andersson 2011
	use lin_sol_gen_int
	use puls_eqt_set_opt03_ex1
	use lo_eqt
use Root_Finding
	!	2 components
	!	Quantity Q = {QQ(i)} [Y] {c(i)} ------ {vec} [2nd Tensor] {vec} inner product
	implicit none
	real(8) :: H0F, KF, cY(1:6,1:5), k2_, Coef(1:5)
	real(8) :: eNu, eLamb, eNu_hf, eLamb_hf, dLamb, dNu, dP, gamma_1, P_, rho_, m_, mu_,  nu_, r_
	real(8) :: pi4, l_1, l_2_hf, ddNu, cs2
	real(8) :: XX(1:6), HH(1:6), KK(1:6), H2_C(1:6), PP(1:6), TT(1:6), PPF
	integer :: mode, i
	real(8) :: YI(1:6, 1:5),  YS(1:6, 1:5), M(1:5, 1:5), Vec(1:5,1), Sol(1:5,1), YS_Comb(1:6), HS, dHS, yR_
real(8) :: YI_Comb(1:6)
		!¡¸ Crust Core interface condition
		P_ = Cri(1)
		rho_ = Cri(2)
		m_ = Cri(3)
		!mu_ = mu_Cr(ii) Andersson's mu seems to be off by 1/2, st it cannot be reduced to RCA case
		mu_ = Cri(4)
		nu_ = Cri(5)
		r_ =  Cri(6)
gamma_1 = lo_gamma(P_)

		pi4 = 4.d0*pi
		l_1 = l_0 * (l_0 + 1.d0)
		l_2_hf = (l_0 + 2.d0)*(l_0 - 1.d0)/2.d0

		eLamb = pes03_eLamb(rho_, m_, r_) ** 2
		eNu = dexp(2.d0*nu_)
		eLamb_hf = eLamb**(0.5d0)
		eNu_hf = eNu**(0.5d0)
		dLamb = pes03_dLamb(rho_, m_, r_) * 2.d0
		dNu = pes03_dnu(P_, rho_, m_, r_) * 2.d0
		ddNu = pes03_ddnu(P_, rho_, m_, r_) * 2.d0
		dP = -dNu/2.d0 * (P_ + rho_)
		cs2 = (P_/(rho_+P_))*gamma_1

		!¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸
		! Solving for WE as given in Andersson (2011) Eqt(53)-(55)
		call lo_solid_Algebraic_Relation(P_, rho_, m_, nu_, mu_, gamma_1, r_, PP, TT)
		PPF = r_**2/2.d0 * (P_ + rho_) * (H0F)		

		YI = 0.d0
		YI(1:6, 1) = cY(1:6, 1)
		YI(1:6, 2) = cY(1:6, 2)
		YI(1:6, 3) = cY(1:6, 3)
		YI(1:6, 4) = cY(1:6, 4)
		YI(1:6, 5) = cY(1:6, 5)


		M(1, 1:5) = 0.d0
		do i = 1, 6
			M(1, 1:5) = M(1, 1:5) + YI(i, 1:5) * (r_**2.d0*PP(i) + TT(i))
		enddo
		Vec(1,1) = PPF						! Fluid, Zero Freq: P_Eu = -0.5 * (rho+P) * H0_LD

! Temparay tools: input the W value right at the Interface (found with Fluid_W code); continuity of T + r^2P replaced by continuity of W
!read(*,*) PPF
!Vec(1,1) = PPF
!M(1, 1:5) = YI(4, 1:5)
		
		M(2, 1:5) = YI(6, 1:5)
		Vec(2,1) = 0.d0

		M(3, 1:5) = YI(1, 1:5)
		Vec(3,1) = H0F
		M(4, 1:5) = YI(3, 1:5)
		Vec(4,1) = KF

		!¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸
		! Solving for Crust Surface BC
		
		YS = 0.d0
		YS(1, 1) = 1.d0
		YS(2, 2) = 1.d0
		YS(3, 3) = 1.d0
		YS(4, 4) = 1.d0
		YS(5, 5) = 1.d0

		!¡¸ T1 + P continuous
		!¡¸ Crust Surface condition
		P_ = Crf(1)
		rho_ = Crf(2)
		m_ = Crf(3)
		!mu_ = mu_Cr(ii) Andersson's mu seems to be off by 1/2, st it cannot be reduced to RCA case
		mu_ = Crf(4)
		nu_ = Crf(5)
		r_ = Crf(6)
gamma_1 = lo_gamma(P_)

		eLamb = pes03_eLamb(rho_, m_, r_) ** 2
		eNu = dexp(2.d0*nu_)
		eLamb_hf = eLamb**(0.5d0)
		eNu_hf = eNu**(0.5d0)
		dLamb = pes03_dLamb(rho_, m_, r_) * 2.d0
		dNu = pes03_dnu(P_, rho_, m_, r_) * 2.d0
		ddNu = pes03_ddnu(P_, rho_, m_, r_) * 2.d0
		dP = -dNu/2.d0 * (P_ + rho_)
		cs2 = (P_/(rho_+P_))*gamma_1

		call lo_solid_Algebraic_Relation(P_, rho_, m_, nu_, mu_, gamma_1, r_, PP, TT)
		XX = 0.d0
		XX(4) = -0.5d0*(P_+rho_)*dNu/r_
		XX(1:6) = XX(1:6) + PP(1:6)

		M(5, 1:5) = 0.d0
		do i = 1, 6
			M(5, 1:5) = M(5, 1:5) + YS(i, 1:5) * (r_**2.d0*XX(i) + TT(i)) 
		enddo

		Vec(5,1) = 0.d0

		call LIN_SOL_GEN(M, Vec, Sol)					!IMSL linear equations solver
		coef(1:5) = Sol(1:5,1)
		
		YS_Comb = 0.d0
		do i = 1,5
			YS_Comb(1:6) = YS_Comb(1:6) + coef(i)*YS(1:6, i)
		enddo

		dHS = YS_Comb(2)
		HS = YS_Comb(1)

		yR_ = r_*dHS/HS
!solR_temp = yR_- 4.d0*pi*rho_*r_**3/m_ + 32.d0*pi*	(dNu* r_)* (mu_*YS_Comb(5)/HS)
!write(*,*) solR_temp

! Confirm that the Linear Equation Solver is accurate
!write(*,*) M(1, 1)*coef(2) + M(1, 2)*coef(4) + M(1, 3)*coef(5)
!write(*,*) Vec(1,1)
!write(*,*) M(2, 1)*coef(2) + M(2, 2)*coef(4) + M(2, 3)*coef(5)
!write(*,*) Vec(2,1)
!write(*,*) M(3, 1)*coef(2) + M(3, 2)*coef(4) + M(3, 3)*coef(5)
!write(*,*) Vec(3,1)
!pause
		call lo_k2(yR_, YS_Comb(5)/HS, P_, rho_, m_, mu_, r_, k2_)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*** Test
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*) "surface W, H0", r_*YS_Comb(4)*dP, 0.5d0*(rho_+P_)*r_**2.d0*YS_Comb(1)	

P_ = Cri(1)
rho_ = Cri(2)
m_ = Cri(3)
mu_ = Cri(4)
nu_ = Cri(5)
r_ =  Cri(6)
gamma_1 = lo_gamma(P_)

pi4 = 4.d0*pi
l_1 = l_0 * (l_0 + 1.d0)
l_2_hf = (l_0 + 2.d0)*(l_0 - 1.d0)/2.d0

eLamb = pes03_eLamb(rho_, m_, r_) ** 2
eNu = dexp(2.d0*nu_)
eLamb_hf = eLamb**(0.5d0)
eNu_hf = eNu**(0.5d0)
dLamb = pes03_dLamb(rho_, m_, r_) * 2.d0
dNu = pes03_dnu(P_, rho_, m_, r_) * 2.d0
ddNu = pes03_ddnu(P_, rho_, m_, r_) * 2.d0
dP = -dNu/2.d0 * (P_ + rho_)
cs2 = (P_/(rho_+P_))*gamma_1

YI_Comb = 0.d0
do i = 1,5
	YI_Comb(1:6) = YI_Comb(1:6) + coef(i)*cY(1:6, i)
enddo

write(*,*) "FE interface W, H0", r_*YI_Comb(4)*dP, 0.5d0*(rho_+P_)*r_**2.d0*YI_Comb(1)
write(*,*) "FE interface T2", YI_Comb(6)

call lo_solid_Algebraic_Relation(P_, rho_, m_, nu_, mu_, gamma_1, r_, PP, TT)
write(*,*) "interface r^2P + T1 continuity", r_**2*vec_dot(YI_Comb(1:6), PP(1:6))+ vec_dot(YI_Comb(1:6), TT(1:6)), r_**2/2.d0 * (P_ + rho_) * (H0F)	
write(*,*) "FE interface r^2 dP + T1, H0 terms", r_**2 * vec_dot(PP(1:6), YI_Comb(1:6)) + vec_dot(TT(1:6), YI_Comb(1:6)), - r_**2 /2.d0*(rho_ + P_) *YI_Comb(1)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*** Test
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	endsubroutine lo_bc_Solve_2Side

!Andersson 2011 Matching---------------------------------------------------
	subroutine lo_bc_Solve_W_2Side(H0F, KF, WF, cY, k2_, Coef)
	! solving with integrating W
	! Method is not correct
	use lin_sol_gen_int
	use puls_eqt_set_opt03_ex1
	use lo_eqt
	!	2 components
	!	Quantity Q = {QQ(i)} [Y] {c(i)} ------ {vec} [2nd Tensor] {vec} inner product
	implicit none
	real(8) :: H0F, KF, WF, cY(1:6, 1:5), k2_, Coef(1:5)
	real(8) :: eNu, eLamb, eNu_hf, eLamb_hf, dLamb, dNu, dP, gamma_1, P_, rho_, m_, mu_,  nu_, r_
	real(8) :: pi4, l_1, l_2_hf, ddNu, cs2
	real(8) :: XX(1:6), HH(1:6), KK(1:6), H2_C(1:6), PP(1:6), TT(1:6), PPF
	integer :: mode, i
	real(8) :: YI(1:6, 1:5),  YS(1:6, 1:5), M(1:5, 1:5), Vec(1:5,1), Sol(1:5,1), YS_Comb(1:6), HS, dHS, yR_
		
		YI(1:6,1:5) = cY(1:6,1:5)
		M(1, 1:5) = YI(4, 1:5)
		Vec(1,1) = WF
		M(2, 1:5) = YI(6, 1:5)
		Vec(2,1) = 0.d0
		M(3, 1:5) = YI(1, 1:5)
		Vec(3,1) = H0F
		M(4, 1:5) = YI(3, 1:5)
		Vec(4,1) = KF

		!¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸
		! Solving for Crust Surface BC
		
		YS = 0.d0
		YS(1, 1) = 1.d0
		YS(2, 2) = 1.d0
		YS(3, 3) = 1.d0
		YS(4, 4) = 1.d0
		YS(5, 5) = 1.d0

		!¡¸ T1 + P continuous
		!¡¸ Crust Surface condition
		P_ = Crf(1)
		rho_ = Crf(2)
		m_ = Crf(3)
		!mu_ = mu_Cr(ii) Andersson's mu seems to be off by 1/2, st it cannot be reduced to RCA case
		mu_ = Crf(4)
		nu_ = Crf(5)
		r_ = Crf(6)
		gamma_1 = lo_gamma(P_)

		eLamb = pes03_eLamb(rho_, m_, r_) ** 2
		eNu = dexp(2.d0*nu_)
		eLamb_hf = eLamb**(0.5d0)
		eNu_hf = eNu**(0.5d0)
		dLamb = pes03_dLamb(rho_, m_, r_) * 2.d0
		dNu = pes03_dnu(P_, rho_, m_, r_) * 2.d0
		ddNu = pes03_ddnu(P_, rho_, m_, r_) * 2.d0
		dP = -dNu/2.d0 * (P_ + rho_)
		cs2 = (P_/(rho_+P_))*gamma_1

		call lo_solid_Algebraic_Relation(P_, rho_, m_, nu_, mu_, gamma_1, r_, PP, TT)
		XX = 0.d0
		XX(4) = -0.5d0*(P_+rho_)*dNu/r_
		XX(1:6) = XX(1:6) + PP(1:6)

		M(5, 1:5) = 0.d0
		do i = 1, 6
			M(5, 1:5) = M(5, 1:5) + YS(i, 1:5) * (r_**2.d0*XX(i) + TT(i)) 
		enddo

		Vec(5,1) = 0.d0

		call LIN_SOL_GEN(M, Vec, Sol)					!IMSL linear equations solver
		coef(1:5) = Sol(1:5,1)
		
		YS_Comb = 0.d0
		do i = 1,5
			YS_Comb(1:6) = YS_Comb(1:6) + coef(i)*YS(1:6, i)
		enddo

		dHS = YS_Comb(2)
		HS = YS_Comb(1)

		yR_ = r_*dHS/HS

		call lo_k2(yR_, YS_Comb(5)/HS, P_, rho_, m_, mu_, r_, k2_)

	endsubroutine lo_bc_Solve_W_2Side

!Andersson 2011 Solid Core Fluid Envelop Interface---------------------------------------------------
	subroutine lo_bc_EF_Solve(cY, Coef, H0F, dH0F, KF)
	! solving with Andersson 2011
	use lin_sol_gen_int
	use puls_eqt_set_opt03_ex1
	use lo_eqt
use Root_Finding
	!	2 components
	!	Quantity Q = {QQ(i)} [Y] {c(i)} ------ {vec} [2nd Tensor] {vec} inner product
	implicit none
	real(8) :: H0F, KF, dH0F, cY(1:6,1:3), Coef(1:3)
	real(8) :: eNu, eLamb, eNu_hf, eLamb_hf, dLamb, dNu, dP, gamma_1, P_, rho_, m_, mu_,  nu_, r_
	real(8) :: pi4, l_1, l_2_hf, ddNu, cs2
	real(8) :: HH(1:6), PP(1:6), TT(1:6)
	integer :: mode, i
	real(8) :: YI(1:6, 1:3),  M(1:2, 1:2), Vec(1:2,1), Sol(1:2,1)
real(8) :: k2_, YI_Comb(1:6), PP_Comb, TT_Comb
		!¡¸ Crust Core interface condition
		P_ = Cri(1)
		rho_ = Cri(2)
		m_ = Cri(3)
		!mu_ = mu_Cr(ii) Andersson's mu seems to be off by 1/2, st it cannot be reduced to RCA case
		mu_ = Cri(4)
		nu_ = Cri(5)
		r_ =  Cri(6)
gamma_1 = lo_gamma(P_)

		pi4 = 4.d0*pi
		l_1 = l_0 * (l_0 + 1.d0)
		l_2_hf = (l_0 + 2.d0)*(l_0 - 1.d0)/2.d0

		eLamb = pes03_eLamb(rho_, m_, r_) ** 2
		eNu = dexp(2.d0*nu_)
		eLamb_hf = eLamb**(0.5d0)
		eNu_hf = eNu**(0.5d0)
		dLamb = pes03_dLamb(rho_, m_, r_) * 2.d0
		dNu = pes03_dnu(P_, rho_, m_, r_) * 2.d0
		ddNu = pes03_ddnu(P_, rho_, m_, r_) * 2.d0
		dP = -dNu/2.d0 * (P_ + rho_)
		cs2 = (P_/(rho_+P_))*gamma_1

		!¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸
		! Solving for WE as given in Andersson (2011) Eqt(53)-(55)
		call lo_solid_Algebraic_Relation(P_, rho_, m_, nu_, mu_, gamma_1, r_, PP, TT)
		HH = 0.d0
		HH(1) = 1.d0

		Coef(1) = 1.d0

		YI = 0.d0
		YI(1:6, 1:3) = cY(1:6,1:3)

		M(1, 1:2) = 0.d0
		Vec(1,1) = 0.d0						! Fluid, Zero Freq: P_Eu = 0.5 * (rho+P) * H0_LD
		do i = 1, 6
			M(1, 1) = M(1, 1) + YI(i, 2) * (r_**2.d0*PP(i) + TT(i) - 0.5d0*(rho_+P_)*r_**2.d0*HH(i))
			M(1, 2) = M(1, 2) + YI(i, 3) * (r_**2.d0*PP(i) + TT(i) - 0.5d0*(rho_+P_)*r_**2.d0*HH(i))
			Vec(1,1) = Vec(1,1) - (Coef(1)*YI(i, 1))* (r_**2.d0*PP(i) + TT(i) - 0.5d0*(rho_+P_)*r_**2.d0*HH(i))		! Continuity of (r^2 P_Eu - T1) 
		enddo

		M(2,1:2) = 0.d0
		Vec(2,1) = 0.d0
			M(2, 1) = YI(6, 2)
			M(2, 2) = YI(6, 3)
			Vec(2,1) = -(Coef(1)*YI(6, 1))		! Continuity of T2

		call LIN_SOL_GEN(M, Vec, Sol)					!IMSL linear equations solver
		coef(2) = Sol(1,1)
		coef(3) = Sol(2,1)

		H0F = coef(1)*YI(1, 1) + coef(2)*YI(1, 2) + coef(3)*YI(1, 3)
		KF = coef(1)*YI(3, 1) + coef(2)*YI(3, 2) + coef(3)*YI(3, 3)
		dH0F =( - (l_1*eLamb-2.d0 + (r_*dnu)**2 - r_*(dnu+dLamb))*H0F + (l_1-2.d0)*eLamb*KF) / (r_**2*dnu)  !## (38) of Andersson (2011)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*** Test
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!write(*,*) "EF Jump across Interface result", dH0F - (coef(1)*YI(2, 1) + coef(2)*YI(2, 2) + coef(3)*YI(2, 3))
!write(*,*) "EF Jump across Interface gap",  32.d0*pi*dnu * (mu_*(coef(1)*YI(5, 1) + coef(2)*YI(5, 2) + coef(3)*YI(5, 3)) )
!pause

call lo_k2(r_*(coef(1)*YI(2, 1) + coef(2)*YI(2, 2) + coef(3)*YI(2, 3))/H0F, (coef(1)*YI(5, 1) + coef(2)*YI(5, 2) + coef(3)*YI(5, 3))/H0F, P_, rho_, m_, mu_, r_, k2_)
write(*,*) k2_
YI_Comb = 0.d0
do i = 1,3
	YI_Comb(1:6) = YI_Comb(1:6) + coef(i)*YI(1:6, i)
enddo
write(*,*) "EF interface W, H0", r_*YI_Comb(4)*dP, 0.5d0*(rho_+P_)*r_**2.d0*YI_Comb(1)
write(*,*) "EF interface T2", YI_Comb(6)
write(*,*) "EF interface r^2P + T1", r_**2*vec_dot(YI_Comb(1:6), PP(1:6)) + vec_dot(TT(1:6),YI_Comb(1:6))
!PP_Comb = 0.d0
!TT_Comb = 0.d0
!do i = 1,6
!	PP_Comb = PP_Comb + YI_Comb(i)*PP(i)
!	TT_Comb = TT_Comb + YI_Comb(i)*TT(i)
!enddo
!write(*,*) r_**2 * PP_Comb + TT_Comb
!write(*,*) r_**2 * PP_Comb + TT_Comb + r_*YI_Comb(4)*dP
!write(*,*) YI_Comb(6)
!pause
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*** Test
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	endsubroutine lo_bc_EF_Solve

!!!!!!!! DELETE THIS IF NOT NEEDED
subroutine test(y, r_)
! test the algebraic equation relating H0`, H0 and K
use puls_eqt_set_opt03_ex1
use lo_eqt
implicit none
real(8) :: y(1:6)
real(8) :: eNu, eLamb, eNu_hf, eLamb_hf, dLamb, dNu, dP, gamma_1, P_, rho_, m_, mu_,  nu_, r_
real(8) :: pi4, l_1, l_2_hf, ddNu, cs2
		P_ = y(4)
		rho_ =dsign(lo_EOS(dabs(y(4)), f = 'rho(p)'), y(4))
		if (rho_ < 0.d0) pause 'err: rho_ < 0'
		m_ = y(5)
		nu_ = y(6)
gamma_1 = lo_gamma(P_)

		pi4 = 4.d0*pi
		l_1 = l_0 * (l_0 + 1.d0)
		l_2_hf = (l_0 + 2.d0)*(l_0 - 1.d0)/2.d0
		eLamb = pes03_eLamb(rho_, m_, r_) ** 2
		eNu = dexp(2.d0*nu_)
		eLamb_hf = eLamb**(0.5d0)
		eNu_hf = eNu**(0.5d0)
		dLamb = pes03_dLamb(rho_, m_, r_) * 2.d0
		dNu = pes03_dnu(P_, rho_, m_, r_) * 2.d0
		ddNu = pes03_ddnu(P_, rho_, m_, r_) * 2.d0
		dP = -dNu/2.d0 * (P_ + rho_)
		cs2 = (P_/(rho_+P_))*gamma_1
write(*,*) "dH0F Eqt", ( - (l_1*eLamb-2.d0 + (r_*dnu)**2 - r_*(dnu+dLamb))*y(1) + (l_1-2.d0)*eLamb*y(3)) / (r_**2*dnu)
write(*,*) "dH0F", y(2)
endsubroutine test

!F1 A9 Solid---------------------------------------------------
	subroutine lo_bc_Solid_Solve(cY, k2_, Coef)
	use lin_sol_gen_int
	use puls_eqt_set_opt03_ex1
	use lo_eqt
	!	1 component: solid
	!	Quantity Q = {QQ(i)} [Y] {c(i)} ------ {vec} [2nd Tensor] {vec} inner product
	implicit none
	real(8) :: H0F, KF, cY(1:6,1:3), k2_, Coef(1:3)
	real(8) :: eNu, eLamb, eNu_hf, eLamb_hf, dLamb, dNu, dP, gamma_1, P_, rho_, m_, mu_,  nu_, r_
	real(8) :: pi4, l_1, l_2_hf, ddNu, cs2
	real(8) :: XX(1:6), HH(1:6), KK(1:6), H2_C(1:6), PP(1:6), TT(1:6), PPF
	integer :: mode, i
	real(8) :: YS(1:6, 1:3), M(1:2, 1:2), Vec(1:2,1), Sol(1:2,1), YS_Comb(1:6), HS, dHS, yR_
	!¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸
		! Solving for Crust Surface BC
		YS = 0.d0
		YS(1:6, 1) = cY(1:6, 1)
		YS(1:6, 2) = cY(1:6, 2)
		YS(1:6, 3) = cY(1:6, 3)
		coef(1) = 1.d0
		!¡¸ T1 + P continuous
		!¡¸ Crust Surface condition
		P_ = Crf(1)
		rho_ = Crf(2)
		m_ = Crf(3)
		mu_ = Crf(4)
		nu_ = Crf(5)
		r_ = Crf(6)
gamma_1 = lo_gamma(P_)

		pi4 = 4.d0*pi
		l_1 = l_0 * (l_0 + 1.d0)
		l_2_hf = (l_0 + 2.d0)*(l_0 - 1.d0)/2.d0

		eLamb = pes03_eLamb(rho_, m_, r_) ** 2
		eNu = dexp(2.d0*nu_)
		eLamb_hf = eLamb**(0.5d0)
		eNu_hf = eNu**(0.5d0)
		dLamb = pes03_dLamb(rho_, m_, r_) * 2.d0
		dNu = pes03_dnu(P_, rho_, m_, r_) * 2.d0
		ddNu = pes03_ddnu(P_, rho_, m_, r_) * 2.d0
		dP = -dNu/2.d0 * (P_ + rho_)
		cs2 = (P_/(rho_+P_))*gamma_1

		call lo_solid_Algebraic_Relation(P_, rho_, m_, nu_, mu_, gamma_1, r_, PP, TT)
		XX = 0.d0
		XX(4) = -0.5d0*(P_+rho_)*dNu/r_
		XX(1:6) = XX(1:6) + PP(1:6)

		Vec(1,1) = 0.d0
		M(1, 1:2) = 0.d0
		do i = 1, 6
			M(1, 1) = M(1, 1) + YS(i, 2) * (r_**2.d0*XX(i) + TT(i))
			M(2, 2) = M(1, 2) + YS(i, 3) * (r_**2.d0*XX(i) + TT(i))
			Vec(1,1) = Vec(1,1) - (Coef(1) * YS(i, 1))* (r_**2.d0*XX(i) + TT(i))
		enddo

		!¡¸ T2 = 0
		M(2, 1) = YS(6, 2)
		M(2, 2) = YS(6, 3)
		Vec(2,1) = -YS(6, 1)

		call LIN_SOL_GEN(M, Vec, Sol)					!IMSL linear equations solver
		coef(2) = Sol(1,1)
		coef(3) = Sol(2,1)
		
		YS_Comb = 0.d0
		do i = 1,3
			YS_Comb(1:6) = YS_Comb(1:6) + coef(i)*YS(1:6, i)
		enddo
		
		dHS = YS_Comb(2)
		HS = YS_Comb(1)

		yR_ = r_*dHS/HS
		call lo_k2(yR_, YS_Comb(5)/HS, P_, rho_, m_, mu_, r_, k2_)

	endsubroutine lo_bc_Solid_Solve

!Andersson 2011 Matching---------------------------------------------------
	subroutine lo_bc_Solid_Solve_2Side(cY, k2_, Coef)
	! solving with Andersson 2011
	use lin_sol_gen_int
	use puls_eqt_set_opt03_ex1
	use lo_eqt
	use Root_Finding
	!	2 components
	!	Quantity Q = {QQ(i)} [Y] {c(i)} ------ {vec} [2nd Tensor] {vec} inner product
	implicit none
	real(8) :: cY(1:6,1:8), k2_, Coef(1:8)
	real(8) :: eNu, eLamb, eNu_hf, eLamb_hf, dLamb, dNu, dP, gamma_1, P_, rho_, m_, mu_,  nu_, r_
	real(8) :: pi4, l_1, l_2_hf, ddNu, cs2
	real(8) :: XX(1:6), HH(1:6), KK(1:6), H2_C(1:6), PP(1:6), TT(1:6), PPF
	integer :: mode, i
	real(8) :: YI(1:6, 1:8),  YS(1:6, 1:8), M(1:7, 1:7), Vec(1:7,1), Sol(1:7,1), YS_Comb(1:6), HS, dHS, yR_
real(8) :: dP_fluid

		Coef(1) = 1.d0

		!¡¸ Crust Core interface condition
		P_ = Cri(1)
		rho_ = Cri(2)
		m_ = Cri(3)
		!mu_ = mu_Cr(ii) Andersson's mu seems to be off by 1/2, st it cannot be reduced to RCA case
		mu_ = Cri(4)
		nu_ = Cri(5)
		r_ =  Cri(6)
		gamma_1 = lo_gamma(P_)

		pi4 = 4.d0*pi
		l_1 = l_0 * (l_0 + 1.d0)
		l_2_hf = (l_0 + 2.d0)*(l_0 - 1.d0)/2.d0

		eLamb = pes03_eLamb(rho_, m_, r_) ** 2
		eNu = dexp(2.d0*nu_)
		eLamb_hf = eLamb**(0.5d0)
		eNu_hf = eNu**(0.5d0)
		dLamb = pes03_dLamb(rho_, m_, r_) * 2.d0
		dNu = pes03_dnu(P_, rho_, m_, r_) * 2.d0
		ddNu = pes03_ddnu(P_, rho_, m_, r_) * 2.d0
		dP = -dNu/2.d0 * (P_ + rho_)
		cs2 = (P_/(rho_+P_))*gamma_1

		!¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸
		! Solving for WE as given in Andersson (2011) Eqt(53)-(55)
		call lo_solid_Algebraic_Relation(P_, rho_, m_, nu_, mu_, gamma_1, r_, PP, TT)
		XX = 0.d0
		XX(4) = -0.5d0*(P_+rho_)*dNu/r_
		XX(1:6) = XX(1:6) + PP(1:6)
		!¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸
		! Continuity at CC Interface
		YI(1:6, 1:8) = cY(1:6,1:8)

			M(1,1) = - vec_dot(YI(1:6, 2),XX(1:6))	! Lagrangian Pressure
			M(1,2) = - vec_dot(YI(1:6, 3),XX(1:6))
			M(1,3) = vec_dot(YI(1:6, 4),XX(1:6))
			M(1,4) = vec_dot(YI(1:6, 5),XX(1:6))
			M(1,5) = vec_dot(YI(1:6, 6),XX(1:6))
			M(1,6) = vec_dot(YI(1:6, 7),XX(1:6))
			M(1,7) = vec_dot(YI(1:6, 8),XX(1:6))
			Vec(1,1) = Coef(1)*vec_dot(YI(1:6, 1),XX(1:6))
			M(2, 1:2) = -YI(1, 2:3)  ! H0
			M(2, 3:7) = YI(1, 4:8)
			Vec(2,1) = (Coef(1)*YI(1, 1))
			M(3:6, 1:2) = -YI(3:6, 2:3)  ! K, W, V, T2
			M(3:6, 3:7) = YI(3:6, 4:8)
			Vec(3:6,1) = (Coef(1)*YI(3:6, 1))
			!M(1:6, 1:2) = -YI(1:6, 2:3)
			!M(1:6, 3:7) = YI(1:6, 4:8)
			!Vec(1:6,1) = (Coef(1)*YI(1:6, 1))

		!¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸
		! Solving for Crust Surface BC
		
		YS = 0.d0
		YS(1, 4) = 1.d0
		YS(2, 5) = 1.d0
		YS(3, 6) = 1.d0
		YS(4, 7) = 1.d0
		YS(5, 8) = 1.d0

		!¡¸ T1 + P continuous
		!¡¸ Crust Surface condition
		P_ = Crf(1)
		rho_ = Crf(2)
		m_ = Crf(3)
		!mu_ = mu_Cr(ii) Andersson's mu seems to be off by 1/2, st it cannot be reduced to RCA case
		mu_ = Crf(4)
		nu_ = Crf(5)
		r_ = Crf(6)
gamma_1 = lo_gamma(P_)

		pi4 = 4.d0*pi
		l_1 = l_0 * (l_0 + 1.d0)
		l_2_hf = (l_0 + 2.d0)*(l_0 - 1.d0)/2.d0

		eLamb = pes03_eLamb(rho_, m_, r_) ** 2
		eNu = dexp(2.d0*nu_)
		eLamb_hf = eLamb**(0.5d0)
		eNu_hf = eNu**(0.5d0)
		dLamb = pes03_dLamb(rho_, m_, r_) * 2.d0
		dNu = pes03_dnu(P_, rho_, m_, r_) * 2.d0
		ddNu = pes03_ddnu(P_, rho_, m_, r_) * 2.d0
		dP = -dNu/2.d0 * (P_ + rho_)
		cs2 = (P_/(rho_+P_))*gamma_1

		call lo_solid_Algebraic_Relation(P_, rho_, m_, nu_, mu_, gamma_1, r_, PP, TT)
		XX = 0.d0
		XX(4) = -0.5d0*(P_+rho_)*dNu/r_
		XX(1:6) = XX(1:6) + PP(1:6)

		Vec(7,1) = 0.d0
		M(7, 1:7) = 0.d0
		do i = 1, 6
			M(7, 3:7) = M(7, 3:7) + YS(i, 4:8) * (r_**2.d0*XX(i) + TT(i))
		enddo
!M(7, 3:7) = M(7, 3:7) - r_/2.d0*(P_+rho_)*(r_*YS(1,4:8) - dNu*YS(4,4:8)) !Fluid surface BC not satisfying XX=0

		call LIN_SOL_GEN(M, Vec, Sol)					!IMSL linear equations solver
		coef(2:3) = Sol(1:2,1)
		coef(4:8) = Sol(3:7,1)

		YS_Comb = 0.d0
		do i = 1,8
			YS_Comb(1:6) = YS_Comb(1:6) + coef(i)*YS(1:6, i)
		enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*** Test
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dP_fluid = 0.d0
do i = 1, 6
	dP_fluid = dP_fluid + YS_Comb(i) * PP(i)
enddo
write(*,*) "dP_Fluid, H0", dP_fluid*r_**2, -0.5d0*(rho_+P_)*r_**2.d0*YS_Comb(1)
write(*,*) "W, H0", r_*YS_Comb(4)*dP, 0.5d0*(rho_+P_)*r_**2.d0*YS_Comb(1)
write(*,*) "dP", r_*YS_Comb(4)*dP+ 0.5d0*(rho_+P_)*r_**2.d0*YS_Comb(1)		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*** Test
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		dHS = YS_Comb(2)
		HS = YS_Comb(1)

		yR_ = r_*dHS/HS

		call lo_k2(yR_, YS_Comb(5)/HS, P_, rho_, m_, mu_, r_, k2_)

	endsubroutine lo_bc_Solid_Solve_2Side

!Andersson 2011 Matching---------------------------------------------------
	subroutine lo_bc_Solid_Solve_2Side_muGap(cY, k2_, Coef)
	! solving with Andersson 2011
	! able to due with 2 Solid models, with density gap/ shear modulus gap
	use lin_sol_gen_int
	use puls_eqt_set_opt03_ex1
	use lo_eqt
use Format_IO
use FWrite
use Root_Finding
	!	2 components
	!	Quantity Q = {QQ(i)} [Y] {c(i)} ------ {vec} [2nd Tensor] {vec} inner product
	implicit none
	real(8) :: cY(1:6,1:8), k2_, Coef(1:8)
	real(8) :: eNu, eLamb, eNu_hf, eLamb_hf, dLamb, dNu, dP, gamma_1, P_, rho_, m_, mu_,  nu_, r_
	real(8) :: pi4, l_1, l_2_hf, ddNu, cs2
	real(8) :: XX(1:6), HH(1:6), KK(1:6), H2_C(1:6), PP(1:6), TT(1:6), PPF
	integer :: mode, i
	real(8) :: YI(1:6, 1:8),  YS(1:6, 1:8), M(1:7, 1:7), Vec(1:7,1), Sol(1:7,1), YS_Comb(1:6), HS, dHS, yR_
real(8) :: YI_Comb_Co(1:6), YI_Comb_Cr(1:6)

		Coef(1) = 1.d0

		!¡¸ Crust Core interface condition
		P_ = Cri(1)
		rho_ = Cri(2)
		m_ = Cri(3)
		!mu_ = mu_Cr(ii) Andersson's mu seems to be off by 1/2, st it cannot be reduced to RCA case
		mu_ = Cri(4)
		nu_ = Cri(5)
		r_ =  Cri(6)
		gamma_1 = lo_gamma(P_)

		pi4 = 4.d0*pi
		l_1 = l_0 * (l_0 + 1.d0)
		l_2_hf = (l_0 + 2.d0)*(l_0 - 1.d0)/2.d0

		eLamb = pes03_eLamb(rho_, m_, r_) ** 2
		eNu = dexp(2.d0*nu_)
		eLamb_hf = eLamb**(0.5d0)
		eNu_hf = eNu**(0.5d0)
		dLamb = pes03_dLamb(rho_, m_, r_) * 2.d0
		dNu = pes03_dnu(P_, rho_, m_, r_) * 2.d0
		ddNu = pes03_ddnu(P_, rho_, m_, r_) * 2.d0
		dP = -dNu/2.d0 * (P_ + rho_)
		cs2 = (P_/(rho_+P_))*gamma_1

		!¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸
		! Solving for WE as given in Andersson (2011) Eqt(53)-(55)
		call lo_solid_Algebraic_Relation(P_, rho_, m_, nu_, mu_, gamma_1, r_, PP, TT)
		XX = 0.d0
		XX(4) = -0.5d0*(P_+rho_)*dNu/r_
		XX(1:6) = XX(1:6) + PP(1:6)
		!¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸
		! Continuity at CC Interface
		YI(1:6, 1:8) = cY(1:6,1:8)

		M(1,1) = - vec_dot(YI(1:6, 2),(r_**2.d0*XX(1:6) + TT(1:6)))	! Lagrangian Pressure
		M(1,2) = - vec_dot(YI(1:6, 3),(r_**2.d0*XX(1:6) + TT(1:6)))
		Vec(1,1) = Coef(1)*vec_dot(YI(1:6, 1),(r_**2.d0*XX(1:6) + TT(1:6)))
		!¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸

		!¡¸ Crust Core interface condition
		P_ = Cri_out(1)
		rho_ = Cri_out(2)
		m_ = Cri_out(3)
		mu_ = Cri_out(4)
		nu_ = Cri_out(5)
		r_ =  Cri_out(6)
		gamma_1 = lo_gamma(P_)

		pi4 = 4.d0*pi
		l_1 = l_0 * (l_0 + 1.d0)
		l_2_hf = (l_0 + 2.d0)*(l_0 - 1.d0)/2.d0

		eLamb = pes03_eLamb(rho_, m_, r_) ** 2
		eNu = dexp(2.d0*nu_)
		eLamb_hf = eLamb**(0.5d0)
		eNu_hf = eNu**(0.5d0)
		dLamb = pes03_dLamb(rho_, m_, r_) * 2.d0
		dNu = pes03_dnu(P_, rho_, m_, r_) * 2.d0
		ddNu = pes03_ddnu(P_, rho_, m_, r_) * 2.d0
		dP = -dNu/2.d0 * (P_ + rho_)
		cs2 = (P_/(rho_+P_))*gamma_1

		!¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸
		! Solving for WE as given in Andersson (2011) Eqt(53)-(55)
		call lo_solid_Algebraic_Relation(P_, rho_, m_, nu_, mu_, gamma_1, r_, PP, TT)
		XX = 0.d0
		XX(4) = -0.5d0*(P_+rho_)*dNu/r_
		XX(1:6) = XX(1:6) + PP(1:6)

			M(1,3) = vec_dot(YI(1:6, 4),(r_**2.d0*XX(1:6) + TT(1:6)))
			M(1,4) = vec_dot(YI(1:6, 5),(r_**2.d0*XX(1:6) + TT(1:6)))
			M(1,5) = vec_dot(YI(1:6, 6),(r_**2.d0*XX(1:6) + TT(1:6)))
			M(1,6) = vec_dot(YI(1:6, 7),(r_**2.d0*XX(1:6) + TT(1:6)))
			M(1,7) = vec_dot(YI(1:6, 8),(r_**2.d0*XX(1:6) + TT(1:6)))
			
			M(2, 1:2) = -YI(1, 2:3)  ! H0
			M(2, 3:7) = YI(1, 4:8)
			Vec(2,1) = (Coef(1)*YI(1, 1))
			M(3:6, 1:2) = -YI(3:6, 2:3)  ! K, W, V, T2
			M(3:6, 3:7) = YI(3:6, 4:8)
			Vec(3:6,1) = (Coef(1)*YI(3:6, 1))
			!M(1:6, 1:2) = -YI(1:6, 2:3)
			!M(1:6, 3:7) = YI(1:6, 4:8)
			!Vec(1:6,1) = (Coef(1)*YI(1:6, 1))
! Un-common to employ Free slip
!M(5,1:2) = 0.d0
!M(5,3:7) = YI(6,4:8)
!Vec(5,1) = 0.d0
!******************************

		!¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸
		! Solving for Crust Surface BC
		
		YS = 0.d0
		YS(1, 4) = 1.d0
		YS(2, 5) = 1.d0
		YS(3, 6) = 1.d0
		YS(4, 7) = 1.d0
		YS(5, 8) = 1.d0

		!¡¸ T1 + P continuous
		!¡¸ Crust Surface condition
		P_ = Crf(1)
		rho_ = Crf(2)
		m_ = Crf(3)
		!mu_ = mu_Cr(ii) Andersson's mu seems to be off by 1/2, st it cannot be reduced to RCA case
		mu_ = Crf(4)
		nu_ = Crf(5)
		r_ = Crf(6)
gamma_1 = lo_gamma(P_)

		eLamb = pes03_eLamb(rho_, m_, r_) ** 2
		eNu = dexp(2.d0*nu_)
		eLamb_hf = eLamb**(0.5d0)
		eNu_hf = eNu**(0.5d0)
		dLamb = pes03_dLamb(rho_, m_, r_) * 2.d0
		dNu = pes03_dnu(P_, rho_, m_, r_) * 2.d0
		ddNu = pes03_ddnu(P_, rho_, m_, r_) * 2.d0
		dP = -dNu/2.d0 * (P_ + rho_)
		cs2 = (P_/(rho_+P_))*gamma_1

		call lo_solid_Algebraic_Relation(P_, rho_, m_, nu_, mu_, gamma_1, r_, PP, TT)
		XX = 0.d0
		XX(4) = -0.5d0*(P_+rho_)*dNu/r_
		XX(1:6) = XX(1:6) + PP(1:6)

		Vec(7,1) = 0.d0
		M(7, 1:7) = 0.d0
		do i = 1, 6
			M(7, 3:7) = M(7, 3:7) + YS(i, 4:8) * (r_**2.d0*XX(i) + TT(i))
		enddo

		call LIN_SOL_GEN(M, Vec, Sol)					!IMSL linear equations solver
		coef(2:3) = Sol(1:2,1)
		coef(4:8) = Sol(3:7,1)

		YS_Comb = 0.d0
		do i = 1,8
			YS_Comb(1:6) = YS_Comb(1:6) + coef(i)*YS(1:6, i)
		enddo
		
		dHS = YS_Comb(2)
		HS = YS_Comb(1)

		yR_ = r_*dHS/HS

		call lo_k2(yR_, YS_Comb(5)/HS, P_, rho_, m_, mu_, r_, k2_)
! for plotting yR_out
solR_temp = yR_- 4.d0*pi*rho_*r_**3/m_ + 32.d0*pi*	(dNu* r_)* (mu_*YS_Comb(5)/HS)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*** Test
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
P_ = Cri(1)
rho_ = Cri(2)
m_ = Cri(3)
mu_ = Cri(4)
nu_ = Cri(5)
r_ =  Cri(6)
gamma_1 = lo_gamma(P_)

pi4 = 4.d0*pi
l_1 = l_0 * (l_0 + 1.d0)
l_2_hf = (l_0 + 2.d0)*(l_0 - 1.d0)/2.d0

eLamb = pes03_eLamb(rho_, m_, r_) ** 2
eNu = dexp(2.d0*nu_)
eLamb_hf = eLamb**(0.5d0)
eNu_hf = eNu**(0.5d0)
dLamb = pes03_dLamb(rho_, m_, r_) * 2.d0
dNu = pes03_dnu(P_, rho_, m_, r_) * 2.d0
ddNu = pes03_ddnu(P_, rho_, m_, r_) * 2.d0
dP = -dNu/2.d0 * (P_ + rho_)
cs2 = (P_/(rho_+P_))*gamma_1

call lo_solid_Algebraic_Relation(P_, rho_, m_, nu_, mu_, gamma_1, r_, PP, TT)
do i = 1,6
YI_Comb_Co(i) = vec_dot(Coef(1:3), YI(i,1:3))
YI_Comb_Cr(i) = vec_dot(Coef(4:8), YI(i,4:8))
enddo
write(*,*) "interface r^2 dP + T1, H0 terms, Co side", r_**2 * vec_dot(PP(1:6), YI_Comb_Co(1:6)) + vec_dot(TT(1:6), YI_Comb_Co(1:6)), - r_**2 /2.d0*(rho_ + P_) * YI_Comb_Co(1)
write(*,*) "interface r^2 dP + T1, H0 terms, Cr side", r_**2 * vec_dot(PP(1:6), YI_Comb_Cr(1:6)) + vec_dot(TT(1:6), YI_Comb_Cr(1:6)), - r_**2 /2.d0*(rho_ + P_) * YI_Comb_Cr(1)
!write(*,*) "interface T1, Co side",vec_dot(TT(1:6), YI_Comb_Co(1:6))
!write(*,*) "interface T1, Cr side",vec_dot(TT(1:6), YI_Comb_Cr(1:6))
!write(*,*) "interface T1 + PP, Co side",(vec_dot(TT(1:6), YI_Comb_Co(1:6)) + Cri(6)**2 * vec_dot(PP(1:6), YI_Comb_Co(1:6)))/YI_Comb_Co(1)
!write(*,*) "interface T1 + PP, Cr side",(vec_dot(TT(1:6), YI_Comb_Cr(1:6)) + Cri_Out(6)**2 * vec_dot(PP(1:6), YI_Comb_Cr(1:6)))/YI_Comb_Cr(1)
write(*,*) "T2, Co side", YI_Comb_Co(6)
write(*,*) "T2, Cr side", YI_Comb_Cr(6)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*** Test
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	endsubroutine lo_bc_Solid_Solve_2Side_muGap

	!Tidal Love Number
	subroutine lo_k2(y, V_H0, P_, rho_, m_, mu_, r_, k2_)
	use puls_eqt_set_opt03_ex1
	implicit none
	real(8) :: y, V_H0, k2_
	real(8) :: m_, r_, P_, rho_, mu_, dNu
	real(8) :: yR_, T_(1:3)
	real(8) :: C_

		dNu = pes03_dnu(P_, rho_, m_, r_) * 2.d0

		yR_ = y - 4.d0*pi*rho_*r_**3/m_ + 32.d0*pi*	(dNu* r_)* (mu_*V_H0) 		! Jump across surface-vaccuum interface ! negative sign: Samson
		! y+ = y- - 4pi rho R^3/M rho(R-) + 32 pi dNu mu (rV/H0); mu is using Andersson 2011's mu, which is 2 times larger than the correct value
		C_ = m_/r_
		T_(1) = 8.d0/5.d0*C_**5*(1.d0-2.d0*C_)**2 * (2.d0*C_*(yR_-1.d0)+2.d0-yR_)
		T_(2) = 2.d0*C_*( 4.d0*(yR_+1.d0)*C_**4 + (6.d0*yR_ -4.d0)*C_**3 + (26.d0 - 22.d0*yR_)*C_**2 + 3.d0*C_*(5.d0*yR_ - 8.d0) -3.d0*yR_ + 6.d0 )
		T_(3) = 3.d0*(1.d0-2.d0*C_)**2 * (2.d0*C_*(yR_ - 1.d0) - yR_ +2.d0) * dlog(1.d0-2.d0*C_)

		k2_ = T_(1) /(T_(2) + T_(3) )
		
	endsubroutine lo_k2


endmodule lo_match