module lo_match_P02
! Love number matching Page 02
use lo_match
use global_var
implicit none
contains


!F1 A9 FEF model---------------------------------------------------
subroutine lo_bc_FEF_Solve(H0_Co, K_Co, cY, Coef, H0F, dH0F, KF)
	! solving with Andersson 2011
	use lin_sol_gen_int
	use puls_eqt_set_opt03_ex1
	use lo_eqt
	!	2 components
	!	Quantity Q = {QQ(i)} [Y] {c(i)} ------ {vec} [2nd Tensor] {vec} inner product
	implicit none
	real(8), intent(in) :: H0_Co, K_Co, cY(1:6,1:5)
	real(8), intent(out) ::H0F, KF, dH0F, Coef(1:5)
	real(8) :: eNu, eLamb, eNu_hf, eLamb_hf, dLamb, dNu, dP, gamma_1, P_, rho_, m_, mu_,  nu_, r_
	real(8) :: pi4, l_1, l_2_hf, ddNu, cs2
	real(8) :: HH(1:6), KK(1:6), H2_C(1:6), PP(1:6), TT(1:6), PP_Co, k2_ 
	integer :: mode, i
	real(8) :: YI(1:6, 1:5),  YF(1:6, 1:5), M(1:3, 1:3), Vec(1:3,1), Sol(1:3,1), YF_Comb(1:6)
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
		PP_Co = r_**2/2.d0 * (P_ + rho_) * (H0_Co)		

		Coef(1) = H0_Co
		Coef(3) = K_Co

		YI = 0.d0
		YI(1, 1) = 1.d0
		YI(2, 2) = 1.d0
		YI(3, 3) = 1.d0
		YI(4, 4) = 1.d0
		YI(5, 5) = 1.d0

		M(1, 1:3) = 0.d0
		Vec(1,1) = PP_Co							! Fluid, Zero Freq: P_Eu = -0.5 * (rho+P) * H0_LD
		do i = 1, 6
			M(1, 1) = M(1, 1) + YI(i, 2) * (r_**2.d0*PP(i) + TT(i))
			M(1, 2) = M(1, 2) + YI(i, 4) * (r_**2.d0*PP(i) + TT(i))
			M(1, 3) = M(1, 3) + YI(i, 5) * (r_**2.d0*PP(i) + TT(i))
			Vec(1,1) = Vec(1,1) - (Coef(1)*YI(i, 1) + Coef(3)*YI(i, 3) )* (r_**2.d0*PP(i) + TT(i))		! Continuity of (r^2 P_Eu - T1) 
		enddo

		!¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸
		! Solving for Crust Surface BC
		
		YF = 0.d0
		YF(1:6, 1) = cY(1:6, 1)
		YF(1:6, 2) = cY(1:6, 2)
		YF(1:6, 3) = cY(1:6, 3)
		YF(1:6, 4) = cY(1:6, 4)
		YF(1:6, 5) = cY(1:6, 5)

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
		HH = 0.d0
		HH(1) = 1.d0

		Vec(2,1) = 0.d0
		M(2, 1:3) = 0.d0
		do i = 1, 6
			M(2, 1) = M(2, 1) + YF(i, 2) * (r_**2.d0*PP(i) + TT(i) - 0.5d0*(rho_+P_)*r_**2.d0*HH(i))
			M(2, 2) = M(2, 2) + YF(i, 4) * (r_**2.d0*PP(i) + TT(i) - 0.5d0*(rho_+P_)*r_**2.d0*HH(i))
			M(2, 3) = M(2, 3) + YF(i, 5) * (r_**2.d0*PP(i) + TT(i) - 0.5d0*(rho_+P_)*r_**2.d0*HH(i))
			Vec(2,1) = Vec(2,1) - (Coef(1) * YF(i, 1) + Coef(3) * YF(i, 3))* (r_**2.d0*PP(i) + TT(i) - 0.5d0*(rho_+P_)*r_**2.d0*HH(i))
		enddo

		!¡¸ T2 = 0
		M(3, 1) = YF(6, 2)
		M(3, 2) = YF(6, 4)
		M(3, 3) = YF(6, 5)
		Vec(3,1) = 0.d0 - (  Coef(1) * YF(6, 1) + Coef(3) * YF(6, 3)  )

		call LIN_SOL_GEN(M, Vec, Sol)					!IMSL linear equations solver
		coef(2) = Sol(1,1)
		coef(4) = Sol(2,1)
		coef(5) = Sol(3,1)
		
		YF_Comb = 0.d0
		do i = 1,5
			YF_Comb(1:6) = YF_Comb(1:6) + coef(i)*YF(1:6, i)
		enddo

		H0F = YF_Comb(1)
		KF = YF_Comb(3)
		dH0F =( - (l_1*eLamb-2.d0 + (r_*dnu)**2 - r_*(dnu+dLamb))*H0F + (l_1-2.d0)*eLamb*KF) / (r_**2*dnu)  !## (38) of Andersson (2011)

	endsubroutine lo_bc_FEF_Solve


!F1 A10 Fluid model W---------------------------------------------------
	subroutine lo_bc_Fluid_W_Solve(cY, k2_, Coef)
	! solving fluid with displacement variable W
	use lin_sol_gen_int
	use puls_eqt_set_opt03_ex1
	use lo_eqt
	!	2 components
	!	Quantity Q = {QQ(i)} [Y] {c(i)} ------ {vec} [2nd Tensor] {vec} inner product
	implicit none
	real(8) :: cY(1:5,1:2), k2_, Coef(1:2)
	real(8) :: eNu, eLamb, eNu_hf, eLamb_hf, dLamb, dNu, dP, gamma_1, P_, rho_, m_, mu_,  nu_, r_
	real(8) :: pi4, l_1, l_2_hf, ddNu, cs2
	real(8) :: XX(1:6), HH(1:6), KK(1:6), H2_C(1:6), PP(1:6), TT(1:6), PPF
	integer :: mode, i
	real(8) :: YI(1:6, 1:8),  YS(1:6, 1:8), M(1:7, 1:7), Vec(1:7,1), Sol(1:7,1), YS_Comb(1:6), HS, dHS, yR_
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

		Coef(1) = 1.d0
		Coef(2) = -(cY(1,1) - dNu/r_*cY(4,1))/(cY(1,2) - dNu/r_*cY(4,2)) * Coef(1)

		dHS = Coef(1)* cY(2,1)+Coef(2)* cY(2,2)
		HS = Coef(1)* cY(1,1)+Coef(2)* cY(1,2)
		

		yR_ = r_*dHS/HS
write(*,*) "Fluid W model H0", Coef(1)*cY(1,1)
write(*,*) "Fluid W model W", (Coef(1)*cY(4,1) + Coef(2)*cY(4,2)) * (-dNu)/r_
		call lo_k2(yR_, 0.d0, P_, rho_, m_, 0.d0, r_, k2_)
	endsubroutine lo_bc_Fluid_W_Solve

	!F1 A10 Fluid model W---------------------------------------------------
	subroutine lo_bc_EF_W_Solve(cZ, cY, k2_, Coef)
	use lin_sol_gen_int
	use puls_eqt_set_opt03_ex1
	use lo_eqt
	implicit none
	real(8),intent(in) :: cZ(1:6,1:3), cY(1:5,4:8)
	real(8),intent(out) :: k2_, Coef(1:8)
	real(8) :: YI(1:6, 1:8), YS(1:5,4:8) 
	real(8) :: eNu, eLamb, eNu_hf, eLamb_hf, dLamb, dNu, dP, gamma_1, P_, rho_, m_, mu_,  nu_, r_
	real(8) :: pi4, l_1, l_2_hf, ddNu, cs2
	real(8) :: M(1:7, 1:7), Vec(1:7,1), Sol(1:7,1), YS_Comb(1:5), XX(1:6), XX_F(1:6), PP(1:6), TT(1:6), dHH_eqt(1:6)
	real(8) :: yR_, HS, dHS
	integer :: i
		
		Coef = 0.d0
		Coef(1) = 1.d0
		YI = 0.d0
		YS = 0.d0
		YI(1:6, 1:3) = cZ(1:6,1:3)
		YI(1:5, 4:8) = cY(1:5,4:8)
		YS(1, 4) = 1.d0
		YS(2, 5) = 1.d0
		YS(3, 6) = 1.d0
		YS(4, 7) = 1.d0
		YS(5, 8) = 1.d0

		!¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸
		! Interface conditions
		P_ = Cri(1)
		rho_ = Cri(2)
		m_ = Cri(3)
		mu_ = Cri(4)
		nu_ = Cri(5)
		r_ = Cri(6)
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
		XX_F = 0.d0
		XX_F(1) = 0.5d0*(P_+rho_)
		XX_F(4) = -0.5d0*(P_+rho_)*dNu/r_

		! H0 Continuity
		M(1, 1:2) = -YI(1, 2:3)
		M(1, 3:7) = YI(1, 4:8)
		Vec(1,1) = Coef(1) * YI(1, 1)
		
		! H0' condition: Andersson (2011) (38) (45)
		dHH_eqt = 0.d0
		dHH_eqt(1) = (l_1*eLamb-2.d0 + (r_*dnu)**2 - r_*(dnu+dLamb))
		dHH_eqt(2) =  (r_**2*dnu)
		dHH_eqt(3) = -(l_1-2.d0)*eLamb
		M(2, 1:7) = 0.d0
		Vec(2,1) = 0.d0
		do i = 1, 5
			M(2,3:7) = M(2,3:7) + YI(i, 4:8) * (dHH_eqt(i))
		enddo

		! K, W Continuity
		M(3:4, 1:2) = -YI(3:4, 2:3)
		M(3:4, 3:7) = YI(3:4, 4:8)
		Vec(3:4,1) = Coef(1) * YI(3:4, 1)

		! T2 = 0
		M(5, 1:2) = -YI(6, 2:3)
		M(5, 3:7) = 0.d0
		Vec(5, 1) = Coef(1) * YI(6, 1)

		!¡¸T1 + P continuous

		Vec(6,1) = 0.d0
		M(6, 1:7) = 0.d0
		do i = 1, 6
			M(6, 1:2) = M(6, 1:2) - YI(i, 2:3) * (r_**2.d0*XX(i) + TT(i))
			Vec(6,1) = Vec(6,1) + Coef(1) * YI(i, 1) * (r_**2.d0*XX(i) + TT(i))
		enddo
		do i = 1,5
			M(6, 3:7) = M(6, 3:7) + YI(i, 4:8) * (r_**2.d0*XX_F(i))
		enddo

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

		XX = 0.d0
		XX(1) = (rho_+P_)/2.d0
		XX(4) = - (rho_+P_)/2.d0* dNu/r_
		
		M(7, 1:7) = 0.d0
		Vec(7,1) = 0.d0
		do i = 1, 5
			M(7,3:7) = M(7,3:7) + YS(i, 4:8) * (r_**2 * XX(i))
		enddo

		call LIN_SOL_GEN(M, Vec, Sol)					!IMSL linear equations solver
		coef(2:8) = Sol(1:7,1)

		YS_Comb = 0.d0
		do i = 4,8
			YS_Comb(1:5) = YS_Comb(1:5) + coef(i)*YS(1:5, i)
		enddo

		dHS = YS_Comb(2)
		HS = YS_Comb(1)

		yR_ = r_*dHS/HS

		call lo_k2(yR_, YS_Comb(5)/HS, P_, rho_, m_, 0.d0, r_, k2_)

	endsubroutine lo_bc_EF_W_Solve

	!F1 A10 Fluid model W---------------------------------------------------
	subroutine lo_bc_EF_W_Solve_No_Slip(cZ, cY, k2_, Coef)
	use lin_sol_gen_int
	use puls_eqt_set_opt03_ex1
	use lo_eqt
	implicit none
	real(8),intent(in) :: cZ(1:6,1:3), cY(1:5,4:8)
	real(8),intent(out) :: k2_, Coef(1:8)
	real(8) :: YI(1:6, 1:8), YS(1:5,4:8) 
	real(8) :: eNu, eLamb, eNu_hf, eLamb_hf, dLamb, dNu, dP, gamma_1, P_, rho_, m_, mu_,  nu_, r_
	real(8) :: pi4, l_1, l_2_hf, ddNu, cs2
	real(8) :: M(1:7, 1:7), Vec(1:7,1), Sol(1:7,1), YS_Comb(1:5), XX(1:6), XX_F(1:6), PP(1:6), TT(1:6), dHH_eqt(1:6)
	real(8) :: yR_, HS, dHS
	integer :: i
		
		Coef = 0.d0
		Coef(1) = 1.d0
		YI = 0.d0
		YS = 0.d0
		YI(1:6, 1:3) = cZ(1:6,1:3)
		YI(1:5, 4:8) = cY(1:5,4:8)
		YS(1, 4) = 1.d0
		YS(2, 5) = 1.d0
		YS(3, 6) = 1.d0
		YS(4, 7) = 1.d0
		YS(5, 8) = 1.d0

		!¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸
		! Interface conditions
		P_ = Cri(1)
		rho_ = Cri(2)
		m_ = Cri(3)
		mu_ = Cri(4)
		nu_ = Cri(5)
		r_ = Cri(6)
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
		XX_F = 0.d0
		XX_F(1) = 0.5d0*(P_+rho_)
		XX_F(4) = -0.5d0*(P_+rho_)*dNu/r_

		! H0 Continuity
		M(1, 1:2) = -YI(1, 2:3)
		M(1, 3:7) = YI(1, 4:8)
		Vec(1,1) = Coef(1) * YI(1, 1)
		
		! H0' condition: Andersson (2011) (38) (45)
		!####################[Bug: T2 in flluid??... need this in No Slip Condition]
		dHH_eqt = 0.d0
		!dHH_eqt(1) = (l_1*eLamb-2.d0 + (r_*dnu)**2 - r_*(dnu+dLamb))
		dHH_eqt(1) = (l_1*eLamb-2.d0 + (r_*dnu)**2)
		dHH_eqt(2) =  (r_**2*dnu)
		dHH_eqt(3) = -(l_1-2.d0)*eLamb
		M(2, 1:7) = 0.d0
		M(2,1:2) = - 8.d0*pi*(r_*dNu + 2.d0)*YI(6,2:3)
		Vec(2,1) = 8.d0*pi*(r_*dNu + 2.d0)*Coef(1)*YI(6,1)
		do i = 1,6
			M(2,1:2) = M(2,1:2) + YI(i,2:3) * ( - 16.d0*pi*eLamb*(r_**2.d0*PP(i) + TT(i)) )
			Vec(2,1) = Vec(2,1) - Coef(1)*YI(2,1)* ( - 16.d0*pi*eLamb*(r_**2.d0*PP(i) + TT(i)) )
		enddo
		do i = 1, 5
			M(2,3:7) = M(2,3:7) + YI(i, 4:8) * (dHH_eqt(i))
		enddo

		! K, W, V Continuity
		M(3:5, 1:2) = -YI(3:5, 2:3)
		M(3:5, 3:7) = YI(3:5, 4:8)
		Vec(3:5,1) = Coef(1) * YI(3:5, 1)

! Un-common for free slip condition		
!M(5, 1:2) = -YI(6, 2:3)
!M(5, 3:7) = 0.d0
!Vec(5,1) = Coef(1) * YI(6, 1)

		!¡¸T1 + P continuous

		Vec(6,1) = 0.d0
		M(6, 1:7) = 0.d0
		do i = 1, 6
			M(6, 1:2) = M(6, 1:2) - YI(i, 2:3) * (r_**2.d0*XX(i) + TT(i))
			Vec(6,1) = Vec(6,1) + Coef(1) * YI(i, 1) * (r_**2.d0*XX(i) + TT(i))
		enddo
		do i = 1,5
			M(6, 3:7) = M(6, 3:7) + YI(i, 4:8) * (r_**2.d0*XX_F(i))
		enddo

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

		XX = 0.d0
		XX(1) = (rho_+P_)/2.d0
		XX(4) = - (rho_+P_)/2.d0* dNu/r_
		
		M(7, 1:7) = 0.d0
		Vec(7,1) = 0.d0
		do i = 1, 5
			M(7,3:7) = M(7,3:7) + YS(i, 4:8) * (r_**2 * XX(i))
		enddo

		call LIN_SOL_GEN(M, Vec, Sol)					!IMSL linear equations solver
		coef(2:8) = Sol(1:7,1)

		YS_Comb = 0.d0
		do i = 4,8
			YS_Comb(1:5) = YS_Comb(1:5) + coef(i)*YS(1:5, i)
		enddo

		dHS = YS_Comb(2)
		HS = YS_Comb(1)

		yR_ = r_*dHS/HS

		call lo_k2(yR_, YS_Comb(5)/HS, P_, rho_, m_, 0.d0, r_, k2_)

	endsubroutine lo_bc_EF_W_Solve_No_Slip

endmodule lo_match_P02