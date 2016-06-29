module puls_eqt_set_opt03_ex1
! contains necessary coefficents for the pulsation equaton set option 03
use global_var
contains
	
!*********************************************************
!	Fcn---------------------------------------------------
!	nu and lambda according to Yoshida's formulism (contains factor of 2); ie g00 = exp(2*nu)
	function pes03_eLamb(rho_, m_, r_)
	implicit none
	real(8) :: pes03_eLamb, rho_, m_, r_
		pes03_eLamb = (1.d0 - rel * 2.d0*m_/r_)**(-0.5d0)
	endfunction pes03_eLamb

	function pes03_dLamb(rho_, m_, r_)
	implicit none
	real(8) :: pes03_dLamb, rho_, m_, r_, c1
		c1 = pes03_eLamb(rho_, m_, r_) **2
		pes03_dLamb = c1 / r_ * rel * (4.d0*pi*r_**2*rho_ - m_/r_)
	endfunction pes03_dLamb

	function pes03_dnu(P_, rho_, m_, r_)
	implicit none
	real(8) :: pes03_dnu, P_, rho_, m_, r_, c1
		c1 = pes03_eLamb(rho_, m_, r_) **2
		pes03_dnu = c1 / r_ * rel * ( rel * 4.d0*pi*r_**(2)*P_ + m_/r_)
	endfunction pes03_dnu

	function pes03_ddnu(P_, rho_, m_, r_)
	implicit none
	real(8) :: pes03_ddnu, P_, rho_, m_, r_, c1
	real(8) :: dLamb, dP, dnu
	real(8) :: B(1:3)

		c1 = pes03_eLamb(rho_, m_, r_) **2
		dLamb = pes03_dLamb(rho_, m_, r_)
		dnu = pes03_dnu(P_, rho_, m_, r_)
		dP = -dnu *(rho_ + P_)
		
		B(1) = 2.d0*dLamb*(4.d0*pi*r_*P_ + m_/r_**2)
		B(2) = 4.d0*pi*(P_+r_*dP)
		B(3) = 4.d0*pi*rho_ - 2.d0*m_/r_**3
		pes03_ddnu = c1 * (B(1) +B(2) +B(3))

	endfunction pes03_ddnu

	function pes03_gamma(mode, P_, rho_, ii, i1, i2)
	!	i1, i2 the range of P_ and rho_; must be consistent with ii
	implicit none
	integer :: ii, mode, i1, i2, w, steps
	real(8), dimension(i1:i2) :: P_, rho_
	real(8) :: c1, pes03_gamma
		
		w = ii
		if	(mode == 1) then 
			if (w == i2) w = i2 -1
			if (P_(w+1)*P_(w)*rho_(w+1)*rho_(w) == 0.d0) then
				write(*,*) "err: gamma_1 pulsation equations"
				pause
			endif
			c1 = (1.d0 + rel * P_(w)/rho_(w))
			pes03_gamma = (P_(w+1)-P_(w))/(rho_(w+1)-rho_(w))*(rho_(w)/P_(w))*c1
			if (rho_(w+1)-rho_(w) == 0) pes03_gamma = (P_(w+2)-P_(w))/(rho_(w+2)-rho_(w))*(rho_(w)/P_(w))*c1 ! avoid infinity
		elseif (mode == 2) then
			if (w == i1) w = i1 +1
			if (P_(w)*P_(w-1)*rho_(w)*rho_(w-1) == 0.d0) then
				write(*,*) "err: gamma_1 pulsation equations"
				pause
			endif
			c1 = (1.d0 + rel * P_(w-1)/rho_(w-1))
			if (rho_(w)-rho_(w-1) /= 0.d0) then 
				pes03_gamma = (P_(w)-P_(w-1))/(rho_(w)-rho_(w-1))*(rho_(w-1)/P_(w-1))*c1
			else
			steps = 1
			do while (rho_(w)-rho_(w-steps) == 0.d0)
				steps = steps + 1
				pes03_gamma = (P_(w)-P_(w-steps))/(rho_(w)-rho_(w-steps))*(rho_(w-steps)/P_(w-steps))*c1
				if (steps >= 200) then
					write(*,*) "err: gamma_1 infinity"
					pause
				endif
			enddo
			endif
		endif
	endfunction pes03_gamma

!*********************************************************
!	Coefficients---------------------------------------------------
!	Coefficients in the basis {H1, K, W, X}
	!	E_(1,x) is the coefficients for H0
	!	E_(2,x) is the coefficients for V
	!	i.e. H0 = E_(1,1)*H1 + E_(1,2)*K + E_(1,3)*W + E_(1,4)*X
!	Fluid LD2---------------------------------------------------
	subroutine pes03_E_Fluid(P_, rho_, m_, Nu_, r_, E_)

	implicit none
	real(8) :: P_, rho_, m_, Nu_, r_, l_1, l_2, pi4, dP
	real(8) :: eLamb, eNu
	complex(8) :: Alpha_(1:3), E_(1:2, 1:4)
		
		pi4 = 4.d0 * pi
		l_1 = l_0 * (l_0 + 1.d0)
		l_2 = (l_0 + 2.d0) * (l_0 - 1.d0)
		eLamb = pes03_eLamb(rho_, m_, r_)**2
		eNu = dexp(Nu_)**2
		dP = -(rho_ + P_)*pes03_dnu(P_, rho_, m_, r_)
		Alpha_(1) = (3.d0 * m_ + 0.5d0 * l_2 * r_ + 4.d0 * pi * r_**3 * P_)	! Coefficient of H0 in eqt (1)
		Alpha_(2) = OmeC_sq * (P_ + rho_) * eNu **  (-0.5d0)	! Coefficient of V in eqt (2)
		Alpha_(3) = 0.5d0 * (P_ + rho_) * eNu **  (0.5d0)		! Coefficient of H0 in eqt (2)

		E_(1,1) = -1.d0/Alpha_(1) * ( 0.5d0 * l_1 * (m_ +pi4*r_**3 * P_) - OmeC_sq*r_**3 *(eLamb*eNu)**(-1.d0))
		E_(1,2) = 1.d0/Alpha_(1) * (0.5d0 * l_2 * r_ - OmeC_sq*r_**3 *(eNu)**(-1.d0) - eLamb/r_ * (m_ + pi4*r_**3 * P_) * (3.d0*m_ - r_ + pi4*r_**3 * P_))
		E_(1,3) = (0.d0, 0.d0)
		E_(1,4) = 1.d0/Alpha_(1) * 2.d0*pi4*r_**3 *(eNu)**(-0.5d0)
if (OmeC_sq /= (0.d0, 0.d0)) then
		E_(2,1) = -Alpha_(3)/Alpha_(2) * E_(1,1)
		E_(2,2) = -Alpha_(3)/Alpha_(2) * E_(1,2)
		E_(2,3) = 1.d0/Alpha_(2) * dP/r_ * (eNu/eLamb) ** (0.5d0)
		E_(2,4) = 1.d0/Alpha_(2) - Alpha_(3)/Alpha_(2)*E_(1,4)
else
		write(*,*) "err: pes03_E_Fluid, frequency 0"
		pause
endif
	endsubroutine pes03_E_Fluid

!	Fluid LD2 r=0---------------------------------------------------
	!	2nd order: UVec = Z; very small corrections (following Andersson 2014)
	!	H0 here is referred to H10
	subroutine pes03_bc_Co_Coef(mode, E_)
	use Eigenvalues_Eigenvectors
	implicit none
	real(8) :: l_1, l_2, pi4_3, pi4_5, pi4, eNu0
	real(8) :: P0, rho0, Nu0, gamma0, P2, rho2, Nu2, P4, Nu4
	real(8) :: Q0, Q1
	complex(8) :: T(1:4), U(1:4, 1:4), E_(1:2, 1:4), Vec(1:4), Z(1:4), H0, K0, W0, X0, OmeC_sqb
	integer :: mode
	real(8) :: err_
		pi4 = 4.d0 * pi
		pi4_3 = 4.d0/3.d0*pi
		pi4_5 = 4.d0/5.d0*pi

		l_1 = l_0 * (l_0 + 1.d0)
		l_2 = (l_0 + 2.d0)*(l_0 - 1.d0)

		P0 = P(0) !P_Co(0)
		rho0 = rho(0) !rho_Co(0)
		Nu0 = nu(0) !nu_Co(0)
		eNu0 = dexp(Nu0)**2.d0
		gamma0 = (P(0) + rho(0))/P(0) * (P(1) - P(0)) /(rho(1) - rho(0))  !dlog(P_Co(1)/P_Co(0))/dlog(rho_Co(1)/rho_Co(0))

		OmeC_sqb = OmeC_sq *eNu0**(-1.d0) 

		P2 = -pi4_3*(P0 + rho0)*(3.d0 *P0 + rho0)
		rho2 = P2*(P0 + rho0)/gamma0/P0
		Nu2 = pi4_3*2.d0*(3.d0*P0 + rho0)

		P4 = -0.5d0*pi4_5 * (P0 + rho0)*(5.d0*P2 + rho2) - 0.5d0*pi4_3*(P2 + rho2)*(3.d0*P0 + rho0) - 32.d0/9.d0*pi**2 * rho0 * (P0+rho0)*(3.d0*P0 + rho0)
		Nu4 = pi4_5 * (5.d0*P2 + rho2) + 4.d0*pi4_3**2 * rho0 * (3.d0*P0 + rho0)

		if (mode == 1) K0 = dcmplx(P0 + rho0, 0.d0)
		if (mode == 2) K0 = -dcmplx(P0 + rho0, 0.d0)
		W0 = (1.d0, 0.d0)
		H0 = 1.d0/l_1 * (2.d0 * l_0 * K0 + 16.d0 * pi * (P0 + rho0) * W0)
			T(1) = (P0 + rho0)*eNu0**(0.5d0)
			T(2) = 0.5d0
			T(3) = 4.d0/3.d0* pi *(3.d0*P0 + rho0) - OmeC_sqb/l_0
		X0 = T(1) *(T(2) * K0 + T(3) * W0)
		E_(1,1) = H0
		E_(1,2) = K0
		E_(1,3) = W0
		E_(1,4) = X0

		Q0 = 4.d0/l_2 *(pi4*2.d0 * eNu0**(-0.5d0) * X0 - (pi4_3*2.d0 * rho0 + OmeC_sqb) * K0 - (pi4_3*0.5d0 * l_1 * (3.d0*P0 + rho0) - OmeC_sqb)*H0 )
		Q1 = 2.d0/l_1 * (eNu0**(-0.5d0)/gamma0/P0* X0 + 1.5d0* K0 + pi4_3*(l_0 + 1.d0)*rho0 *W0)

		U(1,1) = 0.d0
		U(1,2) = -(P0 + rho0)/4.d0
		U(1,3) = 0.5d0 * (P2 + (P0 + rho0)*OmeC_sqb*(l_0+3.d0)/l_1 )
		U(1,4) = 0.5d0* eNu0**(-0.5d0)

		U(2,1) = - l_1/4.d0
		U(2,2) = (l_0 + 2.d0)/2.d0
		U(2,3) = pi4 * (P0 + rho0)
		U(2,4) = 0.d0
		
		U(3,1) = (l_0 + 3.d0)/2.d0
		U(3,2) = -1.d0
		U(3,3) = -pi4*2.d0 * (P0 + rho0)*(l_0+3.d0)/l_1
		U(3,4) = 0.d0

		U(4,1) = -l_1 * (P0 + rho0) * eNu0**(0.5d0) /8.d0
		U(4,2) = 0.d0
		U(4,3) = -(P0 + rho0) * eNu0**(0.5d0) * ((l_0 + 2.d0)/4.d0 * Nu2 - pi4/2.d0 * (P0 + rho0) - OmeC_sqb/2.d0)
		U(4,4) = (l_0 + 2.d0)/2.d0

			T(1) = Nu2*eNu0**(-0.5d0)/4.d0 * X0
			T(2) = (rho2 + P2)/4.d0 * K0
			T(3) = (P0 + rho0)/4.d0* Q0 + OmeC_sqb*(P0 + rho0)/2.d0 * Q1
			T(4) = -(P4 - pi4_3*rho0*P2 + OmeC_sqb/2.d0/l_0*(rho2 + P2 - (P0 + rho0)*Nu2 ) )* W0
		Z(1) = T(1) + T(2) + T(3) + T(4)

			T(1) = pi4_3 * (3.d0*P0 + rho0) * K0
			T(2) = 0.5d0 * Q0
			T(3) = -pi4 * (rho2 + P2 + pi4_3*2.d0 * rho0 * (P0 * rho0))*W0
		Z(2) = T(1) + T(2) + T(3)

			T(1) = pi4*((2.d0*l_0 + 3.d0)*rho0/3.d0 - P0)*H0
			T(2) = pi4*2.d0/l_0 * (P2 + rho2)*W0
			T(3) = -pi4*2.d0 * (P0 + rho0)* Q1 + 0.5d0 * Q0
		Z(3) = T(1) + T(2) + T(3)
					
			T(1) = 0.5d0*(rho2 + P2 + 0.5d0*(P0 + rho0)*Nu2)*l_0/(P0 + rho0) * X0
			T(2) = (P0 + rho0)*eNu0**(0.5d0) * (0.5d0*Nu2* K0 + 0.25d0* Q0 + 0.5d0*OmeC_sqb* H0 -0.25d0*l_1*Nu2* Q1)
			T(3) = (P0 + rho0)*eNu0**(0.5d0) * (0.5d0*(l_0 + 1.d0)*Nu4 - pi4/2.d0*(P2 + rho2) - pi4_3*pi4*rho0*(P0 + rho0) + 0.5d0*(Nu4 - pi4_3*rho0*Nu2) + 0.5d0*OmeC_sqb*(Nu2 - pi4_3*2.d0*rho0))* W0
		Z(4) = T(1) + T(2) + T(3)

		call LinSys_Cramer_c(U, Z, Vec(1:4), 4)

		E_(2,1:4) = Vec(1:4)

	endsubroutine pes03_bc_Co_Coef
!	Crust Andersson---------------------------------------------------

!	Coefficients---------------------------------------------------	
!	Coefficients in the basis {H1, K, H0, W, V, T2}
	subroutine pes03_E_Crust(P_, rho_, m_, Nu_, mu_, gamma_1, r_, E_)

	implicit none
	real(8) :: P_, rho_, m_, Nu_, r_, l_1, l_2_hf, pi4, dP, Q, dNu, mu_, gamma_1
	real(8) :: eLamb, eNu
	complex(8) :: A_(1:4), B_(1:2, 1:6), E_(1:3, 1:6)
	integer :: i
		pi4 = 4.d0 * pi
		l_1 = l_0 * (l_0 + 1.d0)
		l_2_hf = 0.5d0 * (l_0 + 2.d0) * (l_0 - 1.d0)
		eLamb = pes03_eLamb(rho_, m_, r_)**2
		eNu = dexp(Nu_)**2
		dNu = pes03_dnu(P_, rho_, m_, r_)*2.d0
		dP = -(rho_ + P_)*pes03_dnu(P_, rho_, m_, r_)
		Q = 0.5d0 * r_**2 * eLamb**(-1.d0) * dNu

		A_(1) = 2.d0*pi4*r_**3*eNu**(-0.5d0)			! Coefficient of X in eqt (2)
		A_(2) = 2.d0*pi4*r_								! Coefficient of T1 in eqt (2)
		A_(3) = 2.d0/3.d0*eNu**(-0.5d0)*mu_ *r_**2		! Coefficient of X in eqt (3)
		A_(4) = -0.25d0*gamma_1*P_						! Coefficient of T1 in eqt (3)

		B_ = 0.d0
		B_(1,1) = (l_2_hf + 1.d0)*Q - OmeC_sq*r_**3 * (eLamb*eNu)**(-1.d0)						! Coefficient of H1 in eqt (2)
		B_(1,2) = -(l_2_hf*r_ - OmeC_sq*r_**3 *eNu**(-1.d0) - eLamb/r_*Q *(2.d0*m_ + Q - r_))	! Coefficient of K in eqt (2)
		B_(1,3) = 2.d0*m_ + Q + l_2_hf*r_														! Coefficient of H0 in eqt (2)
		B_(1,6) = 4.d0*pi4*r_ * eLamb**(-1.d0)													! Coefficient of T2 in eqt (2)

		B_(2,2) = - mu_ * gamma_1 * P_ * r_**2													! Coefficient of K in eqt (3)
		B_(2,4) = 2.d0 * mu_ * gamma_1 * P_ * eLamb**(-0.5d0)									! Coefficient of W in eqt (3)
		B_(2,5) = l_1 * mu_ * gamma_1 * P_														! Coefficient of V in eqt (3)

		E_ = 0.d0																				! H2
		E_(1,3) = 1.d0
		E_(1,5) = 16.d0*pi4*mu_
		
		E_(2,1:6) = (B_(1,1:6)/A_(2) - B_(2,1:6)/A_(4))/(A_(1)/A_(2) - A_(3)/A_(4))				! X		! Crammer's Rule
		E_(3,1:6) = (B_(1,1:6)/A_(1) - B_(2,1:6)/A_(3))/(A_(2)/A_(1) - A_(4)/A_(3))				! T1	! Crammer's Rule

	endsubroutine pes03_E_Crust

!	Coefficients---------------------------------------------------
!	Coefficients in the basis {H1, K, H0, W, V}
	!	Coefficients for BC in crust surface, 2 component model
	!	At R0, X = 0     ! <---- check equivalent to X - T1 = 0???      !<---- Ans: NO, see Finn 1990, true BC should be X - T1 = 0
	!	Using the algebraic relations from Andersson 2014,
	!	we obtain an equation of [A dot {H1, K, H0, W, V} = 0]
	!
	subroutine pes03_bc_R_Coe(P_, rho_, m_, nu_, mu_, gamma_, r_, A)
	implicit none
	complex(8) :: A(3:7), H1, K, H0, W, V, E_(1:3, 1:6)
	real(8) :: P_, rho_, m_, nu_, r_, mu_, gamma_, l_1, l_2_hf, eLamb, eNu, dNu, Q

		H1 = 1.d0/R0**2			!	From BC_R
		K = 1.d0/R0**2
		H0 = 1.d0/R0**2
		W = 1.d0
		V = 1.d0

		l_1 = l_0 * (l_0 + 1.d0)
		l_2_hf = 0.5d0 * (l_0 + 2.d0) * (l_0 - 1.d0)
		eLamb = pes03_eLamb(rho_, m_, r_)**2
		eNu = dexp(Nu_)**2
		dNu = pes03_dnu(P_, rho_, m_, r_)*2.d0
		Q = 0.5d0 * r_**2 * eLamb**(-1.d0) * dNu

		!¡¸ A dot {H1, K, H0, W, V}
!		A(3) = 1.d0/4.d0 * ((l_2_hf + 1.d0)*Q - OmeC_sq*r_**3 * (eLamb*eNu)**(-1.d0)) * H1
!		A(4) = - ( 8.d0 * pi * r_**3 * mu_ + 1.d0/4.d0 * (l_2_hf*r_ - OmeC_sq*r_**3 *eNu**(-1.d0) - eLamb/r_*Q *(2.d0*m_ + Q - r_)) ) * K
!		A(5) = 1.d0/4.d0 * (2.d0*m_ + Q + l_2_hf*r_) * H0
!		A(6) = 16.d0 * pi * r_ * mu_ * eLamb**(-0.5d0) * W
!		A(7) = 8.d0 * pi * r_ * l_1 * mu_ * V
		
		call pes03_E_Crust(P_, rho_, m_, nu_, mu_, gamma_, r_, E_)
		A(3) = (-r_**2 * eNu**(-0.5d0)*E_(2,1) - E_(3,1) )*H1				! (r^2 P_eu + X = 0); Finn: (P_Lag - T1_And * (r^(l_0-2)))
		A(4) = (-r_**2 * eNu**(-0.5d0)*E_(2,2) - E_(3,2) )*K
		A(5) = (-r_**2 * eNu**(-0.5d0)*E_(2,3) - E_(3,3) )*H0
		A(6) = (-r_**2 * eNu**(-0.5d0)*E_(2,4) - E_(3,4) )*W
		A(7) = (-r_**2 * eNu**(-0.5d0)*E_(2,5) - E_(3,5) )*V
	endsubroutine pes03_bc_R_Coe

!	Outside the star---------------------------------------------------
	function pes03_V_z(m_, r_)
	implicit none
	real(8) :: m_, r_, n, T(1:5)
	real(8) :: pes03_V_z

	n = (l_0 -1.d0)*(l_0+ 2.d0)/2.d0
	T(1) = (1.d0-2.d0*m_/r_)/(r_**3 * (n*r_ + 3.d0 * m_)**2)
	T(2) = 2.d0 *n**2 *(n+1.d0)*r_**3
	T(3) = 6.d0 * n**2 *m_*r_**2
	T(4) = 18.d0 * n *m_**2 *r_
	T(5) = 18.d0 *m_ **3
	
	pes03_V_z = T(1)*(T(2) + T(3) + T(4) + T(5) )
	endfunction pes03_V_z

endmodule puls_eqt_set_opt03_ex1