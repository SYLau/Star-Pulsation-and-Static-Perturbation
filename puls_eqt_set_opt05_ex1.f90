module puls_eqt_set_opt05_ex1
! contains necessary coefficents for the pulsation equaton set option 05
use global_var
contains
	
!*********************************************************
!	Fcn---------------------------------------------------
!	nu and lambda according to Yoshida's formulism (contains factor of 2); ie g00 = exp(2*nu)
	function pes05_eLamb(rho_, m_, r_)
	implicit none
	real(8) :: pes05_eLamb, rho_, m_, r_
		pes05_eLamb = (1.d0 - rel * 2.d0*m_/r_)**(-0.5d0)
	endfunction pes05_eLamb

	function pes05_dLamb(rho_, m_, r_)
	implicit none
	real(8) :: pes05_dLamb, rho_, m_, r_, c1
		c1 = pes05_eLamb(rho_, m_, r_) **2
		pes05_dLamb = c1 / r_ * rel * (4.d0*pi*r_**2*rho_ - m_/r_)
	endfunction pes05_dLamb

	function pes05_dnu(P_, rho_, m_, r_)
	implicit none
	real(8) :: pes05_dnu, P_, rho_, m_, r_, c1
		c1 = pes05_eLamb(rho_, m_, r_) **2
		pes05_dnu = c1 / r_ * rel * ( rel * 4.d0*pi*r_**(2)*P_ + m_/r_)
	endfunction pes05_dnu

	function pes05_ddnu(P_, rho_, m_, r_)
	implicit none
	real(8) :: pes05_ddnu, P_, rho_, m_, r_, c1
	real(8) :: dLamb, dP, dnu
	real(8) :: B(1:3)

		c1 = pes05_eLamb(rho_, m_, r_) **2
		dLamb = pes05_dLamb(rho_, m_, r_)
		dnu = pes05_dnu(P_, rho_, m_, r_)
		dP = -dnu *(rho_ + P_)
		
		B(1) = 2.d0*dLamb*(4.d0*pi*r_*P_ + m_/r_**2)
		B(2) = 4.d0*pi*(P_+r_*dP)
		B(3) = 4.d0*pi*rho_ - 2.d0*m_/r_**3
		pes05_ddnu = c1 * (B(1) +B(2) +B(3))

	endfunction pes05_ddnu

	function pes05_gamma(mode, P_, rho_, ii, i1, i2)
	!	i1, i2 the range of P_ and rho_; must be consistent with ii
	implicit none
	integer :: ii, mode, i1, i2, w, steps
	real(8), dimension(i1:i2) :: P_, rho_
	real(8) :: c1, pes05_gamma
		
		w = ii
		if	(mode == 1) then 
			if (w == i2) w = i2 -1
			if (P_(w+1)*P_(w)*rho_(w+1)*rho_(w) == 0.d0) then
				write(*,*) "err: gamma_1 pulsation equations"
				pause
			endif
			c1 = (1.d0 + rel * P_(w)/rho_(w))
			pes05_gamma = (P_(w+1)-P_(w))/(rho_(w+1)-rho_(w))*(rho_(w)/P_(w))*c1
			if (rho_(w+1)-rho_(w) == 0) pes05_gamma = (P_(w+2)-P_(w))/(rho_(w+2)-rho_(w))*(rho_(w)/P_(w))*c1 ! avoid infinity
		elseif (mode == 2) then
			if (w == i1) w = i1 +1
			if (P_(w)*P_(w-1)*rho_(w)*rho_(w-1) == 0.d0) then
				write(*,*) "err: gamma_1 pulsation equations"
				pause
			endif
			c1 = (1.d0 + rel * P_(w-1)/rho_(w-1))
			if (rho_(w)-rho_(w-1) /= 0.d0) then 
				pes05_gamma = (P_(w)-P_(w-1))/(rho_(w)-rho_(w-1))*(rho_(w-1)/P_(w-1))*c1
			else
			steps = 1
			do while (rho_(w)-rho_(w-steps) == 0.d0)
				steps = steps + 1
				pes05_gamma = (P_(w)-P_(w-steps))/(rho_(w)-rho_(w-steps))*(rho_(w-steps)/P_(w-steps))*c1
				if (steps >= 200) then
					write(*,*) "err: gamma_1 infinity"
					pause
				endif
			enddo
			endif
		endif
	endfunction pes05_gamma


!	Crust LD2 r=0---------------------------------------------------
	!	0th order
	!	H0 here is referred to H00
	!	Using Andersson (2014) Definition of T2; hence mu -> mu/2
	subroutine pes05_bc_Co_Coef(mode,H10, K0, H0, W0, V0, W2, V2)
	implicit none
	real(8) :: l_1, pi4_3, pi4_5, pi4, eNu0
	real(8) :: P0, rho0, Nu0, gamma0, mu0
	complex(8) :: T(1:3), H10, H0, K0, W0, V0, W2, V2
	integer :: mode
	real(8) :: err_
		l_1 = l_0 * (l_0 + 1.d0)

		P0 = P_Co(0)
		rho0 = rho_Co(0)
		Nu0 = nu_Co(0)
		eNu0 = dexp(Nu0)**2.d0
		gamma0 = (P_Co(0) + rho_Co(0))/P_Co(0) * (P_Co(1) - P_Co(0)) /(rho_Co(1) - rho_Co(0))  !dlog(P_Co(1)/P_Co(0))/dlog(rho_Co(1)/rho_Co(0))
		mu0 = mu_Co(0)/2.d0

		K0 = (0.d0, 0.d0)
		V2 = (0.d0, 0.d0)
		V0 = (0.d0, 0.d0)

			if (mode == 1) then
				K0 = dcmplx(P0 + rho0, 0.d0)
			elseif (mode == 2) then
				K0 = -dcmplx(P0 + rho0, 0.d0)
				V2 = dcmplx(1.d0/R0**2, 0.d0)
			elseif (mode == 3) then 
				V0 = dcmplx(1.d0,0.d0)
			endif
			W0 = -l_0 * V0
			H0 = K0 - 64.d0*pi*mu0*V0
			H10 = 1.d0/l_1 * (2.d0 * l_0 * K0 + 16.d0 * pi * (P0 + rho0) * W0)
				
				T(1) = 0.5d0 * ( (3.d0 * gamma0 + 1.d0)*P0 + rho0 )
				T(2) = - (gamma0 * P0 * l_0 + 2.d0 * mu0/3.d0 * (l_0 - 6.d0)) * (l_0 + 1.d0)  ! Finn (1990) missing (l_0 + 1) ???
				T(3) = (gamma0 * P0 * (l_0 + 3.d0) + 2.d0 * mu0/3.d0 * (l_0 + 9.d0))			! Added 2* mu for Andersson 2014
			W2 = (T(1)*H0 + T(2)*V2)/T(3)

	endsubroutine pes05_bc_Co_Coef

!	Crust Andersson---------------------------------------------------

!	Coefficients---------------------------------------------------	
!	Coefficients in the basis {H1, K, H0, W, V, T2}
	subroutine pes05_E_Crust(P_, rho_, m_, Nu_, mu_, gamma_1, r_, E_)

	implicit none
	real(8) :: P_, rho_, m_, Nu_, r_, l_1, l_2_hf, pi4, dP, Q, dNu, mu_, gamma_1
	real(8) :: eLamb, eNu
	complex(8) :: A_(1:4), B_(1:2, 1:6), E_(1:3, 1:6)
	integer :: i
		pi4 = 4.d0 * pi
		l_1 = l_0 * (l_0 + 1.d0)
		l_2_hf = 0.5d0 * (l_0 + 2.d0) * (l_0 - 1.d0)
		eLamb = pes05_eLamb(rho_, m_, r_)**2
		eNu = dexp(Nu_)**2
		dNu = pes05_dnu(P_, rho_, m_, r_)*2.d0
		dP = -(rho_ + P_)*pes05_dnu(P_, rho_, m_, r_)
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

	endsubroutine pes05_E_Crust

!	Coefficients---------------------------------------------------
!	Coefficients in the basis {H1, K, H0, W, V}
	!	Coefficients for BC in crust surface, 2 component model
	!	At R0, X = 0
	!	Using the algebraic relations from Andersson 2014,
	!	we obtain an equation of [A dot {H1, K, H0, W, V} = 0]
	!
	subroutine pes05_bc_R_Coe(P_, rho_, m_, nu_, mu_, gamma_, r_, A)
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
		eLamb = pes05_eLamb(rho_, m_, r_)**2
		eNu = dexp(Nu_)**2
		dNu = pes05_dnu(P_, rho_, m_, r_)*2.d0
		Q = 0.5d0 * r_**2 * eLamb**(-1.d0) * dNu

		!¡¸ A dot {H1, K, H0, W, V}
!		A(3) = 1.d0/4.d0 * ((l_2_hf + 1.d0)*Q - OmeC_sq*r_**3 * (eLamb*eNu)**(-1.d0)) * H1
!		A(4) = - ( 8.d0 * pi * r_**3 * mu_ + 1.d0/4.d0 * (l_2_hf*r_ - OmeC_sq*r_**3 *eNu**(-1.d0) - eLamb/r_*Q *(2.d0*m_ + Q - r_)) ) * K
!		A(5) = 1.d0/4.d0 * (2.d0*m_ + Q + l_2_hf*r_) * H0
!		A(6) = 16.d0 * pi * r_ * mu_ * eLamb**(-0.5d0) * W
!		A(7) = 8.d0 * pi * r_ * l_1 * mu_ * V
		
		call pes05_E_Crust(P_, rho_, m_, nu_, mu_, gamma_, r_, E_)
		A(3) = (-r_**2 * eNu**(-0.5d0)*  E_(2,1) - E_(3,1) )*H1
		A(4) = (-r_**2 * eNu**(-0.5d0)*E_(2,2) - E_(3,2) )*K
		A(5) = (-r_**2 * eNu**(-0.5d0)*E_(2,3) - E_(3,3) )*H0
		A(6) = (-r_**2 * eNu**(-0.5d0)*E_(2,4) - E_(3,4) )*W
		A(7) = (-r_**2 * eNu**(-0.5d0)*E_(2,5) - E_(3,5) )*V
	endsubroutine pes05_bc_R_Coe

!	Outside the star---------------------------------------------------
	function pes05_V_z(m_, r_)
	implicit none
	real(8) :: m_, r_, n, T(1:5)
	real(8) :: pes05_V_z

	n = (l_0 -1.d0)*(l_0+ 2.d0)/2.d0
	T(1) = (1.d0-2.d0*m_/r_)/(r_**3 * (n*r_ + 3.d0 * m_)**2)
	T(2) = 2.d0 *n**2 *(n+1.d0)*r_**3
	T(3) = 6.d0 * n**2 *m_*r_**2
	T(4) = 18.d0 * n *m_**2 *r_
	T(5) = 18.d0 *m_ **3
	
	pes05_V_z = T(1)*(T(2) + T(3) + T(4) + T(5) )
	endfunction pes05_V_z

endmodule puls_eqt_set_opt05_ex1