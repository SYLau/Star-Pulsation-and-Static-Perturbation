module puls_eqt_set_opt05
! nonradial oscillations full relativistic formulism (Andersson 2014)
! Co: xc = {H1, K, W, X} according to LD formalism
! Cr: xc = {H1, K, H0, W, V, T2} according to Andersson 2014
! Input nu defined by g00 = exp(2*nu)
use global_var
use puls_eqt_set_opt05_ex1
contains

!*********************************************************
!	BC---------------------------------------------------
	subroutine pes05_bc_Coi(mode, r_, z)
	implicit none
	real(8) :: r_
	complex(8) :: z(1:6), mu0, H10, K0, H0, W0, V0, W2, V2
	integer :: mode

		call pes05_bc_Co_Coef(mode,H10, K0, H0, W0, V0, W2, V2)		! zeroth order and 2nd order expansion coefficients; ref: LD2
		mu0 = mu_Co(0)/2.d0

			z(1) = H10
			z(2) = K0
			z(3) = H0
			z(4) = W0 + W2 * r_**2
			z(5) = V0 + V2 * r_**2
			z(6) = - 2.d0 * mu0 * (W0 + W2 * r_**2 - (l_0 -2.d0) * (V0 + V2 * r_**2) -2.d0 * (V2 * r_**2))

	endsubroutine pes05_bc_Coi

	subroutine pes05_bc_Cof(z, cY)
	implicit none
	complex(8) :: z(1:6), cY(1:6)
	integer :: i
		cY(1:6) = z(1:6)
	endsubroutine pes05_bc_Cof


	subroutine pes05_bc_Cri(z, cY)
	implicit none
	complex(8) :: z(1:6), cY(1:6)
		cY(1:6) = z(1:6)
	endsubroutine pes05_bc_Cri

	subroutine pes05_bc_Crf(mode, z)
	!	3 components
	implicit none
	complex(8) :: z(1:6)
	integer :: mode

		z = 0.d0
		if (mode == 1) z(1) = dcmplx(1.d0/R0**2, 0.d0)		!!	setting z1-z3 to be 1.d0 gives significant error in z4 in Crammer's rule during matching solution
		if (mode == 2) z(2) = dcmplx(1.d0/R0**2, 0.d0)
		if (mode == 3) z(3) = dcmplx(1.d0/R0**2, 0.d0)
		if (mode == 4) z(4) = (1.d0, 0.d0)
		if (mode == 5) z(5) = (1.d0, 0.d0)
	endsubroutine pes05_bc_Crf
	! Outside Neutron Star
	subroutine pes05_bc_Oui(H0, K, y)
	implicit none
	real(8) :: n, m_, r_, g_, l_, h_, k_, det, tcf
	complex(8) :: a_, b_, Y_R(1:4), E(1:2, 1:4)
	complex(8) :: coef(1:5), y(1:2), H0, K
	integer :: i
		m_ = M0
		r_ = R0
		tcf = (1.d0 - 2.d0*m_/r_)**(-1.d0)

		n = (l_0 -1.d0)*(l_0+ 2.d0)/2.d0
		a_ = -(n*r_ + 3.d0*m_)/(OmeC_sq*r_**2 - (n+1.d0)*m_/r_)
		b_ = (n*r_*(r_-2.d0*m_) - OmeC_sq * r_**4 + m_*(r_-3.d0*m_))/(r_-2.d0*m_)/(OmeC_sq*r_**2 - (n+1.d0)*m_/r_) ! this mistake caused ~15% error in f-mode
		g_ = (n*(n+1.d0)*r_**2 + 3.d0*n*m_*r_ + 6.d0*m_**2)/r_**2/(n*r_ + 3.d0 * m_)
		l_ = 1.d0
		h_ = (-n*r_**2 + 3.d0*n*m_*r_ + 3.d0*m_**2)/(r_-2.d0*m_)/(n*r_ + 3.d0 * m_)
		k_ = -r_**2/(r_ - 2.d0*m_)
		det = g_ * k_ - h_ * l_

		y(1) = (-a_*l_ * H0 + (k_ - b_ * l_) * K )/det
		y(2) = (g_ *a_ * H0 + (-h_ + b_ * g_) * K )/det*tcf

	endsubroutine pes05_bc_Oui

	subroutine pes05_bc_Ouf_LD(r_, y, be, ga)
	implicit none
	real(8) :: n, m_, r_, r_t, tcf
	complex(8) :: alpha(0:2), calpha(0:2), im, OmeC, det, be, ga
	complex(8) :: z_asym(1:2, 1:2), y(1:2)
		m_ = M0
		r_t = r_ + 2.d0*m_* dlog(r_/(2.d0*m_) - 1.d0)
		! wrong r_t = r_ + dlog(r_/(2.d0*m_ - 1.d0))
		OmeC = OmeC_sq ** 0.5d0
		im = (0.d0, 1.d0)
		n = (l_0 -1.d0)*(l_0+ 2.d0)/2.d0
		tcf = (1.d0 - 2.d0*m_/r_)**(-1.d0)

		alpha(0) = 1.d0
		alpha(1) = -im*(n + 1.d0)/OmeC *alpha(0)
		alpha(2) = -0.5d0/OmeC_sq*(n*(n+1.d0) - 1.5d0 * im * m_ * OmeC * (1.d0 + 2.d0/n)) * alpha(0)
		calpha = dconjg(alpha)
		
		z_asym(1,1) = cdexp(-im*OmeC * r_t )*(alpha(0) + alpha(1)*r_**(-1.d0) + alpha(2)*r_**(-2.d0))
		z_asym(1,2) = cdexp(im*OmeC * r_t )*(calpha(0) + calpha(1)*r_**(-1.d0) + calpha(2)*r_**(-2.d0))
		z_asym(2,1) = -im*OmeC * z_asym(1,1) + cdexp(-im*OmeC * r_t)*(-alpha(1)*r_**(-2.d0) - 2.d0 * alpha(2)*r_**(-3.d0)) /tcf
		z_asym(2,2) = im*OmeC * z_asym(1,2) + cdexp(im*OmeC * r_t)*(-calpha(1)*r_**(-2.d0) - 2.d0 * calpha(2)*r_**(-3.d0)) /tcf

		det = z_asym(1,1)*z_asym(2,2) - z_asym(1,2)*z_asym(2,1)

		be = (y(1)*z_asym(2,2) - z_asym(1,2)*y(2)/tcf )/det
		ga = (z_asym(1,1)*y(2)/tcf - y(1)*z_asym(2,1) )/det
!write(*,*) "im*OmeC * r_t", im*OmeC * r_t 
!write(*,*) "cdexp(-im*OmeC * r_t )", cdexp(-im*OmeC * r_t )
!write(*,*) "cdexp(im*OmeC * r_t )", cdexp(im*OmeC * r_t )
!write(*,*) "z_asym(1,1)", z_asym(1,1)
!write(*,*) "z_asym(1,2)", z_asym(1,2)
!write(*,*) "z_asym(2,1)", z_asym(2,1)
!write(*,*) "z_asym(2,2)", z_asym(2,2)
!write(*,*) "dreal(cdexp(-im*OmeC * r_t ))", dreal(cdexp(-im*OmeC * r_t ))
!write(*,*) "dreal(cdexp(im*OmeC * r_t ))", dreal(cdexp(im*OmeC * r_t ))
!write(*,*) "T1", dreal(cdexp(-im*OmeC * r_t ))*alpha(0)
!write(*,*) "T2", dreal(cdexp(-im*OmeC * r_t ))*alpha(1)*r_**(-1.d0)
!write(*,*) "T3", dreal(cdexp(-im*OmeC * r_t ))*alpha(2)*r_**(-2.d0)
!pause
	endsubroutine pes05_bc_Ouf_LD

	subroutine pes05_bc_Ouf_Scatter(r_, y, be, ga)
	implicit none
	real(8) :: n, m_, r_, r_t, tcf
	complex(8) :: T(1:2), dT(1:2), calpha(0:2), im, OmeC, det, be, ga
	complex(8) :: z_asym(1:2, 1:2), y(1:2)
		m_ = M0
		r_t = r_ + 2.d0*m_* dlog(r_/(2.d0*m_) - 1.d0)
		! wrong r_t = r_ + dlog(r_/(2.d0*m_ - 1.d0))
		OmeC = OmeC_sq ** 0.5d0
		im = (0.d0, 1.d0)
		n = (l_0 -1.d0)*(l_0+ 2.d0)/2.d0
		tcf = (1.d0 - 2.d0*m_/r_)**(-1.d0)

		T(1) = 1.d0 - n*(n+1.d0)/2.d0/OmeC_sq/r_**2		!	Coefficients on "Alpha" - Notation refer to Ferrari 1991
		T(2) = -(n+1.d0)/OmeC/r_ + 3.d0 * m_*(1.d0 + 2.d0/n)/2.d0/OmeC/r_**2		!	Coefficients on "Beta" - Notation refer to Ferrari 1991
		dT(1) = n*(n+1.d0)/OmeC_sq/r_**3
		dT(2) = (n+1.d0)/OmeC/r_**2 - 3.d0 * m_*(1.d0 + 2.d0/n)/OmeC/r_**3

		z_asym(1,1) = cdcos(OmeC*r_t)*T(1) + cdsin(OmeC*r_t)*T(2)
		z_asym(1,2) = cdcos(OmeC*r_t)*T(2) + cdsin(OmeC*r_t)*(-T(1))
		z_asym(2,1) = cdcos(OmeC*r_t)*(dT(1) + OmeC*T(2)/tcf) + cdsin(OmeC*r_t)*(dT(2) - OmeC*T(1)/tcf)
		z_asym(2,2) = cdcos(OmeC*r_t)*(dT(2) - OmeC*T(1)/tcf) + cdsin(OmeC*r_t)*(-dT(1) - OmeC*T(2)/tcf)

		det = z_asym(1,1)*z_asym(2,2) - z_asym(1,2)*z_asym(2,1)

		be = (y(1)*z_asym(2,2) - z_asym(1,2)*y(2)/tcf )/det
		ga = (z_asym(1,1)*y(2)/tcf - y(1)*z_asym(2,1) )/det

	endsubroutine pes05_bc_Ouf_Scatter

	subroutine pes05_bc_Ouf_phase(r_, y, be, ga)
	implicit none
	real(8) :: n, m_, r_, r_t, tcf
	complex(8) :: T(1:2), dT(1:2), calpha(0:2), im, OmeC, det, alp1, alp2, be, ga
	complex(8) :: z_asym(1:2, 1:2), y(1:2)
		m_ = M0
		r_t = r_ + 2.d0*m_* dlog(r_/(2.d0*m_) - 1.d0)
		! wrong r_t = r_ + dlog(r_/(2.d0*m_ - 1.d0))
		OmeC = OmeC_sq ** 0.5d0
		im = (0.d0, 1.d0)
		n = (l_0 -1.d0)*(l_0+ 2.d0)/2.d0
		tcf = (1.d0 - 2.d0*m_/r_)**(-1.d0)

		T(1) = 1.d0 - n*(n+1.d0)/2.d0/OmeC_sq/r_**2
		T(2) = -(n+1.d0)/OmeC/r_ + 3.d0 * m_*(1.d0 + 2.d0/n)/2.d0/OmeC/r_**2
		dT(1) = n*(n+1.d0)/OmeC_sq/r_**3
		dT(2) = (n+1.d0)/OmeC/r_**2 - 3.d0 * m_*(1.d0 + 2.d0/n)/OmeC/r_**3

		z_asym(1,1) = cdcos(OmeC*r_t)*T(1) + cdsin(OmeC*r_t)*T(2)
		z_asym(1,2) = cdcos(OmeC*r_t)*T(2) + cdsin(OmeC*r_t)*(-T(1))
		z_asym(2,1) = cdcos(OmeC*r_t)*(dT(1) + OmeC*T(2)/tcf) + cdsin(OmeC*r_t)*(dT(2) - OmeC*T(1)/tcf)
		z_asym(2,2) = cdcos(OmeC*r_t)*(dT(2) - OmeC*T(1)/tcf) + cdsin(OmeC*r_t)*(-dT(1) - OmeC*T(2)/tcf)

!z_asym(1,1) = cdcos(OmeC*r_t)
!z_asym(1,2) = cdsin(OmeC*r_t)
!z_asym(2,1) = cdsin(OmeC*r_t)*(- OmeC/tcf)
!z_asym(2,2) = cdcos(OmeC*r_t)*(- OmeC/tcf)

		det = z_asym(1,1)*z_asym(2,2) - z_asym(1,2)*z_asym(2,1)

		alp1 = (y(1)*z_asym(2,2) - z_asym(1,2)*y(2)/tcf )/det
		alp2 = (z_asym(1,1)*y(2)/tcf - y(1)*z_asym(2,1) )/det

		be = alp1* T(1) + alp2* T(2)
		ga = alp2* T(1) - alp1* T(2)

!be = alp1
!ga = alp2

	endsubroutine pes05_bc_Ouf_phase

function pes05_series(r_)
implicit none
complex(8) :: pes05_series
real(8) :: n, m_, r_, r_t
complex(8) :: alpha(0:2), calpha(0:2), im, OmeC
	m_ = M0
	r_t = r_ + 2.d0*m_* dlog(r_/(2.d0*m_) - 1.d0)

	OmeC = OmeC_sq ** 0.5d0
	im = (0.d0, 1.d0)
	n = (l_0 -1.d0)*(l_0+ 2.d0)/2.d0

	alpha(0) = 1.d0
	alpha(1) = -im*(n + 1.d0)/OmeC *alpha(0)
	alpha(2) = -0.5d0/OmeC_sq*(n*(n+1.d0) - 1.5d0 * im * m_ * OmeC * (1.d0 + 2.d0/n)) * alpha(0)
	calpha = dconjg(alpha)
		
	pes05_series = cdexp(-im*OmeC * r_t )*(alpha(0) + alpha(1)*r_**(-1.d0) + alpha(2)*r_**(-2.d0))

endfunction pes05_series
!*********************************************************
!	Eqts---------------------------------------------------

!	Core ¡¸///////////////////////////////////////
	subroutine pes05_vx_Co(n, t, x, fcn)
	implicit none
	integer :: n
	real(8), dimension(1:n) :: x, fcn
	real(8) :: t
	complex(8), dimension(1:n/2,1:n/2) :: B
	complex(8), dimension(1:3,1:n/2) :: E
	complex(8), dimension(1:n/2) :: xc, fcnc, H2_C, X_C, T1_C
	real(8) :: eNu, eLamb, eNu_hf, eLamb_hf, dLamb, dNu, gamma_1, P_, rho_, m_, nu_, mu_, r_
	real(8) :: pi4, l_1, l_2_hf, V1, ddNu, dP
	integer :: i, j, ii

	if (mod(n,2) /= 0) then
		write(*,*) "err: puls eqt set n not even"
		pause
	endif

	xc(1:n/2) = dcmplx(x(1 :n/2),x(n/2+1 :n))

	ii = rk4_HS(XCo_i, dx_Co,t)
	P_ = P_Co(ii)
	rho_ = rho_Co(ii)
	m_ = m_Co(ii)
	nu_ = nu_Co(ii)
	!mu_ = mu_Co(ii) Andersson's mu seems to be off by 1/2, st it cannot be reduced to RCA case
	mu_ = mu_Co(ii)/2.d0
	r_ = r_Co(ii)
	gamma_1 = pes05_gamma(1, P_Co, rho_Co, ii, 0, pg_N1*2)
	call pes05_E_Crust(P_, rho_, m_, nu_, mu_, gamma_1, r_, E)

	pi4 = 4.d0*pi
	l_1 = l_0 * (l_0 + 1.d0)
	l_2_hf = (l_0 + 2.d0)*(l_0 - 1.d0)/2.d0
	eLamb = pes05_eLamb(rho_, m_, r_) ** 2
	eNu = dexp(2.d0*nu_)
	eLamb_hf = eLamb**(0.5d0)
	eNu_hf = eNu**(0.5d0)
	dLamb = pes05_dLamb(rho_, m_, r_) * 2.d0
	dNu = pes05_dnu(P_, rho_, m_, r_) * 2.d0
	ddNu = pes05_ddnu(P_, rho_, m_, r_) * 2.d0
	dP = -dNu/2.d0 * (P_ + rho_)
	V1 = (1.d0 + rho_/P_)*r_ *dNu/2.d0
	H2_C(1:6) = E(1,1:6)
	X_C(1:6) = E(2,1:6)
	T1_C(1:6) = E(3,1:6)

	B = 0.d0

	B(1,1) = 0.5d0*(dLamb - dNu)*r_ - (l_0 + 1.d0)
	B(1,2) = eLamb
	B(1,5) = -eLamb*4.d0*pi4*(P_ + rho_)
	B(1, 1:6) = B(1, 1:6) + eLamb*H2_C(1:6)	! Adding the algebraic solution of H2 & X & T1

	B(2,1) = l_2_hf + 1.d0
	B(2,2) = 0.5d0* dNu*r_ - (l_0 + 1.d0)
	B(2,4) = -2.d0*pi4*(P_ + rho_)*eLamb_hf
	B(2, 1:6) = B(2, 1:6) + H2_C(1:6)

	B(3,1) = -r_**2 * eNu**(-1.d0) * OmeC_sq
	B(3,2) = l_0
	B(3,3) = -(0.5d0* dNu*r_ + (l_0 - 1.d0))
	B(3,6) = -4.d0*pi4
	B(3, 1:6) = B(3, 1:6) + B(2, 1:6) - (0.5d0* dNu*r_ + 1.d0)*H2_C(1:6)

	B(4,2) = r_**2 * eLamb_hf
	B(4,4) = -(l_0 + 1.d0)
	B(4,5) = -eLamb_hf * l_1
	B(4, 1:6) = B(4, 1:6) + 0.5d0 * r_**2 * eLamb_hf * H2_C(1:6) +  r_**2 * eLamb_hf * eNu_hf**(-1.d0)/ gamma_1 / P_ * X_C(1:6)

	B(5,4) = eLamb_hf
	B(5,5) = (2.d0 - l_0)
	B(5,6) = 0.5d0/mu_

	B(6,3) = -0.5d0* r_**2 *eLamb * (P_ + rho_)
	B(6,4) = eLamb_hf * dP * r_
	B(6,5) = 4.d0 * l_2_hf * eLamb * mu_ - eLamb/eNu * r_**2 * OmeC_sq * (P_ + rho_)
	B(6,6) = 0.5d0 * (dLamb - dNu) * r_ - (l_0 + 1.d0)
	B(6, 1:6) = B(6, 1:6) + r_**2 * eLamb/eNu_hf * (X_C(1:6) - 0.5d0/r_**2 * eNu_hf * T1_C(1:6))

	do i = 1,n/2
		fcnc(i) = 0.d0
		do j = 1,n/2
			fcnc(i) = fcnc(i) + B(i,j)*xc(j)/(1.d0 +V1)
		enddo
	enddo
	fcn(1:n/2) = dreal(fcnc(1:n/2))
	fcn(n/2+1 :n) = dimag(fcnc(1:n/2))

	endsubroutine pes05_vx_Co

!	Crust ¡¸///////////////////////////////////////
	subroutine pes05_vx_Cr(n, t, x, fcn)
	implicit none
	integer :: n
	real(8), dimension(1:n) :: x, fcn
	real(8) :: t
	complex(8), dimension(1:n/2,1:n/2) :: B
	complex(8), dimension(1:3,1:n/2) :: E
	complex(8), dimension(1:n/2) :: xc, fcnc, H2_C, X_C, T1_C
	real(8) :: eNu, eLamb, eNu_hf, eLamb_hf, dLamb, dNu, gamma_1, P_, rho_, m_, nu_, mu_, r_
	real(8) :: pi4, l_1, l_2_hf, V1, ddNu, dP
	integer :: i, j, ii

	if (mod(n,2) /= 0) then
		write(*,*) "err: puls eqt set n not even"
		pause
	endif

	xc(1:n/2) = dcmplx(x(1 :n/2),x(n/2+1 :n))

	ii = rk4_HS(XCr_i, dx_Cr,t)
	P_ = P_Cr(ii)
	rho_ = rho_Cr(ii)
	m_ = m_Cr(ii)
	nu_ = nu_Cr(ii)
	!mu_ = mu_Cr(ii) Andersson's mu seems to be off by 1/2, st it cannot be reduced to RCA case
	mu_ = mu_Cr(ii)/2.d0
	r_ = r_Cr(ii)
	gamma_1 = pes05_gamma(2, P_Cr, rho_Cr, ii, 0, pg_N2*2)
	call pes05_E_Crust(P_, rho_, m_, nu_, mu_, gamma_1, r_, E)
!write(*,*) ii, pg_N2*2, P_, rho_, m_, nu_, mu_, r_, gamma_1
!write(*,*) (P_Cr(pg_N2*2)-P_Cr(pg_N2*2-1))/(rho_Cr(pg_N2*2)-rho_Cr(pg_N2*2-1))*(rho_Cr(pg_N2*2-1)/P_Cr(pg_N2*2-1))*(1.d0 + P_Cr(pg_N2*2-1)/rho_Cr(pg_N2*2-1))
!write(*,*) rho_Cr(pg_N2*2)-rho_Cr(pg_N2*2-1), P_Cr(pg_N2*2-1), rho_Cr(pg_N2*2-1)

	pi4 = 4.d0*pi
	l_1 = l_0 * (l_0 + 1.d0)
	l_2_hf = (l_0 + 2.d0)*(l_0 - 1.d0)/2.d0
	eLamb = pes05_eLamb(rho_, m_, r_) ** 2
	eNu = dexp(2.d0*nu_)
	eLamb_hf = eLamb**(0.5d0)
	eNu_hf = eNu**(0.5d0)
	dLamb = pes05_dLamb(rho_, m_, r_) * 2.d0
	dNu = pes05_dnu(P_, rho_, m_, r_) * 2.d0
	ddNu = pes05_ddnu(P_, rho_, m_, r_) * 2.d0
	dP = -dNu/2.d0 * (P_ + rho_)
	V1 = (1.d0 + rho_/P_)*r_ *dNu/2.d0
	H2_C(1:6) = E(1,1:6)
	X_C(1:6) = E(2,1:6)
	T1_C(1:6) = E(3,1:6)

	B = 0.d0

	B(1,1) = 0.5d0*(dLamb - dNu)*r_ - (l_0 + 1.d0)
	B(1,2) = eLamb
	B(1,5) = -eLamb*4.d0*pi4*(P_ + rho_)
	B(1, 1:6) = B(1, 1:6) + eLamb*H2_C(1:6)	! Adding the algebraic solution of H2 & X & T1

	B(2,1) = l_2_hf + 1.d0
	B(2,2) = 0.5d0* dNu*r_ - (l_0 + 1.d0)
	B(2,4) = -2.d0*pi4*(P_ + rho_)*eLamb_hf
	B(2, 1:6) = B(2, 1:6) + H2_C(1:6)

	B(3,1) = -r_**2 * eNu**(-1.d0) * OmeC_sq
	B(3,2) = l_0
	B(3,3) = -(0.5d0* dNu*r_ + (l_0 - 1.d0))
	B(3,6) = -4.d0*pi4
	B(3, 1:6) = B(3, 1:6) + B(2, 1:6) - (0.5d0* dNu*r_ + 1.d0)*H2_C(1:6)

	B(4,2) = r_**2 * eLamb_hf
	B(4,4) = -(l_0 + 1.d0)
	B(4,5) = -eLamb_hf * l_1
	B(4, 1:6) = B(4, 1:6) + 0.5d0 * r_**2 * eLamb_hf * H2_C(1:6) +  r_**2 * eLamb_hf * eNu_hf**(-1.d0)/ gamma_1 / P_ * X_C(1:6)

	B(5,4) = eLamb_hf
	B(5,5) = (2.d0 - l_0)
	B(5,6) = 0.5d0/mu_

	B(6,3) = -0.5d0* r_**2 *eLamb * (P_ + rho_)
	B(6,4) = eLamb_hf * dP * r_
	B(6,5) = 4.d0 * l_2_hf * eLamb * mu_ - eLamb/eNu * r_**2 * OmeC_sq * (P_ + rho_)
	B(6,6) = 0.5d0 * (dLamb - dNu) * r_ - (l_0 + 1.d0)
	B(6, 1:6) = B(6, 1:6) + r_**2 * eLamb/eNu_hf * (X_C(1:6) - 0.5d0/r_**2 * eNu_hf * T1_C(1:6))

!!!!!!!! T1_C, X_C etc
	do i = 1,n/2
		fcnc(i) = 0.d0
		do j = 1,n/2
			fcnc(i) = fcnc(i) + B(i,j)*xc(j)/(1.d0 +V1)
		enddo
	enddo
	fcn(1:n/2) = dreal(fcnc(1:n/2))
	fcn(n/2+1 :n) = dimag(fcnc(1:n/2))

	endsubroutine pes05_vx_Cr

!	Vaccuum ¡¸///////////////////////////////////////
	subroutine pes05_vr_Ou(n, t, x, fcn)
	implicit none
	integer :: n
	real(8), dimension(1:n) :: x, fcn
	real(8) :: t
	complex(8), dimension(1:n/2,1:n/2) :: B
	complex(8), dimension(1:n/2) :: xc, fcnc
	real(8) :: V_z, m_, r_, tcf
	integer :: i, j

	if (mod(n,2) /= 0) then
		write(*,*) "err: puls eqt set n not even"
		pause
	endif
	xc(1:n/2) = dcmplx(x(1 :n/2),x(n/2+1 :n))

	m_ = M0
	r_ = t
	tcf = (1.d0 - 2.d0*m_/r_)**(-1.d0)		! tortoise coordinate factor: exp[(nu-lambda)/2]
	V_z = pes05_V_z(m_, r_)
	
	B(1,1) = (0.d0, 0.d0)
	B(1,2) = (1.d0, 0.d0)
	!B(2,1) = V_z - OmeC_sq
	!B(2,2) = (0.d0, 0.d0)
	B(2,1) = (V_z - OmeC_sq)* tcf **2
	B(2,2) = -2.d0*m_/r_**2 * tcf

	do i = 1,n/2
		fcnc(i) = 0.d0
		do j = 1,n/2
			fcnc(i) = fcnc(i) + B(i,j)*xc(j)
		enddo
	enddo
	fcn(1:n/2) = dreal(fcnc(1:n/2))
	fcn(n/2+1 :n) = dimag(fcnc(1:n/2))

	endsubroutine pes05_vr_Ou

	function rk4_HS(xi, dx,x)
	implicit none
	real(8) :: xi, dx, x
	real(8) :: dx_hf
	integer :: rk4_HS
		dx_hf = dx/2.d0
		rk4_HS = nint((x-xi)/dx_hf)
	endfunction rk4_HS
endmodule puls_eqt_set_opt05