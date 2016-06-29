module puls_eqt_set_opt01
! nonradial oscillations newtonian cowling approximation
! FM is used for fluid model
!	Part1: functions
!	Part2: boundary conditions
!	Part3: differential equations f(x)=dy_dx(x)
use global_var
contains

!*********************************************************
!	Fcn---------------------------------------------------
	function pes01_V_temp(P_, rho_, m_, r_)
	implicit none
	real(8) :: pes01_V_temp, P_, rho_, m_, r_
	real(8) :: c1, c2, c3
		c1 = (1.d0 + P_/rho_/c**2)
		c2 = (1.d0 + 4.d0*pi*r_**3*P_/m_/c**2)
		c3 = (1.d0 - 2.d0*Grav_Const*m_/r_/c**2)**(-1.d0)

		pes01_V_temp = Grav_Const * rho_ * m_/ r_/P_*c1*c2*c3
	endfunction pes01_V_temp

	function pes01_V_tilde(P_, rho_, m_, r_)
	implicit none
	real(8) :: pes01_V_tilde, P_, rho_, m_, r_

		pes01_V_tilde = Grav_Const * rho_ * m_/ r_/P_
	endfunction pes01_V_tilde

	function pes01_gamma(mode, P_, rho_, ii, i1, i2)
	!	i1, i2 the range of P_ and rho_; must be consistent with ii
	implicit none
	integer :: ii, mode, i1, i2, w
	real(8), dimension(i1:i2) :: P_, rho_
	real(8) :: pes01_gamma
		
		w = ii
		if	(mode == 1) then 
			if (w == i2) w = i2 -1
			if (P_(w+1)*P_(w)*rho_(w+1)*rho_(w) == 0.d0) then
				write(*,*) "err: gamma_1 pulsation equations"
				pause
			endif
			pes01_gamma = (P_(w+1)-P_(w))/(rho_(w+1)-rho_(w))*(rho_(w)/P_(w))
		elseif (mode == 2) then
			if (w == i1) w = i1 +1
			if (P_(w)*P_(w-1)*rho_(w)*rho_(w-1) == 0.d0) then
				write(*,*) "err: gamma_1 pulsation equations"
				pause
			endif
			pes01_gamma = (P_(w)-P_(w-1))/(rho_(w)-rho_(w-1))*(rho_(w-1)/P_(w-1))
		endif
	endfunction pes01_gamma

	function pes01_U_tilde(rho_, m_, r_)
	implicit none
	real(8) :: pes01_U_Tilde, rho_, m_, r_
		pes01_U_Tilde = 4.d0*pi*r_**3*rho_/m_
	endfunction pes01_U_tilde

!*********************************************************
!	BC---------------------------------------------------
	subroutine pes01_bc_Coi(y)
	implicit none
	real(8) :: c1
	real(8) :: y(1:2)
		c1 = (r_Co(0)/R0)**3 *M0/m_Co(0)
		y(1) = -1.d-4
		y(2) = c1*Omega_sq/l_0*y(1)
	endsubroutine pes01_bc_Coi

	subroutine pes01_bc_Cof(y, z_M)
	implicit none
	real(8), dimension(1:NMat, 1:NMat) :: z_M
	real(8) :: y(1:2), V_tilde
	integer :: ii
	ii= 2*pg_N1
	if (newt_V_opt == 1) then
		V_tilde = pes01_V_tilde(P_Co(ii), rho_Co(ii), m_Co(ii), r_Co(ii))
	elseif (newt_V_opt == 2 .or. newt_V_opt == 3) then
V_tilde = pes01_V_temp(P_Co(ii), rho_Co(ii), m_Co(ii), r_Co(ii))
	else
		write(*,*) "err: newt_V_opt"
		pause
	endif

	z_M(1,1) = y(1)
	z_M(2,1) = V_tilde*(y(1) - y(2))
	z_M(3,1) = 0.d0
	endsubroutine pes01_bc_Cof

	subroutine pes01_bc_Cri(mode,z, z_M)
	implicit none
	real(8), dimension(1:NMat, 1:NMat) :: z_M
	real(8) :: z(1:4)
	integer :: mode, i
	if (mode == 1) then
		i = 2
	elseif (mode == 2) then
		i = 3
	else
		write(*,*) "err: pulsation equations crust core interface bc"
		pause
	endif

	z_M(1,i) = z(1)
	z_M(2,i) = z(2)
	z_M(3,i) = z(4)
	endsubroutine pes01_bc_Cri

	subroutine pes01_bc_Crf(mode,z)
	!	2 components
	implicit none
	real(8) :: z(1:4)
	integer :: mode
		z = 0.d0
		if (mode == 1) then
			z(1) = 1.d0
		elseif (mode == 2) then
			z(3) = 1.d0
		else
			write(*,*) "err: pulsation equations crust surface"
			pause
		endif
	endsubroutine pes01_bc_Crf

	subroutine pes01_bc_Oci(mode, y, z)
	implicit none
	real(8) :: y(1:2), z(1:4), V_tilde
	integer :: mode
	V_tilde = pes01_V_tilde(P_Oc(0), rho_Oc(0), m_Oc(0), r_Oc(0))
		z = 0.d0
		if (mode == 1) then
			z(1) = y(1)
			z(2) = V_tilde*(y(1) - y(2))
		elseif (mode == 2) then
			z(3) = 1.d0
		else
			write(*,*) "err: pulsation equations crust surface"
			pause
		endif
	endsubroutine pes01_bc_Oci

	subroutine pes01_bc_Ocf(y)
	implicit none
	real(8) :: y(1:2)
		y(1) = 1.d0
		y(2) = 1.d0
	endsubroutine pes01_bc_Ocf

	subroutine pes01_bc_Co_FM_f(y, z_M)
	implicit none
	real(8), dimension(1:NMat, 1:NMat) :: z_M
	real(8) :: y(1:2), V_tilde
	integer :: ii
	
	ii= 2*pg_N1
	if (newt_V_opt == 1) then
		V_tilde = pes01_V_tilde(P_Co(ii), rho_Co(ii), m_Co(ii), r_Co(ii))
	elseif (newt_V_opt == 2 .or. newt_V_opt == 3) then
V_tilde = pes01_V_temp(P_Co(ii), rho_Co(ii), m_Co(ii), r_Co(ii))
	else
		write(*,*) "err: newt_V_opt"
		pause
	endif

		z_M(1,1) = y(1)
		z_M(2,1) = V_tilde*(y(1) - y(2))
	endsubroutine pes01_bc_Co_FM_f

	! Fluid Model
	subroutine pes01_bc_Cr_FM_i(y, z_M)
	implicit none
	real(8), dimension(1:NMat, 1:NMat) :: z_M
	real(8) :: y(1:2), V_tilde

	if (newt_V_opt == 1) then
		V_tilde = pes01_V_tilde(P_Cr(0), rho_Cr(0), m_Cr(0), r_Cr(0))
	elseif (newt_V_opt == 2 .or. newt_V_opt == 3) then
V_tilde = pes01_V_temp(P_Cr(0), rho_Cr(0), m_Cr(0), r_Cr(0))
	else
		write(*,*) "err: newt_V_opt"
		pause
	endif
		z_M(1,2) = y(1)
		z_M(2,2) = V_tilde*(y(1) - y(2))
	endsubroutine pes01_bc_Cr_FM_i

	subroutine pes01_bc_Cr_FM_f(y)
	!	2 components
	implicit none
	real(8) :: y(1:2)
	integer :: mode
		y(1) = 1.d0
		y(2) = 1.d0
	endsubroutine pes01_bc_Cr_FM_f

!*********************************************************
!	Eqts---------------------------------------------------
	subroutine pes01_vx_Co(n, t, x, fcn)
	implicit none
	integer :: n
	real(8), dimension(1:n) :: x, fcn
	real(8) :: t
	real(8), dimension(1:n,1:n) :: B
	real(8) :: V_tilde, gamma_1, c1, A, U_tilde, l_
real(8) :: e2Lamb	! test
				real(8) :: factor, A_, B_, C_
	integer :: i, j, ii

	ii = rk4_HS(XCo_i, dx_Co,t)

	if (newt_V_opt == 1) then
	V_tilde = pes01_V_tilde(P_Co(ii), rho_Co(ii), m_Co(ii), r_Co(ii))
	elseif (newt_V_opt == 2 .or. newt_V_opt == 3) then
V_tilde = pes01_V_temp(P_Co(ii), rho_Co(ii), m_Co(ii), r_Co(ii))
!V_tilde = pes01_V_tilde(P_Co(ii), rho_Co(ii), m_Co(ii), r_Co(ii))		!M4
	else
		write(*,*) "err: newt_V_opt"
		pause
	endif
e2Lamb = (1.d0 - 2.d0*Grav_Const*m_Co(ii)/c**2/r_Co(ii))**(-1.d0)
	gamma_1 = pes01_gamma(1, P_Co, rho_Co, ii, 0, pg_N1*2)
	c1 = r_Co(ii)**3/R0**3*M0/m_Co(ii)
	A = 0.d0
				factor = (1.d0 +pes01_V_temp(P_Co(ii), rho_Co(ii), m_Co(ii), r_Co(ii)))/(1.d0 + pes01_V_tilde(P_Co(ii), rho_Co(ii), m_Co(ii), r_Co(ii)))
!				A_ = (1.d0 + P_Co(ii)/rho_Co(ii)/c**2)
!				B_ = (1.d0 + 4.d0 * pi* r_Co(ii) **3 * P_Co(ii)/m_Co(ii)/c**2)
!				C_ = 1.d0 / (1.d0 - 2.d0 * Grav_Const*m_Co(ii)/c**2/r_Co(ii))
	U_tilde = pes01_U_tilde(rho_Co(ii), m_Co(ii), r_Co(ii))
	l_ = l_0*(l_0+1.d0)

	B(1,1) = V_tilde/gamma_1 -3.d0
!				B(1,2) = l_/c1/Omega_sq/dexp(-2.d0*nu_Co(ii)) -V_tilde/gamma_1	 ! freqeucy redshift
!				B(1,2) = l_/c1/Omega_sq*pes01_V_temp(P_Co(ii), rho_Co(ii), m_Co(ii), r_Co(ii))/V_tilde -V_tilde/gamma_1	 ! M8
!				B(1,2) = l_/c1/Omega_sq/dexp(-2.d0*nu_Co(ii))*C_*B_ -V_tilde/gamma_1	 !M7: M2 to RCA
	B(1,2) = l_/c1/Omega_sq -V_tilde/gamma_1

!				B(2,1) = c1*Omega_sq*dexp(-2.d0*nu_Co(ii)) +A*r_Co(ii) ! freqeucy redshift
!				B(2,1) = c1*Omega_sq*factor +A*r_Co(ii)		! M3 reduce to M1
!				B(2,1) = c1*Omega_sq*dexp(-2.d0*nu_Co(ii))/B_ +A*r_Co(ii)		!M7: M2 to RCA
	B(2,1) = c1*Omega_sq +A*r_Co(ii)
	B(2,2) = 1.d0 -U_tilde -A*r_Co(ii)
!				B(2,1) = B(2,1)/(1.d0 - Grav_Const*m_Co(ii)/c**2/r_Co(ii))		! M5: M3 correct for f-mode

	do i = 1,n
		fcn(i) = 0.d0
		do j = 1,n
			fcn(i) = fcn(i) + B(i,j)*x(j)/(1.d0 +V_tilde)
		enddo
	enddo
	endsubroutine pes01_vx_Co

	subroutine pes01_vx_Cr(n, t, x, fcn)
	implicit none
	integer :: n
	real(8), dimension(1:n) :: x, fcn
	real(8) :: t
	real(8), dimension(1:n,1:n) :: B
	real(8) :: V_tilde, gamma_1, c1, A, U_tilde, l_, a1, a2, a3, U1, U2 
	real(8) :: e2Lamb	! test
				real(8) :: factor
	integer :: i, j, ii

	ii = rk4_HS(XCr_i, dx_Cr,t)

	if (newt_V_opt == 1) then
	V_tilde = pes01_V_tilde(P_Cr(ii), rho_Cr(ii), m_Cr(ii), r_Cr(ii))
	elseif (newt_V_opt == 2 .or. newt_V_opt == 3) then
V_tilde = pes01_V_temp(P_Cr(ii), rho_Cr(ii), m_Cr(ii), r_Cr(ii))
	else
		write(*,*) "err: newt_V_opt"
		pause
	endif
e2Lamb = (1.d0 - 2.d0*Grav_Const*m_Cr(ii)/c**2/r_Cr(ii))**(-1.d0)
gamma_1 = pes01_gamma(2, P_Cr, rho_Cr, ii, 0, pg_N2*2)
	c1 = r_Cr(ii)**3/R0**3*M0/m_Cr(ii)
	A = 0.d0
				factor = (1.d0 +pes01_V_temp(P_Cr(ii), rho_Cr(ii), m_Cr(ii), r_Cr(ii)))/(1.d0 + pes01_V_tilde(P_Cr(ii), rho_Cr(ii), m_Cr(ii), r_Cr(ii)))
	U_tilde = pes01_U_tilde(rho_Cr(ii), m_Cr(ii), r_Cr(ii))
U1 = - Grav_Const /c**2 *2.d0 * m_Cr(ii)/r_Cr(ii) /(1.d0 - Grav_Const /c**2 *2.d0 * m_Cr(ii)/r_Cr(ii))
U2 = - Grav_Const /c**2 *m_Cr(ii)/r_Cr(ii) /(1.d0 - Grav_Const /c**2 *2.d0 * m_Cr(ii)/r_Cr(ii))
	l_ = l_0*(l_0+1.d0)
!a1 = mu_Cr(ii)/P_Cr(ii) /dexp(-2.d0 * nu_Cr(ii)) ! freqeucy redshift
	a1 = mu_Cr(ii)/P_Cr(ii)
	a2 = gamma_1 - 2.d0/3.d0*a1
	a3 = gamma_1 + 4.d0/3.d0*a1

	B(1,1) = -(1.d0 +2.d0*a2/a3)
	B(1,2) = 1.d0/a3
	B(1,3) = l_*a2/a3
	B(1,4) = 0.d0

!				B(2,1) = -c1*V_tilde*Omega_sq*dexp(-2.d0*nu_Cr(ii)) -4.d0*V_tilde +U_tilde*V_tilde + 12.d0*gamma_1*a1/a3 ! freqeucy redshift	! M6
	B(2,1) = -c1*V_tilde*Omega_sq -4.d0*V_tilde +U_tilde*V_tilde + 12.d0*gamma_1*a1/a3
	B(2,2) = V_tilde -4.d0*a1/a3
	B(2,3) = l_*(V_tilde -6.d0*gamma_1*a1/a3)
	B(2,4) = l_

	B(3,1) = -1.d0
	B(3,2) = 0.d0
	B(3,3) = 0.d0
!				B(3,4) = 1.d0/a1 *e2Lamb		! gives correct s modes by keeping dPdr = TOV	!M3 M4 M5 M6
	B(3,4) = 1.d0/a1

	B(4,1) = V_tilde - 6.d0*gamma_1*a1/a3
	B(4,2) = -a2/a3
!				B(4,3) = -c1*V_tilde*Omega_sq*dexp(-2.d0*nu_Cr(ii)) + 2.d0/a3*((2.d0*l_ -1.d0)*a1*a2 + 2.d0*(l_ - 1.d0)*a1**2)		! freqeucy redshift
	B(4,3) = -c1*V_tilde*Omega_sq + 2.d0/a3*((2.d0*l_ -1.d0)*a1*a2 + 2.d0*(l_ - 1.d0)*a1**2)
	B(4,4) = V_tilde -3.d0
!				B(1:2,3) = B(1:2,3)/dexp(-2.d0*nu_Cr(ii))			! M8	fixing i-mode
	
	do i = 1,n
		fcn(i) = 0.d0
		do j = 1,n
			fcn(i) = fcn(i) + B(i,j)*x(j)/(1.d0 +V_tilde)
		enddo
	enddo

!if (ii == nint(0.2d0 * 2* pg_N2)) then
!write(*,*) "B(2,1), B(2,2), B(2,3), B(2,4)", B(2,1), B(2,2), B(2,3), B(2,4)
!write(*,*) "B(4,1), B(4,2), B(4,3), B(4,4)", B(4,1), B(4,2), B(4,3), B(4,4)
!pause
!endif

	endsubroutine pes01_vx_Cr

	subroutine pes01_vx_Oc(n, t, x, fcn)
	implicit none
	integer :: n
	real(8), dimension(1:n) :: x, fcn
	real(8) :: t
	real(8), dimension(1:n,1:n) :: B
	real(8) :: V_tilde, gamma_1, c1, A, U_tilde, l_
	integer :: i, j, ii

	ii = rk4_HS(XOc_i, dx_Oc,t)

	if (newt_V_opt == 1) then
	V_tilde = pes01_V_tilde(P_Oc(ii), rho_Oc(ii), m_Oc(ii), r_Oc(ii))
	elseif (newt_V_opt == 2 .or. newt_V_opt == 3) then
V_tilde = pes01_V_temp(P_Oc(ii), rho_Oc(ii), m_Oc(ii), r_Oc(ii))
	else
		write(*,*) "err: newt_V_opt"
		pause
	endif

	gamma_1 = pes01_gamma(1, P_Oc, rho_Oc, ii, 0, pg_N3*2)
	c1 = r_Oc(ii)**3/R0**3*M0/m_Oc(ii)
	A = 0.d0
	U_tilde = pes01_U_tilde(rho_Oc(ii), m_Oc(ii), r_Oc(ii))
	l_ = l_0*(l_0+1.d0)

	B(1,1) = V_tilde/gamma_1 -3.d0
	B(1,2) = l_/c1/Omega_sq -V_tilde/gamma_1
	B(2,1) = c1*Omega_sq +A*r_Oc(ii)
	B(2,2) = 1.d0 -U_tilde -A*r_Oc(ii)

	do i = 1,n
		fcn(i) = 0.d0
		do j = 1,n
			fcn(i) = fcn(i) + B(i,j)*x(j)/(1.d0 +V_tilde)
		enddo
	enddo
	endsubroutine pes01_vx_Oc

	subroutine pes01_vx_Cr_FM(n, t, x, fcn)
	implicit none
	integer :: n
	real(8), dimension(1:n) :: x, fcn
	real(8) :: t
	real(8), dimension(1:n,1:n) :: B
	real(8) :: V_tilde, gamma_1, c1, A, U_tilde, l_
				real(8) :: factor, A_, B_, C_
	integer :: i, j, ii
!real(8) :: V2, ep

	ii = rk4_HS(XCr_i, dx_Cr,t)

	if (newt_V_opt == 1) then
	V_tilde = pes01_V_tilde(P_Cr(ii), rho_Cr(ii), m_Cr(ii), r_Cr(ii))
	elseif (newt_V_opt == 2 .or. newt_V_opt == 3) then
V_tilde = pes01_V_temp(P_Cr(ii), rho_Cr(ii), m_Cr(ii), r_Cr(ii))
!V_tilde = pes01_V_tilde(P_Cr(ii), rho_Cr(ii), m_Cr(ii), r_Cr(ii))		!M4
	else
		write(*,*) "err: newt_V_opt"
		pause
	endif
!V2 = pes01_V_temp(P_Cr(ii), rho_Cr(ii), m_Cr(ii), r_Cr(ii))
!ep = (1.d0 + V_tilde)/(1.d0 + V2)
	gamma_1 = pes01_gamma(2, P_Cr, rho_Cr, ii, 0, pg_N2*2)
	c1 = r_Cr(ii)**3/R0**3*M0/m_Cr(ii)
	A = 0.d0
				factor = (1.d0 +pes01_V_temp(P_Cr(ii), rho_Cr(ii), m_Cr(ii), r_Cr(ii)))/(1.d0 + pes01_V_tilde(P_Cr(ii), rho_Cr(ii), m_Cr(ii), r_Cr(ii)))
!				A_ = (1.d0 + P_Cr(ii)/rho_Cr(ii)/c**2)
!				B_ = (1.d0 + 4.d0 * pi* r_Cr(ii) **3 * P_Cr(ii)/m_Cr(ii)/c**2)
!				C_ = 1.d0 / (1.d0 - 2.d0 * Grav_Const*m_Cr(ii)/c**2/r_Cr(ii))
	U_tilde = pes01_U_tilde(rho_Cr(ii), m_Cr(ii), r_Cr(ii))
	l_ = l_0*(l_0+1.d0)

	B(1,1) = V_tilde/gamma_1 -3.d0
!				B(1,2) = l_/c1/Omega_sq/dexp(-2.d0*nu_Cr(ii)) -V_tilde/gamma_1					!M8: M1 to RCA
!				B(1,2) = l_/c1/Omega_sq/dexp(-2.d0*nu_Cr(ii))*C_*B_ -V_tilde/gamma_1			!M7: M2 to RCA
	B(1,2) = l_/c1/Omega_sq -V_tilde/gamma_1

!				B(2,1) = c1*Omega_sq*dexp(-2.d0*nu_Cr(ii)) +A*r_Cr(ii)
!				B(2,1) = c1*Omega_sq*factor +A*r_Cr(ii)			! M3 reduce to M1
!				B(2,1) = c1*Omega_sq*dexp(-2.d0*nu_Cr(ii))/B_ +A*r_Cr(ii)			!M7: M2 to RCA
	B(2,1) = c1*Omega_sq +A*r_Cr(ii)
	B(2,2) = 1.d0 -U_tilde -A*r_Cr(ii)

!B(2,1) = c1*Omega_sq +ep*A*r_Cr(ii) + (1.d0 - ep) - (1.d0 - ep) *U_tilde  - V_tilde + ep *V2 + V_tilde/gamma_1 - ep*V2/gamma_1
!B(2,2) = ep - ep * U_tilde -ep*A*r_Cr(ii)	+ V_tilde - ep *V2 - V_tilde/gamma_1 + ep*V2/gamma_1

	do i = 1,n
		fcn(i) = 0.d0
		do j = 1,n
			fcn(i) = fcn(i) + B(i,j)*x(j)/(1.d0 +V_tilde)
		enddo
	enddo

!if (ii == nint(0.1d0 * 2* pg_N2)) then
!write(*,*) "B(1,1), B(1,2)", B(1,1), B(1,2)
!write(*,*) "B(2,1), B(2,2)", B(2,1), B(2,2)
!write(*,*) " opt 1 B(1,1)", B(1,1)/ (1.d0 - 2.d0 * Grav_Const * m_Cr(ii)/c**2/r_Cr(ii))
!write(*,*) " opt 1 B(1,2)", B(1,2)/ (1.d0 - 2.d0 * Grav_Const * m_Cr(ii)/c**2/r_Cr(ii))
!write(*,*) " opt 1 B(2,1)", B(2,1)/ (1.d0 - 2.d0 * Grav_Const * m_Cr(ii)/c**2/r_Cr(ii))
!write(*,*) " opt 1 B(2,2)", B(2,2)/ (1.d0 - 2.d0 * Grav_Const * m_Cr(ii)/c**2/r_Cr(ii))
!pause
!endif
	endsubroutine pes01_vx_Cr_FM

!*********************************************************
!	Half Steps for RK4---------------------------------------------------
	function rk4_HS(xi, dx,x)
	implicit none
	real(8) :: xi, dx, x
	real(8) :: dx_hf
	integer :: rk4_HS
		dx_hf = dx/2.d0
		rk4_HS = nint((x-xi)/dx_hf)
	endfunction rk4_HS

endmodule puls_eqt_set_opt01