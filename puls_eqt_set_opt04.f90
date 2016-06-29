module puls_eqt_set_opt04
! nonradial oscillations newtonian full
! FM is used for fluid model
!	Part1: functions
!	Part2: boundary conditions
!	Part3: differential equations f(x)=dy_dx(x)
use global_var
contains

!*********************************************************
!	Fcn---------------------------------------------------
	function pes04_V_temp(P_, rho_, m_, r_)
	implicit none
	real(8) :: pes04_V_temp, P_, rho_, m_, r_
	real(8) :: c1, c2, c3
		c1 = (1.d0 + P_/rho_/c**2)
		c2 = (1.d0 + 4.d0*pi*r_**3*P_/m_/c**2)
		c3 = (1.d0 - 2.d0*Grav_Const*m_/r_/c**2)**(-1.d0)

		pes04_V_temp = Grav_Const * rho_ * m_/ r_/P_*c1*c2*c3
	endfunction pes04_V_temp

	function pes04_V_tilde(P_, rho_, m_, r_)
	implicit none
	real(8) :: pes04_V_tilde, P_, rho_, m_, r_

		pes04_V_tilde = Grav_Const * rho_ * m_/ r_/P_
	endfunction pes04_V_tilde

	function pes04_gamma(mode, P_, rho_, ii, i1, i2)
	!	i1, i2 the range of P_ and rho_; must be consistent with ii
	implicit none
	integer :: ii, mode, i1, i2, w
	real(8), dimension(i1:i2) :: P_, rho_
	real(8) :: pes04_gamma
		
		w = ii
		if	(mode == 1) then 
			if (w == i2) w = i2 -1
			if (P_(w+1)*P_(w)*rho_(w+1)*rho_(w) == 0.d0) then
				write(*,*) "err: gamma_1 pulsation equations"
				pause
			endif
			pes04_gamma = (P_(w+1)-P_(w))/(rho_(w+1)-rho_(w))*(rho_(w)/P_(w))
		elseif (mode == 2) then
			if (w == i1) w = i1 +1
			if (P_(w)*P_(w-1)*rho_(w)*rho_(w-1) == 0.d0) then
				write(*,*) "err: gamma_1 pulsation equations"
				pause
			endif
			pes04_gamma = (P_(w)-P_(w-1))/(rho_(w)-rho_(w-1))*(rho_(w-1)/P_(w-1))
		endif
	endfunction pes04_gamma

	function pes04_U_tilde(rho_, m_, r_)
	implicit none
	real(8) :: pes04_U_Tilde, rho_, m_, r_
		pes04_U_Tilde = 4.d0*pi*r_**3*rho_/m_
	endfunction pes04_U_tilde

!*********************************************************
!	BC---------------------------------------------------
	subroutine pes04_bc_Coi(mode, y)
	implicit none
	real(8) :: c1
	real(8) :: y(1:4)
	real(8) :: r_, m_, U0, Y0, Phi0, g0
	integer :: mode
		
		r_ = r_Co(0)
		m_ = m_Co(0)
		g0 = Grav_Const*m_/r_**3
		c1 = (r_/R0)**3 *M0/m_

		if (mode == 1) then
			U0 = -1.d0
			Y0 = c1*Omega_sq/l_0 * U0
			Phi0 = 0.d0
		elseif (mode == 2) then
			U0 = -1.d0
			Y0 = 0.d0
			Phi0 = c1*Omega_sq/l_0 * U0
		endif

		
		y(1) = U0
		y(2) = r_**(l_0-1.d0)/g0*(Y0 + Phi0)
		y(3) = r_**(l_0-1.d0)/g0* Phi0
		y(4) = l_0 * r_**(l_0-1.d0)/g0* Phi0

	endsubroutine pes04_bc_Coi

	subroutine pes04_bc_Cof(y, z_M)
	implicit none
	real(8), dimension(1:5) :: z_M
	real(8) :: y(1:4), V_tilde
	integer :: ii

	ii = 2*pg_N1
	V_tilde = pes04_V_tilde(P_Co(ii), rho_Co(ii), m_Co(ii), r_Co(ii))

	z_M(1) = y(1)
	z_M(2) = V_tilde*(y(1) - y(2) + y(3))
	z_M(3) = 0.d0
	z_M(4) = y(3)
	z_M(5) = y(4)
	endsubroutine pes04_bc_Cof

	subroutine pes04_bc_Cri(z, z_M)
	implicit none
	real(8), dimension(1:5) :: z_M
	real(8) :: z(1:6)
	integer :: i

	z_M(1) = z(1)
	z_M(2) = z(2)
	z_M(3) = z(4)
	z_M(4) = z(5)
	z_M(5) = z(6)
	endsubroutine pes04_bc_Cri

	subroutine pes04_bc_Crf(mode,z)
	!	2 components
	implicit none
	real(8) :: z(1:6)
	integer :: mode
		z = 0.d0
		if (mode == 1) then
			z(1) = 1.d0
		elseif (mode == 2) then
			z(2) = 1.d0
		elseif (mode == 3) then
			z(3) = 1.d0
		elseif (mode == 4) then
			z(5) = 1.d0
		elseif (mode == 5) then
			z(6) = 1.d0
		else
			write(*,*) "err: pulsation equations crust surface"
			pause
		endif
	endsubroutine pes04_bc_Crf

	subroutine pes04_bc_Oci(mode, y, z)
	implicit none
	real(8) :: y(1:4), z(1:6), V_tilde
	integer :: mode

	endsubroutine pes04_bc_Oci

	subroutine pes04_bc_Ocf(y)
	implicit none
	real(8) :: y(1:4)

	endsubroutine pes04_bc_Ocf

	subroutine pes04_bc_Co_FM_f(y, z_M)
	implicit none
	real(8), dimension(1:4) :: z_M
	real(8) :: y(1:4), V_tilde
	integer :: 	ii = 2*pg_N1

	V_tilde = pes04_V_tilde(P_Co(ii), rho_Co(ii), m_Co(ii), r_Co(ii))

	z_M(1) = y(1)
	z_M(2) = y(2)
	z_M(3) = y(3)
	z_M(4) = y(4)
	endsubroutine pes04_bc_Co_FM_f

	! Fluid Model
	subroutine pes04_bc_Cr_FM_i(y, z_M)
	implicit none
	real(8), dimension(1:4) :: z_M
	real(8) :: y(1:4), V_tilde
	integer :: ii= 0

	V_tilde = pes04_V_tilde(P_Cr(ii), rho_Cr(ii), m_Cr(ii), r_Cr(ii))

	z_M(1) = y(1)
	z_M(2) = y(2)
	z_M(3) = y(3)
	z_M(4) = y(4)
	endsubroutine pes04_bc_Cr_FM_i

	subroutine pes04_bc_Cr_FM_f(mode, y)
	!	2 components
	implicit none
	real(8) :: rho_, m_, r_, y(1:4)
	integer :: mode, ii
		ii = 2*pg_N2
		y = 0.d0
		if (mode == 1) then
			y(2) = 1.d0
			y(3) = 1.d0
			y(4) = (-l_0-1.d0)
		elseif (mode == 2) then
			rho_ = rho_Cr(ii) 
			m_ = m_Cr(ii)
			r_ = r_Cr(ii)

			y(1)=1.D0
			y(3)=-1.D0
			y(4)=l_0+1.d0-4.D0*pi*rho_*r_**3/m_
		else
			write(*,*) "err: pulsation equations crust surface"
			pause
		endif
	endsubroutine pes04_bc_Cr_FM_f

!*********************************************************
!	Eqts---------------------------------------------------
	subroutine pes04_vx_Co(n, t, x, fcn)
	implicit none
	integer :: n
	real(8), dimension(1:n) :: x, fcn
	real(8) :: t
	real(8), dimension(1:n,1:n) :: B
	real(8) :: V_tilde, gamma_1, c1, A, U_tilde, l_, P_, rho_, m_ , r_
	integer :: i, j, ii

	if (n /= 4) then
		write(*,*) "err: puls eqt set: pes04 n /= 4"
		pause
	endif

	ii = rk4_HS(XCo_i, dx_Co,t)
	P_ = P_Co(ii)
	rho_ = rho_Co(ii)
	m_ = m_Co(ii)
	r_ = r_Co(ii)
	V_tilde = pes04_V_tilde(P_, rho_, m_, r_)
	gamma_1 = pes04_gamma(1, P_Co, rho_Co, ii, 0, pg_N1*2)
	c1 = r_**3/R0**3*M0/m_
	A = 0.d0
	U_tilde = pes04_U_tilde(rho_, m_, r_)
	l_ = l_0*(l_0+1.d0)

	B = 0.d0
	B(1,1) = V_tilde/gamma_1 -3.d0
	B(1,2) = l_/c1/Omega_sq -V_tilde/gamma_1
	B(1,3) = V_tilde/gamma_1

	B(2,1) = c1*Omega_sq +A*r_
	B(2,2) = 1.d0 -U_tilde -A*r_
	B(2,3) = A*r_

	B(3,3) = 1.d0 -U_tilde
	B(3,4) = 1.d0

	B(4,1) = -U_tilde* A *r_
	B(4,2) = U_tilde*V_tilde/gamma_1
	B(4,3) = l_ - U_tilde*V_tilde/gamma_1
	B(4,4) = -U_tilde

	do i = 1,n
		fcn(i) = 0.d0
		do j = 1,n
			fcn(i) = fcn(i) + B(i,j)*x(j)/(1.d0 +V_tilde)
		enddo
	enddo
	endsubroutine pes04_vx_Co

	subroutine pes04_vx_Cr(n, t, x, fcn)
	implicit none
	integer :: n
	real(8), dimension(1:n) :: x, fcn
	real(8) :: t
	real(8), dimension(1:n,1:n) :: B
	real(8) :: V_tilde, gamma_1, c1, A, U_tilde, l_, a1, a2, a3, P_, rho_, m_, mu_, r_
	integer :: i, j, ii

	if (n /= 6) then
		write(*,*) "err: puls eqt set: pes04 n /= 6"
		pause
	endif

	ii = rk4_HS(XCr_i, dx_Cr,t)
	P_ = P_Cr(ii)
	rho_ = rho_Cr(ii)
	m_ = m_Cr(ii)
	mu_ = mu_Cr(ii)
	r_ = r_Cr(ii)
	V_tilde = pes04_V_tilde(P_, rho_, m_, r_)
	gamma_1 = pes04_gamma(2, P_Cr, rho_Cr, ii, 0, pg_N2*2)
	c1 = r_**3/R0**3*M0/m_
	A = 0.d0
	U_tilde = pes04_U_tilde(rho_, m_, r_)
	l_ = l_0*(l_0+1.d0)
	a1 = mu_/P_
	a2 = gamma_1 - 2.d0/3.d0*a1
	a3 = gamma_1 + 4.d0/3.d0*a1

	B = 0.d0
	B(1,1) = -(1.d0 +2.d0*a2/a3)
	B(1,2) = 1.d0/a3
	B(1,3) = l_*a2/a3

	B(2,1) = -c1*V_tilde*Omega_sq -4.d0*V_tilde +U_tilde*V_tilde + 12.d0*gamma_1*a1/a3
	B(2,2) = V_tilde -4.d0*a1/a3
	B(2,3) = l_*(V_tilde -6.d0*gamma_1*a1/a3)
	B(2,4) = l_
	B(2,6) = V_tilde

	B(3,1) = -1.d0
	B(3,4) = 1.d0/a1

	B(4,1) = V_tilde - 6.d0*gamma_1*a1/a3
	B(4,2) = -a2/a3
	B(4,3) = -c1*V_tilde*Omega_sq + 2.d0/a3*((2.d0*l_ -1.d0)*a1*a2 + 2.d0*(l_ - 1.d0)*a1**2)
	B(4,4) = V_tilde -3.d0
	B(4,5) = V_tilde

	B(5,5) = 1.d0 - U_tilde
	B(5,6) = 1.d0

	B(6,1) = U_tilde * (- A*r_ + V_tilde/gamma_1 - 2.d0 + 2.d0*a2/a3)
	B(6,2) = -U_tilde/a3
	B(6,3) = l_ * U_tilde * (1.d0 - a2/a3)
	B(6,5) = l_
	B(6,6) = -U_tilde

	do i = 1,n
		fcn(i) = 0.d0
		do j = 1,n
			fcn(i) = fcn(i) + B(i,j)*x(j)/(1.d0 +V_tilde)
		enddo
	enddo

!if (r_ > R_mid .and. r_ <= R_mid*1.1d0) then
!write(*,*) pes04_V_temp(P_, rho_, m_, r_)
!write(*,*) pes04_V_tilde(P_, rho_, m_, r_)
!pause
!endif

	endsubroutine pes04_vx_Cr

	subroutine pes04_vx_Oc(n, t, x, fcn)
	implicit none
	integer :: n
	real(8), dimension(1:n) :: x, fcn
	real(8) :: t
	real(8), dimension(1:n,1:n) :: B
	real(8) :: V_tilde, gamma_1, c1, A, U_tilde, l_, P_, rho_, m_ , r_
	integer :: i, j, ii

	if (n /= 4) then
		write(*,*) "err: puls eqt set: pes04 n /= 4"
		pause
	endif

	ii = rk4_HS(XOc_i, dx_Oc,t)
	P_ = P_Oc(ii)
	rho_ = rho_Oc(ii)
	m_ = m_Oc(ii)
	r_ = r_Oc(ii)
	V_tilde = pes04_V_tilde(P_, rho_, m_, r_)
	gamma_1 = pes04_gamma(1, P_Oc, rho_Oc, ii, 0, pg_N3*2)
	c1 = r_**3/R0**3*M0/m_
	A = 0.d0
	U_tilde = pes04_U_tilde(rho_, m_, r_)
	l_ = l_0*(l_0+1.d0)

	B = 0.d0
	B(1,1) = V_tilde/gamma_1 -3.d0
	B(1,2) = l_/c1/Omega_sq -V_tilde/gamma_1
	B(1,3) = V_tilde/gamma_1

	B(2,1) = c1*Omega_sq +A*r_
	B(2,2) = 1.d0 -U_tilde -A*r_
	B(2,3) = A*r_

	B(3,3) = 1.d0 -U_tilde
	B(3,4) = 1.d0

	B(4,1) = -U_tilde* A *r_
	B(4,2) = U_tilde*V_tilde/gamma_1
	B(4,3) = l_ - U_tilde*V_tilde/gamma_1
	B(4,4) = -U_tilde

	do i = 1,n
		fcn(i) = 0.d0
		do j = 1,n
			fcn(i) = fcn(i) + B(i,j)*x(j)/(1.d0 +V_tilde)
		enddo
	enddo
	endsubroutine pes04_vx_Oc

	subroutine pes04_vx_Cr_FM(n, t, x, fcn)
	implicit none
	integer :: n
	real(8), dimension(1:n) :: x, fcn
	real(8) :: t
	real(8), dimension(1:n,1:n) :: B
	real(8) :: V_tilde, gamma_1, c1, A, U_tilde, l_, P_, rho_, m_ , r_
	integer :: i, j, ii

	if (n /= 4) then
		write(*,*) "err: puls eqt set: pes04 n /= 4"
		pause
	endif

	ii = rk4_HS(XCr_i, dx_Cr,t)
	P_ = P_Cr(ii)
	rho_ = rho_Cr(ii)
	m_ = m_Cr(ii)
	r_ = r_Cr(ii)
	V_tilde = pes04_V_tilde(P_, rho_, m_, r_)
	gamma_1 = pes04_gamma(1, P_Cr, rho_Cr, ii, 0, pg_N2*2)
	c1 = r_**3/R0**3*M0/m_
	A = 0.d0
	U_tilde = pes04_U_tilde(rho_, m_, r_)
	l_ = l_0*(l_0+1.d0)

	B = 0.d0
	B(1,1) = V_tilde/gamma_1 -3.d0
	B(1,2) = l_/c1/Omega_sq -V_tilde/gamma_1
	B(1,3) = V_tilde/gamma_1

	B(2,1) = c1*Omega_sq +A*r_
	B(2,2) = 1.d0 -U_tilde -A*r_
	B(2,3) = A*r_

	B(3,3) = 1.d0 -U_tilde
	B(3,4) = 1.d0

	B(4,1) = -U_tilde* A *r_
	B(4,2) = U_tilde*V_tilde/gamma_1
	B(4,3) = l_ - U_tilde*V_tilde/gamma_1
	B(4,4) = -U_tilde

	do i = 1,n
		fcn(i) = 0.d0
		do j = 1,n
			fcn(i) = fcn(i) + B(i,j)*x(j)/(1.d0 +V_tilde)
		enddo
	enddo
	endsubroutine pes04_vx_Cr_FM

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

endmodule puls_eqt_set_opt04