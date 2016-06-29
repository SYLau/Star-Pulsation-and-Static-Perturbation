module puls_eqt_set_opt02
! nonradial oscillations relativistic cowling approximation (Yoshida 2002)
! FM is used for fluid model
!	Part1: functions
!	Part2: boundary conditions
!	Part3: differential equations f(x)=dy_dx(x)
use global_var
contains

!*********************************************************
!	Fcn---------------------------------------------------
	function pes02_eLamb(rho_, m_, r_)
	implicit none
	real(8) :: pes02_eLamb, rho_, m_, r_
		pes02_eLamb = (1.d0 - rel * 2.d0*Grav_Const*m_/r_/c**2)**(-0.5d0)
	endfunction pes02_eLamb

	function pes02_dLamb(rho_, m_, r_)
	implicit none
	real(8) :: pes02_dLamb, rho_, m_, r_, c1
		c1 = pes02_eLamb(rho_, m_, r_) **2
		pes02_dLamb = c1 / r_ * rel * (Grav_Const/c**2)*(4.d0*pi*r_**2*rho_ - m_/r_)
	endfunction pes02_dLamb

	function pes02_dnu(P_, rho_, m_, r_)
	implicit none
	real(8) :: pes02_dnu, P_, rho_, m_, r_, c1
		c1 = pes02_eLamb(rho_, m_, r_) **2
		pes02_dnu = c1 / r_ * rel * (Grav_Const/c**2)*( rel * 4.d0*pi*r_**(2)*P_/c**(2) + m_/r_)
	endfunction pes02_dnu

	function pes02_V1(P_, rho_, m_, r_)
	implicit none
	real(8) :: pes02_V1, P_, rho_, m_, r_
	real(8) :: c1
		c1 = (rho_*c**2 /P_/ rel + 1.d0)
		pes02_V1 = pes02_dnu(P_, rho_, m_, r_)*r_*c1
	endfunction pes02_V1

	function pes02_V2(P_, rho_, m_, r_)
	implicit none
	real(8) :: pes02_V2, P_, rho_, m_, r_
	real(8) :: c1
		c1 = (rho_*c**2/P_) / rel
		pes02_V2 = pes02_dnu(P_, rho_, m_, r_)*r_*c1
	endfunction pes02_V2

	function pes02_gamma(mode, P_, rho_, ii, i1, i2)
	!	i1, i2 the range of P_ and rho_; must be consistent with ii
	implicit none
	integer :: ii, mode, i1, i2, w
	real(8), dimension(i1:i2) :: P_, rho_
	real(8) :: c1, pes02_gamma
		
		w = ii
		if	(mode == 1) then 
			if (w == i2) w = i2 -1
			if (P_(w+1)*P_(w)*rho_(w+1)*rho_(w) == 0.d0) then
				write(*,*) "err: gamma_1 pulsation equations"
				pause
			endif
			c1 = (1.d0 + rel * P_(w)/rho_(w)/c**2)
			pes02_gamma = (P_(w+1)-P_(w))/(rho_(w+1)-rho_(w))*(rho_(w)/P_(w))*c1
		elseif (mode == 2) then
			if (w == i1) w = i1 +1
			if (P_(w)*P_(w-1)*rho_(w)*rho_(w-1) == 0.d0) then
				write(*,*) "err: gamma_1 pulsation equations"
				pause
			endif
			c1 = (1.d0 + rel * P_(w-1)/rho_(w-1)/c**2)
			pes02_gamma = (P_(w)-P_(w-1))/(rho_(w)-rho_(w-1))*(rho_(w-1)/P_(w-1))*c1
		endif
	endfunction pes02_gamma

	function pes02_U1(P_, rho_, m_, r_)
	implicit none
	real(8) :: pes02_U1, P_, rho_, m_, r_, c1, dLamb, dP, dnu
	real(8) :: B(1:3)
		c1 = pes02_eLamb(rho_, m_, r_) **2
		dLamb = pes02_dLamb(rho_, m_, r_)
		dnu = pes02_dnu(P_, rho_, m_, r_)
		dP = -dnu *(rho_* c**2 + P_)
		
		B(1) = 2.d0*dLamb*(rel * 4.d0*pi*r_**2*P_/c**2 + m_/r_)
		B(2) = rel * 4.d0*pi/c**2*(2.d0*r_*P_+r_**(2)*dP)
		B(3) = 4.d0*pi*r_*rho_ - m_/r_**2
		pes02_U1 = rel * (Grav_Const/c**2) * c1 * (B(1) +B(2) +B(3))/dnu
	endfunction pes02_U1

	function pes02_U2(rho_, m_, r_)
	implicit none
	real(8) :: pes02_U2, rho_, m_, r_
		pes02_U2 = pes02_dLamb(rho_, m_, r_) * r_
	endfunction pes02_U2

	function pes02_c1(P_, rho_, m_, eNu_, r_)
	implicit none
	real(8) :: pes02_c1, P_, rho_, m_, eNu_, r_, dNu
		dNu = pes02_dNu(P_, rho_, m_, r_)
		pes02_c1 = rel * (Grav_Const/c**2)*M0/R0**(3)*r_* eNu_**(-2.d0 * rel) /dNu
	endfunction pes02_c1

!*********************************************************
!	BC---------------------------------------------------
	subroutine pes02_bc_Coi(y)
	implicit none
	real(8) :: c1, eNu, dNu
	real(8) :: y(1:2)
		eNu = dexp(nu_Co(0))
		c1 = pes02_c1(P_Co(0), rho_Co(0), m_Co(0), eNu, r_Co(0))
		y(1) = -1.d-4
		y(2) = c1*Omega_sq/l_0*y(1)
	endsubroutine pes02_bc_Coi

	subroutine pes02_bc_Cof(y, z_M)
	implicit none
	real(8), dimension(1:NMat, 1:NMat) :: z_M
	real(8) :: y(1:2), V1
	integer :: ii
	ii= 2*pg_N1
	V1 = pes02_V1(P_Co(ii), rho_Co(ii), m_Co(ii), r_Co(ii))

	z_M(1,1) = y(1)
	z_M(2,1) = V1*(y(1) - y(2))
	z_M(3,1) = 0.d0
	endsubroutine pes02_bc_Cof

	subroutine pes02_bc_Cri(mode,z, z_M)
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
	endsubroutine pes02_bc_Cri

	subroutine pes02_bc_Crf(mode,z)
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
	endsubroutine pes02_bc_Crf

	subroutine pes02_bc_Oci(mode, y, z)
	implicit none
	real(8) :: y(1:2), z(1:4), V1
	integer :: mode
	V1 = pes02_V1(P_Oc(0), rho_Oc(0), m_Oc(0), r_Oc(0))
		z = 0.d0
		if (mode == 1) then
			z(1) = y(1)
			z(2) = V1*(y(1) - y(2))
		elseif (mode == 2) then
			z(3) = 1.d0
		else
			write(*,*) "err: pulsation equations crust surface"
			pause
		endif
	endsubroutine pes02_bc_Oci

	subroutine pes02_bc_Ocf(y)
	implicit none
	real(8) :: y(1:2)
		y(1) = 1.d0
		y(2) = 1.d0
	endsubroutine pes02_bc_Ocf

	! Fluid Model
	subroutine pes02_bc_Co_FM_f(y, z_M)
	implicit none
	real(8), dimension(1:NMat, 1:NMat) :: z_M
	real(8) :: y(1:2), V1
	integer :: ii = 2*pg_N1
	V1 = pes02_V1(P_Co(ii), rho_Co(ii), m_Co(ii), r_Co(ii))
		z_M(1,1) = y(1)
		z_M(2,1) = V1*(y(1) - y(2))
	endsubroutine pes02_bc_Co_FM_f

	subroutine pes02_bc_Cr_FM_i(y, z_M)
	implicit none
	real(8), dimension(1:NMat, 1:NMat) :: z_M
	real(8) :: y(1:2), V1
	V1 = pes02_V1(P_Cr(0), rho_Cr(0), m_Cr(0), r_Cr(0))
		z_M(1,2) = y(1)
		z_M(2,2) = V1*(y(1) - y(2))
	endsubroutine pes02_bc_Cr_FM_i
	
	subroutine pes02_bc_Cr_FM_f(y)
	!	2 components
	implicit none
	real(8) :: y(1:2)
		y(1) = 1.d0
		y(2) = 1.d0
	endsubroutine pes02_bc_Cr_FM_f

!*********************************************************
!	Eqts---------------------------------------------------
	subroutine pes02_vx_Co(n, t, x, fcn)
	implicit none
	integer :: n
	real(8), dimension(1:n) :: x, fcn
	real(8) :: t
	real(8), dimension(1:n,1:n) :: B
	real(8) :: V1, V2, gamma_1, c1, A, U1, U2, l_, eNu, eLamb
	integer :: i, j, ii

	ii = rk4_HS(XCo_i, dx_Co,t)

	eNu = dexp(nu_Co(ii))
	eLamb = pes02_eLamb(rho_Co(ii), m_Co(ii), r_Co(ii))
	V1 = pes02_V1(P_Co(ii), rho_Co(ii), m_Co(ii), r_Co(ii))
	V2 = pes02_V2(P_Co(ii), rho_Co(ii), m_Co(ii), r_Co(ii))
	gamma_1 = pes02_gamma(1, P_Co, rho_Co, ii, 0, pg_N1*2)
	c1 = pes02_c1(P_Co(ii), rho_Co(ii), m_Co(ii), eNu, r_Co(ii))
	A = 0.d0
	U1 = pes02_U1(P_Co(ii), rho_Co(ii), m_Co(ii), r_Co(ii))
	U2 =pes02_U2(rho_Co(ii), m_Co(ii), r_Co(ii))
	l_ = l_0*(l_0+1.d0)

	B(1,1) = -(3.d0 - V1/gamma_1 + U2)
!B(1,1) = -(3.d0 - V1/gamma_1)
	B(1,2) = -(V1/gamma_1 - l_/c1/Omega_sq)
!B(1,2) = -(V1/gamma_1 - l_/c1/Omega_sq/(1.d0 - 2.d0*Grav_Const * m_Co(ii)/r_Co(ii)/c**2))	!trying f mode

	B(2,1) = (eLamb**2)*c1*Omega_sq +A*r_Co(ii)
!B(2,1) = c1*Omega_sq +A*r_Co(ii)
	B(2,2) = -(U1 + A*r_Co(ii))
!B(2,2) = -(r_Co(ii)/m_Co(ii)*4.d0*pi*r_Co(ii)**2*rho_Co(ii) -1.d0  + A*r_Co(ii))

	do i = 1,n
		fcn(i) = 0.d0
		do j = 1,n
			fcn(i) = fcn(i) + B(i,j)*x(j)/(1.d0 +V1)
		enddo
	enddo
	endsubroutine pes02_vx_Co

	subroutine pes02_vx_Cr(n, t, x, fcn)
	implicit none
	integer :: n
	real(8), dimension(1:n) :: x, fcn
	real(8) :: t
	real(8), dimension(1:n,1:n) :: B
	real(8) :: V1, V2, gamma_1, c1, A, U1, U2, l_, a1, a2, a3, eNu, eLamb
	integer :: i, j, ii

	ii = rk4_HS(XCr_i, dx_Cr,t)

	eNu = dexp(nu_Cr(ii))
	eLamb = pes02_eLamb(rho_Cr(ii), m_Cr(ii), r_Cr(ii))
	V1 = pes02_V1(P_Cr(ii), rho_Cr(ii), m_Cr(ii), r_Cr(ii))
	V2 = pes02_V2(P_Cr(ii), rho_Cr(ii), m_Cr(ii), r_Cr(ii))
	gamma_1 = pes02_gamma(2, P_Cr, rho_Cr, ii, 0, pg_N2*2)
	c1 = pes02_c1(P_Cr(ii), rho_Cr(ii), m_Cr(ii), eNu, r_Cr(ii))
	A = 0.d0
	U1 = pes02_U1(P_Cr(ii), rho_Cr(ii), m_Cr(ii), r_Cr(ii))
	U2 = pes02_U2(rho_Cr(ii), m_Cr(ii), r_Cr(ii))
	l_ = l_0*(l_0+1.d0)
	a1 = mu_Cr(ii)/P_Cr(ii)
	a2 = gamma_1 - 2.d0/3.d0*a1
	a3 = gamma_1 + 4.d0/3.d0*a1

	B(1,1) = -(1.d0 +2.d0*a2/a3 + U2)
	B(1,2) = 1.d0/a3
	B(1,3) = l_*a2/a3
	B(1,4) = 0.d0

	B(2,1) = (-3.d0 - U2 + U1 - eLamb**(2) *c1 *Omega_sq )*V1 + 4.d0*a1/a3*(3.d0*a2 +2.d0*a1)
	B(2,2) = V2 -4.d0*a1/a3
	B(2,3) = l_*(V1 -2.d0*a1*(1.d0 +2.d0*a2/a3))
	B(2,4) = l_*eLamb**2

	B(3,1) = -eLamb**2
	B(3,2) = 0.d0
	B(3,3) = 0.d0
	B(3,4) = eLamb**(2)/a1

	B(4,1) = -(-V1 + 6.d0*gamma_1*a1/a3)
	B(4,2) = -a2/a3
	B(4,3) = -(c1*Omega_sq*V1 + 2.d0*a1-2.d0*a1/a3*(a2+a3)*l_)
	B(4,4) = -(3.d0 +U2 -V2)

! Consider m constant
!if (r_Cr(ii) <= r_Cr(0) + (R0-r_Cr(0))*0.99d0) then
!write(*,*) "- U2 + U1", - U2 + U1
!write(*,*) "Newtonian", -1.d0 + r_Cr(ii)**3/m_Cr(ii)*4.d0 * pi * rho_Cr(ii)
!write(*,*) "predict 1", (-1.d0 + Grav_Const *m_Cr(ii)/r_Cr(ii)/c**2)/(1.d0 - 2.d0 * Grav_Const *m_Cr(ii)/r_Cr(ii)/c**2)
!write(*,*) "predict 2", (-1.d0 + Grav_Const *m_Cr(ii)/r_Cr(ii)/c**2 + 4.d0*pi*r_Cr(ii)**2*Grav_Const *rho_Cr(ii)/c**2)/(1.d0 - 2.d0 * Grav_Const *m_Cr(ii)/r_Cr(ii)/c**2)
!write(*,*) "predict 3", Grav_Const /c**2 * (-1.d0*m_Cr(ii)/r_Cr(ii) - 4.d0 * pi * r_Cr(ii)**2 * rho_Cr(ii))/(1.d0 - 2.d0 * Grav_Const *m_Cr(ii)/r_Cr(ii)/c**2) - 1.d0
!write(*,*) "m/r, rho r^3/P", m_Cr(ii)/r_Cr(ii), r_Cr(ii)**3 * rho_Cr(ii)/P_Cr(ii)
!write(*,*) "m", m_Cr(ii), m_Cr(ii) + 4.d0*pi*r_Cr(ii)**3*P_Cr(ii)/c**2
!write(*,*) "rho", -m_Cr(ii), -m_Cr(ii) + 4.d0*pi*r_Cr(ii)**3*rho_Cr(ii)
!write(*,*) "dmdr", m_Cr(ii)/r_Cr(ii), 4.d0*pi*r_Cr(ii)**2*rho_Cr(ii)
!pause
!endif

	do i = 1,n
		fcn(i) = 0.d0
		do j = 1,n
			fcn(i) = fcn(i) + B(i,j)*x(j)/(1.d0 + V1)
		enddo
	enddo
	endsubroutine pes02_vx_Cr

	subroutine pes02_vx_Oc(n, t, x, fcn)
	implicit none
	integer :: n
	real(8), dimension(1:n) :: x, fcn
	real(8) :: t
	real(8), dimension(1:n,1:n) :: B
	real(8) :: V1, V2, gamma_1, c1, A, U1, U2, l_, eNu, eLamb
	integer :: i, j, ii

	ii = rk4_HS(XOc_i, dx_Oc,t)

	eNu = dexp(nu_Oc(ii))
	eLamb = pes02_eLamb(rho_Oc(ii), m_Oc(ii), r_Oc(ii))
	V1 = pes02_V1(P_Oc(ii), rho_Oc(ii), m_Oc(ii), r_Oc(ii))
	V2 = pes02_V2(P_Oc(ii), rho_Oc(ii), m_Oc(ii), r_Oc(ii))
	gamma_1 = pes02_gamma(1, P_Oc, rho_Oc, ii, 0, pg_N3*2)
	c1 = pes02_c1(P_Oc(ii), rho_Oc(ii), m_Oc(ii), eNu, r_Oc(ii))
	A = 0.d0
	U1 = pes02_U1(P_Oc(ii), rho_Oc(ii), m_Oc(ii), r_Oc(ii))
	U2 = pes02_U2(rho_Oc(ii), m_Oc(ii), r_Oc(ii))
	l_ = l_0*(l_0+1.d0)

	B(1,1) = -(3.d0 - V1/gamma_1 + U2)
	B(1,2) = -(V1/gamma_1 - l_/c1/Omega_sq)
	B(2,1) = (eLamb**2)*c1*Omega_sq +A*r_Oc(ii)
	B(2,2) = -(U1 + A*r_Oc(ii))

	do i = 1,n
		fcn(i) = 0.d0
		do j = 1,n
			fcn(i) = fcn(i) + B(i,j)*x(j)/(1.d0 +V1)
		enddo
	enddo
	endsubroutine pes02_vx_Oc

	subroutine pes02_vx_Cr_FM(n, t, x, fcn)
	implicit none
	integer :: n
	real(8), dimension(1:n) :: x, fcn
	real(8) :: t
	real(8), dimension(1:n,1:n) :: B
	real(8) :: V1, V2, gamma_1, c1, A, U1, U2, l_, eNu, eLamb
	integer :: i, j, ii

	ii = rk4_HS(XCr_i, dx_Cr,t)

	eNu = dexp(nu_Cr(ii))
	eLamb = pes02_eLamb(rho_Cr(ii), m_Cr(ii), r_Cr(ii))
	V1 = pes02_V1(P_Cr(ii), rho_Cr(ii), m_Cr(ii), r_Cr(ii))
	V2 = pes02_V2(P_Cr(ii), rho_Cr(ii), m_Cr(ii), r_Cr(ii))
	gamma_1 = pes02_gamma(1, P_Cr, rho_Cr, ii, 0, pg_N2*2)
	c1 = pes02_c1(P_Cr(ii), rho_Cr(ii), m_Cr(ii), eNu, r_Cr(ii))
	A = 0.d0
	U1 = pes02_U1(P_Cr(ii), rho_Cr(ii), m_Cr(ii), r_Cr(ii))
	U2 = pes02_U2(rho_Cr(ii), m_Cr(ii), r_Cr(ii))
	l_ = l_0*(l_0+1.d0)

	B(1,1) = -(3.d0 - V1/gamma_1 + U2)
	B(1,2) = -(V1/gamma_1 - l_/c1/Omega_sq)
!B(1,1) = -(3.d0 - V1/gamma_1)
	B(2,1) = (eLamb**2)*c1*Omega_sq +A*r_Cr(ii)
	B(2,2) = -(U1 + A*r_Cr(ii))
!B(2,2) = -(r_Cr(ii)/m_Cr(ii)*4.d0*pi*r_Cr(ii)**2*rho_Cr(ii) -1.d0  + A*r_Cr(ii))

!if (ii <= nint(pg_N2*2 * 0.15d0)) then
!write(*,*) "U2, - V1/gamma_1", U2, - V1/gamma_1
!write(*,*) "V1/gamma_1, - l_/c1/Omega_sq", V1/gamma_1, - l_/c1/Omega_sq
!pause
!endif

!if (ii == nint(0.1d0 * 2* pg_N2)) then
!write(*,*) "B(1,1), B(1,2)", B(1,1), B(1,2)
!write(*,*) "B(2,1), B(2,2)", B(2,1), B(2,2)
!pause
!endif
	do i = 1,n
		fcn(i) = 0.d0
		do j = 1,n
			fcn(i) = fcn(i) + B(i,j)*x(j)/(1.d0 +V1)
		enddo
	enddo
	endsubroutine pes02_vx_Cr_FM

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
endmodule puls_eqt_set_opt02