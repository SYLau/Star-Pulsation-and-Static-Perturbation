module love_opt01_P06
! Coupled with TOV equation
! Use the background grid to reach lowest density
! To compare with Andersson 2011, with Polytropic EOS and analytic expression of Shear modulus
! Modified the algebraic expressions of Andersson 2011 (Since it is inconsistent with Andersson 2011 and Finn 1990)
! Implement Adaptive Mesh
use global_var
implicit none
character(10) :: region
real(8) :: mu_factor
real(8), parameter :: acc = 1.d-6

contains

	!F1 A1---------------------------------------------------
	subroutine lo01_iterate_fluid
	use Format_IO
	use FWrite
	use RK4_Set
	use lo_match
	implicit none
	real(8), dimension(1:6) :: x
	real(8) :: r_, dr_, k2_
	integer :: i, nok, nbad

write(*,*) "Multiply the shear modulus by mu_factor"
mu_factor = 1.d0
write(*,*) "mu_factor = "
read(*,*) mu_factor
!mu_factor = 1.d-3
		open(01, file = pef_LD_Co1, status = 'replace')
		open(02, file = 'data\sol_stat profile.txt', status = 'replace')
		call lo01_bc_Coi(x(1:6))

		region = 'core'
		do i = 1, sp_N1-1
				r_ = sp_r(i)
				dr_ = delr(i)
call odeint(x(1:6),6,r_,r_ + dr_,acc, dr_, 0.d0 ,nok,nbad,lo01_fluid_eqt, rkqs)
				call Write1R81AR8(01, FR8,r_ + dr_, (x(1:3)), m=3)
				call Write4R8(02, FR8,r_ + dr_, x(4)/(Grav_Const/c**4), x(5)/(Grav_Const/c**2), x(6))
		enddo
		close(01)
write(*,*) sp_r(sp_N1-1), sp_r(sp_N1)
pause
		open(01, file = pef_LD_Cr1, status = 'replace')

		region = 'crust'
		call Write1R81AR8(01, FR8,r_ + dr_, (x(1:3)), m=3)
		do i = sp_N1, sp_N2-1
				r_ = sp_r(i)
				dr_ = delr(i)
call odeint(x(1:6),6,r_,r_ + dr_,acc, dr_, 0.d0 ,nok,nbad,lo01_fluid_eqt, rkqs)
				call Write1R81AR8(01, FR8,r_ + dr_, (x(1:3)), m=3)
				call Write4R8(02, FR8,r_ + dr_, x(4)/(Grav_Const/c**4), x(5)/(Grav_Const/c**2), x(6))
		enddo

call lo_k2(sp_r(sp_N2)*x(2)/x(1), x(5), sp_r(sp_N2-1), k2_)
		write(*,*) "k2", k2_
		write(*,*) "=========================================================="
		close(01)
		close(02)
pause
	endsubroutine lo01_iterate_fluid

	!F1 A2---------------------------------------------------
	subroutine lo01_iterate_FE
	! FE - Fluid Core Elastic Crust
	use Format_IO
	use FWrite
	use RK4_Set
	use lo_match
	implicit none
	real(8), dimension(1:9) :: x
	real(8), allocatable :: r_out(:)
	real(8), allocatable:: z_out(:,:,:)
	real(8) :: Coef(1:5)
	real(8) :: H0_Co, K_Co
	real(8), dimension(1:6, 1:5) :: cY
	real(8) :: r_, r0_, dr_, Pc_, mc_, nuc_, k2_
	integer :: i, it, mode , nok, nbad
	
	allocate(r_out(sp_N1:sp_N2), z_out(1:5, 1:6, sp_N1:sp_N2))

		open(01, file = pef_LD_Co1, status = 'replace')
		cY = 0.d0
		call lo01_bc_Coi(x(1:6))
		region = 'core'
		do i = 1, sp_N1-1
				r_ = sp_r(i)
				dr_ = delr(i)
			call odeint(x(1:6),6,r_,r_ + dr_,acc, dr_, 0.d0 ,nok,nbad,lo01_fluid_eqt, rkqs)
			call Write1R81AR8(01, FR8,r_ + dr_, (x(1:3)), m=3)
		enddo
		close(01)

		H0_Co = x(1)
		K_Co = x(3)

		Pc_ = x(4)
		mc_ = x(5)
		nuc_ = x(6)
		r0_ = sp_r(sp_N1)

! Transfer to lo01_bc_Solve
Cri(1) = Pc_
Cri(2) = dsign((dabs(Pc_)/poly_K/1.d10)**(poly_n/(poly_n + 1.d0)), Pc_)
if (Cri(2) < 0.d0) pause 'err: rho_ < 0'
Cri(3) = mc_
Cri(4) = (Ki_poly*Pc_ + mu_poly/1.d10)/2.d0  *mu_factor
Cri(5) = nuc_
Cri(6) = r0_

		region = 'crust'
		do it = 2, 6
			mode = it - 1
			call lo01_bc_Cri(mode, x(1:6))
			x(7) = Pc_
			x(8) = mc_
			x(9) = nuc_
			z_out(mode, 1:6, sp_N1) = x(1:6)
			r_out(sp_N1) = r0_

			do i = sp_N1, sp_N2-1
				r_ = sp_r(i)
				dr_ = delr(i)
				call odeint(x(1:9),9,r_,r_ + dr_,acc, dr_, 0.d0 ,nok,nbad,lo01_solid_eqt, rkqs)
				r_out(i+1) = r_ + dr_
				z_out(mode, 1:6, i+1) = x(1:6)
			enddo
			cY(1:6,it - 1) = x(1:6)
		enddo

! Transfer to lo01_bc_Solve
Crf(1) = x(7)
Crf(2) = dsign((dabs(x(7))/poly_K/1.d10)**(poly_n/(poly_n + 1.d0)), x(7))
if (Crf(2) < 0.d0) pause 'err: rho_ < 0'
Crf(3) = x(8)
Crf(4) = (Ki_poly*x(7) + mu_poly/1.d10)/2.d0 *mu_factor
Crf(5) = x(9)
Crf(6) = sp_r(sp_N2)
		call lo_bc_Solve(H0_Co, K_Co, cY(1:6, 1:5), k2_, Coef(1:5))


		open(01, file = pef_LD_Cr1, status = 'replace')
		do i = sp_N1, sp_N2-1
			call Write1R81AR8(01, FR8,r_out(i), Coef(1)*z_out(1,1:6, i)+Coef(2)*z_out(2,1:6, i)+Coef(3)*z_out(3,1:6, i)+Coef(4)*z_out(4,1:6, i)+Coef(5)*z_out(5,1:6, i), m=6)
		enddo
		close(01)

		write(*,*) "k2", k2_
		write(*,*) "=========================================================="
		deallocate(r_out, z_out)

	endsubroutine lo01_iterate_FE

	!F1 A3---------------------------------------------------
	subroutine lo01_bc_Coi(y)
	use EOS
	use puls_eqt_set_opt03_ex1
	implicit none
	real(8) :: gamma_1, P_, rho_, m_, nu_, r_, dr_, a0_, a2_
	real(8) :: y(1:6)
		P_ = P(1)
		rho_ = rho(1)
		m_ = m(1)
		nu_ = nu(1)
		r_ = sp_r(1)
gamma_1 = 2.d0 * (1.d0 + P_/rho_)
		a0_ = 1.d0
		a2_ = -2.d0/(2.d0*l_0+3.d0) * pi * (5.d0*rho_+9.d0*P_+(rho_+P_)/(gamma_1*P_/(rho_+P_)) )

		y(1) = a0_ * r_**l_0 * (1.d0 + r_ ** 2 *  a2_ )
		y(2) = a0_ * r_**(l_0-1.d0) * (l_0 + r_ ** 2 * (l_0+2.d0) * a2_  )
		y(3) = a0_ * r_**l_0  * ( 1.d0 +  ((l_0+2.d0)*a2_ + 8.d0*pi/3.d0*(3.d0*P_+rho_))/(l_0+2.d0) )   ! need second order expansion !!

		y(4) = P_
		y(5) = m_
		y(6) = nu_
	endsubroutine lo01_bc_Coi

	!F1 A4---------------------------------------------------
	subroutine lo01_bc_Cri(mode, z)
	use puls_eqt_set_opt03_ex1
	!	2 components
	!	Quantity Q = {QQ(i)} [Y] {c(i)} ------ {vec} [2nd Tensor] {vec} inner product
	implicit none
	real(8) :: z(1:6)
	integer :: mode
		
		z = 0.d0
		if (mode == 1) z(1) = 1.d0 !/r_Cr(0)**2
		if (mode == 2) z(2) = 1.d0 !/r_Cr(0)**2
		if (mode == 3) z(3) = 1.d0 !/r_Cr(0)**2
		if (mode == 4) z(4) = 1.d0
		if (mode == 5) z(5) = 1.d0
z(5) = 1.d0
! Gives better result for Adaptive mesh: first point at CC Interface dy(5)dx = y(5) = 0, causing the error scale in odeint = TINY = 1.d-30 giving huge error
	endsubroutine lo01_bc_Cri


	!F1 A5---------------------------------------------------
	subroutine lo01_bc_Ocf(y)
	use puls_eqt_set_opt03_ex1
	implicit none
	real(8) :: gamma_1, P_, rho_, m_, nu_, r_
	real(8) :: y(1:3)

	endsubroutine lo01_bc_Ocf

	!F1 A6---------------------------------------------------
	subroutine lo01_fluid_eqt(n, t, x, fcn)
	! variables {H0, H0`, K}, defining metric H0
	use puls_eqt_set_opt03_ex1
	implicit none
	integer :: n
	real(8), dimension(1:n) :: x, fcn
	real(8) :: t
	real(8) :: eNu, eLamb, eNu_hf, eLamb_hf, dLamb, dNu, gamma_1, P_, rho_, m_, nu_, r_, Q_
	real(8) :: pi4, l_1, ddNu
	integer :: i, j
	real(8), dimension(1:n,1:n) :: B
	
		if (n /= 6) then
			write(*,*) "err: love fluid eqt set n not 6"
			pause
		endif
		
		if ( region == 'core') then
			P_ = x(4)
			rho_ = dsign((dabs(P_)/poly_K/1.d10)**(poly_n/(poly_n + 1.d0)), P_)
			if (rho_ < 0.d0) pause 'err: rho_ < 0'
			m_ = x(5)
			nu_ = x(6)
			r_ = t
gamma_1 = 2.d0 * (1.d0 + P_/rho_)
		elseif ( region == 'crust') then
			P_ = x(4)
			rho_ = dsign((P_/poly_K/1.d10)**(poly_n/(poly_n + 1.d0)), P_)
			if (rho_ < 0.d0) pause 'err: rho_ < 0'
			m_ = x(5)
			nu_ = x(6)
			r_ = t
gamma_1 = 2.d0 * (1.d0 + P_/rho_)
		endif

		pi4 = 4.d0*pi
		l_1 = l_0 * (l_0 + 1.d0)
		eLamb = pes03_eLamb(rho_, m_, r_) ** 2
		eNu = dexp(2.d0*nu_)
		eLamb_hf = eLamb**(0.5d0)
		eNu_hf = eNu**(0.5d0)
		dLamb = pes03_dLamb(rho_, m_, r_) * 2.d0
		dNu = pes03_dnu(P_, rho_, m_, r_) * 2.d0
		ddNu = pes03_ddnu(P_, rho_, m_, r_) * 2.d0
		
		Q_ = (-l_1*eLamb/r_**2 + pi4*eLamb*(5.d0*rho_+9.d0*P_+(rho_+P_)/(gamma_1*P_/(rho_+P_)) ) - dnu**2)

		B = 0.d0

		B(1,2) = r_
		B(2,1) = - Q_ * r_
		B(2,2) = - (2.d0/r_+ 0.5d0*(dnu-dLamb)) * r_
		B(3,1) = dnu * r_
		B(3,2) = r_

		do i = 1,3
			fcn(i) = 0.d0
			do j = 1,3
				fcn(i) = fcn(i) + B(i,j)*x(j)/r_
			enddo
		enddo

		fcn(4) = - (rho_+P_)*dNu/2.d0
		fcn(5) = pi4*r_**2*rho_
		fcn(6) = dNu/2.d0

	endsubroutine lo01_fluid_eqt

	!F1 A7---------------------------------------------------
	subroutine lo01_solid_eqt(n, t, x, fcn)
	! variables {H1, K, H0, W, V, T2}
	use puls_eqt_set_opt03_ex1
	implicit none
	integer :: n
	real(8), dimension(1:n) :: x, fcn
	real(8) :: t
	real(8), dimension(1:n,1:n) :: B
	real(8), dimension(1:6) :: PP, TT
	real(8) :: eNu, eLamb, eNu_hf, eLamb_hf, dLamb, dNu, gamma_1, P_, rho_, m_, nu_, mu_, r_
	real(8) :: pi4, l_1, l_2_hf, cs2, ddNu, dP, dmu_
	integer :: i, j

	if (n /= 9) then
		write(*,*) "err: love solid eqt set n not 9"
		pause
	endif

	P_ = x(7)
	rho_ = dsign((dabs(P_)/poly_K/1.d10)**(poly_n/(poly_n + 1.d0)), P_)
	if (rho_ < 0.d0) pause 'err: rho_ < 0'
	m_ = x(8)
	nu_ = x(9)
	mu_ = (Ki_poly*P_ + mu_poly/1.d10)/2.d0  *mu_factor
	r_ = t
gamma_1 = 2.d0 * (1.d0 + P_/rho_)

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

	call lo01_solid_Algebraic_Relation(P_, rho_, m_, nu_, mu_, gamma_1, r_, PP, TT)
	call lo01_dmu(dP, dmu_)

	B = 0.d0

	B(1,2) = 1.d0

	B(2,1) = 1.d0/r_**2 * (l_1*eLamb + 2.d0*(eLamb-1.d0) - r_*(dLamb+3.d0*dNu) + (r_*dNu)**2)
	B(2,2) = 1.d0/r_ * (0.5d0*r_*(dLamb-dNu) -2.d0)
	B(2,5) = (-128.d0*pi*mu_/r_**2) * (1.d0-eLamb+r_*(dNu+0.5d0*dLamb) - 0.25d0*(r_*dNu)**2) + (-32.d0*pi*dNu)*dmu_
	B(2,6) = (-8.d0*pi*dNu/r_)
	B(2,1:6) = B(2,1:6) + (-8.d0*pi*eLamb)*(3.d0+1.d0/cs2)*PP(1:6)

	B(3,1) = dNu
	B(3,2) = 1.d0
	B(3,5) = 32.d0*pi*mu_/r_*(r_*dNu + 2.d0)
	B(3,6) = -8.d0*pi/r_

	B(4,1) = -r_/2.d0
	B(4,3) = r_/2.d0
	B(4,4) = 2.d0/r_ - dLamb/2.d0
	B(4,5) = -(32.d0*pi*mu_*r_ + l_1/2.d0/r_)
	B(4,1:6) = B(4,1:6) + (-3.d0/8.d0/mu_/r_)*TT(1:6)

	B(5,4) = -eLamb/r_
	B(5,5) = 2.d0/r_
	B(5,6) = -1.d0/4.d0/mu_/r_

	B(6,1) = 1.d0/8.d0/pi*(dNu+dLamb)
	B(6,5) = 4.d0*mu_/r_*eLamb*(-2.d0*l_2_hf)
	B(6,6) = -(0.5d0*(dNu-dLamb) + 1.d0/r_)
	B(6,1:6) = B(6,1:6) + (-2.d0*eLamb*r_)*PP(1:6) + eLamb/r_*TT(1:6)

	B(2,1:6) = B(2,1:6) + (-32.d0*pi*dNu)*mu_*B(5,1:6)
	B = B * r_
	do i = 1,6
		fcn(i) = 0.d0
		do j = 1,6
			fcn(i) = fcn(i) + B(i,j)*x(j)/r_
		enddo
	enddo
	fcn(7) = - (rho_+P_)*dNu/2.d0
	fcn(8) = pi4*r_**2*rho_
	fcn(9) = dNu/2.d0

	endsubroutine lo01_solid_eqt


	!F1 A7 Add 1---------------------------------------------------
	subroutine lo01_solid_Algebraic_Relation(P_, rho_, m_, nu_, mu_, gamma_1, r_, PP, TT)
	use puls_eqt_set_opt03_ex1
	implicit none
	real(8), dimension(1:6) :: PP, TT
	real(8) :: eNu, eLamb, eNu_hf, eLamb_hf, dLamb, dNu, gamma_1, P_, rho_, m_, nu_, mu_, r_
	real(8) :: pi4, l_1, l_2_hf, cs2, ddNu, dP
	real(8) :: Coef_, mu_8_3, pi_8
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
	
	mu_8_3 = 8.d0*mu_/3.d0
	pi_8 =  -16.d0*pi*eLamb               !8.d0*pi*(eLamb-3.d0)			! consistent with Andersson 2014
	Coef_ = 16.d0*pi*r_**2*eLamb - mu_8_3 * r_**2/(rho_+P_)/cs2 * pi_8

	PP = 0.d0
	PP(1) = (l_1*eLamb-2.d0+(r_*dNu)**2)
	PP(2) = r_**2 *dNu
	PP(3) = -2.d0*l_2_hf*eLamb + pi_8 * mu_8_3 * 1.5d0*r_**2
	PP(4) = pi_8*mu_8_3* (3.d0-r_*dNu/2.d0/cs2)
	PP(5) = 32.d0*pi*mu_* (r_*dNu)**2 + pi_8*mu_8_3 * (-1.5d0*l_1) !32.d0*pi*mu_* ((r_*dNu)**2 + 2.d0*l_2_hf*(1.d0-eLamb)) + pi_8*mu_8_3 * (-1.5d0*l_1)
	PP(6) = -8.d0*pi*(r_*dNu + 2.d0)
	PP(1:6) = PP(1:6) / Coef_

	TT = 0.d0
	TT(3) = 1.5d0*r_**2
	TT(4) = 3.d0-r_*dNu/2.d0/cs2
	TT(5) = -1.5d0*l_1
	TT(1:6) = TT(1:6) + r_**2/(rho_+P_)/cs2*PP(1:6)
	TT(1:6) = TT(1:6) * mu_8_3
	endsubroutine lo01_solid_Algebraic_Relation

	!F1 A7 Add 2---------------------------------------------------
	subroutine lo01_dmu(dP, dmu_)
	use puls_eqt_set_opt03_ex1
	implicit none
	real(8) :: dmu_, dP

	dmu_ = Ki_poly * dP
dmu_ = dmu_/2.d0   * mu_factor
	endsubroutine lo01_dmu

	!F1 B1---------------------------------------------------
	function rk4_HS(xi, dx,x)
	implicit none
	real(8) :: xi, dx, x
	real(8) :: dx_hf
	integer :: rk4_HS
		dx_hf = dx/2.d0  ! dx/2; half step size
		rk4_HS = nint((x-xi)/dx_hf)
	endfunction rk4_HS

endmodule love_opt01_P06