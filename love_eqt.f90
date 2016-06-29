module lo_eqt
use global_var
implicit none
real(8) :: mu_factor
character(10) :: region, eos_choice

contains
	!Andersson 2011 Fluid Part eqt---------------------------------------------------
	subroutine lo_fluid_eqt_TOV_dr(n, t, x, fcn)
	! variables {H0, H0`, K}, defining metric H0
	! independent variable r
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
			rho_ = dsign(lo_EOS(dabs(P_), f = 'rho(p)'), P_)
			if (rho_ < 0.d0) pause 'err: rho_ < 0'
			m_ = x(5)
			nu_ = x(6)
			r_ = t
gamma_1 = lo_gamma(P_)
		elseif ( region == 'crust') then
			P_ = x(4)
			rho_ = dsign(lo_EOS(dabs(P_), f = 'rho(p)'), P_)
			if (rho_ < 0.d0) pause 'err: rho_ < 0'
			m_ = x(5)
			nu_ = x(6)
			r_ = t
gamma_1 = lo_gamma(P_)
		elseif ( region == 'ocean') then
			P_ = x(4)
			rho_ = dsign(lo_EOS(dabs(P_), f = 'rho(p)'), P_)
			if (rho_ < 0.d0) pause 'err: rho_ < 0'
			m_ = x(5)
			nu_ = x(6)
			r_ = t
gamma_1 = lo_gamma(P_)
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

	endsubroutine lo_fluid_eqt_TOV_dr


	!Andersson 2011 Solid Part eqt---------------------------------------------------
	subroutine lo_solid_eqt_TOV_dr(n, t, x, fcn)
	! variables {H1, K, H0, W, V, T2}
	! independent variable r
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
	rho_ = dsign(lo_EOS(dabs(P_), f = 'rho(p)'), P_)
	if (rho_ < 0.d0) pause 'err: rho_ < 0'
	m_ = x(8)
	nu_ = x(9)
	mu_ = lo_mu(P_)
	r_ = t
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
	call lo_dmu(P_, dP, dmu_)

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

	endsubroutine lo_solid_eqt_TOV_dr

!Andersson 2011 Fluid Part eqt with Fluid Displacement W---------------------------------------------------
	subroutine lo_fluid_W_eqt_TOV_dr(n, t, x, fcn)
	! variables {H0, H0`, K}, defining metric H0
	! independent variable r
	use puls_eqt_set_opt03_ex1
	implicit none
	integer :: n
	real(8), dimension(1:n) :: x, fcn
	real(8) :: t
	real(8) :: eNu, eLamb, eNu_hf, eLamb_hf, dLamb, dNu, gamma_1, P_, rho_, m_, nu_, r_, Q_
	real(8) :: pi4, l_1, l_2_hf, ddNu
	real(8) :: DV_Coef_H0_Nom, DV_Coef_H1_DeNom, DV_Coef_K_Nom
	integer :: i, j
	real(8), dimension(1:n,1:n) :: B
	
		if (n /= 8) then
			write(*,*) "err: love fluid eqt set n not 8"
			pause
		endif

		P_ = x(6)
		rho_ = dsign(lo_EOS(dabs(P_), f = 'rho(p)'), P_)
		if (rho_ < 0.d0) pause 'err: rho_ < 0'
		m_ = x(7)
		nu_ = x(8)
		r_ = t
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
		
		Q_ = (-l_1*eLamb/r_**2 + pi4*eLamb*(5.d0*rho_+9.d0*P_+(rho_+P_)/(gamma_1*P_/(rho_+P_)) ) - dnu**2)
		
		DV_Coef_H0_Nom = 4.d0*(-3.d0*m_ - l_2_hf*r_ + pi4*r_**3*rho_)
		DV_Coef_H1_DeNom = l_1 * r_/eLamb * dNu
		DV_Coef_K_Nom = 4.d0*(l_2_hf*r_ - 0.5d0*dNu*r_*(3.d0*m_-r_+pi4*r_**3*P_))

		B = 0.d0

		B(1,2) = r_
		B(2,1) = - Q_ * r_
		B(2,2) = - (2.d0/r_+ 0.5d0*(dnu-dLamb)) * r_
		B(3,1) = dnu * r_
		B(3,2) = r_

		B(4,1) = - 0.5d0*r_*((rho_+P_)/gamma_1/P_ + 1.d0) 
		B(4,3) = -r_
		B(4,4) = 0.5d0*(rho_+P_)/gamma_1/P_*dNu - 0.5d0*dLamb - 1.d0/r_
		B(4,5) = l_1/r_

		B(5,1) = DV_Coef_H0_Nom/DV_Coef_H1_DeNom
		B(5,3) = DV_Coef_K_Nom/DV_Coef_H1_DeNom
		B(5,4) = -eLamb* (16.d0*pi/l_1*(rho_ + P_)*r_ -1.d0/r_) 
		B(5,5) = dNu
		
		B(4:5, 1:5) = B(4:5, 1:5) * r_

		do i = 1,5
			fcn(i) = 0.d0
			do j = 1,5
				fcn(i) = fcn(i) + B(i,j)*x(j)/r_
			enddo
		enddo

		fcn(6) = - (rho_+P_)*dNu/2.d0
		fcn(7) = pi4*r_**2*rho_
		fcn(8) = dNu/2.d0

	endsubroutine lo_fluid_W_eqt_TOV_dr

	!Andersson 2011 Fluid Part eqt with Fluid Displacement W---------------------------------------------------
	subroutine lo_fluid_W_eqt_TOV_dr_COPY(n, t, x, fcn)
	! Backup: this set of equations is wrong, didn't account for the existence of H1
	use puls_eqt_set_opt03_ex1
	implicit none
	integer :: n
	real(8), dimension(1:n) :: x, fcn
	real(8) :: t
	real(8) :: eNu, eLamb, eNu_hf, eLamb_hf, dLamb, dNu, gamma_1, P_, rho_, m_, nu_, r_, Q_
	real(8) :: pi4, l_1, ddNu
	integer :: i, j
	real(8), dimension(1:n,1:n) :: B
	
		if (n /= 8) then
			write(*,*) "err: love fluid eqt set n not 8"
			pause
		endif
		
		if ( region == 'core') then
			P_ = x(6)
			rho_ = dsign(lo_EOS(dabs(P_), f = 'rho(p)'), P_)
			if (rho_ < 0.d0) pause 'err: rho_ < 0'
			m_ = x(7)
			nu_ = x(8)
			r_ = t
gamma_1 = lo_gamma(P_)
		elseif ( region == 'crust') then
			P_ = x(6)
			rho_ = dsign(lo_EOS(dabs(P_), f = 'rho(p)'), P_)
			if (rho_ < 0.d0) pause 'err: rho_ < 0'
			m_ = x(7)
			nu_ = x(8)
			r_ = t
gamma_1 = lo_gamma(P_)
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

		B(4,1) = - 0.5d0*r_*((rho_+P_)/gamma_1/P_ + 1.d0) 
		B(4,3) = -r_
		B(4,4) = 0.5d0*(rho_+P_)/gamma_1/P_*dNu - 0.5d0*dLamb - 1.d0/r_
		B(4,5) = l_1/r_
		B(5,4) = eLamb/r_
		B(5,5) = dNu
		
		B(4:5, 1:5) = B(4:5, 1:5) * r_

		do i = 1,5
			fcn(i) = 0.d0
			do j = 1,5
				fcn(i) = fcn(i) + B(i,j)*x(j)/r_
			enddo
		enddo

		fcn(6) = - (rho_+P_)*dNu/2.d0
		fcn(7) = pi4*r_**2*rho_
		fcn(8) = dNu/2.d0

	endsubroutine lo_fluid_W_eqt_TOV_dr_COPY

	!F1 A7 Add 1---------------------------------------------------
	subroutine lo_solid_Algebraic_Relation(P_, rho_, m_, nu_, mu_, gamma_1, r_, PP, TT)
	! CORRECT VERSION
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
	pi_8 = -16.d0*pi*eLamb               !8.d0*pi*(eLamb-3.d0)			! consistent with Andersson 2014
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
	endsubroutine lo_solid_Algebraic_Relation


	subroutine lo_solid_Algebraic_Relation_Wrong(P_, rho_, m_, nu_, mu_, gamma_1, r_, PP, TT)
	! WRONG VERSION (Andersson 2011)
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
	pi_8 = 8.d0*pi*(eLamb-3.d0)			! consistent with Andersson 2014
	Coef_ = 16.d0*pi*r_**2*eLamb - mu_8_3 * r_**2/(rho_+P_)/cs2 * pi_8

	PP = 0.d0
	PP(1) = (l_1*eLamb-2.d0+(r_*dNu)**2)
	PP(2) = r_**2 *dNu
	PP(3) = -2.d0*l_2_hf*eLamb + pi_8 * mu_8_3 * 1.5d0*r_**2
	PP(4) = pi_8*mu_8_3* (3.d0-r_*dNu/2.d0/cs2)
	PP(5) = 32.d0*pi*mu_* ((r_*dNu)**2 + 2.d0*l_2_hf*(1.d0-eLamb)) + pi_8*mu_8_3 * (-1.5d0*l_1)
	PP(6) = -8.d0*pi*(r_*dNu + 2.d0)
	PP(1:6) = PP(1:6) / Coef_

	TT = 0.d0
	TT(3) = 1.5d0*r_**2
	TT(4) = 3.d0-r_*dNu/2.d0/cs2
	TT(5) = -1.5d0*l_1
	TT(1:6) = TT(1:6) + r_**2/(rho_+P_)/cs2*PP(1:6)
	TT(1:6) = TT(1:6) * mu_8_3
	endsubroutine lo_solid_Algebraic_Relation_Wrong
    !EOS
	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	function lo_EOS(x_, f)
	! in rho or P in km unit, return P or rho in km unit
	use EOS
	implicit none
	real(8) :: lo_EOS, x_
	character(len=*), optional :: f
		if (eos_choice == 'NS_A11') then
			lo_EOS = lo_A11_NS_EOS(x_, f)
		elseif (eos_choice == 'QS') then 
			lo_EOS = lo_qs_EOS(x_, f)
		else
			pause 'err: love eos'
		endif
	endfunction lo_EOS
	
	function lo_gamma(P_)
	! input P in km unit and output gamma
	implicit none
	real(8) :: lo_gamma, P_
		if (eos_choice == 'NS_A11') lo_gamma = lo_A11_NS_gamma(P_)
		if (eos_choice == 'QS') lo_gamma = lo_qs_gamma(P_)
	endfunction lo_gamma

	function lo_mu(P_)
	implicit none
	real(8) :: lo_mu, P_, y_, Conv_P
		if (eos_choice == 'NS_A11') then
			y_ = lo_A11_NS_mu(P_)
		elseif (eos_choice == 'QS') then
		Conv_P = Grav_Const/c**4
			y_ = lo_qs_mu(P_/Conv_P)
			y_ = y_ * Conv_P
		endif
		if (region == 'core') then
			lo_mu = y_  /2.d0 * mu_factor  ! divide by 2 to cancel out the typo in Andersson (2010)
		elseif (region == 'crust') then
			lo_mu = y_  /2.d0 * mu_factor
		else
			lo_mu = y_  /2.d0 * mu_factor 
		endif
	endfunction lo_mu

	subroutine lo_dmu(P_, dP_, dmu_)
	implicit none
	real(8) :: dmu_, P_, dP_, Conv_P, y_
		if (eos_choice == 'NS_A11') then
			y_ = lo_A11_NS_dmu(dP_)
			y_ = y_
		elseif (eos_choice == 'QS') then
			Conv_P = Grav_Const/c**4
			y_ = lo_qs_dmu(P_/Conv_P, dP_/Conv_P)
			y_ = y_ * Conv_P  ! divide by 2 to cancel out the typo in Andersson (2010)
		endif
		if (region == 'core') then
			dmu_ = y_  /2.d0 * mu_factor ! divide by 2 to cancel out the typo in Andersson (2010)
		elseif (region == 'crust') then
			dmu_ = y_  /2.d0 * mu_factor
		else
			dmu_ = y_  /2.d0 * mu_factor
		endif

	endsubroutine lo_dmu
	! QS EOS
	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	function lo_qs_EOS(x_, f)
	! input rho or P in km unit, return P or rho in km unit
	use EOS
	implicit none
	real(8) :: lo_qs_EOS, x_, y_, Conv_rho, Conv_P
	character(len=*), optional :: f
		Conv_rho = Grav_Const/c**2
		Conv_P = Grav_Const/c**4
		if (present(f) /= .true.) then
			y_ = sp_qs_rho_P(x_/Conv_P)
			y_ = y_ * Conv_rho
		elseif (f == 'rho(p)') then
			y_ = sp_qs_rho_P(x_/Conv_P)
			y_ = y_ * Conv_rho
		elseif (f == 'p(rho)') then
			y_ = sp_qs_P_rho(x_/Conv_rho)
			y_ = y_ * Conv_P
		endif
		lo_qs_EOS = y_
	endfunction lo_qs_EOS
	
	function lo_qs_gamma(x_)
	! input P in km unit and output gamma
	implicit none
	real(8) :: lo_qs_gamma, P_, x_, y_
	real(8) :: C1_, A_, B_, C_, mu_2, Conv_P, Conv_MeV_P	! mu_ and mu_2 are chemical potential
	real(8), parameter :: e_mks = 1.60217657d-19
		!¡¸Quark Star EOS with ref to Pagliaroli 2014; 
		!¡¸Convert to natural units: h_bar = 1, c = 1, P in dimension [E/L^3] = [E^4/((h_bar*c)^3)]
		Conv_MeV_P = ((1.d13)*e_mks)**4 /(h_bar*c)**3
		Conv_P = Grav_Const/c**4
		C1_ = 3.d0/4.d0/pi**(2.d0)
		A_ = C1_ * QS_a4
		B_ = - C1_ * QS_a2
		C_ = -B_eff
		P_ = x_ / Conv_P / Conv_MeV_P
		mu_2 = -B_/2.d0/A_ + dsqrt(B_**2 - 4.d0*A_*(C_-P_))/2.d0/A_
		y_ = mu_2 * (4.d0*A_*mu_2 +2.d0*B_)**2 / (12.d0*A_*mu_2 + 2.d0*B_) / (A_*mu_2**2 + B_ * mu_2 + C_)
		lo_qs_gamma = y_
	endfunction lo_qs_gamma
	! QS Shear Modulus
	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	function lo_qs_mu(x_)
	! input P in cgs unit, returns mu in cgs unit
	implicit none
	real(8) :: lo_qs_mu, x_, P_, y_
	real(8) :: C1_, A_, B_, C_, mu_2, mu_, rho_, Conv_MeV_P	! mu_ and mu_2 are chemical potential
	real(8), parameter :: e_mks = 1.60217657d-19
		!¡¸Quark Star EOS with ref to Pagliaroli 2014; 
		!¡¸Convert to natural units: h_bar = 1, c = 1, P in dimension [E/L^3] = [E^4/((h_bar*c)^3)]
		Conv_MeV_P = ((1.d13)*e_mks)**4 /(h_bar*c)**3
		C1_ = 3.d0/4.d0/pi**(2.d0)
		A_ = C1_ * QS_a4
		B_ = - C1_ * QS_a2
		C_ = -B_eff
		P_ = x_ / Conv_MeV_P
			mu_2 = -B_/2.d0/A_ + dsqrt(B_**2 - 4.d0*A_*(C_-P_))/2.d0/A_
		y_ = QS_mu0 *(QS_Gap/10.d0)**2 * mu_2/(400.d0)**2
		y_ = y_ * ((1.d13)*e_mks) * (1.d39)	! convert unit from MeV/fm^3 to cgs
		lo_qs_mu = y_
	endfunction lo_qs_mu
	! QS Shear Modulus Derivative
	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	function lo_qs_dmu(x_, dx_)
	! input P and dP in cgs unit, return dmu in cgs unit
	implicit none
	real(8) :: lo_qs_dmu, P_, dP_, x_, dx_, y_
	real(8) :: C1_, A_, B_, C_, Conv_MeV_P	! mu_ and mu_2 are chemical potential
	real(8), parameter :: e_mks = 1.60217657d-19
		!¡¸Quark Star EOS with ref to Pagliaroli 2014; 
		!¡¸Convert to natural units: h_bar = 1, c = 1, P in dimension [E/L^3] = [E^4/((h_bar*c)^3)]
		Conv_MeV_P = ((1.d13)*e_mks)**4 /(h_bar*c)**3
		C1_ = 3.d0/4.d0/pi**(2.d0)
		A_ = C1_ * QS_a4
		B_ = - C1_ * QS_a2
		C_ = -B_eff
			P_ = x_ / Conv_MeV_P
			dP_ = dx_ / Conv_MeV_P
		y_ = dP_ / dsqrt(B_**2 - 4.d0*A_*(C_-P_))
		y_ = y_ * QS_mu0 *(QS_Gap/10.d0)**2 /(400.d0)**2 
		y_ = y_ * ((1.d13)*e_mks) * (1.d39)	! convert unit from MeV/fm^3 to cgs
		lo_qs_dmu = y_
	endfunction lo_qs_dmu
	


	!Andersson 2011 NS EOS
	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	function lo_A11_NS_EOS(x_, f)
	! in rho or P in km unit, return P or rho in km unit
	use global_var
	implicit none
	real(8) :: lo_A11_NS_EOS, x_, y_, Conv_rho, Conv_P
	character(len=*), optional :: f
		if (present(f) /= .true.) then
			y_ = ((x_)/poly_K/1.d10)**(poly_n/(poly_n + 1.d0))
		elseif (f == 'rho(p)') then
			y_ = ((x_)/poly_K/1.d10)**(poly_n/(poly_n + 1.d0))
		elseif (f == 'p(rho)') then
			y_ = ((x_))**((poly_n + 1.d0)/poly_n) *poly_K*1.d10
		endif
		lo_A11_NS_EOS = y_
	endfunction lo_A11_NS_EOS

	function lo_A11_NS_gamma(x_)
	implicit none
	real(8) :: lo_A11_NS_gamma, x_ , P_, rho_
		P_ = x_
		rho_ =  lo_A11_NS_EOS(x_, f='rho(p)')
		lo_A11_NS_gamma = 2.d0 * (1.d0 + P_/rho_)
	endfunction lo_A11_NS_gamma
	!Andersson 2011 NS mu
	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	function lo_A11_NS_mu(x_)
	! input P in cgs unit, returns mu in cgs unit
	implicit none
	real(8) :: lo_A11_NS_mu, x_, P_, y_
	real(8) :: C1_, A_, B_, C_, mu_2, mu_, rho_, Conv_MeV_P	! mu_ and mu_2 are chemical potential
	real(8), parameter :: e_mks = 1.60217657d-19
		y_ = (Ki_poly*x_ + mu_poly/1.d10)
		lo_A11_NS_mu = y_
	endfunction lo_A11_NS_mu
	!Andersson 2011 NS dmu
	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	function lo_A11_NS_dmu(dP)
	implicit none
	real(8) :: lo_A11_NS_dmu, dP
		lo_A11_NS_dmu = Ki_poly * dP
	endfunction lo_A11_NS_dmu
	
endmodule lo_eqt