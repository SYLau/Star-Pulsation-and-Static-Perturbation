module love_opt01_P01
use global_var
implicit none
character(10) :: region

contains

	!F1 A1---------------------------------------------------
	subroutine lo01_iterate_fluid
	use Format_IO
	use FWrite
	use RK4_Set
	use lo_match
	implicit none
	real(8), dimension(1:3) :: x, y, fcn(1:3)
	real(8) :: k2_
	integer :: i, ii
write(*,*) "##Love gamma_1 analytic for accuracy issue. To compare with Andersson 2011"
write(*,*) "Changes the value of first point in Cr to be the same as last point in Co"
P_Cr(0) = P_Co(2*pg_N1)
rho_Cr(0) = rho_Co(2*pg_N1)
m_Cr(0) = m_Co(2*pg_N1)
mu_Cr(0) = mu_Co(2*pg_N1)
nu_Cr(0) = nu_Co(2*pg_N1)
r_Cr(0) = r_Co(2*pg_N1)

		open(01, file = pef_LD_Co1, status = 'replace')
		call lo01_bc_Coi(y(1:3))
		region = 'core'
		do i = 0, pg_N1-1
				ii = i*2
				x(1:3) = y(1:3)
				call RK4(lo01_fluid_eqt, 3, x_Co(ii), x_Co(ii+2), x(1:3), y(1:3))
				call Write1R81AR8(01, FR8,r_Co(ii), (y(1:3)), m=3)
		enddo
		close(01)

		open(01, file = pef_LD_Cr1, status = 'replace')

		region = 'crust'
		do i = 0, pg_N2-1
				ii = i*2
				x(1:3) = y(1:3)
				call RK4(lo01_fluid_eqt, 3, x_Cr(ii), x_Cr(ii+2), x(1:3), y(1:3))
				call Write1R81AR8(01, FR8,r_Cr(ii), (y(1:3)), m=3)
		enddo

call lo_k2(r_Cr(pg_N2*2)*y(2)/y(1), 0.d0, P_Cr(pg_N2*2), rho_Cr(pg_N2*2), m_Cr(pg_N2*2), 0.d0, r_Cr(pg_N2*2), k2_)
		write(*,*) "k2", k2_
		write(*,*) "=========================================================="
		close(01)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Euler's method
!call lo01_bc_Coi(y(1:3))
!region = 'core'
!do i = 0, pg_N1-1
!		ii = 2*i
!		x(1:3) = y(1:3)
!		call lo01_fluid_eqt(3, x_Co(ii), x(1:3), fcn(1:3))
!		y(1:3) = x(1:3) + fcn(1:3)*(x_Co(ii+2)- x_Co(ii))
!enddo
!region = 'crust'
!do i = 0, pg_N2-1
!		ii = 2*i
!		x(1:3) = y(1:3)
!		call lo01_fluid_eqt(3, x_Cr(ii), x(1:3), fcn(1:3))
!		y(1:3) = x(1:3) + fcn(1:3)*(x_Cr(ii+2)- x_Cr(ii))
!enddo
!	call lo01_k2(r_Cr(pg_N2*2)*y(2)/y(1), m_Cr(pg_N2*2), r_Cr(pg_N2*2), k2_)
!write(*,*) "Euler k2", k2_
!write(*,*) "=========================================================="
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call lo01_bc_Coi2(y(1))
		region = 'core'
		do i = 0, pg_N1-1
				ii = i*2
				x(1) = y(1)
				call RK4(lo01_fluid_eqt2, 1, x_Co(ii), x_Co(ii+2), x(1), y(1))
		enddo
		region = 'crust'
		do i = 0, pg_N2-1
				ii = i*2
				x(1) = y(1)
				call RK4(lo01_fluid_eqt2, 1, x_Cr(ii), x_Cr(ii+2), x(1), y(1))
		enddo

call lo_k2(y(1), 0.d0, P_Cr(pg_N2*2), rho_Cr(pg_N2*2), m_Cr(pg_N2*2), 0.d0, r_Cr(pg_N2*2), k2_)
		write(*,*) "k2, non linear", k2_
		write(*,*) "=========================================================="
!		write(*,*) "Surface Jump", rho_Cr(2*pg_N2)/rho_Co(0), 4.d0*pi*rho_Cr(2*pg_N2)*r_Cr(2*pg_N2)**3/m_Cr(2*pg_N2), y(1)


		call lo01_bc_Coi3(y(1:3))
		y(1:3) = y(1:3)
		region = 'core'
		do i = 0, pg_N1-1
				ii = i*2
				x(1:3) = y(1:3)
				call RK4(lo01_fluid_eqt3, 3, x_Co(ii), x_Co(ii+2), x(1:3), y(1:3))
		enddo

		region = 'crust'
		do i = 0, pg_N2-1
				ii = i*2
				x(1:3) = y(1:3)
				call RK4(lo01_fluid_eqt3, 3, x_Cr(ii), x_Cr(ii+2), x(1:3), y(1:3))
		enddo
call lo_k2(r_Cr(pg_N2*2)*y(2)/y(1) + l_0, 0.d0, P_Cr(pg_N2*2), rho_Cr(pg_N2*2), m_Cr(pg_N2*2), 0.d0, r_Cr(pg_N2*2), k2_)
		!call lo_k2(r_Cr(pg_N2*2)*y(2)/y(1) + l_0, m_Cr(pg_N2*2), r_Cr(pg_N2*2), k2_)
		write(*,*) "k2, LD definition", k2_
		write(*,*) "=========================================================="

	endsubroutine lo01_iterate_fluid
	
	!F1 A2---------------------------------------------------
	subroutine lo01_iterate_FE
	! FE - Fluid Core Elastic Crust
	use Format_IO
	use FWrite
	use RK4_Set
	implicit none
	real(8), dimension(1:6) :: x, y, z
	real(8), dimension(0:pg_N2) :: r_out
	real(8), dimension(1:5, 1:6, 0:pg_N2):: z_out
	real(8) :: Coef(1:5)
	real(8) :: H0_Co, K_Co
	real(8), dimension(1:6, 1:5) :: cY
	real(8) :: k2_
	integer :: i, ii, it, mode 
		
		open(01, file = pef_LD_Co1, status = 'replace')
		cY = 0.d0
		call lo01_bc_Coi(y(1:3))
		region = 'core'
		do i = 0, pg_N1-1
				ii = i*2
				x(1:3) = y(1:3)
				call RK4(lo01_fluid_eqt, 3, x_Co(ii), x_Co(ii+2), x(1:3), y(1:3))
				call Write1R81AR8(01, FR8,r_Co(ii), (y(1:3)), m=3)
		enddo
		close(01)


		H0_Co = y(1)
		K_Co = y(3)
!mu_Cr = mu_Cr*1.d3
!mu_Cr = mu_Cr*2.d0
!mu_Cr = mu_Cr*1.d-2

		region = 'crust'
		do it = 2, 6
			mode = it - 1
			call lo01_bc_Cri(mode, z(1:6))

			do i = 0, pg_N2-1
				ii = i*2
				x(1:6) = z(1:6)
				call RK4(lo01_solid_eqt, 6, x_Cr(ii), x_Cr(ii+2), x(1:6), z(1:6))
				r_out(i) = r_Cr(ii)
				z_out(mode, 1:6, i) = z(1:6)
			enddo
			cY(1:6,it - 1) = z(1:6)
		enddo
		

		call lo01_bc_Solve(H0_Co, K_Co, cY(1:6, 1:5), k2_, Coef(1:5))


		open(01, file = pef_LD_Cr1, status = 'replace')
		do i = 0, pg_N2-1
			call Write1R81AR8(01, FR8,r_out(i), Coef(1)*z_out(1,1:6, i)+Coef(2)*z_out(2,1:6, i)+Coef(3)*z_out(3,1:6, i)+Coef(4)*z_out(4,1:6, i)+Coef(5)*z_out(5,1:6, i), m=6)
		enddo
		close(01)

		write(*,*) "k2", k2_
		write(*,*) "=========================================================="

	endsubroutine lo01_iterate_FE


!F1 A2---------------------------------------------------
	subroutine lo01_iterate_FEF
	! FE - Fluid Core Elastic Crust
	use Format_IO
	use FWrite
	use RK4_Set
	implicit none
	real(8), dimension(1:6) :: x, y, z
	real(8), dimension(0:pg_N2) :: r_out
	real(8), dimension(1:5, 1:6, 0:pg_N2):: z_out
	real(8) :: Coef(1:5)
	real(8) :: H0_Co, K_Co, H0_Oc, K_Oc
	real(8), dimension(1:6, 1:5) :: cY
	real(8) :: k2_
	integer :: i, ii, it, mode 
		
		open(01, file = pef_LD_Co1, status = 'replace')
		cY = 0.d0
		call lo01_bc_Coi(y(1:3))
		region = 'core'
		do i = 0, pg_N1-1
				ii = i*2
				x(1:3) = y(1:3)
				call RK4(lo01_fluid_eqt, 3, x_Co(ii), x_Co(ii+2), x(1:3), y(1:3))
				call Write1R81AR8(01, FR8,r_Co(ii), (y(1:3)), m=3)
		enddo
		close(01)


		H0_Co = y(1)
		K_Co = y(3)
!mu_Cr = mu_Cr*8.d0
!mu_Cr = mu_Cr*2.d0
!mu_Cr = mu_Cr*1.d-2

		region = 'crust'
		do it = 2, 6
			mode = it - 1
			call lo01_bc_Cri(mode, z(1:6))

			do i = 0, pg_N2-1
				ii = i*2
				x(1:6) = z(1:6)
				call RK4(lo01_solid_eqt, 6, x_Cr(ii), x_Cr(ii+2), x(1:6), z(1:6))
				r_out(i) = r_Cr(ii)
				z_out(mode, 1:6, i) = z(1:6)
			enddo
			cY(1:6,it - 1) = z(1:6)
		enddo

		open(01, file = pef_LD_Oc1, status = 'replace')
		call lo01_bc_Ocf(y(1:3))
		region = 'ocean'
		do i = pg_N3, 1, -1
				ii = i*2
				x(1:3) = y(1:3)
				call RK4(lo01_fluid_eqt, 3, x_Oc(ii), x_Oc(ii-2), x(1:3), y(1:3))
				call Write1R81AR8(01, FR8,r_Oc(ii), (y(1:3)), m=3)
		enddo
		close(01)
		H0_Oc = y(1)
		K_Oc = y(3)


		call lo01_bc_SolveFEF(H0_Co, K_Co, H0_Oc, K_Oc, cY(1:6, 1:5), k2_, Coef(1:5))


		open(01, file = pef_LD_Cr1, status = 'replace')
		do i = 0, pg_N2-1
			call Write1R81AR8(01, FR8,r_out(i), Coef(1)*z_out(1,1:6, i)+Coef(2)*z_out(2,1:6, i)+Coef(3)*z_out(3,1:6, i)+Coef(4)*z_out(4,1:6, i)+Coef(5)*z_out(5,1:6, i), m=6)
		enddo
		close(01)

		write(*,*) "k2", k2_
		write(*,*) "=========================================================="

	endsubroutine lo01_iterate_FEF

	!F1 A3---------------------------------------------------
	subroutine lo01_bc_Coi(y)
	use puls_eqt_set_opt03_ex1
	implicit none
	real(8) :: gamma_1, P_, rho_, m_, nu_, r_, a0_, a2_
	real(8) :: y(1:3)
		P_ = P_Co(0)
		rho_ = rho_Co(0)
		m_ = m_Co(0)
		nu_ = nu_Co(0)
		r_ = r_Co(0)
		gamma_1 = pes03_gamma(1, P_Co, rho_Co, 0, 0, pg_N1*2)
gamma_1 = 2.d0 * (1.d0 + P_/rho_)
		a0_ = 1.d0
		a2_ = -2.d0/(2.d0*l_0+3.d0) * pi * (5.d0*rho_+9.d0*P_+(rho_+P_)/(gamma_1*P_/(rho_+P_)) )

		y(1) = a0_ * r_**l_0 * (1.d0 + r_ ** 2 *  a2_ )
		y(2) = a0_ * r_**(l_0-1.d0) * (l_0 + r_ ** 2 * (l_0+2.d0) * a2_  )
		y(3) = a0_ * r_**l_0  * ( 1.d0 +  ((l_0+2.d0)*a2_ + 8.d0*pi/3.d0*(3.d0*P_+rho_))/(l_0+2.d0) )   ! need second order expansion !!
	endsubroutine lo01_bc_Coi

	subroutine lo01_bc_Coi2(y)
	use puls_eqt_set_opt03_ex1
	implicit none
	real(8) :: gamma_1, P_, rho_, m_, nu_, r_
	real(8) :: y(1)
		P_ = P_Co(0)
		rho_ = rho_Co(0)
		m_ = m_Co(0)
		nu_ = nu_Co(0)
		r_ = r_Co(0)
		gamma_1 = pes03_gamma(1, P_Co, rho_Co, 0, 0, pg_N1*2)
gamma_1 = 2.d0 * (1.d0 + P_/rho_)

		y(1) = l_0 - 4.d0*pi/5.d0*((6.d0-2.d0*l_0*l_0/3.d0)*rho_+10.d0*P_+(P_+rho_)/(gamma_1*P_/(rho_+P_)) )*r_*r_
	endsubroutine lo01_bc_Coi2

	subroutine lo01_bc_Coi3(y)
	use puls_eqt_set_opt03_ex1
	implicit none
	real(8) :: gamma_1, P_, rho_, m_, nu_, r_
	real(8) :: y(1:3)
		P_ = P_Co(0)
		rho_ = rho_Co(0)
		m_ = m_Co(0)
		nu_ = nu_Co(0)
		r_ = r_Co(0)
		gamma_1 = pes03_gamma(1, P_Co, rho_Co, 0, 0, pg_N1*2)
gamma_1 = 2.d0 * (1.d0 + P_/rho_)

		y(1) = (1.d0 - r_ ** 2 * 2.d0/(2.d0*l_0+3.d0) * pi * (5.d0*rho_+9.d0*P_+(rho_+P_)/(gamma_1*P_/(rho_+P_)) )  )
		y(2) = - r_ * (2.d0) *2.d0/(2.d0*l_0+3.d0) * pi * (5.d0*rho_+9.d0*P_+(rho_+P_)/(gamma_1*P_/(rho_+P_)) )  
		y(3) = (1.d0 - r_ ** 2 * 2.d0/(2.d0*l_0+3.d0) * pi * (5.d0*rho_+9.d0*P_+(rho_+P_)/(gamma_1*P_/(rho_+P_)) )  + 8.d0*pi*r_**2*(rho_/3.d0 + P_)/ (l_0 + 2.d0)  ) ! need second order expansion !!

	endsubroutine lo01_bc_Coi3
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
	real(8) :: pi4, l_1, V1, ddNu
	integer :: i, j, ii
	real(8), dimension(1:n,1:n) :: B
	
		if (n /= 3) then
			write(*,*) "err: love fluid eqt set n not 3"
			pause
		endif
		
		if ( region == 'core') then

			ii = rk4_HS(XCo_i, dx_Co,t)
			P_ = P_Co(ii)
			rho_ = rho_Co(ii)
			m_ = m_Co(ii)
			nu_ = nu_Co(ii)
			r_ = r_Co(ii)
			gamma_1 = pes03_gamma(1, P_Co, rho_Co, ii, 0, pg_N1*2)
gamma_1 = 2.d0 * (1.d0 + P_/rho_)
		elseif ( region == 'crust') then
			ii = rk4_HS(XCr_i, dx_Cr,t)
			P_ = P_Cr(ii)
			rho_ = rho_Cr(ii)
			m_ = m_Cr(ii)
			nu_ = nu_Cr(ii)
			r_ = r_Cr(ii)
			gamma_1 = pes03_gamma(1, P_Cr, rho_Cr, ii, 0, pg_N2*2)
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
		V1 = (1.d0 + rho_/P_)*r_ *dNu/2.d0
		
		Q_ = (-l_1*eLamb/r_**2 + pi4*eLamb*(5.d0*rho_+9.d0*P_+(rho_+P_)/(gamma_1*P_/(rho_+P_)) ) - dnu**2)

		B = 0.d0

		B(1,2) = r_
		B(2,1) = - Q_ * r_
		B(2,2) = - (2.d0/r_+ 0.5d0*(dnu-dLamb)) * r_
		B(3,1) = dnu * r_
		B(3,2) = r_

!B(1,2) = r_
!B(2,1) = - (Q_ + l_0*(l_0-1.d0)/r_**2 + l_0/r_ *(2.d0/r_+ 0.5d0*(dnu-dLamb)) ) * r_
!B(2,2) = - (2.d0/r_+ 0.5d0*(dnu-dLamb) + 2.d0*l_0/r_) * r_
!B(3,1) = ( l_0/r_ + dnu) * r_
!B(3,2) = r_
!B(3,3) = -l_0

		do i = 1,n
			fcn(i) = 0.d0
			do j = 1,n
				fcn(i) = fcn(i) + B(i,j)*x(j)/(1.d0 +V1)
			enddo
		enddo

	endsubroutine lo01_fluid_eqt

	subroutine lo01_fluid_eqt2(n, t, x, fcn)
	! variables y = rH0`/H0
	use puls_eqt_set_opt03_ex1
	implicit none
	integer :: n
	real(8), dimension(1:n) :: x, fcn
	real(8) :: t
	real(8) :: eNu, eLamb, eNu_hf, eLamb_hf, dLamb, dNu, gamma_1, P_, rho_, m_, nu_, r_, Q_
	real(8) :: pi4, l_1, V1, ddNu
	integer :: i, j, ii
	real(8), dimension(1:1,0:2) :: B
	
		if (n /= 1) then
			write(*,*) "err: love fluid eqt set n not 1"
			pause
		endif
		
		if ( region == 'core') then
			ii = rk4_HS(XCo_i, dx_Co,t)
			P_ = P_Co(ii)
			rho_ = rho_Co(ii)
			m_ = m_Co(ii)
			nu_ = nu_Co(ii)
			r_ = r_Co(ii)
			gamma_1 = pes03_gamma(1, P_Co, rho_Co, ii, 0, pg_N1*2)
gamma_1 = 2.d0 * (1.d0 + P_/rho_)
		elseif ( region == 'crust') then
			ii = rk4_HS(XCr_i, dx_Cr,t)
			P_ = P_Cr(ii)
			rho_ = rho_Cr(ii)
			m_ = m_Cr(ii)
			nu_ = nu_Cr(ii)
			r_ = r_Cr(ii)
			gamma_1 = pes03_gamma(1, P_Cr, rho_Cr, ii, 0, pg_N2*2)
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
		V1 = (1.d0 + rho_/P_)*r_ *dNu/2.d0

		Q_ = (-l_1*eLamb/r_**2 + pi4*eLamb*(5.d0*rho_+9.d0*P_+(rho_+P_)/(gamma_1*P_/(rho_+P_)) ) - dnu**2)

		B(1,0) = - Q_ * r_**2
		B(1,1) = -eLamb*(1.d0 + pi4*r_**2 *(P_ - rho_))
		B(1,2) = -1.d0
		
		fcn(1) = (B(1,0) + B(1,1)*x(1)+ B(1,2)*x(1)**2 )/(1.d0 +V1)

	endsubroutine lo01_fluid_eqt2


	!F1 A6---------------------------------------------------
	subroutine lo01_fluid_eqt3(n, t, x, fcn)
	! variables {H0, H0`, K}, defining metric r^l_0 H0
	use puls_eqt_set_opt03_ex1
	implicit none
	integer :: n
	real(8), dimension(1:n) :: x, fcn
	real(8) :: t
	real(8) :: eNu, eLamb, eNu_hf, eLamb_hf, dLamb, dNu, gamma_1, P_, rho_, m_, nu_, r_, Q_
	real(8) :: pi4, l_1, V1, ddNu
	integer :: i, j, ii
	real(8), dimension(1:n,1:n) :: B
	
		if (n /= 3) then
			write(*,*) "err: love fluid eqt set n not 3"
			pause
		endif
		
		if ( region == 'core') then

			ii = rk4_HS(XCo_i, dx_Co,t)
			P_ = P_Co(ii)
			rho_ = rho_Co(ii)
			m_ = m_Co(ii)
			nu_ = nu_Co(ii)
			r_ = r_Co(ii)
			gamma_1 = pes03_gamma(1, P_Co, rho_Co, ii, 0, pg_N1*2)
gamma_1 = 2.d0 * (1.d0 + P_/rho_)
		elseif ( region == 'crust') then
			ii = rk4_HS(XCr_i, dx_Cr,t)
			P_ = P_Cr(ii)
			rho_ = rho_Cr(ii)
			m_ = m_Cr(ii)
			nu_ = nu_Cr(ii)
			r_ = r_Cr(ii)
			gamma_1 = pes03_gamma(1, P_Cr, rho_Cr, ii, 0, pg_N2*2)
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
		V1 = (1.d0 + rho_/P_)*r_ *dNu/2.d0
		
		Q_ = (-l_1*eLamb/r_**2 + pi4*eLamb*(5.d0*rho_+9.d0*P_+(rho_+P_)/(gamma_1*P_/(rho_+P_)) ) - dnu**2)

		B = 0.d0

		B(1,2) = r_
		B(2,1) = - (Q_ + l_0*(l_0-1.d0)/r_**2 + l_0/r_ *(2.d0/r_+ 0.5d0*(dnu-dLamb)) ) * r_
		B(2,2) = - (2.d0/r_+ 0.5d0*(dnu-dLamb) + 2.d0*l_0/r_) * r_
		B(3,1) = ( l_0/r_ + dnu) * r_
		B(3,2) = r_
		B(3,3) = -l_0

		do i = 1,n
			fcn(i) = 0.d0
			do j = 1,n
				fcn(i) = fcn(i) + B(i,j)*x(j)/(1.d0 +V1)
			enddo
		enddo

	endsubroutine lo01_fluid_eqt3
	!F1 A7---------------------------------------------------
	subroutine lo01_solid_eqt(n, t, x, fcn)
	! variables {H1, K, H0, W, V, T2}
	use puls_eqt_set_opt03_ex1
	implicit none
	integer :: n
	real(8), dimension(1:n) :: x, fcn
	real(8) :: t
	complex(8), dimension(1:n,1:n) :: B
	complex(8), dimension(1:3,1:n) :: E
	complex(8), dimension(1:n) :: H2_C, X_C, T1_C
	real(8) :: eNu, eLamb, eNu_hf, eLamb_hf, dLamb, dNu, gamma_1, P_, rho_, m_, nu_, mu_, r_
	real(8) :: pi4, l_1, l_2_hf, V1, ddNu, dP
	integer :: i, j, ii

	if (n /= 6) then
		write(*,*) "err: love solid eqt set n not 6"
		pause
	endif

	ii = rk4_HS(XCr_i, dx_Cr,t)
	P_ = P_Cr(ii)
	rho_ = rho_Cr(ii)
	m_ = m_Cr(ii)
	nu_ = nu_Cr(ii)
	!mu_ = mu_Cr(ii) Andersson's mu seems to be off by 1/2, st it cannot be reduced to RCA case
	mu_ = mu_Cr(ii)/2.d0
	r_ = r_Cr(ii)
	gamma_1 = pes03_gamma(1, P_Cr, rho_Cr, ii, 0, pg_N2*2)
gamma_1 = 2.d0 * (1.d0 + P_/rho_)
	call pes03_E_Crust(P_, rho_, m_, nu_, mu_, gamma_1, r_, E)

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

	B(3,1) = -r_**2 * eNu**(-1.d0) * (0.d0)
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
	B(6,5) = 4.d0 * l_2_hf * eLamb * mu_ - eLamb/eNu * r_**2 * (0.d0) * (P_ + rho_)
	B(6,6) = 0.5d0 * (dLamb - dNu) * r_ - (l_0 + 1.d0)
	B(6, 1:6) = B(6, 1:6) + r_**2 * eLamb/eNu_hf * (X_C(1:6) - 0.5d0/r_**2 * eNu_hf * T1_C(1:6))

!!!!!!!! T1_C, X_C etc
	do i = 1,n
		fcn(i) = 0.d0
		do j = 1,n
			fcn(i) = fcn(i) + dreal(B(i,j)*x(j)/(1.d0 +V1))
		enddo
	enddo

	endsubroutine lo01_solid_eqt

	!F1 A8---------------------------------------------------
	subroutine lo01_bc_Solve(H0F, KF, cY, k2_, Coef)
	use lin_sol_gen_int
	use puls_eqt_set_opt03_ex1
	use lo_match
	!	2 components
	!	Quantity Q = {QQ(i)} [Y] {c(i)} ------ {vec} [2nd Tensor] {vec} inner product
	implicit none
	real(8) :: H0F, KF, cY(1:6,1:5), k2_, Coef(1:5)
	real(8) :: eNu, eLamb, eNu_hf, eLamb_hf, dLamb, dNu, dP, gamma_1, P_, rho_, m_, mu_,  nu_, r_
	real(8) :: pi4, l_1, l_2_hf, ddNu, cs2
	real(8) :: XX(1:6), HH(1:6), KK(1:6), H2_C(1:6), PP(1:6), TT(1:6), TTF
	complex(8) :: E_(1:3,1:6), AA(1:6)
	integer :: mode, i
	real(8) :: YI(1:5, 1:5),  YS(1:6, 1:5), M(1:3, 1:3), Vec(1:3,1), Sol(1:3,1), YS_Comb(1:6), HS, dHS, yR_
real(8) :: dP_Eu, T1_E, dH0_E
		!¡¸ Crust Core interface conditiion
		P_ = P_Cr(0)
		rho_ = rho_Cr(0)
		m_ = m_Cr(0)
		!mu_ = mu_Cr(ii) Andersson's mu seems to be off by 1/2, st it cannot be reduced to RCA case
		mu_ = mu_Cr(0)/2.d0
		nu_ = nu_Cr(0)
		r_ = r_Cr(0)
		gamma_1 = pes03_gamma(1, P_Cr, rho_Cr, 0, 0, pg_N2*2)
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


		!¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸
		! Solving for WE as given in Andersson (2011) Eqt(53)-(55)
		call pes03_E_Crust(P_, rho_, m_, nu_, mu_, gamma_1, r_, E_)
		XX(1:6) = dreal(E_(2,1:6))
		TT(1:6) = dreal(E_(3,1:6))

		PP(1:6) = -r_**l_0 /eNu_hf * XX(1:6)
		PP(4) = PP(4) - r_**l_0 /r_/eLamb_hf * dP
		PP = PP/r_**l_0							! Account for the difference in Metric definition H0 vs r^l_0 H0
		TTF = -r_**2/2.d0 * (P_ + rho_) * (-H0F/r_**l_0)		! Account for the difference in Metric definition H0 vs -r^l_0 H0

!!!!?SDSF>SGKFKJAFJASKF
		Coef(2) = -KF/r_**l_0						! Account for the difference in Metric definition H0 vs r^l_0 H0
		Coef(3) = -H0F/r_**l_0

		YI = 0.d0
		YI(1, 1) = 1.d0 !/r_**2
		YI(2, 2) = 1.d0
		YI(3, 3) = 1.d0
		YI(4, 4) = 1.d0
		YI(5, 5) = 1.d0

		M(1, 1:3) = 0.d0
		Vec(1,1) = TTF							! Fluid, Zero Freq: P_Eu = -0.5 * (rho+P) * H0_LD
		do i = 1, 5
			M(1, 1) = M(1, 1) + YI(i, 1) * (r_**2.d0*PP(i) - TT(i))
			M(1, 2) = M(1, 2) + YI(i, 4) * (r_**2.d0*PP(i) - TT(i))
			M(1, 3) = M(1, 3) + YI(i, 5) * (r_**2.d0*PP(i) - TT(i))
			Vec(1,1) = Vec(1,1) - (Coef(2)*YI(i, 2) + Coef(3)*YI(i, 3) )* (r_**2.d0*PP(i) - TT(i))		! Continuity of (r^2 P_Eu - T1) 
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
		P_ = P_Cr(2*pg_N2)
		rho_ = rho_Cr(2*pg_N2)
		m_ = m_Cr(2*pg_N2)
		!mu_ = mu_Cr(ii) Andersson's mu seems to be off by 1/2, st it cannot be reduced to RCA case
		mu_ = mu_Cr(2*pg_N2)/2.d0
		nu_ = nu_Cr(2*pg_N2)
		r_ = r_Cr(2*pg_N2)
		gamma_1 = pes03_gamma(1, P_Cr, rho_Cr, 2*pg_N2, 0, pg_N2*2)
gamma_1 = 2.d0 * (1.d0 + P_/rho_)

		eLamb = pes03_eLamb(rho_, m_, r_) ** 2
		eNu = dexp(2.d0*nu_)
		eLamb_hf = eLamb**(0.5d0)
		eNu_hf = eNu**(0.5d0)
		dLamb = pes03_dLamb(rho_, m_, r_) * 2.d0
		dNu = pes03_dnu(P_, rho_, m_, r_) * 2.d0
		ddNu = pes03_ddnu(P_, rho_, m_, r_) * 2.d0
		dP = -dNu/2.d0 * (P_ + rho_)
		cs2 = (P_/(rho_+P_))*gamma_1

!		call pes03_bc_R_Coe(P_, rho_, m_, nu_, mu_, gamma_1, r_, AA(1:5))
!		AA(1:3) = AA(1:3)*R0**2											! pes03_bc_R_Coe normalization constants
!		XXS(1:5) = dreal(AA(1:5))
!		XXS(6) = 0.d0

		call pes03_E_Crust(P_, rho_, m_, nu_, mu_, gamma_1, r_, E_)
		XX(1:6) = dreal(E_(2,1:6))	
		TT(1:6) = dreal(E_(3,1:6))

		Vec(2,1) = 0.d0
		M(2, 1:3) = 0.d0
		do i = 1, 6
			M(2, 1) = M(2, 1) + YS(i, 1) * (r_**2.d0/eNu_hf*XX(i) + TT(i))
			M(2, 2) = M(2, 2) + YS(i, 4) * (r_**2.d0/eNu_hf*XX(i) + TT(i))
			M(2, 3) = M(2, 3) + YS(i, 5) * (r_**2.d0/eNu_hf*XX(i) + TT(i))
			Vec(2,1) = Vec(2,1) - (Coef(2) * YS(i, 2) + Coef(3) * YS(i, 3))* (r_**2.d0/eNu_hf*XX(i) + TT(i))
		enddo
write(*,*) "rho jump:"
write(*,*) rho_Co(2*pg_N1)
write(*,*) rho_Cr(0)
pause
		!¡¸ T2 = 0
		M(3, 1) = YS(6, 1)
		M(3, 2) = YS(6, 4)
		M(3, 3) = YS(6, 5)
		Vec(3,1) = 0.d0 - (     Coef(2) * YS(6, 2) + Coef(3) * YS(6, 3)    )

		call LIN_SOL_GEN(M, Vec, Sol)					!IMSL linear ODE solver
		coef(1) = Sol(1,1)
		coef(4) = Sol(2,1)
		coef(5) = Sol(3,1)
		
		YS_Comb = 0.d0
		do i = 1,5
			YS_Comb(1:6) = YS_Comb(1:6) + coef(i)*YS(1:6, i)
		enddo
		call pes03_E_Crust(P_, rho_, m_, nu_, mu_, gamma_1, r_, E_)

		H2_C(1:6) = E_(1,1:6)
		KK = 0.d0
		HH = 0.d0

		KK(1) = l_2_hf + 1.d0
		KK(2) = 0.5d0* dNu*r_ - (l_0 + 1.d0)
		KK(4) = -2.d0*pi4*(P_ + rho_)*eLamb_hf
		KK(1:6) = KK(1:6) + H2_C(1:6)

		HH(2) = l_0
		HH(3) = -(0.5d0* dNu*r_ + (l_0 - 1.d0))
		HH(6) = -4.d0*pi4
		HH(1:6) = (HH(1:6) + KK(1:6) - (0.5d0* dNu*r_ + 1.d0)*H2_C(1:6))/r_

		dHS = YS_Comb(1)*HH(1)+YS_Comb(2)*HH(2)+YS_Comb(3)*HH(3)+YS_Comb(4)*HH(4)+YS_Comb(5)*HH(5)+YS_Comb(6)*HH(6)
		HS = YS_Comb(3)

		yR_ = r_*dHS/HS + l_0
		
! Confirm that the Linear Equation Solver is accurate
!write(*,*) M(1, 1)*coef(1) + M(1, 2)*coef(4) + M(1, 3)*coef(5)
!write(*,*) Vec(1,1)
!write(*,*) M(2, 1)*coef(1) + M(2, 2)*coef(4) + M(2, 3)*coef(5)
!write(*,*) Vec(2,1)
!write(*,*) M(3, 1)*coef(1) + M(3, 2)*coef(4) + M(3, 3)*coef(5)
!write(*,*) Vec(3,1)
!pause
call lo_k2(yR_, YS_Comb(5)/HS, P_, rho_, m_, mu_, r_, k2_)
		!call lo_k2(yR_, m_, r_, k2_)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TEST TEST TEST TEST TEST TEST TEST
!¡¸ Crust Core interface conditiion
		P_ = P_Cr(0)
		rho_ = rho_Cr(0)
		m_ = m_Cr(0)
		!mu_ = mu_Cr(ii) Andersson's mu seems to be off by 1/2, st it cannot be reduced to RCA case
		mu_ = mu_Cr(0)/2.d0
		nu_ = nu_Cr(0)
		r_ = r_Cr(0)
		gamma_1 = pes03_gamma(1, P_Cr, rho_Cr, 0, 0, pg_N2*2)
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


		!¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸
		! Solving for WE as given in Andersson (2011) Eqt(53)-(55)
		call pes03_E_Crust(P_, rho_, m_, nu_, mu_, gamma_1, r_, E_)
		XX(1:6) = dreal(E_(2,1:6))
		TT(1:6) = dreal(E_(3,1:6))
		H2_C(1:6) = E_(1,1:6)
		KK = 0.d0
		HH = 0.d0

		KK(1) = l_2_hf + 1.d0
		KK(2) = 0.5d0* dNu*r_ - (l_0 + 1.d0)
		KK(4) = -2.d0*pi4*(P_ + rho_)*eLamb_hf
		KK(1:6) = KK(1:6) + H2_C(1:6)

		HH(2) = l_0
		HH(3) = -(0.5d0* dNu*r_ + (l_0 - 1.d0))
		HH(6) = -4.d0*pi4
		HH(1:6) = (HH(1:6) + KK(1:6) - (0.5d0* dNu*r_ + 1.d0)*H2_C(1:6))/r_

		PP(1:6) = -r_**l_0 /eNu_hf * XX(1:6)
		PP(4) = PP(4) - r_**l_0 /r_/eLamb_hf * dP
		PP = PP/r_**l_0							! Account for the difference in Metric definition H0 vs r^l_0 H0
		TTF = -r_**2/2.d0 * (P_ + rho_) * (-H0F/r_**l_0)		! Account for the difference in Metric definition H0 vs -r^l_0 H0

!!!!?SDSF>SGKFKJAFJASKF
		Coef(2) = -KF/r_**l_0				! Account for the difference in Metric definition H0 vs r^l_0 H0
		Coef(3) = -H0F/r_**l_0

		YI = 0.d0
		YI(1, 1) = 1.d0 !/r_**2
		YI(2, 2) = 1.d0
		YI(3, 3) = 1.d0
		YI(4, 4) = 1.d0
		YI(5, 5) = 1.d0

		dP_Eu = 0.d0
		T1_E = 0.d0
		dH0_E = 0.d0
		do i = 1, 5
			dP_Eu = dP_Eu + Coef(i)*YI(i, i) * PP(i)
			T1_E = T1_E + TT(i)*coef(i)
			dH0_E = dH0_E + HH(i)*coef(i)
		enddo
write(*,*) "LHS (54)", -3.d0/8.d0/mu_ * T1_E
write(*,*) "RHS (54)", r_**2 *dP_Eu/(rho_+P_)/cs2 - 1.5d0 * r_**2 * coef(2) +1.5d0*l_1 * coef(5)  + (3.d0 + 2.d0* r_*dP/(rho_+P_)/2.d0/cs2) * coef(4) /eLamb_hf
pause

write(*,*) "LHS (55)", -T1_E
write(*,*) "RHS (55)", r_**2 *(  0.5d0 * (rho_+P_)*(H0F/r_**l_0) -  dP_Eu)
pause

write(*,*) "LHS (53)", 16.d0*pi*r_**2 * eLamb * dP_Eu
write(*,*) "RHS (53)", r_**2 * dNu * ((-dH0_E) + l_0 /r_ * (-coef(3))) + (-2.d0*l_2_hf)*eLamb*(-coef(2)) + (l_1*eLamb - 2.d0 + r_**2 * dNu**2) * (-coef(3)) + 32.d0*pi*mu_*( r_**2*dNu**2 + (2.d0*l_2_hf)*(1.d0-eLamb) )*(-Coef(5)) + 8.d0*pi*(eLamb - 3.d0)*(-T1_E)
write(*,*) "Delta (53)",(16.d0*pi*r_**2 * eLamb * dP_Eu - (r_**2 * dNu * ((-dH0_E) + l_0 /r_ * (-coef(3))) + (-2.d0*l_2_hf)*eLamb*(-coef(2)) + (l_1*eLamb - 2.d0 + r_**2 * dNu**2) * (-coef(3)) + 32.d0*pi*mu_*( r_**2*dNu**2 + (2.d0*l_2_hf)*(1.d0-eLamb) )*(-Coef(5)) + 8.d0*pi*(eLamb - 3.d0)*(-T1_E))   )/(16.d0*pi*r_**2 * eLamb * dP_Eu)
pause
		
	endsubroutine lo01_bc_Solve

!F1 A8---------------------------------------------------
	subroutine lo01_bc_SolveFEF(H0F, KF, H0Oc, KOc, cY, k2_, Coef)
	use lin_sol_gen_int
	use puls_eqt_set_opt03_ex1
	use lo_match
	!	2 components
	!	Quantity Q = {QQ(i)} [Y] {c(i)} ------ {vec} [2nd Tensor] {vec} inner product
	implicit none
	real(8) :: H0F, KF, H0Oc, KOc, cY(1:6,1:5), k2_, Coef(1:5)
	real(8) :: eNu, eLamb, eNu_hf, eLamb_hf, dLamb, dNu, dP, gamma_1, P_, rho_, m_, mu_,  nu_, r_
	real(8) :: pi4, l_1, l_2_hf, ddNu, cs2
	real(8) :: XX(1:6), XXS(1:6), HH(1:6), KK(1:6), H2_C(1:6), PP(1:6), TT(1:6), TTF
	complex(8) :: E_(1:3,1:6), AA(1:6)
	integer :: mode, i
	real(8) :: YI(1:5, 1:5),  YS(1:6, 1:5), M(1:3, 1:3), Vec(1:3,1), Sol(1:3,1), YS_Comb(1:6), HS, dHS, yR_
		
		!¡¸ Crust Core interface conditiion
		P_ = P_Cr(0)
		rho_ = rho_Cr(0)
		m_ = m_Cr(0)
		!mu_ = mu_Cr(ii) Andersson's mu seems to be off by 1/2, st it cannot be reduced to RCA case
		mu_ = mu_Cr(0)/2.d0
		nu_ = nu_Cr(0)
		r_ = r_Cr(0)
		gamma_1 = pes03_gamma(1, P_Cr, rho_Cr, 0, 0, pg_N2*2)
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


		!¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸¡¸
		! Solving for WE as given in Andersson (2011) Eqt(53)-(55)
		call pes03_E_Crust(P_, rho_, m_, nu_, mu_, gamma_1, r_, E_)
		XX(1:6) = dreal(E_(2,1:6))
		TT(1:6) = dreal(E_(3,1:6))

		PP(1:6) = -r_**l_0 /eNu_hf * XX(1:6)
		PP(4) = PP(4) - r_**l_0 /r_/eLamb_hf * dP
		PP = PP/r_**l_0							! Account for the difference in Metric definition H0 vs r^l_0 H0
		TTF = -r_**2/2.d0 * (P_ + rho_) * (-H0F/r_**l_0)		! Account for the difference in Metric definition H0 vs -r^l_0 H0

!!!!?SDSF>SGKFKJAFJASKF
		Coef(2) = -KF/r_**l_0				! Account for the difference in Metric definition H0 vs r^l_0 H0
		Coef(3) = -H0F/r_**l_0

		YI = 0.d0
		YI(1, 1) = 1.d0 !/r_**2
		YI(2, 2) = 1.d0
		YI(3, 3) = 1.d0
		YI(4, 4) = 1.d0
		YI(5, 5) = 1.d0

		M(1, 1:3) = 0.d0
		Vec(1,1) = TTF							! Fluid, Zero Freq: P_Eu = -0.5 * (rho+P) * H0_LD
		do i = 1, 5
			M(1, 1) = M(1, 1) + YI(i, 1) * (r_**2.d0*PP(i) - TT(i))
			M(1, 2) = M(1, 2) + YI(i, 4) * (r_**2.d0*PP(i) - TT(i))
			M(1, 3) = M(1, 3) + YI(i, 5) * (r_**2.d0*PP(i) - TT(i))
			Vec(1,1) = Vec(1,1) - (Coef(2)*YI(i, 2) + Coef(3)*YI(i, 3) )* (r_**2.d0*PP(i) - TT(i))		! Continuity of (r^2 P_Eu - T1) 
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
		P_ = P_Cr(2*pg_N2)
		rho_ = rho_Cr(2*pg_N2)
		m_ = m_Cr(2*pg_N2)
		!mu_ = mu_Cr(ii) Andersson's mu seems to be off by 1/2, st it cannot be reduced to RCA case
		mu_ = mu_Cr(2*pg_N2)/2.d0
		nu_ = nu_Cr(2*pg_N2)
		r_ = r_Cr(2*pg_N2)
		gamma_1 = pes03_gamma(1, P_Cr, rho_Cr, 2*pg_N2, 0, pg_N2*2)
gamma_1 = 2.d0 * (1.d0 + P_/rho_)

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
		call pes03_E_Crust(P_, rho_, m_, nu_, mu_, gamma_1, r_, E_)
		XX(1:6) = dreal(E_(2,1:6))
		TT(1:6) = dreal(E_(3,1:6))

		PP(1:6) = -r_**l_0 /eNu_hf * XX(1:6)
		PP(4) = PP(4) - r_**l_0 /r_/eLamb_hf * dP
		PP = PP/r_**l_0							! Account for the difference in Metric definition H0 vs r^l_0 H0
		TTF = -r_**2/2.d0 * (P_ + rho_) * (-H0Oc/r_**l_0)		! Account for the difference in Metric definition H0 vs -r^l_0 H0
		
		Vec(2,1) = 0.d0
		M(2, 1:3) = 0.d0
		do i = 1, 6
			M(2, 1) = M(2, 1) + YS(i, 1) * (r_**2.d0*PP(i) - TT(i))
			M(2, 2) = M(2, 2) + YS(i, 4) * (r_**2.d0*PP(i) - TT(i))
			M(2, 3) = M(2, 3) + YS(i, 5) * (r_**2.d0*PP(i) - TT(i))
			Vec(2,1) = Vec(2,1) - (Coef(2) * YS(i, 2) + Coef(3) * YS(i, 3))* (r_**2.d0*PP(i) - TT(i))
		enddo

		!¡¸ T2 = 0
		M(3, 1) = YS(6, 1)
		M(3, 2) = YS(6, 4)
		M(3, 3) = YS(6, 5)
		Vec(3,1) = 0.d0 - (     Coef(2) * YS(6, 2) + Coef(3) * YS(6, 3)    )

		call LIN_SOL_GEN(M, Vec, Sol)					!IMSL linear equations solver
		coef(1) = Sol(1,1)
		coef(4) = Sol(2,1)
		coef(5) = Sol(3,1)
		
		YS_Comb = 0.d0
		do i = 1,5
			YS_Comb(1:6) = YS_Comb(1:6) + coef(i)*YS(1:6, i)
		enddo
		call pes03_E_Crust(P_, rho_, m_, nu_, mu_, gamma_1, r_, E_)

		H2_C(1:6) = E_(1,1:6)
		KK = 0.d0
		HH = 0.d0

		KK(1) = l_2_hf + 1.d0
		KK(2) = 0.5d0* dNu*r_ - (l_0 + 1.d0)
		KK(4) = -2.d0*pi4*(P_ + rho_)*eLamb_hf
		KK(1:6) = KK(1:6) + H2_C(1:6)

		HH(2) = l_0
		HH(3) = -(0.5d0* dNu*r_ + (l_0 - 1.d0))
		HH(6) = -4.d0*pi4
		HH(1:6) = (HH(1:6) + KK(1:6) - (0.5d0* dNu*r_ + 1.d0)*H2_C(1:6))/r_

		dHS = YS_Comb(1)*HH(1)+YS_Comb(2)*HH(2)+YS_Comb(3)*HH(3)+YS_Comb(4)*HH(4)+YS_Comb(5)*HH(5)+YS_Comb(6)*HH(6)
		HS = YS_Comb(3)

		yR_ = r_*dHS/HS + l_0
		
! Confirm that the Linear Equation Solver is accurate
!write(*,*) M(1, 1)*coef(1) + M(1, 2)*coef(4) + M(1, 3)*coef(5)
!write(*,*) Vec(1,1)
!write(*,*) M(2, 1)*coef(1) + M(2, 2)*coef(4) + M(2, 3)*coef(5)
!write(*,*) Vec(2,1)
!write(*,*) M(3, 1)*coef(1) + M(3, 2)*coef(4) + M(3, 3)*coef(5)
!write(*,*) Vec(3,1)
!pause
call lo_k2(yR_, YS_Comb(5)/HS, P_, rho_, m_, mu_, r_, k2_)		
		!call lo_k2(yR_, m_, r_, k2_)
		
	endsubroutine lo01_bc_SolveFEF

	!F1 B1---------------------------------------------------
	function rk4_HS(xi, dx,x)
	implicit none
	real(8) :: xi, dx, x
	real(8) :: dx_hf
	integer :: rk4_HS
		dx_hf = dx/2.d0  ! dx/2; half step size
		rk4_HS = nint((x-xi)/dx_hf)
	endfunction rk4_HS

endmodule love_opt01_P01