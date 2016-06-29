module moI_opt01

! Finding the GR moment of inertia


use global_var
contains

subroutine mi01_Ctrl(I_, I_bar, Eta)
implicit none
real(8) :: I_, I_bar, Eta
	
	if (pe_fcall == .true.) then
		write(100,*) "moI_opt01: linear coupled ode solver"
			
		if (pes_opt == 3) then
			write(*,*) "Fluid Star"
			write(*,*) "=========================================================="
		elseif (pes_opt == 5) then
			write(*,*) "Star with Solid Crust"
			write(*,*) "=========================================================="
		else
			write(*,*) "err: pes_opt"
			pause
		endif
	endif

		call mi01_iterate(I_, I_bar, Eta)

	if (pe_fcall == .true.) then
		write(100,*) "mi01_: passed"
		write(100,*) "=========================================================="
		pe_fcall = .false.
	endif

endsubroutine mi01_Ctrl


subroutine mi01_iterate(I_, I_bar, Eta)
use Format_IO
use FWrite
use RK4_Set
implicit none
real(8), dimension(1:2) :: x, y
real(8) :: I_, I_bar, Eta
integer :: i, ii
	
	!* Integrate drag angular velocity
	call mi01_bc_Coi(y(1:2))
	open(01, file = 'data/moI_y_Co.txt', status = 'replace')
	call Write4R8(01, FR8,r_Co(0),x_Co(0), y(1), y(2))
		do i = 0, pg_N1-1
			ii = i*2
			x(1:2) = y(1:2)
			call RK4(mi01_vx_Co, 2, x_Co(ii), x_Co(ii+2), x(1:2), y(1:2))
			call Write4R8(01, FR8,r_Co(ii+2),x_Co(ii+2), y(1), y(2))
		enddo
	close(01)
	
	x(1:2) = y(1:2)

	open(01, file = 'data/moI_y_Cr.txt', status = 'replace')
	call Write4R8(01, FR8,r_Cr(0),x_Cr(0), x(1), x(2))
		do i = 0, pg_N2-1
			ii = i*2
			x(1:2) = y(1:2)
			call RK4(mi01_vx_Cr, 2, x_Cr(ii), x_Cr(ii+2), x(1:2), y(1:2))
			call Write4R8(01, FR8,r_Cr(ii+2),x_Cr(ii+2), y(1), y(2))
		enddo
	close(01)

	!* Calculate and output moment of inertia
	call mi01_I(y(1:2), I_, I_bar, Eta)

endsubroutine mi01_iterate


!*********************************************************
!	Fcn---------------------------------------------------

function mi01_eLamb(rho_, m_, r_)
implicit none
real(8) :: mi01_eLamb, rho_, m_, r_
	mi01_eLamb = (1.d0 - rel * 2.d0*m_/r_)**(-0.5d0)
endfunction mi01_eLamb

function mi01_dnu(P_, rho_, m_, r_)
	implicit none
	real(8) :: mi01_dnu, P_, rho_, m_, r_, c1
		c1 = mi01_eLamb(rho_, m_, r_) **2
		mi01_dnu = c1 / r_ * rel * ( rel * 4.d0*pi*r_**(2)*P_ + m_/r_)
endfunction mi01_dnu

function rk4_HS(xi, dx,x)
implicit none
real(8) :: xi, dx, x
real(8) :: dx_hf
integer :: rk4_HS
		dx_hf = dx/2.d0
		rk4_HS = nint((x-xi)/dx_hf)
endfunction rk4_HS

function mi01_I_uni_R(Eta_)
! ref: H K Lau, P T Leung, L M Lin (2010)
implicit none
real(8) :: Eta_
real(8) :: mi01_I_uni_R
		mi01_I_uni_R = -0.0047d0 + 0.133d0 * Eta_ + 0.575d0 * Eta_ **2
endfunction mi01_I_uni_R

function mi01_I_uni_I(Eta_)
implicit none
real(8) :: Eta_
real(8) :: mi01_I_uni_I
		mi01_I_uni_I = 0.00694d0 - 0.0256d0 * Eta_ **2
endfunction mi01_I_uni_I

!*********************************************************
!	BC---------------------------------------------------

subroutine mi01_bc_Coi(y)
implicit none
real(8) :: y(1:)
real(8) :: C0, P_, rho_, m_, r_, eLamb
	C0 = 1.d0
	P_ = P_Co(0)
	rho_ = rho_Co(0)
	m_ = m_Co(0)
	r_ = r_Co(0)
	eLamb = mi01_eLamb(rho_, m_, r_) ** 2

	y(1) = C0 + 8.d0/5.d0 * pi * (P_ + rho_) * eLamb * C0 * r_ **2
	y(2) = 16.d0/5.d0 * pi * (P_ + rho_) * eLamb * C0 * r_

endsubroutine mi01_bc_Coi


subroutine mi01_I(y, I_, I_bar, Eta)
use Format_IO
use FWrite
implicit none
real(8) :: y(1:)
real(8) :: m_, r_, I_, I_bar, Eta

	r_ = R0
	m_ = M0

	I_ = r_**4/ (6.d0 * y(1)/y(2) + 2.d0 * r_)
	I_bar = I_/m_**3
	Eta = dsqrt(1.d0/I_bar)
	
	write(*,*) "GR Moment of Inertia:"
	write(*,*) "Moment of inertia = ", I_
	write(*,*) "Moment of inertia normalized = ", I_bar
	write(*,*) "Eta = ", Eta
	write(*,*) "freq real part (normalized)", m_ *dreal(cdsqrt(OmeC_sq))
	write(*,*) "freq imaginary part (normalized)", I_**2/m_**5  * dimag(cdsqrt(OmeC_sq))
	write(*,*) "Predicted f-freq real part", mi01_I_uni_R(Eta)
	write(*,*) "Predicted f-freq imaginary part", mi01_I_uni_I(Eta)
	write(*,*) "=========================================================="

	open(01, file = 'data/moI.txt', status = 'replace')
		!! Choose Output
!		call Write5R8(01, FR8,I_,I_bar, Eta, m_ *dreal(cdsqrt(OmeC_sq)), I_**2/m_**5  * dimag(cdsqrt(OmeC_sq)) )
!		call Write11R8(01, FR8,rho_0, m_/r_, I_,I_bar, Eta, 	dreal(cdsqrt(OmeC_sq))*c/2.d0/pi, dimag(cdsqrt(OmeC_sq))*c/2.d0/pi, m_ *dreal(cdsqrt(OmeC_sq)), I_**2/m_**5  * dimag(cdsqrt(OmeC_sq)), mi01_I_uni_R(Eta), mi01_I_uni_I(Eta) )
!		call Write12R8(01, FR8,rho_0, m_/r_, m_Co(2*pg_N1)/M0, I_,I_bar, Eta, 	dreal(cdsqrt(OmeC_sq))*c/2.d0/pi, dimag(cdsqrt(OmeC_sq))*c/2.d0/pi, m_ *dreal(cdsqrt(OmeC_sq)), I_**2/m_**5  * dimag(cdsqrt(OmeC_sq)), mi01_I_uni_R(Eta), mi01_I_uni_I(Eta) )
		call Write12R8(01, FR8,rho_0, r_Co(2*pg_N1)/R0, m_Co(2*pg_N1)/M0, I_,I_bar, Eta, 	dreal(cdsqrt(OmeC_sq))*c/2.d0/pi, dimag(cdsqrt(OmeC_sq))*c/2.d0/pi, m_ *dreal(cdsqrt(OmeC_sq)), I_**2/m_**5  * dimag(cdsqrt(OmeC_sq)), mi01_I_uni_R(Eta), mi01_I_uni_I(Eta) )
!		call Write11R8(01, FR8,mu(sp_N1), mu(sp_N1)/P(sp_N1), I_,I_bar, Eta, 	dreal(cdsqrt(OmeC_sq))*c/2.d0/pi, dimag(cdsqrt(OmeC_sq))*c/2.d0/pi, m_ *dreal(cdsqrt(OmeC_sq)), I_**2/m_**5  * dimag(cdsqrt(OmeC_sq)), mi01_I_uni_R(Eta), mi01_I_uni_I(Eta) )
	close(01)


		write(*,*) "=========================================================="

endsubroutine mi01_I


!*********************************************************
!	Eqts---------------------------------------------------
	subroutine mi01_vx_Co(n, t, x, fcn)
	implicit none
	integer :: n
	real(8), dimension(1:n) :: x, fcn
	real(8) :: t
	real(8), dimension(1:2,1:2) :: B
	real(8) :: eLamb, dNu, V1, P_, rho_, m_ , r_
	integer :: i, j, ii

	if (n /= 2) then
		write(*,*) "err: mom I 01:  n /= 2"
		pause
	endif

	ii = rk4_HS(XCo_i, dx_Co,t)
	P_ = P_Co(ii)
	rho_ = rho_Co(ii)
	m_ = m_Co(ii)
	r_ = r_Co(ii)

	eLamb = mi01_eLamb(rho_, m_, r_) ** 2
	dNu = mi01_dnu(P_, rho_, m_, r_) * 2.d0
	V1 = (1.d0 + rho_/P_)*r_ *dNu/2.d0

	B(1,1) = 0.d0
	B(1,2) = 1.d0

	B(2,1) = 16.d0 * pi * (P_ + rho_) * eLamb
	B(2,2) = -4.d0 / r_ + 4.d0 * pi * r_ * (P_ + rho_) * eLamb

	do i = 1,n
		fcn(i) = 0.d0
		do j = 1,n
			fcn(i) = fcn(i) + B(i,j)*x(j)*r_/(1.d0 +V1)
		enddo
	enddo
	endsubroutine mi01_vx_Co

	subroutine mi01_vx_Cr(n, t, x, fcn)
	implicit none
	integer :: n
	real(8), dimension(1:n) :: x, fcn
	real(8) :: t
	real(8), dimension(1:2,1:2) :: B
	real(8) :: eLamb, dNu, V1, P_, rho_, m_ , r_
	integer :: i, j, ii

	if (n /= 2) then
		write(*,*) "err: mom I 01:  n /= 2"
		pause
	endif

	ii = rk4_HS(XCr_i, dx_Cr,t)
	P_ = P_Cr(ii)
	rho_ = rho_Cr(ii)
	m_ = m_Cr(ii)
	r_ = r_Cr(ii)

	eLamb = mi01_eLamb(rho_, m_, r_) ** 2
	dNu = mi01_dnu(P_, rho_, m_, r_) * 2.d0
	V1 = (1.d0 + rho_/P_)*r_ *dNu/2.d0

	B(1,1) = 0.d0
	B(1,2) = 1.d0

	B(2,1) = 16.d0 * pi * (P_ + rho_) * eLamb
	B(2,2) = -4.d0 / r_ + 4.d0 * pi * r_ * (P_ + rho_) * eLamb

	do i = 1,n
		fcn(i) = 0.d0
		do j = 1,n
			fcn(i) = fcn(i) + B(i,j)*x(j)*r_/(1.d0 +V1)
		enddo
	enddo
	endsubroutine mi01_vx_Cr


endmodule moI_opt01