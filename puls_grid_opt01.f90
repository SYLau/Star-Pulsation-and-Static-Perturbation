module puls_grid_opt01
use global_var

contains

!C2---------------------------------------------------
subroutine pg01_Ctrl
implicit none

	write(100,*) "puls_grid_opt01: interpolate for the new variable x"
	write(100,*) "- compute the grid for pulsation equations"

	call pg01_gridx
	write(100,*) "pg01_gridx: passed"
	
	call pg01_display
	write(100,*) "=========================================================="
endsubroutine pg01_Ctrl

!C2 A1---------------------------------------------------
subroutine pg01_gridx
implicit none
real(8), dimension(1:sp_Nt) :: x_j

	call pg01_xj(x_j)

	call pg01_xi_iniCo
	call pg01_xi_finCo
	call pg01_xi_Co(x_j)


	call pg01_xi_iniCr
	call pg01_xi_finCr
	call pg01_xi_Cr(x_j)

	call pg01_ri_Ou
endsubroutine pg01_gridx

!C2 B1---------------------------------------------------
subroutine pg01_xj(x_j)
use Format_IO
use FWrite
implicit none
real(8) :: x_j(1:sp_Nt)
real(8) :: Ri, Rs
integer :: j

	Ri = R0 * pg_ri_frac
	Rs = R0 - (R0 - R_mid)* pg_rs_frac
R0 = Rs			!!!!!!!!!!!!******************** Rs /= R0 Causes a weird mode to Arise in LD formalism (Conversion from H0 K to Z Z`; a singular point is involved)

	do j = 1, sp_N1-2
		if ((sp_r(j) - Ri)*(sp_r(j+1) - Ri) <= 0.d0) pg_ji = j
		x_j(j) = dlog(sp_r(j)/P(j))
	enddo
	if (pg_ji < 1) then
		write(*,*) "err: starting r too small, plz adjust pg_ri_frac"
		pause
	endif
	x_j(sp_N1-1) = dlog(sp_r(sp_N1-1)/P(sp_N1-1))
	XCo_i = x_j(pg_ji)
	XCo_f = x_j(sp_N1-1)

	do j = sp_N1, sp_N2-1
		if ((sp_r(j) - Rs)*(sp_r(j+1) - Rs) <= 0.d0) pg_js = j+1       ! taking the point closer to star surface (a bug for LD program, r_Cr /= R0, a false mode found)
		x_j(j) = dlog(sp_r(j)/P(j))
	enddo
	x_j(sp_N2) = dlog(sp_r(sp_N2)/P(sp_N2))
	XCr_i = x_j(sp_N1)
	XCr_f = x_j(pg_js)

endsubroutine pg01_xj

!C2 B2---------------------------------------------------
subroutine pg01_xi_iniCo
implicit none

	dx_Co = (XCo_f - XCo_i)/pg_N1
	r_Co(0) = sp_r(pg_ji)
	x_Co(0) = XCo_i
	P_Co(0) = P(pg_ji)
	rho_Co(0) = rho(pg_ji)
	m_Co(0) = m(pg_ji)
	mu_Co(0) = mu(pg_ji)

endsubroutine pg01_xi_iniCo

!C2 B3---------------------------------------------------
subroutine pg01_xi_Co(x_j)
use interpolation
use Format_IO
use FWrite
implicit none
real(8), dimension(1:sp_Nt) :: x_j
integer :: ii, iiN1, mode
	mode = 4
	iiN1 = pg_N1 *2
	open(41, file='data/PG_Co.txt',status='replace')
	call Write6R8(41, FR8, r_Co(0), x_Co(0), P_Co(0), rho_Co(0), m_Co(0), mu_Co(0))
	do ii = 1, iiN1-1
		x_Co(ii) = dx_Co/2.d0 + x_Co(ii-1)
		r_Co(ii) = Int_ada(mode, pg_ji, sp_N1-1, x_j, sp_r(1:sp_Nt), x_Co(ii), "pulsation grid interpolation; core radius")
		P_Co(ii) = Int_ada(mode, pg_ji, sp_N1-1, x_j, P(1:sp_Nt), x_Co(ii), "pulsation grid interpolation; core pressure")
		rho_Co(ii) = Int_ada(mode, pg_ji, sp_N1-1, x_j, rho(1:sp_Nt), x_Co(ii), "pulsation grid interpolation; core density")
		m_Co(ii) = Int_ada(mode, pg_ji, sp_N1-1, x_j, m(1:sp_Nt), x_Co(ii), "pulsation grid interpolation; core mass")
		mu_Co(ii) = Int_ada(mode, pg_ji, sp_N1-1, x_j, mu(1:sp_Nt), x_Co(ii), "pulsation grid interpolation; core shear modulus")
		call Write6R8(41, FR8, r_Co(ii), x_Co(ii), P_Co(ii), rho_Co(ii), m_Co(ii), mu_Co(ii))
	enddo
	call Write6R8(41, FR8, r_Co(pg_N1 *2), x_Co(pg_N1 *2), P_Co(pg_N1 *2), rho_Co(pg_N1 *2), m_Co(pg_N1 *2), mu_Co(pg_N1 *2))
	close(41)
endsubroutine pg01_xi_Co
 
!C2 B4---------------------------------------------------
subroutine pg01_xi_finCo
implicit none

	r_Co(pg_N1 *2) = sp_r(sp_N1-1)
	x_Co(pg_N1 *2) = XCo_f
	P_Co(pg_N1 *2) = P(sp_N1-1)
	rho_Co(pg_N1 *2) = rho(sp_N1-1)
	m_Co(pg_N1 *2) = m(sp_N1-1)
	mu_Co(pg_N1 *2) = mu(sp_N1-1)

endsubroutine pg01_xi_finCo

!C2 B5---------------------------------------------------
subroutine pg01_xi_iniCr
implicit none

	dx_Cr = (XCr_f - XCr_i)/pg_N2
	r_Cr(0) = sp_r(sp_N1)
	x_Cr(0) = XCr_i
	P_Cr(0) = P(sp_N1)
	rho_Cr(0) = rho(sp_N1)
	m_Cr(0) = m(sp_N1)
	mu_Cr(0) = mu(sp_N1)

endsubroutine pg01_xi_iniCr
 
!C2 B6---------------------------------------------------
subroutine pg01_xi_Cr(x_j)
use interpolation
use Format_IO
use FWrite
implicit none
real(8), dimension(1:sp_Nt) :: x_j
integer :: ii, iiN2, mode
	mode = 4
	iiN2 = pg_N2*2
	open(42, file='data/PG_Cr.txt',status='replace')
	call Write6R8(42, FR8, r_Cr(0), x_Cr(0), P_Cr(0), rho_Cr(0), m_Cr(0), mu_Cr(0))

	do ii = 1, iiN2-1
		x_Cr(ii) = dx_Cr/2.d0 + x_Cr(ii-1)
		r_Cr(ii) = Int_ada(mode, sp_N1, pg_js, x_j, sp_r(1:sp_Nt), x_Cr(ii), "pulsation grid interpolation; crust radius")
		P_Cr(ii) = Int_ada(mode, sp_N1, pg_js, x_j, P(1:sp_Nt), x_Cr(ii), "pulsation grid interpolation; crust pressure")
		rho_Cr(ii) = Int_ada(mode, sp_N1, pg_js, x_j, rho(1:sp_Nt), x_Cr(ii), "pulsation grid interpolation; crust density")
		m_Cr(ii) = Int_ada(mode, sp_N1, pg_js, x_j, m(1:sp_Nt), x_Cr(ii), "pulsation grid interpolation; crust mass")
		mu_Cr(ii) = Int_ada(mode, sp_N1, pg_js, x_j, mu(1:sp_Nt), x_Cr(ii), "pulsation grid interpolation; crust shear modulus")
		call Write6R8(42, FR8, r_Cr(ii), x_Cr(ii), P_Cr(ii), rho_Cr(ii), m_Cr(ii), mu_Cr(ii))
	enddo
	call Write6R8(42, FR8, r_Cr(pg_js), x_Cr(pg_js), P_Cr(pg_js), rho_Cr(pg_js), m_Cr(pg_js), mu_Cr(pg_js))
	close(42)
endsubroutine pg01_xi_Cr

!C2 B7---------------------------------------------------
subroutine pg01_xi_finCr
implicit none

	r_Cr(pg_N2*2) = sp_r(pg_js)
	x_Cr(pg_N2*2) = XCr_f
	P_Cr(pg_N2*2) = P(pg_js)
	rho_Cr(pg_N2*2) = rho(pg_js)
	m_Cr(pg_N2*2) = m(pg_js)
	mu_Cr(pg_N2*2) = mu(pg_js)

endsubroutine pg01_xi_finCr

!C2 B8---------------------------------------------------
subroutine pg01_ri_Ou
use Format_IO
use FWrite
implicit none
integer :: i
	P_Ou = P(pg_js)
	rho_Ou = rho(pg_js)
	m_Ou = m(pg_js)
	nu_Ou = nu(pg_js)
	ROu_i = R0
endsubroutine pg01_ri_Ou

!C2 A2---------------------------------------------------
subroutine pg01_display
implicit none
	write(*,*) "Pulsation Grid Parameters:"
	write(*,*) "XCo_i", XCo_i
	write(*,*) "XCo_f", XCo_f
	write(*,*) "XCr_i", XCr_i
	write(*,*) "XCr_f", XCr_f
	write(*,*) "dx_Co", dx_Co
	write(*,*) "dx_Cr", dx_Cr
	write(*,*) "=========================================================="
endsubroutine pg01_display

endmodule puls_grid_opt01