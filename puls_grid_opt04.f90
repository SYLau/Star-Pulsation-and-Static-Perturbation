module puls_grid_opt04
!	pulsation grid for relativistic cowling approximation, compute nu(r)
use global_var

contains

!C2---------------------------------------------------
subroutine pg04_Ctrl
implicit none

	write(100,*) "puls_grid_opt04: interpolate for the new variable x"
	write(100,*) "- compute the grid for pulsation equations"
	write(100,*) "- includes also nu(r) for relativistic cowling puls eqt"
	write(100,*) "- multi-interface"

	call pg04_gridx
	write(100,*) "pg04_gridx: passed"

	call pg04_display
	write(100,*) "=========================================================="
endsubroutine pg04_Ctrl

!C2 A1---------------------------------------------------
subroutine pg04_gridx
implicit none
real(8), dimension(1:sp_Nt) :: x_j
integer :: n
character(len=2) :: no

	call pg04_xj(x_j)

	call pg04_xi_iniCo
	call pg04_xi_finCo
	call pg04_xi_Co(x_j)

	do n = 1, NIface-1
		write(no,'(i2)') n
		pef_ziC(n) = pef_zC(1: len(pef_zC)-4)//trim(no)//'.txt'
		pef_ziD(n) = pef_zD(1: len(pef_zD)-4)//trim(no)//'.txt'
		pef_LD_Cri1(n) = pef_LD_Cr1(1: len(pef_LD_Cr1)-4)//trim(no)//'.txt'
		pef_LD_Cri2(n) = pef_LD_Cr2(1: len(pef_LD_Cr2)-4)//trim(no)//'.txt'
		pef_LD_Cri3(n) = pef_LD_Cr3(1: len(pef_LD_Cr3)-4)//trim(no)//'.txt'
		pef_LD_Cri4(n) = pef_LD_Cr4(1: len(pef_LD_Cr4)-4)//trim(no)//'.txt'
		pef_LD_Cri5(n) = pef_LD_Cr5(1: len(pef_LD_Cr5)-4)//trim(no)//'.txt'
		call pg04_xi_iniCri(n)
		call pg04_xi_finCri(n)
		call pg04_xi_Cri(n,x_j)
	enddo

	call pg04_xi_iniOc
	call pg04_xi_finOc
	call pg04_xi_Oc(x_j)

	call pg04_ri_Ou

endsubroutine pg04_gridx

!C2 B1---------------------------------------------------
subroutine pg04_xj(x_j)
use Format_IO
use FWrite
implicit none
real(8) :: x_j(1:sp_Nt)
real(8) :: Ri, Rs
integer :: j, n

	Ri = R0 * pg_ri_frac
	Rs = R0 - (R0 - R_g)* pg_rs_frac
R0 = Rs			!!!!!!!!!!!!******************** Rs /= R0 Causes a weird mode to Arise in LD formalism (Conversion from H0 K to Z Z`; a singular point is involved)

	do j = 1, sp_Ni(1)-2
		if ((sp_r(j) - Ri)*(sp_r(j+1) - Ri) <= 0.d0) pg_ji = j
		x_j(j) = dlog(sp_r(j)/P(j))
	enddo
	if (pg_ji < 1) then
		write(*,*) "err: starting r too small, plz adjust pg_ri_frac"
		pause
	endif
	x_j(sp_Ni(1)-1) = dlog(sp_r(sp_Ni(1)-1)/P(sp_Ni(1)-1))
	XCo_i = x_j(pg_ji)
	XCo_f = x_j(sp_Ni(1)-1)

	do j = sp_Ni(1), sp_Ni(NIface)-2
		x_j(j) = dlog(sp_r(j)/P(j))
	enddo
	x_j(sp_Ni(NIface)-1) = dlog(sp_r(sp_Ni(NIface)-1)/P(sp_Ni(NIface)-1))
	do n = 1,NIface-1
		XCri_i(n) = x_j(sp_Ni(n))
		XCri_f(n) = x_j(sp_Ni(n+1)-1)
	enddo
	do j = sp_Ni(NIface), sp_Nt-1
		if ((sp_r(j) - Rs)*(sp_r(j+1) - Rs) <= 0.d0) pg_js = j+1       ! taking the point closer to star surface (a bug for LD program, r_Cr /= R0, a false mode found)
		x_j(j) = dlog(sp_r(j)/P(j))
	enddo

	x_j(sp_Nt) = dlog(sp_r(sp_Nt)/P(sp_Nt))
	XOc_i = x_j(sp_Ni(NIface))
	XOc_f = x_j(pg_js)

endsubroutine pg04_xj

!C2 B2---------------------------------------------------
subroutine pg04_xi_iniCo
implicit none

	dx_Co = (XCo_f - XCo_i)/pg_N1
	r_Co(0) = sp_r(pg_ji)
	x_Co(0) = XCo_i
	P_Co(0) = P(pg_ji)
	rho_Co(0) = rho(pg_ji)
	m_Co(0) = m(pg_ji)
	mu_Co(0) = mu(pg_ji)
	nu_Co(0) = nu(pg_ji)

endsubroutine pg04_xi_iniCo

!C2 B3---------------------------------------------------
subroutine pg04_xi_Co(x_j)
use interpolation
use Format_IO
use FWrite
implicit none
real(8), dimension(1:sp_Nt) :: x_j
integer :: ii, iiN1, mode
	mode = 4
	iiN1 = pg_N1 *2
	open(41, file='data/PG_Co.txt',status='replace')
	call Write7R8(41, FR8, r_Co(0), x_Co(0), P_Co(0), rho_Co(0), m_Co(0), mu_Co(0), nu_Co(0))
	do ii = 1, iiN1-1
		x_Co(ii) = dx_Co/2.d0 + x_Co(ii-1)
		r_Co(ii) = Int_ada(mode, pg_ji, sp_N1-1, x_j, sp_r(1:sp_Nt), x_Co(ii), "pulsation grid interpolation; core radius")
		P_Co(ii) = Int_ada(mode, pg_ji, sp_N1-1, x_j, P(1:sp_Nt), x_Co(ii), "pulsation grid interpolation; core pressure")
		rho_Co(ii) = Int_ada(mode, pg_ji, sp_N1-1, x_j, rho(1:sp_Nt), x_Co(ii), "pulsation grid interpolation; core density")
		m_Co(ii) = Int_ada(mode, pg_ji, sp_N1-1, x_j, m(1:sp_Nt), x_Co(ii), "pulsation grid interpolation; core mass")
		mu_Co(ii) = Int_ada(mode, pg_ji, sp_N1-1, x_j, mu(1:sp_Nt), x_Co(ii), "pulsation grid interpolation; core shear modulus")
		nu_Co(ii) = Int_ada(mode, pg_ji, sp_N1-1, x_j, nu(1:sp_Nt), x_Co(ii), "pulsation grid interpolation; core nu(r)")
		call Write6R8(41, FR8, r_Co(ii), x_Co(ii), P_Co(ii), rho_Co(ii), m_Co(ii), nu_Co(ii))
	enddo
	call Write7R8(41, FR8, r_Co(pg_N1 *2), x_Co(pg_N1 *2), P_Co(pg_N1 *2), rho_Co(pg_N1 *2), m_Co(pg_N1 *2), mu_Co(pg_N1 *2), nu_Co(pg_N1 *2))
	close(41)
endsubroutine pg04_xi_Co
 
!C2 B4---------------------------------------------------
subroutine pg04_xi_finCo
implicit none

	r_Co(pg_N1 *2) = sp_r(sp_Ni(1)-1)
	x_Co(pg_N1 *2) = XCo_f
	P_Co(pg_N1 *2) = P(sp_Ni(1)-1)
	rho_Co(pg_N1 *2) = rho(sp_Ni(1)-1)
	m_Co(pg_N1 *2) = m(sp_Ni(1)-1)
	mu_Co(pg_N1 *2) = mu(sp_Ni(1)-1)
	nu_Co(pg_N1 *2) = nu(sp_Ni(1)-1)

endsubroutine pg04_xi_finCo

!C2 B5---------------------------------------------------
subroutine pg04_xi_iniCri(i)
implicit none
integer :: i
	dx_Cri(i) = (XCri_f(i) - XCri_i(i))/pg_Ni
	r_Cri(i,0) = sp_r(sp_Ni(i))
	x_Cri(i,0) = XCri_i(i)
	P_Cri(i,0) = P(sp_Ni(i))
	rho_Cri(i,0) = rho(sp_Ni(i))
	m_Cri(i,0) = m(sp_Ni(i))
	mu_Cri(i,0) = mu(sp_Ni(i))
	nu_Cri(i,0) = nu(sp_Ni(i))
endsubroutine pg04_xi_iniCri

!C2 B6---------------------------------------------------
subroutine pg04_xi_Cri(i, x_j)
use interpolation
use Format_IO
use FWrite
implicit none
real(8), dimension(1:sp_Nt) :: x_j
integer :: ii, iiN2, i, mode
character(len=2) :: no
	
	mode = 4
	iiN2 = pg_Ni*2
	write( no, '(i2)' ) i

	open(42, file='data/PG_Cr'//trim(no)//'.txt',status='replace')
	call Write7R8(42, FR8, r_Cri(i,0), x_Cri(i,0), P_Cri(i,0), rho_Cri(i,0), m_Cri(i,0), mu_Cri(i,0), nu_Cri(i,0))

	do ii = 1, iiN2-1
		x_Cri(i,ii) = dx_Cri(i)/2.d0 + x_Cri(i,ii-1)
		r_Cri(i,ii) = Int_ada(mode, sp_Ni(i), sp_Ni(i+1)-1, x_j, sp_r(1:sp_Nt), x_Cri(i,ii), "pulsation grid interpolation; crust radius")
		P_Cri(i,ii) = Int_ada(mode, sp_Ni(i), sp_Ni(i+1)-1, x_j, P(1:sp_Nt), x_Cri(i,ii), "pulsation grid interpolation; crust pressure")
		rho_Cri(i,ii) = Int_ada(mode, sp_Ni(i), sp_Ni(i+1)-1, x_j, rho(1:sp_Nt), x_Cri(i,ii), "pulsation grid interpolation; crust density")
		m_Cri(i,ii) = Int_ada(mode, sp_Ni(i), sp_Ni(i+1)-1, x_j, m(1:sp_Nt), x_Cri(i,ii), "pulsation grid interpolation; crust mass")
		mu_Cri(i,ii) = Int_ada(mode, sp_Ni(i), sp_Ni(i+1)-1, x_j, mu(1:sp_Nt), x_Cri(i,ii), "pulsation grid interpolation; crust shear modulus")
		nu_Cri(i,ii) = Int_ada(mode, sp_Ni(i), sp_Ni(i+1)-1, x_j, nu(1:sp_Nt), x_Cri(i,ii), "pulsation grid interpolation; crust nu(r)")
		call Write7R8(42, FR8, r_Cri(i,ii), x_Cri(i,ii), P_Cri(i,ii), rho_Cri(i,ii), m_Cri(i,ii), mu_Cri(i,ii), nu_Cri(i,ii))
	enddo
	call Write7R8(42, FR8, r_Cri(i,pg_Ni*2), x_Cri(i,pg_Ni*2), P_Cri(i,pg_Ni*2), rho_Cri(i,pg_Ni*2), m_Cri(i,pg_Ni*2), mu_Cri(i,pg_Ni*2), nu_Cri(i,pg_Ni*2))
	close(42)
endsubroutine pg04_xi_Cri

!C2 B7---------------------------------------------------
subroutine pg04_xi_finCri(i)
implicit none
integer :: i
	r_Cri(i,pg_Ni*2) = sp_r(sp_Ni(i+1)-1)
	x_Cri(i,pg_Ni*2) = XCri_f(i)
	P_Cri(i,pg_Ni*2) = P(sp_Ni(i+1)-1)
	rho_Cri(i,pg_Ni*2) = rho(sp_Ni(i+1)-1)
	m_Cri(i,pg_Ni*2) = m(sp_Ni(i+1)-1)
	mu_Cri(i,pg_Ni*2) = mu(sp_Ni(i+1)-1)
	nu_Cri(i,pg_Ni*2) = nu(sp_Ni(i+1)-1)
endsubroutine pg04_xi_finCri
!C2 B8---------------------------------------------------
subroutine pg04_xi_iniOc
implicit none

	dx_Oc = (XOc_f - XOc_i)/pg_N3
	r_Oc(0) = sp_r(sp_N2)
	x_Oc(0) = XOc_i
	P_Oc(0) = P(sp_N2)
	rho_Oc(0) = rho(sp_N2)
	m_Oc(0) = m(sp_N2)
	nu_Oc(0) = nu(sp_N2)

endsubroutine pg04_xi_iniOc
 
!C2 B9---------------------------------------------------
subroutine pg04_xi_Oc(x_j)
use interpolation
use Format_IO
use FWrite
implicit none
real(8), dimension(1:sp_Nt) :: x_j
integer :: ii, iiN3, mode
	mode = 4
	iiN3 = pg_N3*2
	open(43, file='data/PG_Oc.txt',status='replace')
	call Write6R8(43, FR8, r_Oc(0), x_Oc(0), P_Oc(0), rho_Oc(0), m_Oc(0), nu_Oc(0))

	do ii = 1, iiN3-1
		x_Oc(ii) = dx_Oc/2.d0 + x_Oc(ii-1)
		r_Oc(ii) = Int_ada(mode, sp_N2, pg_js, x_j, sp_r(1:sp_Nt), x_Oc(ii), "pulsation grid interpolation; ocean radius")
		P_Oc(ii) = Int_ada(mode, sp_N2, pg_js, x_j, P(1:sp_Nt), x_Oc(ii), "pulsation grid interpolation; ocean pressure")
		rho_Oc(ii) = Int_ada(mode, sp_N2, pg_js, x_j, rho(1:sp_Nt), x_Oc(ii), "pulsation grid interpolation; ocean density")
		m_Oc(ii) = Int_ada(mode, sp_N2, pg_js, x_j, m(1:sp_Nt), x_Oc(ii), "pulsation grid interpolation; ocean mass")
		nu_Oc(ii) = Int_ada(mode, sp_N2, pg_js, x_j, nu(1:sp_Nt), x_Oc(ii), "pulsation grid interpolation; ocean nu(r)")
		call Write6R8(43, FR8, r_Oc(ii), x_Oc(ii), P_Oc(ii), rho_Oc(ii), m_Oc(ii), nu_Oc(ii))
	enddo
	call Write6R8(43, FR8, r_Oc(pg_N3*2), x_Oc(pg_N3*2), P_Oc(pg_N3*2), rho_Oc(pg_N3*2), m_Oc(pg_N3*2), nu_Oc(pg_N3*2))
	close(43)
endsubroutine pg04_xi_Oc

!C2 B10---------------------------------------------------
subroutine pg04_xi_finOc
implicit none

	r_Oc(pg_N3*2) = sp_r(pg_js)
	x_Oc(pg_N3*2) = XOc_f
	P_Oc(pg_N3*2) = P(pg_js)
	rho_Oc(pg_N3*2) = rho(pg_js)
	m_Oc(pg_N3*2) = m(pg_js)
	nu_Oc(pg_N3*2) = nu(pg_js)

endsubroutine pg04_xi_finOc

!C2 B11---------------------------------------------------
subroutine pg04_ri_Ou
use Format_IO
use FWrite
implicit none
integer :: i
	P_Ou = P(pg_js)
	rho_Ou = rho(pg_js)
	m_Ou = m(pg_js)
	nu_Ou = nu(pg_js)
	ROu_i = R0
endsubroutine pg04_ri_Ou

!C2 A2---------------------------------------------------
subroutine pg04_display
implicit none
	write(*,*) "Pulsation Grid Parameters:"
	write(*,*) "XCo_i", XCo_i
	write(*,*) "XCo_f", XCo_f
	write(*,*) "XCr_i", XCr_i
	write(*,*) "XCr_f", XCr_f
	write(*,*) "XOc_i", XOc_i
	write(*,*) "XOc_f", XOc_f
	write(*,*) "dx_Co", dx_Co
	write(*,*) "dx_Cr", dx_Cr
	write(*,*) "dx_Oc", dx_Oc
	write(*,*) "=========================================================="
endsubroutine pg04_display

endmodule puls_grid_opt04