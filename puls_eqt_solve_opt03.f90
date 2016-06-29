module puls_eqt_solve_opt03
!	3 components model
use global_var

contains
!D1---------------------------------------------------
subroutine pe03_Ctrl(beta, gamma, wFile)
implicit none
complex(8) :: beta, gamma
logical :: wFile

	if (NMat /= 4) then
		write(*,*) "err: pe03_Ctrl; NMat mismatch"
		pause
	endif
	if (pe_fcall == .true.) then
		write(100,*) "puls_eqt_solve_opt03: linear coupled ode solver"
		write(100,*) "- 1 component neutron star model"
		write(100,*) "- Fully relativistic formalism; LD2 formalism"
		if (pes_opt == 3) then
			write(*,*) "LD2 pulsation equations"
			write(*,*) "=========================================================="
		else
			write(*,*) "err: pes_opt"
			pause
		endif
	endif

		call pe03_iterate(beta, gamma, wFile)

	if (pe_fcall == .true.) then
		write(100,*) "pe03_iterate: passed"
		write(100,*) "=========================================================="
		pe_fcall = .false.
	endif

endsubroutine pe03_Ctrl

!D1 A1---------------------------------------------------
subroutine pe03_iterate(beta, gamma, wFile)
use Format_IO
use FWrite
use puls_eqt_set_opt03
use puls_eqt_set_opt03_ex1
implicit none
complex(8), dimension(1: Neqt) :: x, y, y_R
complex(8), dimension(1:Neqt+1) :: coef
complex(8) ::  H0, K, W, W_mid(1:5), W_norm,  beta, gamma, V(1:5), E(1:2,1:4), E1(1:2,1:4), E2(1:2,1:4)
complex(8) :: h0_Love_R0, h0_Love_R1, dh0_Love
real(8) :: r_Ou
integer :: i, j, ii, n, it, mode
logical :: wFile
complex(8), allocatable, dimension(:, :) :: cY
allocate(cY(1:NMat, 1:NMat+1))

	h0_Love_R0 = (0.d0, 0.d0)
	h0_Love_R1 = (0.d0, 0.d0)

	!¡¸Within the star, solutions represented by vectors: core cY1, cY2; fluid crust cY3-cY5
	do it = 1, 5
		if (it == 1 .or. it == 2) then
			if (it == 1) then 
				mode = 1
				if (wFile == .true.) open(01, file = pef_LD_Co1, status = 'replace')
			elseif (it == 2) then
				mode = 2
				if (wFile == .true.) open(01, file = pef_LD_Co2, status = 'replace')
			endif

			call pe03_bc_Coi(mode, r_Co(0), y)

			if (wFile == .true.) call Write1R82AR8(01, FR8,r_Co(0), dreal(CoeC(it)*y(1:4)), dimag(CoeC(it)*y(1:4)), m=4)
			do i = 0, pg_N1-1
				ii = i*2
				x(1:4) = y(1:4)
				call pe03_it_Co(x_Co(ii), x_Co(ii+2), x(1:4), y(1:4))
				if (wFile == .true.) call Write1R82AR8(01, FR8,r_Co(ii+2), dreal(CoeC(it)*y(1:4)), dimag(CoeC(it)*y(1:4)), m=4)

if (  (r_Co(ii)- R0/2.d0)*(r_Co(ii+2) -R0/2.d0) <= 0.d0) then
W_mid(it) = y(3)
endif
			enddo

			call pe03_bc_Cof(y, cY(1:4 ,it))



		elseif (it == 3 .or. it == 4 .or. it == 5) then
			if (it == 3) then 
				mode = 1
				if (wFile == .true.) open(01, file = pef_LD_Cr1, status = 'replace')
			elseif (it == 4) then
				mode = 2
				if (wFile == .true.) open(01, file = pef_LD_Cr2, status = 'replace')
			elseif (it == 5) then
				mode = 3
				if (wFile == .true.) open(01, file = pef_LD_Cr3, status = 'replace')
			endif
			call pe03_bc_Crf(mode, y)
			if (wFile == .true.) call Write1R82AR8(01, FR8,r_Cr(2*pg_N2), dreal(CoeC(it)*y(1:4)), dimag(CoeC(it)*y(1:4)), m=4)
			do i = pg_N2, 1, -1
				ii = i*2
				x(1:4) = y(1:4)
				call pe03_it_Cr(x_Cr(ii), x_Cr(ii-2), x(1:4), y(1:4))
				if (wFile == .true.) call Write1R82AR8(01, FR8,r_Cr(ii-2), dreal(CoeC(it)*y(1:4)), dimag(CoeC(it)*y(1:4)), m=4)
if (  (r_Cr(ii-2)- R0/2.d0)*(r_Cr(ii) -R0/2.d0) <= 0.d0) then
W_mid(it) = y(3)
endif
			
			if (i == pg_N2) then 
				call pes03_E_Fluid(P_Cr(2*pg_N2), rho_Cr(2*pg_N2), m_Cr(2*pg_N2), nu_Cr(2*pg_N2), r_Cr(2*pg_N2), E)
				h0_Love_R0 = h0_Love_R0+  CoeC(it)* (E(1,1)*y(1) + E(1,2)*y(2) + E(1,3)*y(3) + E(1,4)*y(4))
			elseif (i == pg_N2-1) then
				call pes03_E_Fluid(P_Cr(2*pg_N2-2), rho_Cr(2*pg_N2-2), m_Cr(2*pg_N2-2), nu_Cr(2*pg_N2-2), r_Cr(2*pg_N2-2), E)
				h0_Love_R1 = h0_Love_R1+  CoeC(it)* (E(1,1)*y(1) + E(1,2)*y(2) + E(1,3)*y(3) + E(1,4)*y(4))
			endif

			enddo
			call pe03_bc_Cri(y, cY(1:4 ,it))
		endif
		if (wFile == .true.) close(01)
	enddo

	call pe03_cY_solve(cY(1:4 ,1:5), coef(1:5))
	CoeC(1:5) = coef(1:5)
	Y_R(1:3) = coef(3:5)
	Y_R(4) = 0.d0
	

	!¡¸ Calculating Love Number
	dh0_Love = (h0_Love_R0-h0_Love_R1)/(r_Cr(2*pg_N2)-r_Cr(2*pg_N2-2))
	call pe03_k2(h0_Love_R0, dh0_Love, M0, R0, k2)


	!¡¸Solve algebraric equations from BCs to obtain H0 and K
	call pes03_E_Fluid(P_Ou, rho_Ou, m_Ou, nu_Ou, ROu_i, E)
	H0 = (0.d0, 0.d0)
	do i = 1,4
		H0 = H0 + E(1,i) * Y_R(i) 
	enddo
	K = Y_R(2)
	W = coef(5)

W_norm = W_mid(1)*coef(1) + W_mid(2)*coef(2) + W_mid(3)*coef(3) + W_mid(4)*coef(4) + W_mid(5)*coef(5)


!¡¸Normalize H0, K by W at surface
H0 = H0/W_norm
K = K/W_norm

	if (zero_freq == .true.) return ! Skip the outter integration for zero frequency


!write(*,*) "coef(1)", coef(1)
!write(*,*) "coef(2)", coef(2)
!write(*,*) "coef(3)", coef(3)
!write(*,*) "coef(4)", coef(4)
!write(*,*) "coef(5)", coef(5)
!V(1:5) = 0.d0
!call pes03_E_Fluid(P_Co(2*pg_N1), rho_Co(2*pg_N1), m_Co(2*pg_N1), nu_Co(2*pg_N1), r_Co(2*pg_N1), E1)
!call pes03_E_Fluid(P_Cr(0), rho_Cr(0), m_Cr(0), nu_Cr(0), r_Cr(0), E2)
!do i = 1,4
!	V(1:2) = V(1:2) + E1(2,i) * cY(i, 1:2) 
!	V(3:5) = V(3:5) + E2(2,i) * cY(i, 3:5) 
!enddo
!write(*,*) "Left", coef(1) * real(V(1)) + coef(2) * real(V(2))
!write(*,*) "Right", coef(3) * real(V(3)) + coef(4) * real(V(4)) + coef(5) * real(V(5))
!write(*,*) "Left", coef(1) * real(cY(1:4, 1)) + coef(2) * real(cY(1:4, 2))
!write(*,*) "Right", coef(3) * real(cY(1:4, 3)) + coef(4) * real(cY(1:4, 4)) + coef(5) * real(cY(1:4, 5))
!pause

!open(02, file = 'data/debug_Ze_Series.txt')

	!	Outside the Neutron Star
	if (wFile == .true.) open(01, file = pef_LD_Ou, status = 'replace')

	call pe03_bc_Oui(H0, K, y(1:2))
!call Write2R8(1234, FR8, dreal(cdsqrt(OmeC_sq*c**2)/2.d0/pi), dreal(y(1)/y(2))) ! Open in eigen_freq_opt02
	ROu_f = R_inf_f/cdabs(cdsqrt(OmeC_sq))
	ROu_i = R0
	r_Ou = ROu_i
	dr_Ou = (ROu_f - ROu_i)/pg_NOu
	if (wFile == .true.) call Write1R82AR8(01, FR8,ROu_i*dreal(cdsqrt(OmeC_sq))/pi, dreal(y(1:2)), dimag(y(1:2)), m=2)
	!r normalized by pi/omega, to visualize phase angle

	do i = 0, pg_NOu -1
		r_Ou = i * dr_Ou + ROu_i
		x(1:2) = y(1:2)
		call pe03_it_Ou(r_Ou, r_Ou + dr_Ou, x(1:2), y(1:2))
		if (wFile == .true.) call Write1R82AR8(01, FR8,(r_Ou + dr_Ou)*dreal(cdsqrt(OmeC_sq))/pi, dreal(y(1:2)), dimag(y(1:2)), m=2)
!call Write3R8(02, FR8,r_Ou + dr_Ou, dreal(pes03_series(r_Ou + dr_Ou)), dimag(pes03_series(r_Ou + dr_Ou)))
	enddo
	r_Ou = ROu_f

	call pe03_bc_Ouf(ROu_f, y(1:2), beta, gamma)

	if (wFile == .true.) close(01)

!write(*,*) OmeC_sq** 0.5d0
!write(*,*) dreal(-dcmplx(0.d0, 1.d0)*OmeC_sq** 0.5d0 * ROu_f)
!write(*,*) cdabs(cdexp(-dcmplx(0.d0, 1.d0)*OmeC_sq** 0.5d0 * ROu_f ))
!write(*,*) beta, gamma
!pause
!close(02)

endsubroutine pe03_iterate

!D1 B1---------------------------------------------------
subroutine pe03_it_Co(xi, xf, xc, yc)
use puls_eqt_set_opt03
use RK4_Set
implicit none
complex(8), dimension(1:) :: xc, yc
real(8), dimension(1:2*size(yc)) :: x, y
real(8) :: xi, xf
integer :: dim
	dim = size(yc)
	x(1:dim) = dreal(xc(1:dim))
	x(dim+1: 2*dim) = dimag(xc(1:dim))
	
	call RK4(pes03_vx_Co, 8, xi, xf, x, y)
	yc(1:dim) = dcmplx(y(1:dim), y(dim+1: 2*dim))
endsubroutine pe03_it_Co

!D1 B2---------------------------------------------------
subroutine pe03_it_Cr(xi, xf, xc, yc)
use puls_eqt_set_opt03
use RK4_Set
implicit none
complex(8), dimension(1:) :: xc, yc
real(8), dimension(1:2*size(yc)) :: x, y
real(8) :: xi, xf
integer :: dim
	dim = size(yc)
	x(1:dim) = dreal(xc(1:dim))
	x(dim+1: 2*dim) = dimag(xc(1:dim))

	call RK4(pes03_vx_Cr_FM, 8, xi, xf, x, y)
	yc(1:dim) = dcmplx(y(1:dim), y(dim+1: 2*dim))
endsubroutine pe03_it_Cr

!D1 B3---------------------------------------------------
subroutine pe03_it_Ou(xi, xf, xc, yc)
use puls_eqt_set_opt03
use RK4_Set
implicit none
complex(8), dimension(1:) :: xc, yc
real(8), dimension(1:2*size(yc)) :: x, y
real(8) :: xi, xf
integer :: dim
	dim = size(yc)
	x(1:dim) = dreal(xc(1:dim))
	x(dim+1: 2*dim) = dimag(xc(1:dim))
	
	call RK4(pes03_vr_Ou, 4, xi, xf, x, y)
	yc(1:dim) = dcmplx(y(1:dim), y(dim+1: 2*dim))
endsubroutine pe03_it_Ou

!D1 B4---------------------------------------------------
subroutine pe03_bc_Coi(mode, r_, yc)
use puls_eqt_set_opt03
implicit none
real(8) :: r_
complex(8) :: yc(1:4)
integer :: mode
	call pes03_bc_Coi(mode, r_, yc)
endsubroutine pe03_bc_Coi

!D1 B5---------------------------------------------------
subroutine pe03_bc_Cof(yc, cY)
use puls_eqt_set_opt03
implicit none
complex(8), dimension(1:4) :: cY
complex(8) :: yc(1:4)
	call pes03_bc_Co_FM_f(yc, cY)
endsubroutine pe03_bc_Cof

!D1 B6---------------------------------------------------
subroutine pe03_bc_Cri(yc, cY)
use puls_eqt_set_opt03
implicit none
complex(8), dimension(1:4) :: cY
complex(8) :: yc(1:4)
	call pes03_bc_Cr_FM_i(yc, cY)
endsubroutine pe03_bc_Cri

!D1 B7---------------------------------------------------
subroutine pe03_bc_Crf(mode, yc)
use puls_eqt_set_opt03
implicit none
complex(8) :: yc(1:4)
integer :: mode
	call pes03_bc_Cr_FM_f(mode,yc)
endsubroutine pe03_bc_Crf

!D1 B8---------------------------------------------------
subroutine pe03_cY_solve(cY, coef)
!use imsl
use lin_sol_gen_int
use Eigenvalues_Eigenvectors
implicit none
complex(8) :: cY(1:4 ,1:5), coef(1:5), M(1:4, 1:4), Vec(1:4,1), X(1:4,1)
real(8) :: err_
	coef(1) = (1.d0, 0.d0)
	M(1:4,1) = -cY(1:4, 2)
	M(1:4,2:4) = cY(1:4, 3:5)

	call LinSys_Cramer_c(M(1:4, 1:4), cY(1:4,1), coef(2:5), 4)
!	Check accuracy of Crammer Rule
!Vec(1:4,1) = cY(1:4,1)
!call LIN_SOL_GEN(M, Vec, X)
!coef(2:5) = X(1:4,1)
!write(*,*) 'Coef',coef(2:5)
!write(*,*) 'X',X(1:4,1)
!pause
!	Nearly singular for some modes
!write(*,*) -coef(2)*cY(1:4,2) + coef(3)*cY(1:4,3) + coef(4)*cY(1:4,4) + coef(5)*cY(1:4,5)
!write(*,*) "&&&&&&&&&&&&&&&&&&&&&&&&&"
!write(*,*) cY(1:4,1)
!pause
endsubroutine pe03_cY_solve

!D1 B9---------------------------------------------------
subroutine pe03_bc_Oui(H0, K, yc)
use puls_eqt_set_opt03
implicit none
complex(8) :: H0, K, yc(1:2)
	call pes03_bc_Oui(H0, K, yc(1:2))
endsubroutine pe03_bc_Oui

!D1 B10---------------------------------------------------
subroutine pe03_bc_Ouf(r_, yc, be, ga)
use puls_eqt_set_opt03
implicit none
complex(8) :: yc(1:2), be, ga
real(8) :: r_
	if (gr_opt == 1) then
		call pes03_bc_Ouf_LD(r_, yc, be, ga)
	elseif (gr_opt == 2) then
		call pes03_bc_Ouf_Scatter(r_, yc, be, ga)
	elseif (gr_opt == 3) then
		call pes03_bc_Ouf_phase(r_, yc, be, ga)
	else
		write(*,*) "err: pe gr_opt"
		pause
	endif
endsubroutine pe03_bc_Ouf

!Tidal Love Number
subroutine pe03_k2(h0_, dh0_, m_, r_, k2_)
implicit none
complex(8) :: h0_, dh0_, k2_
real(8) :: m_, r_
complex(8) :: yR_, T_(1:3)
real(8) :: C_

	yR_ = r_*dh0_/h0_
	C_ = m_/r_
	T_(1) = 8.d0/5.d0*C_**5*(1.d0-2.d0*C_)**2 * (2.d0*C_*(yR_-1.d0)+2.d0-yR_)
	T_(2) = 2.d0*C_*( 4.d0*(yR_+1.d0)*C_**4 + (6.d0*yR_ -4.d0)*C_**3 + (26.d0 - 22.d0*yR_)*C_**2 + 3.d0*C_*(5.d0*yR_ - 8.d0) -3.d0*yR_ + 6.d0 )
	T_(3) = 3.d0*(1.d0-2.d0*C_)**2 * (2.d0*C_*(yR_ - 1.d0) - yR_ +2.d0) * dlog(1.d0-2.d0*C_)

	k2_ = T_(1) /(T_(2) + T_(3) )
	
endsubroutine pe03_k2
endmodule puls_eqt_solve_opt03