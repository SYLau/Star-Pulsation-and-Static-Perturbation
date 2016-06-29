module puls_eqt_solve_opt09
!	1 components model
use global_var

contains
!D1---------------------------------------------------
subroutine pe09_Ctrl(beta, gamma, wFile)
implicit none
complex(8) :: beta, gamma
logical :: wFile

	if (NMat /= 6) then
		write(*,*) "err: pe09_Ctrl; NMat mismatch"
		pause
	endif
	if (pe_fcall == .true.) then
		write(100,*) "puls_eqt_solve_opt09: linear coupled ode solver"
		write(100,*) "- 1 components solid quark star model"
		write(100,*) "- Fully relativistic formalism; LD2 formalism"
		if (pes_opt == 5) then
			write(*,*) "LD2/ Andersson pulsation equations"
			write(*,*) "=========================================================="
		else
			write(*,*) "err: pes_opt"
			pause
		endif
	endif

		call pe09_iterate(beta, gamma, wFile)

	if (pe_fcall == .true.) then
		write(100,*) "pe09_iterate: passed"
		write(100,*) "=========================================================="
		pe_fcall = .false.
	endif

endsubroutine pe09_Ctrl

!D1 A1---------------------------------------------------
subroutine pe09_iterate(beta, gamma, wFile)
use Format_IO
use FWrite
use puls_eqt_set_opt05
use puls_eqt_set_opt05_ex1
implicit none
complex(8) :: E(1:3, 1:6), X_R
real(8) :: gamma_1
complex(8), dimension(1:2) :: y
complex(8), dimension(1:6) :: x, z
complex(8), dimension(1:8) :: coef
complex(8) :: W, V, H0, K, beta, gamma, E1(1:2,1:4), E2(1:2,1:4)
complex(8) :: W_mid(1:2)
real(8) :: r_Ou
integer :: i, j, ii, n, it, mode
logical :: wFile
complex(8), dimension(1:6, 1:8) :: cY_
!complex(8), allocatable, dimension(:, :) :: cY_
!allocate(cY_(1:NMat, 1:7))

	!��Within the star, solutions represented by vectors: core cY1, cY2; crust cY3-cY7
	do it = 1, 8
		if (it == 1 .or. it == 2 .or. it == 3) then
			if (it == 1) then 
				mode = 1
				if (wFile == .true.) open(01, file = pef_LD_Co1, status = 'replace')
			elseif (it == 2) then
				mode = 2
				if (wFile == .true.) open(01, file = pef_LD_Co2, status = 'replace')
			elseif (it == 3) then
				mode = 3
				if (wFile == .true.) open(01, file = pef_LD_Co3, status = 'replace')
			endif

			call pe09_bc_Coi(mode, r_Co(0), z)

			if (wFile == .true.) call Write1R82AR8(01, FR8,r_Co(0), dreal(CoeC(it)*z(1:6)), dimag(CoeC(it)*z(1:6)), m=6)
			do i = 0, pg_N1-1
				ii = i*2
				x(1:6) = z(1:6)
				call pe09_it_Co(x_Co(ii), x_Co(ii+2), x(1:6), z(1:6))
				if (wFile == .true.) call Write1R82AR8(01, FR8,r_Co(ii+2), dreal(CoeC(it)*z(1:6)), dimag(CoeC(it)*z(1:6)), m=6)
			enddo

			call pe09_bc_Cof(z(1:6), cY_(1:6 ,it))

		elseif ( it == 4 .or. it == 5 .or. it == 6 .or. it == 7 .or. it == 8) then
			if (it == 4) then 
				mode = 1
				if (wFile == .true.) open(01, file = pef_LD_Cr1, status = 'replace')
			elseif (it == 5) then
				mode = 2
				if (wFile == .true.) open(01, file = pef_LD_Cr2, status = 'replace')
			elseif (it == 6) then
				mode = 3
				if (wFile == .true.) open(01, file = pef_LD_Cr3, status = 'replace')
			elseif (it == 7) then
				mode = 4
				if (wFile == .true.) open(01, file = pef_LD_Cr4, status = 'replace')
			elseif (it == 8) then
				mode = 5
				if (wFile == .true.) open(01, file = pef_LD_Cr5, status = 'replace')
			endif
			call pe09_bc_Crf(mode, z)

			if (wFile == .true.) call Write1R82AR8(01, FR8,r_Cr(2*pg_N2), dreal(CoeC(it)*z(1:6)), dimag(CoeC(it)*z(1:6)), m=6)
			do i = pg_N2, 1, -1
				ii = i*2
				x(1:6) = z(1:6)
				call pe09_it_Cr(x_Cr(ii), x_Cr(ii-2), x(1:6), z(1:6))
				if (wFile == .true.) call Write1R82AR8(01, FR8,r_Cr(ii-2), dreal(CoeC(it)*z(1:6)), dimag(CoeC(it)*z(1:6)), m=6)
			enddo
			call pe09_bc_Cri(z(1:6), cY_(1:6 ,it))
		endif
		if (wFile == .true.) close(01)
	enddo

	!��Solve algebraric equations from BCs to obtain H0 and K
	call pe09_cY_solve(cY_(1:6 ,1:8), coef(1:8))
	CoeC(1:8) = coef(1:8)
	H0 = coef(6)/R0**2
	K = coef(5)/R0**2
	W = coef(7)
	V = coef(8)

	!��Normalize H0, K by W at surface
	H0 = H0/W
	K = K/W

	!	Outside the Neutron Star
	if (wFile == .true.) open(01, file = pef_LD_Ou, status = 'replace')

	call pe09_bc_Oui(H0, K, y(1:2))
	ROu_f = R_inf_f/cdabs(cdsqrt(OmeC_sq))
	ROu_i = R0
	r_Ou = ROu_i
	dr_Ou = (ROu_f - ROu_i)/pg_NOu
	
	if (wFile == .true.) call Write1R82AR8(01, FR8,ROu_i*dreal(cdsqrt(OmeC_sq))/pi, dreal(y(1:2)), dimag(y(1:2)), m=2)	
	!r normalized by pi/omega, to visualize phase angle

	do i = 0, pg_NOu -1
		r_Ou = i * dr_Ou + ROu_i
		x(1:2) = y(1:2)
		call pe09_it_Ou(r_Ou, r_Ou + dr_Ou, x(1:2), y(1:2))
		if (wFile == .true.) call Write1R82AR8(01, FR8,(r_Ou + dr_Ou)*dreal(cdsqrt(OmeC_sq))/pi, dreal(y(1:2)), dimag(y(1:2)), m=2)
	enddo

	r_Ou = ROu_f

	call pe09_bc_Ouf(ROu_f, y(1:2), beta, gamma)
	
	if (wFile == .true.) close(01)
	!deallocate(cY_)

endsubroutine pe09_iterate


subroutine pe09_iterate2(beta, gamma, wFile)
use Format_IO
use FWrite
use puls_eqt_set_opt05
use puls_eqt_set_opt05_ex1
implicit none
complex(8) :: E(1:3, 1:6), X_R
real(8) :: gamma_1
complex(8), dimension(1:2) :: y
complex(8), dimension(1:6) :: x, z
complex(8), dimension(1:7) :: coef
complex(8) :: W, V, H0, K, beta, gamma, E1(1:2,1:4), E2(1:2,1:4)
complex(8) :: W_mid(1:2)
real(8) :: r_Ou
integer :: i, j, ii, n, it, mode
logical :: wFile
complex(8), dimension(1:6, 1:7) :: cY_
!complex(8), allocatable, dimension(:, :) :: cY_
!allocate(cY_(1:NMat, 1:7))

	!��Within the star, solutions represented by vectors: core cY1, cY2; crust cY3-cY7
	do it = 1, 7
		if (it == 1 .or. it == 2) then
			if (it == 1) then 
				mode = 1
				if (wFile == .true.) open(01, file = pef_LD_Co1, status = 'replace')
			elseif (it == 2) then
				mode = 2
				if (wFile == .true.) open(01, file = pef_LD_Co2, status = 'replace')
			endif

			call pe09_bc_Coi(mode, r_Co(0), z)

			if (wFile == .true.) call Write1R82AR8(01, FR8,r_Co(0), dreal(CoeC(it)*z(1:6)), dimag(CoeC(it)*z(1:6)), m=6)
			do i = 0, pg_N1-1
				ii = i*2
				x(1:6) = z(1:6)
				call pe09_it_Co(x_Co(ii), x_Co(ii+2), x(1:6), z(1:6))
				if (wFile == .true.) call Write1R82AR8(01, FR8,r_Co(ii+2), dreal(CoeC(it)*z(1:6)), dimag(CoeC(it)*z(1:6)), m=6)
			enddo

			call pe09_bc_Cof(z(1:6), cY_(1:6 ,it))

		elseif ( it == 3 .or. it == 4 .or. it == 5 .or. it == 6 .or. it == 7) then
			if (it == 3) then 
				mode = 1
				if (wFile == .true.) open(01, file = pef_LD_Cr1, status = 'replace')
			elseif (it == 4) then
				mode = 2
				if (wFile == .true.) open(01, file = pef_LD_Cr2, status = 'replace')
			elseif (it == 5) then
				mode = 3
				if (wFile == .true.) open(01, file = pef_LD_Cr3, status = 'replace')
			elseif (it == 6) then
				mode = 4
				if (wFile == .true.) open(01, file = pef_LD_Cr4, status = 'replace')
			elseif (it == 7) then
				mode = 5
				if (wFile == .true.) open(01, file = pef_LD_Cr5, status = 'replace')
			endif
			call pe09_bc_Crf(mode, z)

			if (wFile == .true.) call Write1R82AR8(01, FR8,r_Cr(2*pg_N2), dreal(CoeC(it)*z(1:6)), dimag(CoeC(it)*z(1:6)), m=6)
			do i = pg_N2, 1, -1
				ii = i*2
				x(1:6) = z(1:6)
				call pe09_it_Cr(x_Cr(ii), x_Cr(ii-2), x(1:6), z(1:6))
				if (wFile == .true.) call Write1R82AR8(01, FR8,r_Cr(ii-2), dreal(CoeC(it)*z(1:6)), dimag(CoeC(it)*z(1:6)), m=6)
			enddo
			call pe09_bc_Cri(z(1:6), cY_(1:6 ,it))
		endif
		if (wFile == .true.) close(01)
	enddo

	!��Solve algebraric equations from BCs to obtain H0 and K
	call pe09_cY_solve2(cY_(1:6 ,1:7), coef(1:7))
	CoeC(1:7) = coef(1:7)
	H0 = coef(5)/R0**2
	K = coef(4)/R0**2
	W = coef(6)
	V = coef(7)

	!��Normalize H0, K by W at surface
	H0 = H0/W
	K = K/W

	!	Outside the Neutron Star
	if (wFile == .true.) open(01, file = pef_LD_Ou, status = 'replace')

	call pe09_bc_Oui(H0, K, y(1:2))
	ROu_f = R_inf_f/cdabs(cdsqrt(OmeC_sq))
	ROu_i = R0
	r_Ou = ROu_i
	dr_Ou = (ROu_f - ROu_i)/pg_NOu
	
	if (wFile == .true.) call Write1R82AR8(01, FR8,ROu_i*dreal(cdsqrt(OmeC_sq))/pi, dreal(y(1:2)), dimag(y(1:2)), m=2)	
	!r normalized by pi/omega, to visualize phase angle

	do i = 0, pg_NOu -1
		r_Ou = i * dr_Ou + ROu_i
		x(1:2) = y(1:2)
		call pe09_it_Ou(r_Ou, r_Ou + dr_Ou, x(1:2), y(1:2))
		if (wFile == .true.) call Write1R82AR8(01, FR8,(r_Ou + dr_Ou)*dreal(cdsqrt(OmeC_sq))/pi, dreal(y(1:2)), dimag(y(1:2)), m=2)
	enddo

	r_Ou = ROu_f

	call pe09_bc_Ouf(ROu_f, y(1:2), beta, gamma)
	
	if (wFile == .true.) close(01)
	!deallocate(cY_)

endsubroutine pe09_iterate2
!D1 B1---------------------------------------------------
subroutine pe09_it_Co(xi, xf, xc, yc)
use puls_eqt_set_opt05
use RK4_Set
implicit none
complex(8), dimension(1:) :: xc, yc
real(8), dimension(1:2*size(yc)) :: x, y
real(8) :: xi, xf
integer :: dim
	dim = size(yc)
	x(1:dim) = dreal(xc(1:dim))
	x(dim+1: 2*dim) = dimag(xc(1:dim))
	
	call RK4(pes05_vx_Co, 12, xi, xf, x, y)
	yc(1:dim) = dcmplx(y(1:dim), y(dim+1: 2*dim))
endsubroutine pe09_it_Co

!D1 B2---------------------------------------------------
subroutine pe09_it_Cr(xi, xf, xc, yc)
use puls_eqt_set_opt05
use RK4_Set
implicit none
complex(8), dimension(1:) :: xc, yc
real(8), dimension(1:2*size(yc)) :: x, y
real(8) :: xi, xf
integer :: dim
	dim = size(yc)
	x(1:dim) = dreal(xc(1:dim))
	x(dim+1: 2*dim) = dimag(xc(1:dim))

	call RK4(pes05_vx_Cr, 12, xi, xf, x, y)
	yc(1:dim) = dcmplx(y(1:dim), y(dim+1: 2*dim))
endsubroutine pe09_it_Cr

!D1 B3---------------------------------------------------
subroutine pe09_it_Ou(xi, xf, xc, yc)
use puls_eqt_set_opt05
use RK4_Set
implicit none
complex(8), dimension(1:) :: xc, yc
real(8), dimension(1:2*size(yc)) :: x, y
real(8) :: xi, xf
integer :: dim
	dim = size(yc)
	x(1:dim) = dreal(xc(1:dim))
	x(dim+1: 2*dim) = dimag(xc(1:dim))
	
	call RK4(pes05_vr_Ou, 4, xi, xf, x, y)
	yc(1:dim) = dcmplx(y(1:dim), y(dim+1: 2*dim))
endsubroutine pe09_it_Ou

!D1 B4---------------------------------------------------
subroutine pe09_bc_Coi(mode, r_, yc)
use puls_eqt_set_opt05
implicit none
real(8) :: r_
complex(8) :: yc(1:)
integer :: mode
	call pes05_bc_Coi(mode, r_, yc)
endsubroutine pe09_bc_Coi

!D1 B5---------------------------------------------------
subroutine pe09_bc_Cof(zc, cY)
use puls_eqt_set_opt05
implicit none
complex(8), dimension(1:) :: cY
complex(8) :: zc(1:)
	call pes05_bc_Cof(zc, cY)
endsubroutine pe09_bc_Cof

!D1 B6---------------------------------------------------
subroutine pe09_bc_Cri(zc, cY)
use puls_eqt_set_opt05
implicit none
complex(8), dimension(1:) :: cY
complex(8) :: zc(1:)
	call pes05_bc_Cri(zc, cY)
endsubroutine pe09_bc_Cri

!D1 B7---------------------------------------------------
subroutine pe09_bc_Crf(mode, zc)
use puls_eqt_set_opt05
implicit none
complex(8) :: zc(1:)
integer :: mode
	call pes05_bc_Crf(mode,zc)
endsubroutine pe09_bc_Crf

!D1 B8---------------------------------------------------
subroutine pe09_cY_solve(cY, coef)
!	M x coef = Vec
!use numerical_libraries
use lin_sol_gen_int
use puls_eqt_set_opt05_ex1
use Eigenvalues_Eigenvectors
implicit none
complex(8) :: cY(1: ,1:), coef(1:8), M(1:7, 1:7), X(1:7, 1), Vec(1:7, 1), A(4:8), det
real(8) :: gamma_, err_(1:6)

	gamma_ = pes05_gamma(2, P_Cr(2*pg_N2), rho_Cr(2*pg_N2), pg_N2*2, 0, pg_N2*2)
	call pes05_bc_R_Coe(P_Cr(2*pg_N2), rho_Cr(2*pg_N2), m_Cr(2*pg_N2), nu_Cr(2*pg_N2), mu_Cr(2*pg_N2), gamma_, r_Cr(2*pg_N2), A(4:8))
	coef(1) = (1.d0, 0.d0)
	M(1:6,1:2) = -cY(1:6, 2:3)
	M(1:6,3:7) = cY(1:6, 4:8)
	M(7,1:2) = (0.d0, 0.d0)
	M(7,3:7) = A(4:8)			! Summation A = 0; condition eqv to X(R) = 0
	Vec(1:6, 1) = cY(1:6,1)
	Vec(7, 1) = (0.d0, 0.d0)

	call LIN_SOL_GEN(M, Vec, X)
	coef(2:8) = X(1:7,1)

call Determinant_N_c(M(1:7,1:7), det, 7)
write(*,*) "det", det
write(*,*) "*******************************************"

!write(*,*) -coef(2)*cY(1:6,2) - coef(3)*cY(1:6,3) + coef(4)*cY(1:6,4) + coef(5)*cY(1:6,5) + coef(6)*cY(1:6,6) + coef(7)*cY(1:6,7) + coef(8)*cY(1:6,8)
!write(*,*) "&&&&&&&&&&&&&&&&&&&&&&&&&"
!write(*,*) cY(1:6,1)
!pause
endsubroutine pe09_cY_solve

subroutine pe09_cY_solve2(cY, coef)
!	M x coef = Vec
!use numerical_libraries
use lin_sol_gen_int
use puls_eqt_set_opt05_ex1
use Eigenvalues_Eigenvectors
implicit none
complex(8) :: cY(1: ,1:), coef(1:7), M(1:6, 1:6), X(1:6, 1), Vec(1:6, 1), A(3:7)
real(8) :: gamma_, err_(1:6)

	gamma_ = pes05_gamma(2, P_Cr(2*pg_N2), rho_Cr(2*pg_N2), pg_N2*2, 0, pg_N2*2)
	call pes05_bc_R_Coe(P_Cr(2*pg_N2), rho_Cr(2*pg_N2), m_Cr(2*pg_N2), nu_Cr(2*pg_N2), mu_Cr(2*pg_N2), gamma_, r_Cr(2*pg_N2), A(3:7))
	coef(1) = (1.d0, 0.d0)
	M(1:4,1) = -cY(1:4, 2)
	M(5,1) = -cY(6, 2)
	M(1:4,2:6) = cY(1:4, 3:7)
	M(5,2:6) = cY(6, 3:7)
	M(6,1) = (0.d0, 0.d0)
	M(6,2:6) = A(3:7)			! Summation A = 0; condition eqv to X(R) = 0
	Vec(1:4, 1) = cY(1:4,1)
	Vec(5,1) = cY(6,1)
	Vec(6, 1) = (0.d0, 0.d0)

!call LinSys_Cramer_c(M(1:6, 1:6), Vec(1:6,1), coef(2:7), 6, err_p = err_)
!write(*,*) err_
!pause
!write(*,*) -coef(2)*cY(1:5,2) + coef(3)*cY(1:5,3) + coef(4)*cY(1:5,4) + coef(5)*cY(1:5,5) + coef(6)*cY(1:5,6) + coef(7)*cY(1:5,7)
!write(*,*) "&&&&&&&&&&&&&&&&&&&&&&&&&"
!write(*,*) cY(1:5,1)
!pause

	call LIN_SOL_GEN(M, Vec, X)
	coef(2:7) = X(1:6,1)

!write(*,*) -coef(2)*cY(1:6,2) + coef(3)*cY(1:6,3) + coef(4)*cY(1:6,4) + coef(5)*cY(1:6,5) + coef(6)*cY(1:6,6) + coef(7)*cY(1:6,7)
!write(*,*) "&&&&&&&&&&&&&&&&&&&&&&&&&"
!write(*,*) cY(1:6,1)
!pause

endsubroutine pe09_cY_solve2

!D1 B9---------------------------------------------------
subroutine pe09_bc_Oui(H0, K, yc)
use puls_eqt_set_opt05
implicit none
complex(8) :: H0, K, yc(1:2)
	call pes05_bc_Oui(H0, K, yc)
endsubroutine pe09_bc_Oui

!D1 B10---------------------------------------------------
subroutine pe09_bc_Ouf(r_, yc, be, ga)
use puls_eqt_set_opt05
implicit none
complex(8) :: yc(1:2), be, ga
real(8) :: r_
	if (gr_opt == 1) then
		call pes05_bc_Ouf_LD(r_, yc, be, ga)
	elseif (gr_opt == 2) then
		call pes05_bc_Ouf_Scatter(r_, yc, be, ga)
	elseif (gr_opt == 3) then
		call pes05_bc_Ouf_phase(r_, yc, be, ga)
	else
		write(*,*) "err: pe gr_opt"
		pause
	endif
endsubroutine pe09_bc_Ouf

endmodule puls_eqt_solve_opt09