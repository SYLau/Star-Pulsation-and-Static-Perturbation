module puls_eqt_solve_opt06
!	3 components model
use global_var

contains
!D1---------------------------------------------------
subroutine pe06_Ctrl(beta, gamma, wFile)
implicit none
complex(8) :: beta, gamma
logical :: wFile

	if (NMat /= 5) then
		write(*,*) "err: pe06_Ctrl; NMat mismatch"
		pause
	endif
	if (pe_fcall == .true.) then
		write(100,*) "puls_eqt_solve_opt05: linear coupled ode solver"
		write(100,*) "- 3 components neutron star model"
		write(100,*) "- Fully relativistic formalism; LD2 formalism"
		if (pes_opt == 3) then
			write(*,*) "LD2/ Andersson pulsation equations"
			write(*,*) "=========================================================="
		else
			write(*,*) "err: pes_opt"
			pause
		endif
	endif

		call pe06_iterate(beta, gamma, wFile)

	if (pe_fcall == .true.) then
		write(100,*) "pe06_iterate: passed"
		write(100,*) "=========================================================="
		pe_fcall = .false.
	endif

endsubroutine pe06_Ctrl

!D1 A1---------------------------------------------------
subroutine pe06_iterate(beta, gamma, wFile)
use Format_IO
use FWrite
use puls_eqt_set_opt03
use puls_eqt_set_opt03_ex1
implicit none
complex(8) :: E(1:2, 1:4), X_R
real(8) :: gamma_1
complex(8), dimension(1:4) :: y, y_R
complex(8), dimension(1:6) :: x, z
complex(8), dimension(1:10) :: coef
complex(8) :: W, V, H0, K, beta, gamma
complex(8) :: W_mid(1:2)
real(8) :: r_Ou
integer :: i, j, ii, n, it, mode
logical :: wFile
complex(8), allocatable, dimension(:, :) :: cY
allocate(cY(1:NMat, 1:10))

	!¡¸Within the star, solutions represented by vectors: core cY1, cY2; ocean cY3, cY4;crust cY5-cY9
	do it = 1, 10
		if (it == 1 .or. it == 2) then
			if (it == 1) then 
				mode = 1
				if (wFile == .true.) open(01, file = pef_LD_Co1, status = 'replace')
			elseif (it == 2) then
				mode = 2
				if (wFile == .true.) open(01, file = pef_LD_Co2, status = 'replace')
			endif
	
			call pe06_bc_Coi(mode, r_Co(0), y)

			if (wFile == .true.) call Write1R82AR8(01, FR8,r_Co(0), dreal(CoeC(it)*y(1:4)), dimag(CoeC(it)*y(1:4)), m=4)
			do i = 0, pg_N1-1
				ii = i*2
				x(1:4) = y(1:4)
				call pe06_it_Co(x_Co(ii), x_Co(ii+2), x(1:4), y(1:4))
				if (wFile == .true.) call Write1R82AR8(01, FR8,r_Co(ii+2), dreal(CoeC(it)*y(1:4)), dimag(CoeC(it)*y(1:4)), m=4)
			enddo

			call pe06_bc_Cof(y, cY(1:5 ,it))

		elseif (it == 3 .or. it == 4 .or. it == 5) then
			if (it == 3) then 
				mode = 1
				if (wFile == .true.) open(01, file = pef_LD_Oc1, status = 'replace')
			elseif (it == 4) then
				mode = 2
				if (wFile == .true.) open(01, file = pef_LD_Oc2, status = 'replace')
			elseif (it == 5) then
				mode = 3
				if (wFile == .true.) open(01, file = pef_LD_Oc3, status = 'replace')
			endif
			call pe06_bc_Ocf(mode, y)

			if (wFile == .true.) call Write1R82AR8(01, FR8,r_Oc(2*pg_N3), dreal(CoeC(it)*y(1:4)), dimag(CoeC(it)*y(1:4)), m=4)
			do i = pg_N3, 1, -1
				ii = i*2
				x(1:4) = y(1:4)
				call pe06_it_Oc(x_Oc(ii), x_Oc(ii-2), x(1:4), y(1:4))
				if (wFile == .true.) call Write1R82AR8(01, FR8,r_Oc(ii-2), dreal(CoeC(it)*y(1:4)), dimag(CoeC(it)*y(1:4)), m=4)
			enddo
			call pe06_bc_Oci(y, cY(1:5 ,it))

		elseif (it == 6 .or. it == 7 .or. it == 8 .or. it == 9 .or. it == 10) then
			if (it == 6) then 
				mode = 1
				if (wFile == .true.) open(01, file = pef_LD_Cr1, status = 'replace')
			elseif (it == 7) then
				mode = 2
				if (wFile == .true.) open(01, file = pef_LD_Cr2, status = 'replace')
			elseif (it == 8) then
				mode = 3
				if (wFile == .true.) open(01, file = pef_LD_Cr3, status = 'replace')
			elseif (it == 9) then
				mode = 4
				if (wFile == .true.) open(01, file = pef_LD_Cr4, status = 'replace')
			elseif (it == 10) then
				mode = 5
				if (wFile == .true.) open(01, file = pef_LD_Cr5, status = 'replace')
			endif
			call pe06_bc_Crf(mode, z)
			if (wFile == .true.) call Write1R82AR8(01, FR8,r_Cr(2*pg_N2), dreal(CoeC(it)*z(1:6)), dimag(CoeC(it)*z(1:6)), m=6)
			do i = pg_N2, 1, -1
				ii = i*2
				x(1:6) = z(1:6)
				call pe06_it_Cr(x_Cr(ii), x_Cr(ii-2), x(1:6), z(1:6))
				if (wFile == .true.) call Write1R82AR8(01, FR8,r_Cr(ii-2), dreal(CoeC(it)*z(1:6)), dimag(CoeC(it)*z(1:6)), m=6)
			enddo
			call pe06_bc_Cri(z, cY(1:5 ,it))
		endif
		if (wFile == .true.) close(01)
	enddo

	!¡¸Solve algebraric equations from BCs to obtain H0 and K
	call pe06_cY_solve(cY(1:5 ,1:10), coef(1:10))
	CoeC(1:10) = coef(1:10)
	
	Y_R(1) = coef(3)
	Y_R(2) = coef(4)
	Y_R(3) = coef(5)
	Y_R(4) = (0.d0, 0.d0)
	call pes03_E_Fluid(P_Ou, rho_Ou, m_Ou, nu_Ou, ROu_i, E)
	H0 = (0.d0, 0.d0)
	do i = 1,4
		H0 = H0 + E(1,i) * Y_R(i) 
	enddo
	K = coef(4)
	W = coef(5)

	!¡¸Normalize H0, K by W at surface
	H0 = H0/W
	K = K/W

	if (zero_freq == .true.) return ! Skip the outter integration for zero frequency

!gamma_1 = pes03_gamma(2, P_Cr, rho_Cr, 2*pg_N2, 0, pg_N2*2)
!call pes03_E_Crust(P_Cr(2*pg_N2), rho_Cr(2*pg_N2), m_Cr(2*pg_N2), nu_Cr(2*pg_N2), mu_Cr(2*pg_N2), gamma_1, r_Cr(2*pg_N2), E)
!X_R = E(2,1)*coef(3)/R0**2 + E(2,2)*coef(4)/R0**2 + E(2,3)*coef(5)/R0**2 + E(2,4)*coef(6)  + E(2,5)*coef(7)
!write(*,*) P_Ou, rho_Ou, m_Ou, nu_Ou, ROu_i
!write(*,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
!write(*,*) E(1,1:4)
!write(*,*) coef(1)*cY(1,1) + coef(2)*cY(1,2), coef(6)*cY(1,6) +coef(7)*cY(1,7) +coef(8)*cY(1,8) +coef(9)*cY(1,9) +coef(10)*cY(1,10)
!pause

	!	Outside the Neutron Star
	if (wFile == .true.) open(01, file = pef_LD_Ou, status = 'replace')

	call pe06_bc_Oui(H0, K, y(1:2))
	ROu_f = R_inf_f/cdabs(cdsqrt(OmeC_sq))
	ROu_i = R0
	r_Ou = ROu_i
	dr_Ou = (ROu_f - ROu_i)/pg_NOu
	
	if (wFile == .true.) call Write1R82AR8(01, FR8,ROu_i*dreal(cdsqrt(OmeC_sq))/pi, dreal(y(1:2)), dimag(y(1:2)), m=2)	
	!r normalized by pi/omega, to visualize phase angle

	do i = 0, pg_NOu -1
		r_Ou = i * dr_Ou + ROu_i
		x(1:2) = y(1:2)
		call pe06_it_Ou(r_Ou, r_Ou + dr_Ou, x(1:2), y(1:2))
		if (wFile == .true.) call Write1R82AR8(01, FR8,(r_Ou + dr_Ou)*dreal(cdsqrt(OmeC_sq))/pi, dreal(y(1:2)), dimag(y(1:2)), m=2)
	enddo

	r_Ou = ROu_f

	call pe06_bc_Ouf(ROu_f, y(1:2), beta, gamma)

	if (wFile == .true.) close(01)

endsubroutine pe06_iterate

!D1 B1---------------------------------------------------
subroutine pe06_it_Co(xi, xf, xc, yc)
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
endsubroutine pe06_it_Co

!D1 B2---------------------------------------------------
subroutine pe06_it_Cr(xi, xf, xc, yc)
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

	call RK4(pes03_vx_Cr, 12, xi, xf, x, y)
	yc(1:dim) = dcmplx(y(1:dim), y(dim+1: 2*dim))
endsubroutine pe06_it_Cr

!D1 B3---------------------------------------------------
subroutine pe06_it_Oc(xi, xf, xc, yc)
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
	
	call RK4(pes03_vx_Oc, 8, xi, xf, x, y)
	yc(1:dim) = dcmplx(y(1:dim), y(dim+1: 2*dim))
endsubroutine pe06_it_Oc

!D1 B3---------------------------------------------------
subroutine pe06_it_Ou(xi, xf, xc, yc)
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
endsubroutine pe06_it_Ou

!D1 B4---------------------------------------------------
subroutine pe06_bc_Coi(mode, r_, yc)
use puls_eqt_set_opt03
implicit none
real(8) :: r_
complex(8) :: yc(1:)
integer :: mode
	call pes03_bc_Coi(mode, r_, yc)
endsubroutine pe06_bc_Coi

!D1 B5---------------------------------------------------
subroutine pe06_bc_Cof(yc, cY)
use puls_eqt_set_opt03
implicit none
complex(8), dimension(1:) :: cY
complex(8) :: yc(1:)
	call pes03_bc_Cof(yc, cY)
endsubroutine pe06_bc_Cof

!D1 B6---------------------------------------------------
subroutine pe06_bc_Cri(zc, cY)
use puls_eqt_set_opt03
implicit none
complex(8), dimension(1:) :: cY
complex(8) :: zc(1:)
	call pes03_bc_Cri(zc, cY)
endsubroutine pe06_bc_Cri

!D1 B7---------------------------------------------------
subroutine pe06_bc_Crf(mode, zc)
use puls_eqt_set_opt03
implicit none
complex(8) :: zc(1:)
integer :: mode
	call pes03_bc_Crf(mode,zc)
endsubroutine pe06_bc_Crf

!D1 B9---------------------------------------------------
subroutine pe06_bc_Oci(yc, cY)
use puls_eqt_set_opt03
implicit none
complex(8), dimension(1:) :: cY
complex(8) :: yc(1:)
	call pes03_bc_Oci(yc, cY)
endsubroutine pe06_bc_Oci

!D1 B9---------------------------------------------------
subroutine pe06_bc_Ocf(mode, yc)
use puls_eqt_set_opt03
implicit none
complex(8) :: yc(1:)
integer :: mode
	call pes03_bc_Ocf(mode,yc)
endsubroutine pe06_bc_Ocf

!D1 B8---------------------------------------------------
subroutine pe06_cY_solve(cY, coef)
!	M x coef = Vec
!use numerical_libraries
!	Core: coef(1) cY(1) + coef(2) cY(2) ||| Crust: Sum[coef(6-10) cY(6-10)] ||| Ocean: Sum[coef(3-5) cY(3-5)]
!	M1= Transfer matrix C/O Int; M2= Transfer matrix C/C Int; M3= M1 x M2
use lin_sol_gen_int
use puls_eqt_set_opt03_ex1
use Eigenvalues_Eigenvectors
implicit none
complex(8) :: cY(1: ,1:), coef(1:10), M(1:5, 1:5), M1(1:4, 1:3), M2(1:5, 1:4), M3(1:5, 1:3), X(1:5, 1), Vec(1:5, 1)
real(8) :: err_(1:6)
	
	coef(1) = 1.d0
	M1(1:4, 1:3) = cY(1:4, 3:5)
	M1(1:3, 1:3) = M1(1:3, 1:3) * R0**2	! Metric Variables have a scaling factor 1/R0**2

	M2(1:4, 1:4) = cY(1:4, 6:9)
	M2(5, 1:4) = cY(5, 6:9)

	M3 = matmul(M2, M1)

	M(1:5,1) = -cY(1:5, 2)
	M(1:5,2:4) = M3(1:5, 1:3)
	M(1:5,5) = cY(1:5, 10)
	Vec(1:5, 1) = cY(1:5,1)

!call LinSys_Cramer_c(M(1:6, 1:6), Vec(1:6,1), coef(2:7), 6, err_p = err_)
!write(*,*) err_
!pause
!write(*,*) -coef(2)*cY(1:5,2) + coef(3)*cY(1:5,3) + coef(4)*cY(1:5,4) + coef(5)*cY(1:5,5) + coef(6)*cY(1:5,6) + coef(7)*cY(1:5,7)
!write(*,*) "&&&&&&&&&&&&&&&&&&&&&&&&&"
!write(*,*) cY(1:5,1)
!pause
	call LIN_SOL_GEN(M, Vec, X)
	coef(2) = X(1,1)
	coef(3:5) = X(2:4,1)
	coef(10) = X(5,1)

	coef(6:9) = matmul(M1, coef(3:5))

!write(*,*) X(1,1)*(-cY(1:5, 2)) + X(2,1)*M3(1:5, 1) + X(3,1)*M3(1:5, 2) + X(4,1)*M3(1:5, 3) + X(5,1)*cY(1:5, 10)
!write(*,*) "&&&&&&&&&&&&&&&&&&&&&&&&&"
!write(*,*) cY(1:5,1)
!write(*,*) "&&&&&&&&&&&&&&&&&&&&&&&&&"
!write(*,*) coef(3)*M1(1:4, 1) + coef(4)*M1(1:4, 2) + coef(5)*M1(1:4, 3)
!write(*,*) coef(6), coef(7), coef(8), coef(9)
!write(*,*) "&&&&&&&&&&&&&&&&&&&&&&&&&"
!write(*,*) coef(6)*M2(1:5, 1) + coef(7)*M2(1:5, 2) +  coef(8)*M2(1:5, 3) +  coef(9)*M2(1:5, 4) + coef(10)*cY(1:5, 10)
!write(*,*) coef(1)*cY(1:5,1) + coef(2)*cY(1:5,2)
!write(*,*) "&&&&&&&&&&&&&&&&&&&&&&&&&"
!write(*,*) matmul(M3,coef(3:5))
!write(*,*) matmul(M2,coef(6:9))
!pause
endsubroutine pe06_cY_solve

!D1 B9---------------------------------------------------
subroutine pe06_bc_Oui(H0, K, yc)
use puls_eqt_set_opt03
implicit none
complex(8) :: H0, K, yc(1:2)
	call pes03_bc_Oui(H0, K, yc)
endsubroutine pe06_bc_Oui

!D1 B10---------------------------------------------------
subroutine pe06_bc_Ouf(r_, yc, be, ga)
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
endsubroutine pe06_bc_Ouf

endmodule puls_eqt_solve_opt06