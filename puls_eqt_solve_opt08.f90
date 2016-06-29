module puls_eqt_solve_opt08
!	1 components model; Full Newtonian
use global_var

contains
!D1---------------------------------------------------
subroutine pe08_Ctrl(det, wFile)
implicit none
logical :: wFile
real(8) :: det
	
	if (NMat /= 4) then
		write(*,*) "err: pe08_Ctrl; NMat mismatch"
		pause
	endif
	if (pe_fcall == .true.) then
		write(100,*) "puls_eqt_solve_opt08: linear coupled ode solver"
		write(100,*) "- 1 component neutron star model"
		write(100,*) "- no Cowling Approximation used"
		write(100,*) "- independent variable x"
		write(100,*) "- matching at crust core interface; with continuity conditions"
		write(100,*) "- solve by shooting frequency; with 2 independent solutions in the crust"
		write(100,*) "- matching condition determined by the continuity of the COMBINATION of the solns"
		write(100,*) "- up to an arbitrary constant"
			
		if (pes_opt == 4) then
			write(*,*) "Full Newtonian pulsation equations"
			write(*,*) "=========================================================="
		else
			write(*,*) "err: pes_opt"
			pause
		endif
	endif

		call pe08_iterate(det, wFile)

	if (pe_fcall == .true.) then
		write(100,*) "pe08_iterate: passed"
		write(100,*) "=========================================================="
		pe_fcall = .false.
	endif

endsubroutine pe08_Ctrl

!D1 A1---------------------------------------------------
subroutine pe08_iterate(det, wFile)
use Format_IO
use FWrite
implicit none
real(8), dimension(1:Neqt) :: x, y, z
real(8) :: det
integer :: i, j, ii, it, mode
logical :: wFile
real(8), allocatable, dimension(:, :) :: z_M

allocate(z_M(1:NMat, 1:NMat))
	
	!¡¸
	do it = 1, 4
		if (it == 1 .or. it == 2) then
			if (it == 1) then 
				mode = 1
				if (wFile == .true.) open(01, file = pef_FN_Co1, status = 'replace')
			elseif (it == 2) then
				mode = 2
				if (wFile == .true.) open(01, file = pef_FN_Co2, status = 'replace')
			endif

		! core
			call pe08_bc_Coi(mode, y(1:4))
			if (wFile == .true.) call Write5R8(01, FR8, r_Co(0), Coe(it)*y(1), Coe(it)*y(2), Coe(it)*y(3), Coe(it)*y(4))
			do i = 0, pg_N1-1
				ii = i*2
				x(1:4) = y(1:4)
				call pe08_it_Co(x_Co(ii), x_Co(ii+2), x(1:4), y(1:4))
				if (wFile == .true.) call Write5R8(01, FR8, r_Co(ii+2), Coe(it)*y(1), Coe(it)*y(2), Coe(it)*y(3), Coe(it)*y(4))
			enddo
			if (wFile == .true.) close(44)
			call pe08_bc_Cof(y(1:4), z_M(1:5, it))

		elseif (it == 3 .or. it == 4) then

			if (it == 3) then 
				mode = 1
				if (wFile == .true.) open(01, file = pef_FN_Cr1, status = 'replace')
			elseif (it == 4) then
				mode = 2
				if (wFile == .true.) open(01, file = pef_FN_Cr2, status = 'replace')
			endif
		
		! crust backward integration
			call pe08_bc_Crf(mode, y(1:4))
			if (wFile == .true.) call Write5R8(01, FR8, r_Cr(pg_N2*2), Coe(it)*y(1), Coe(it)*y(2), Coe(it)*y(3), Coe(it)*y(4))
			do i = pg_N2, 1, -1
				ii = i*2
				x(1:4) = y(1:4)
				call pe08_it_Cr(x_Cr(ii), x_Cr(ii-2), x(1:4), y(1:4))
				if (wFile == .true.) call Write5R8(01, FR8, r_Cr(ii-2), Coe(it)*y(1), Coe(it)*y(2), Coe(it)*y(3), Coe(it)*y(4))
			enddo
			if (wFile == .true.) close(45)
			call pe08_bc_Cri(y(1:4), z_M(1:4, it))
		endif
		if (wFile == .true.) close(01)
	enddo

	if (wFile /= .true.) then
		call pe08_zM_solve(1, z_M(1:NMat ,1:NMat), det)
	else 
		call pe08_zM_solve(2, z_M(1:NMat ,1:NMat), det)
	endif

endsubroutine pe08_iterate

!D1 B1---------------------------------------------------
subroutine pe08_it_Co(xi, xf, x, y)
use puls_eqt_set_opt04
use RK4_Set
implicit none
real(8), dimension(:) :: x, y
real(8) :: xi, xf
	
	if (pes_opt == 4) then
		call RK4(pes04_vx_Co, 4, xi, xf, x, y)
	else
		write(*,*) "err: pe08_it_Co"
		pause
	endif

endsubroutine pe08_it_Co

!D1 B2---------------------------------------------------
subroutine pe08_it_Cr(xi, xf, x, y)
use puls_eqt_set_opt04
use RK4_Set
implicit none
real(8), dimension(:) :: x, y
real(8) :: xi, xf

	if (pes_opt == 4) then
		call RK4(pes04_vx_Cr_FM, 4, xi, xf, x, y)
	else
		write(*,*) "err: pe08_it_Cr"
		pause
	endif

endsubroutine pe08_it_Cr

!D1 B3---------------------------------------------------
subroutine pe08_bc_Coi(mode, y)
use puls_eqt_set_opt04
implicit none
real(8) :: y(1:4)
integer :: mode
	if (pes_opt == 4) then
		call pes04_bc_Coi(mode, y)
	else
		write(*,*) "err: pe08_bc_Coi"
		pause
	endif

endsubroutine pe08_bc_Coi

!D1 B4---------------------------------------------------
subroutine pe08_bc_Cof(y, z_M)
use puls_eqt_set_opt04
implicit none
real(8), dimension(1:) :: z_M
real(8) :: y(1:4)
	if (pes_opt == 4) then
		call pes04_bc_Co_FM_f(y, z_M)
	else
		write(*,*) "err: pe08_bc_Cof"
		pause
	endif
endsubroutine pe08_bc_Cof

!D1 B5---------------------------------------------------
subroutine pe08_bc_Cri(y, z_M)
use puls_eqt_set_opt04
implicit none
real(8), dimension(1:) :: z_M
real(8) :: y(1:4)
integer :: mode
	if (pes_opt == 4) then
		call pes04_bc_Cr_FM_i(y, z_M)
	else
		write(*,*) "err: pe08_bc_Cri"
		pause
	endif
endsubroutine pe08_bc_Cri

!D1 B6---------------------------------------------------
subroutine pe08_bc_Crf(mode,y)
use puls_eqt_set_opt04
implicit none
real(8) :: y(1:4)
integer :: mode
	if (pes_opt == 4) then
		call pes04_bc_Cr_FM_f(mode, y)
	else
		write(*,*) "err: pe08_bc_Crf"
		pause
	endif
endsubroutine pe08_bc_Crf

!D1 B8---------------------------------------------------
subroutine pe08_zM_solve(mode, z_M, det)
use Eigenvalues_Eigenvectors
implicit none
real(8) :: z_M(1:4 ,1:4), M(1:4, 1:4), det
real(8) :: r_, rho_, g_, z1_f, z5_f, z6_f
integer :: mode
	
	r_ = R0
	rho_ = rho_Cr(2*pg_N2)
	g_ = Grav_Const * M0/R0**2

	M(1:4, 1:4) = z_M(1:4 ,1:4)

	if (mode == 1) then
		call Determinant_N(M(1:4 ,1:4), det, 4)
	elseif (mode == 2) then
		call Eigenvector(M(1:4 ,1:4), Coe, 4)
		Coe(3:4) = -Coe(3:4)
		write(*,*) "Coefficients of linear combination:"
		write(*,*) Coe
		write(*,*) "=========================================================="
		write(*,*) "matching matrix:"
		write(*,*) "Co: ", z_M(1:4,1)*Coe(1) + z_M(1:4,2)*Coe(2)
		write(*,*) "Cr: ", z_M(1:4,3)*Coe(3) +z_M(1:4,4)*Coe(4)
		write(*,*) "=========================================================="
	else
		write(*,*) "err: pe08_zM_solve"
	endif

endsubroutine pe08_zM_solve

endmodule puls_eqt_solve_opt08