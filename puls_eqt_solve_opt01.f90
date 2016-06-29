module puls_eqt_solve_opt01
!	2 components model
use global_var

contains
!D1---------------------------------------------------
subroutine pe01_Ctrl(det, wFile)
implicit none
logical :: wFile
real(8) :: det
	
	if (NMat /= 3) then
		write(*,*) "err: pe01_Ctrl; NMat mismatch"
		pause
	endif
	if (pe_fcall == .true.) then
		write(100,*) "puls_eqt_solve_opt01: linear coupled ode solver"
		write(100,*) "- 2 components neutron star model; fluid core, solid crust"
		write(100,*) "- independent variable x"
		write(100,*) "- matching at crust core interface; with 3 continuity conditions"
		write(100,*) "- solve by shooting frequency; with 2 independent solutions in the crust"
		write(100,*) "- matching condition determined by the continuity of the COMBINATION of the solns"
		write(100,*) "- up to an arbitrary constant"
			
		if (pes_opt == 1) then
			write(*,*) "Newtonian Cowling pulsation equations"
			write(*,*) "=========================================================="
		elseif (pes_opt == 2) then
			write(*,*) "Relativistic Cowling pulsation equations"
			write(*,*) "=========================================================="
		else
			write(*,*) "err: pes_opt"
			pause
		endif
	endif

		call pe01_iterate(det, wFile)

	if (pe_fcall == .true.) then
		write(100,*) "pe01_iterate: passed"
		write(100,*) "=========================================================="
		pe_fcall = .false.
	endif

endsubroutine pe01_Ctrl

!D1 A1---------------------------------------------------
subroutine pe01_iterate(det, wFile)
use Format_IO
use FWrite
implicit none
real(8), dimension(1:Neqt) :: x, y, z
real(8) :: det
integer :: i, j, ii
logical :: wFile
real(8), allocatable, dimension(:, :) :: z_M

allocate(z_M(1:NMat, 1:NMat))

! core
	call pe01_bc_Coi(y(1:2))
	if (wFile == .true.) then
		open(44, file= pef_yA, status = 'replace')
		call Write3R8(44, FR8, r_Co(0), Coe(1)*y(1), Coe(1)*y(2))
	endif
	do i = 0, pg_N1-1
		ii = i*2
		x(1:2) = y(1:2)
		call pe01_it_Co(x_Co(ii), x_Co(ii+2), x(1:2), y(1:2))
		if (wFile == .true.) call Write3R8(44, FR8, r_Co(ii+2), Coe(1)*y(1), Coe(1)*y(2))
	enddo
	if (wFile == .true.) close(44)
	call pe01_bc_Cof(y(1:2), z_M)

! crust backward integration 01
	call pe01_bc_Crf(1, z(1:4))
	if (wFile == .true.) then
		open(45, file= pef_zC, status = 'replace')
		call Write5R8(45, FR8, r_Cr(pg_N2*2), Coe(2)*z(1), Coe(2)*z(2), Coe(2)*z(3), Coe(2)*z(4))
	endif
	do i = pg_N2, 1, -1
		ii = i*2
		x(1:4) = z(1:4)
		call pe01_it_Cr(x_Cr(ii), x_Cr(ii-2), x(1:4), z(1:4))
		if (wFile == .true.) call Write5R8(45, FR8, r_Cr(ii-2), Coe(2)*z(1), Coe(2)*z(2), Coe(2)*z(3), Coe(2)*z(4))
	enddo
	if (wFile == .true.) close(45)
	call pe01_bc_Cri(1,z(1:4), z_M)

! crust backward integration 02
	call pe01_bc_Crf(2, z(1:4))
	if (wFile == .true.) then
		open(46, file= pef_zD, status = 'replace')
		call Write5R8(46, FR8, r_Cr(pg_N2*2), Coe(3)*z(1), Coe(3)*z(2), Coe(3)*z(3), Coe(3)*z(4))
	endif
	do i = pg_N2, 1, -1
		ii = i*2
		x(1:4) = z(1:4)
		call pe01_it_Cr(x_Cr(ii), x_Cr(ii-2), x(1:4), z(1:4))
		if (wFile == .true.) call Write5R8(46, FR8, r_Cr(ii-2), Coe(3)*z(1), Coe(3)*z(2), Coe(3)*z(3), Coe(3)*z(4))
	enddo
	if (wFile == .true.) close(46)
	call pe01_bc_Cri(2,z(1:4), z_M)

	if (wFile /= .true.) then
		call pe01_zM_solve(1, z_M(1:NMat ,1:NMat), det)
	else 
		call pe01_zM_solve(2, z_M(1:NMat ,1:NMat), det)
	endif

endsubroutine pe01_iterate

!D1 B1---------------------------------------------------
subroutine pe01_it_Co(xi, xf, x, y)
use puls_eqt_set_opt01
use puls_eqt_set_opt02
use RK4_Set
implicit none
real(8), dimension(:) :: x, y
real(8) :: xi, xf
	
	if (pes_opt == 1) then
		call RK4(pes01_vx_Co, 2, xi, xf, x, y)
	elseif (pes_opt == 2) then
		call RK4(pes02_vx_Co, 2, xi, xf, x, y)
	endif

endsubroutine pe01_it_Co

!D1 B2---------------------------------------------------
subroutine pe01_it_Cr(xi, xf, x, y)
use puls_eqt_set_opt01
use puls_eqt_set_opt02
use RK4_Set
implicit none
real(8), dimension(:) :: x, y
real(8) :: xi, xf
	if (pes_opt == 1) then
		call RK4(pes01_vx_Cr, 4, xi, xf, x, y)
	elseif (pes_opt == 2) then
		call RK4(pes02_vx_Cr, 4, xi, xf, x, y)
	endif
endsubroutine pe01_it_Cr

!D1 B3---------------------------------------------------
subroutine pe01_bc_Coi(y)
use puls_eqt_set_opt01
use puls_eqt_set_opt02
implicit none
real(8) :: y(1:2)
	if (pes_opt == 1) then
		call pes01_bc_Coi(y)
	elseif (pes_opt == 2) then
		call pes02_bc_Coi(y)
	endif

endsubroutine pe01_bc_Coi

!D1 B4---------------------------------------------------
subroutine pe01_bc_Cof(y, z_M)
use puls_eqt_set_opt01
use puls_eqt_set_opt02
implicit none
real(8), dimension(1:, 1:) :: z_M
real(8) :: y(1:2)
	if (pes_opt == 1) then
		call pes01_bc_Cof(y, z_M)
	elseif (pes_opt == 2) then
		call pes02_bc_Cof(y, z_M)
	endif
endsubroutine pe01_bc_Cof

!D1 B5---------------------------------------------------
subroutine pe01_bc_Cri(mode,z, z_M)
use puls_eqt_set_opt01
use puls_eqt_set_opt02
implicit none
real(8), dimension(1:, 1:) :: z_M
real(8) :: z(1:4)
integer :: mode
	if (pes_opt == 1) then
		call pes01_bc_Cri(mode,z, z_M)
	elseif (pes_opt == 2) then
		call pes02_bc_Cri(mode,z, z_M)
	endif
endsubroutine pe01_bc_Cri

!D1 B6---------------------------------------------------
subroutine pe01_bc_Crf(mode,z)
use puls_eqt_set_opt01
use puls_eqt_set_opt02
implicit none
real(8) :: z(1:4)
integer :: mode
	if (pes_opt == 1) then
		call pes01_bc_Crf(mode,z)
	elseif (pes_opt == 2) then
		call pes02_bc_Crf(mode,z)
	endif
endsubroutine pe01_bc_Crf

!D1 B8---------------------------------------------------
subroutine pe01_zM_solve(mode, z_M, det)
use Eigenvalues_Eigenvectors
implicit none
real(8) :: z_M(1:3 ,1:3), det
integer :: mode
	if (mode == 1) then
		call Determinant_N(z_M(1:NMat ,1:NMat), det, NMat)
	elseif (mode == 2) then
		call Eigenvector(z_M(1:NMat,1:NMat), Coe, NMat)
		Coe(2:3) = -Coe(2:3)
		write(*,*) "Coefficients of linear combination:"
		write(*,*) Coe
		write(*,*) "=========================================================="
		write(*,*) "matching matrix:"
		write(*,*) "Co: ", z_M(1:3,1)*Coe(1)
		write(*,*) "Cr01: ", z_M(1:3,2)*Coe(2)
		write(*,*) "Cr02: ", z_M(1:3,3)*Coe(3)
		write(*,*) "=========================================================="
	else
		write(*,*) "err: pe01_zM_solve"
	endif

endsubroutine pe01_zM_solve

endmodule puls_eqt_solve_opt01