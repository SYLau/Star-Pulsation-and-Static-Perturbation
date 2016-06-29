module main

contains

!A1---------------------------------------------------
subroutine main_body
implicit none

	call settings

	call initialization

	call stat_profile
	call puls_grid

	call eigen_freq

	call o_puls_soln
	call verify_soln

endsubroutine main_body


!S1---------------------------------------------------
subroutine settings
use settings_opt01
implicit none
	call se01_Ctrl
endsubroutine settings

!B1---------------------------------------------------
subroutine initialization
use global_var
implicit none
	
	if (allocated(Coe)) deallocate(Coe)
	if (allocated(CoeC)) deallocate(CoeC)

	if (sp_opt == 4) then
		P_i(1) = P_t
		P_i(2) = P_g
		!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		sp_Ni = 500
		!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	endif
	
	if (pe_opt == 1) then
		allocate(Coe(1:3))
		NMat = 3
	elseif (pe_opt == 2) then
		allocate(Coe(1:3))
		NMat = 3
	elseif (pe_opt == 3) then
		allocate(CoeC(1:5))
		NMat = 4
	elseif (pe_opt == 4) then
		allocate(Coe(1:2))
		NMat = 2
	elseif (pe_opt == 5) then
		allocate(CoeC(1:7))
		NMat = 5
	elseif (pe_opt == 6) then
		allocate(CoeC(1:10))
		NMat = 5
	elseif (pe_opt == 7) then
		allocate(Coe(1:7))
		NMat = 5
	elseif (pe_opt == 8) then
		allocate(Coe(1:4))
		NMat = 4
	elseif (pe_opt == 9) then
		allocate(CoeC(1:8))
		NMat = 6
	else
		write(*,*) "err: pe_opt; initialization"
		pause
	endif
	
endsubroutine initialization 

!B2---------------------------------------------------
subroutine stat_profile
use stat_profile_opt01
use stat_profile_opt02
use stat_profile_opt03
use stat_profile_opt04
use stat_profile_opt05
use stat_profile_opt06
implicit none
	if (sp_opt == 1) then
		call sp01_Ctrl
	elseif (sp_opt == 2) then
		call sp02_Ctrl
	elseif (sp_opt == 3) then
		call sp03_Ctrl
	elseif (sp_opt == 4) then
		call sp04_Ctrl
	elseif (sp_opt == 5) then
		call sp05_Ctrl
	elseif (sp_opt == 6) then
		call sp06_Ctrl
	else
		write(*,*) "err: sp_opt"
		pause
	endif
endsubroutine stat_profile

!B3---------------------------------------------------
subroutine puls_grid
use puls_grid_opt01
use puls_grid_opt02
use puls_grid_opt03
use puls_grid_opt04
implicit none
	if (pg_opt == 1) then
		call pg01_Ctrl
	elseif (pg_opt == 2) then
		call pg02_Ctrl
	elseif (pg_opt == 3) then
		call pg03_Ctrl
	elseif (pg_opt == 4) then
		call pg04_Ctrl
	else
		write(*,*) "err: pg_opt"
		pause
	endif
endsubroutine puls_grid

!B4---------------------------------------------------
subroutine eigen_freq
use eigen_freq_opt01
use eigen_freq_opt02
implicit none
	if (ei_opt == 1) then
		call ei01_Ctrl
	elseif (ei_opt == 2) then
		call ei02_Ctrl
	endif
endsubroutine eigen_freq

!B5---------------------------------------------------
subroutine o_puls_soln
implicit none
	! Null
endsubroutine o_puls_soln

!B6---------------------------------------------------
subroutine verify_soln
use veri_soln_opt01
implicit none
	if (ei_opt == 1) then
		call vs01_Ctrl
	elseif (ei_opt == 2) then

	endif
endsubroutine verify_soln



! OTHER STUFF
! Test the convergence of RK4
subroutine test_rk4
use RK4_Set
implicit none
integer :: i, n_max
real(8) :: x_final, x, dx, yi, yf
	
do while (.true.)
	read(*,*) n_max
	write(*,*) n_max
	x_final = 100.d0
	dx = (x_final - 0.d0)/n_max
	x = 0.d0
	yi = 0.d0
	do i = 1, n_max
		call RK4(yy, 1, x, x+dx, yi, yf)
!		call oderk(x,x + dx,yi,1,yy1) 
		x = x + dx
		yi = yf
	enddo
	yf = yi
	write(*,*) yf
enddo
endsubroutine test_rk4

subroutine yy(n, t, x, fcn)
implicit none
integer :: n
real(8) :: t, x(1:n), fcn(1:n)
	fcn(1) = t**4
endsubroutine yy

subroutine yy1(t, x, fcn)
implicit none
real(8) :: t, x, fcn
	fcn = t**4
endsubroutine yy1

endmodule main