module eigen_freq_opt01
use global_var
use moI_opt01
contains

!C3---------------------------------------------------
subroutine ei01_Ctrl
!	- Lengthy parts serve to provide interactive IO
!	- Important parts are highlighted with " !**** "
use Root_Finding
use eigen_freq_opt01_ex1
implicit none
real(8) :: input, bisacc = 1.D-8
integer :: N, ans
real(8) :: I_, I_bar, Eta
real(8), allocatable, dimension(:, :) :: z_M
allocate(z_M(1:NMat, 1:NMat))

	write(100,*) "eigen_freq_opt01: solve eigenvalue problem by shooting"
	pe_fcall = .true.

34		write(*,*) "Options:"
		write(*,*) "1. Scan frequency"
		write(*,*) "2. Find eigenfunctions"
		write(*,*) "0. Back"
		write(*,*) "input: "
		read(*,*) ans
		write(*,*) "=========================================================="
		
		if (ans /= 1 .and. ans /= 2 .and. ans /= 0 ) then
			write(*,*) "invalid input"
			goto 34
		endif

		if (ans == 1) then
35			write(*,*) "Choose your input"
			write(*,*) "1. angular frequency"
			write(*,*) "2. frequency"
			write(*,*) "3. Omega^2"
			write(*,*) "0. back"
			read(*,*) ans
			write(*,*) "=========================================================="
			
			if (ans == 0) then
				goto 34
			elseif (ans /= 1 .and. ans /= 2 .and. ans /= 3 ) then
				write(*,*) "invalid input"
				goto 35
			endif

			write(*,*) "input lower bound"
			read(*,*) input
			afreqL = freq_conv(ans, input)

			write(*,*) "input upper bound"
			read(*,*) input
			afreqH = freq_conv(ans, input)
			
			write(*,*) "input no of steps"
			read(*,*) N

			write(*,*) "=========================================================="
		!********** Shooting **********
		call ei01_shoot(N)
		write(100,*) "ei01_shoot: passed"
		write(*,*) "Shooting complete!"
		print *, char(7)	! "BEEP"
		write(*,*) "=========================================================="
			goto 34

		elseif (ans == 2) then
36			write(*,*) "Choose your method"
			write(*,*) "1. input frequency"
			write(*,*) "2. bisection"
			write(*,*) "0. back"
			read(*,*) ans
			write(*,*) "=========================================================="
			
			if (ans == 0) then
				goto 34
			elseif (ans /= 1 .and. ans /= 2 .and. ans /= 3 ) then
				write(*,*) "invalid input"
				goto 36
			endif
			
			if (ans == 1) then
			
37				write(*,*) "Choose your input"
				write(*,*) "1. angular frequency"
				write(*,*) "2. frequency"
				write(*,*) "3. Omega^2"
				write(*,*) "0. back"
				read(*,*) ans
				write(*,*) "=========================================================="
		
				if (ans == 0) then
					goto 36
				elseif (ans /= 1 .and. ans /= 2 .and. ans /= 3 ) then
					write(*,*) "invalid input"
					goto 37
				endif
				
				write(*,*) "input the frequency"
				read(*,*) input
				afreq = freq_conv(ans, input)
				write(*,*) "=========================================================="

			!********** Eigen function by direct input of frequency**********
			call ei01_W_soln
			write(100,*) "ei01_soln: passed"

			elseif (ans == 2) then
38				write(*,*) "Choose your input"
				write(*,*) "1. angular frequency"
				write(*,*) "2. frequency"
				write(*,*) "3. Omega^2"
				write(*,*) "0. back"
				read(*,*) ans
				write(*,*) "=========================================================="
			
				if (ans == 0) then
					goto 34
				elseif (ans /= 1 .and. ans /= 2 .and. ans /= 3 ) then
					write(*,*) "invalid input"
					goto 38
				endif

				write(*,*) "input lower bound"
				read(*,*) input
				afreqL = freq_conv(ans, input)

				write(*,*) "input upper bound"
				read(*,*) input
				afreqH = freq_conv(ans, input)
				write(*,*) "=========================================================="

			!********** Eigen function by bisection**********
			afreq = rtbis(ei01_bis,afreqL,afreqH,bisacc)
			call ei01_w_soln

			write(100,*) "ei01_soln: passed"
			
			!********** Calculate Moment of Inertia **********
			OmeC_sq = dcmplx(afreq**2/c**2, 0.d0)
			call ei01_unit(1)
			call mi01_Ctrl(I_, I_bar, Eta)
			call ei01_unit(2)

			endif
		
		! Output solution; eigen_freq_opt01_ex1
		if (pe_opt == 1) then
			call ei01_o_soln_C2L1R2
		elseif (pe_opt == 2) then
			call ei01_o_soln_C3L1R2
		elseif (pe_opt == 4) then
			call ei01_o_soln_C1L1R1
		elseif (pe_opt == 7) then
			call ei01_o_soln_C2L2R5
		elseif (pe_opt == 8) then
			call ei01_o_soln_C1L2R2
		endif
		write(100,*) "ei01_o_soln: passed"

		elseif (ans == 0) then
			write(100,*) "=========================================================="
			return
		endif
	
	write(*,*) "Solution Written!"
	print *, char(7)	! "BEEP"
	write(*,*) "=========================================================="
	goto 34
	
endsubroutine ei01_Ctrl

!C3 A1---------------------------------------------------
function freq_conv(i, input)
implicit none
real(8) :: freq_conv, input
integer :: i
	if (i == 1) freq_conv = input
	if (i == 2) freq_conv = input*2.d0*pi
	if (i == 3) freq_conv = dsqrt(Grav_Const* M0/R0**3)*input
endfunction freq_conv

!C3 A2---------------------------------------------------
subroutine ei01_shoot(N)
use Eigenvalues_Eigenvectors
use global_fcn
use Format_IO
use FWrite
implicit none
real(8) :: d_afreq, det, freq, d
integer :: i, N
	
	open(61, file = 'data/shooting.txt', status='replace')
	d_afreq = (afreqH-afreqL)/N
	do i = 1, N
		afreq = i*d_afreq + afreqL
		Omega_sq = afreq**2 * R0**3/ Grav_Const/ M0

		call puls_eqt_solve(det, .false.)

		freq = afreq/2.d0/pi
		call write5R8(61, FR8, afreq, freq, Omega_sq, det, d)
	enddo
	close(61)

endsubroutine ei01_shoot

!C3 A3---------------------------------------------------
subroutine ei01_w_soln
implicit none
real(8):: det
	
	Coe = 1.d0
	call puls_eqt_solve(det, .true.)	!	Find Coe to write in later puls_eqt_solve

	Omega_sq = afreq**2 * R0**3/ Grav_Const/ M0
	write(*,*) "angular freq = ", afreq
	write(*,*) "freq = ", afreq/2.d0/pi
	write(*,*) "Omega_sq", Omega_sq
	write(*,*) "=========================================================="

	call puls_eqt_solve(det, .true.)

endsubroutine ei01_w_soln

!C3 A4---------------------------------------------------
function ei01_bis(x)
use Eigenvalues_Eigenvectors
implicit none
real(8) :: x, ei01_bis, det
	afreq = x
	Omega_sq = afreq**2 * R0**3/ Grav_Const/ M0

	call puls_eqt_solve(det, .false.)
	ei01_bis = det
endfunction ei01_bis

!Unit----------------------------------------------------
subroutine ei01_unit(mode)
implicit none
integer :: mode
	if (mode == 1) then
		P = P * Grav_Const/c**4
		rho = rho * Grav_Const/c**2
		M0 = M0 * Grav_Const/c**2
		mu = mu * Grav_Const/c**4
		XCo_i = XCo_i + dlog(c**4/Grav_Const) 
		XCo_f = XCo_f + dlog(c**4/Grav_Const)  
		XCr_i = XCr_i + dlog(c**4/Grav_Const)  
		XCr_f = XCr_f + dlog(c**4/Grav_Const)  
		XOc_i = XOc_i + dlog(c**4/Grav_Const)  
		XOc_f = XOc_f + dlog(c**4/Grav_Const)  
		x_Co = x_Co + dlog(c**4/Grav_Const)
		P_Co = P_Co * Grav_Const/c**4
		rho_Co = rho_Co * Grav_Const/c**2
		m_Co = m_Co * Grav_Const/c**2
		mu_Co = mu_Co * Grav_Const/c**4
		x_Cr = x_Cr + dlog(c**4/Grav_Const)
		P_Cr = P_Cr * Grav_Const/c**4
		rho_Cr = rho_Cr * Grav_Const/c**2
		m_Cr = m_Cr * Grav_Const/c**2
		mu_Cr = mu_Cr * Grav_Const/c**4
		x_Oc = x_Oc + dlog(c**4/Grav_Const)
		P_Oc = P_Oc * Grav_Const/c**4
		rho_Oc = rho_Oc * Grav_Const/c**2
		m_Oc = m_Oc * Grav_Const/c**2
		P_Ou = P_Ou * Grav_Const/c**4
		rho_Ou = rho_Ou * Grav_Const/c**2
		m_Ou = m_Ou * Grav_Const/c**2
	elseif (mode == 2) then
		P = P / (Grav_Const/c**4)
		rho = rho / (Grav_Const/c**2)
		M0 = M0 / (Grav_Const/c**2)
		mu = mu / (Grav_Const/c**4)
		XCo_i = XCo_i - dlog(c**4/Grav_Const) 
		XCo_f = XCo_f - dlog(c**4/Grav_Const)  
		XCr_i = XCr_i - dlog(c**4/Grav_Const)  
		XCr_f = XCr_f - dlog(c**4/Grav_Const) 
		XOc_i = XOc_i - dlog(c**4/Grav_Const)  
		XOc_f = XOc_f - dlog(c**4/Grav_Const) 
		x_Co = x_Co - dlog(c**4/Grav_Const)
		P_Co = P_Co / (Grav_Const/c**4)
		rho_Co = rho_Co / (Grav_Const/c**2)
		m_Co = m_Co / (Grav_Const/c**2)
		mu_Co = mu_Co / (Grav_Const/c**4)
		x_Cr = x_Cr - dlog(c**4/Grav_Const)
		P_Cr = P_Cr / (Grav_Const/c**4)
		rho_Cr = rho_Cr / (Grav_Const/c**2)
		m_Cr = m_Cr / (Grav_Const/c**2)
		mu_Cr = mu_Cr / (Grav_Const/c**4)
		x_Oc = x_Oc - dlog(c**4/Grav_Const)
		P_Oc = P_Oc / (Grav_Const/c**4)
		rho_Oc = rho_Oc / (Grav_Const/c**2)
		m_Oc = m_Oc / (Grav_Const/c**2)
		P_Ou = P_Ou / (Grav_Const/c**4)
		rho_Ou = rho_Ou / (Grav_Const/c**2)
		m_Ou = m_Ou / (Grav_Const/c**2)
	endif
endsubroutine ei01_unit

!C3 A5---------------------------------------------------
subroutine puls_eqt_solve(det, wFile)
use puls_eqt_solve_opt01
use puls_eqt_solve_opt02
use puls_eqt_solve_opt04
use puls_eqt_solve_opt07
use puls_eqt_solve_opt08
implicit none
real(8):: det
logical :: wFile
	if (pe_opt == 1) then
		call pe01_Ctrl(det, wFile)
	elseif (pe_opt == 2) then
		call pe02_Ctrl(det, wFile)
	elseif (pe_opt == 4) then
		call pe04_Ctrl(det, wFile)
	elseif (pe_opt == 7) then
		call pe07_Ctrl(det, wFile)
	elseif (pe_opt == 8) then
		call pe08_Ctrl(det, wFile)
	else
		write(*,*) "err: pe_opt"
		pause
	endif
endsubroutine puls_eqt_solve
 
endmodule eigen_freq_opt01