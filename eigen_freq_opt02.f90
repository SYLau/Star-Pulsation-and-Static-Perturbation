module eigen_freq_opt02
use global_var

contains

!C3---------------------------------------------------
subroutine ei02_Ctrl
!	- Lengthy parts serve to provide interactive IO
!	- Important parts are highlighted with " !**** "
use Root_Finding
use Minimum_Maximum
use eigen_freq_opt02_ex1
use moI_opt01
use love_opt01_P05
use love_opt01_P08
use lo_eqt
! P06 adaptive mesh, remember to change rk4 in stat_profile_hyd to odeint if using P06
implicit none
complex(8) :: OmeC_, OmeC1, OmeC2, OmeC3, gamma
complex(8) :: input
real(8) :: OmeRL, OmeRM, OmeRH, OmeIL, OmeIH, OmeRShift, fmin
real(8) :: I_, I_bar, Eta									! Dimensionless factor of moment of inertia
real(8) :: secacc = 1.D-12 
integer :: Ninput, ans, NR, NI, err_n

	write(100,*) "eigen_freq_opt02: solve eigenvalue problem by shooting, muller's method"
	pe_fcall = .true.
	zero_freq = .false.
	call ei02_unit(1)

33		write(*,*) "Methods:"
		write(*,*) "<<LD2 Root Finding>>"
		write(*,*) "1. Muller's Method"
		write(*,*) "2. 1D scanning"
		write(*,*) "3. 2D scanning"
		write(*,*) "<<Ferrari Minimum>>"
		write(*,*) "4. Scattering fitting"
		write(*,*) "5. Scattering bracketing min"
		write(*,*) "6. Scattering scanning"
		write(*,*) "7. Scattering scanning phase angle"
		write(*,*) "8. Scattering input frequency"
		write(*,*) "<<Static Perturbation>>"
		write(*,*) "9. Love Number"
		write(*,*) "0. Back"
		write(*,*) "input: "
		read(*,*) ans
		write(*,*) "=========================================================="
		
		if (ans == 1) then
			gr_opt = 1
34			write(*,*) "Muller's Method Options: (note: you need to input '< >' for c no)"
			write(*,*) "1. Input first trial values"
			write(*,*) "2. Input first & second trial values"
			write(*,*) "3. Input all trial values"
			write(*,*) "0. Back"
			write(*,*) "input: "
			read(*,*) ans
			write(*,*) "=========================================================="
			
			if (ans /= 1 .and. ans /= 2 .and. ans /= 3 .and. ans /= 0 ) then
				write(*,*) "invalid input"
				goto 34
			endif
			if (ans == 1) then
				Ninput = 1
			elseif (ans == 2) then
				Ninput = 2
			elseif (ans == 3) then
				Ninput = 3
			elseif (ans == 0) then
				goto 33
			endif

35				write(*,*) "Choose your input"
				write(*,*) "1. angular frequency"
				write(*,*) "2. frequency"
				write(*,*) "0. back"
				read(*,*) ans
				write(*,*) "=========================================================="

				if (ans == 0) then
					goto 34
				elseif (ans /= 1 .and. ans /= 2 ) then
					write(*,*) "invalid input"
					goto 35
				endif
				if (Ninput == 1) then
					write(*,*) "input first trial"
					read(*,*) input
					OmeC1 = freq_conv(ans, input)

					!********** Shooting **********
					call ei02_shoot(OmeC1, secacc, err_n)
				elseif (Ninput == 2) then
					write(*,*) "input first trial"
					read(*,*) input
					OmeC1 = freq_conv(ans, input)

					write(*,*) "input second trial"
					read(*,*) input
					OmeC2 = freq_conv(ans, input)

					!********** Shooting **********
					call ei02_shoot(OmeC1, secacc, err_n, optx2 = OmeC2)
				elseif (Ninput == 3) then
					write(*,*) "input first trial"
					read(*,*) input
					OmeC1 = freq_conv(ans, input)

					write(*,*) "input second trial"
					read(*,*) input
					OmeC2 = freq_conv(ans, input)

					write(*,*) "input third trial"
					read(*,*) input
					OmeC3 = freq_conv(ans, input)

					!********** Shooting **********
					call ei02_shoot(OmeC1, secacc, err_n, optx2 = OmeC2, optx3 = OmeC3)
				endif
				if (err_n /= 0) then 
					goto 33
				else
					write(100,*) "ei02_shoot: passed"
				endif
		elseif (ans == 2) then
				gr_opt = 1
36				write(*,*) "Choose your input"
				write(*,*) "1. angular frequency"
				write(*,*) "2. frequency"
				write(*,*) "0. back"
				read(*,*) ans
				write(*,*) "=========================================================="
				
				if (ans == 0) then
					goto 33
				elseif (ans /= 1 .and. ans /= 2 ) then
					write(*,*) "invalid input"
					goto 36
				endif
				write(*,*) "input lower bound"
				read(*,*) input
				OmeRL = cdsqrt(freq_conv(ans, input))

				write(*,*) "input upper bound"
				read(*,*) input
				OmeRH = cdsqrt(freq_conv(ans, input))
				
				write(*,*) "input no of steps"
				read(*,*) NR

				write(*,*) "input imaginary part"
				read(*,*) input
				OmeIL = cdsqrt(freq_conv(ans, input))

				write(*,*) "=========================================================="
			!********** Scanning **********
			call ei02_scan1D(NR, OmeRL, OmeRH, OmeIL)
			write(100,*) "ei02_scan1D: passed"
			write(*,*) "Scanning complete!"
			print *, char(7)	! "BEEP"
			write(*,*) "=========================================================="
			goto 33
		elseif (ans == 3) then
				gr_opt = 1
37				write(*,*) "Choose your input"
				write(*,*) "1. angular frequency"
				write(*,*) "2. frequency"
				write(*,*) "0. back"
				read(*,*) ans
				write(*,*) "=========================================================="
				
				if (ans == 0) then
					goto 33
				elseif (ans /= 1 .and. ans /= 2 ) then
					write(*,*) "invalid input"
					goto 37
				endif
				write(*,*) "input real part lower bound"
				read(*,*) input
				OmeRL = cdsqrt(freq_conv(ans, input))

				write(*,*) "input real part upper bound"
				read(*,*) input
				OmeRH = cdsqrt(freq_conv(ans, input))
	
				write(*,*) "input imaginary part lower bound"
				read(*,*) input
				OmeIL = cdsqrt(freq_conv(ans, input))

				write(*,*) "input imaginary part upper bound"
				read(*,*) input
				OmeIH = cdsqrt(freq_conv(ans, input))
							
				write(*,*) "input real part no of steps"
				read(*,*) NR

				write(*,*) "input imaginary part no of steps"
				read(*,*) NI

				write(*,*) "=========================================================="
			!********** Scanning **********
			call ei02_scan2D(NR, NI, OmeRL, OmeRH, OmeIL, OmeIH)
			write(100,*) "ei02_scan2D: passed"
			write(*,*) "Scanning complete!"
			print *, char(7)	! "BEEP"
			write(*,*) "=========================================================="
			goto 33
		elseif (ans == 4) then
			gr_opt = 2
38				write(*,*) "Choose your input, input real values"
				write(*,*) "1. angular frequency"
				write(*,*) "2. frequency"
				write(*,*) "0. back"
				read(*,*) ans
				write(*,*) "=========================================================="

				if (ans == 0) then
					goto 33
				elseif (ans /= 1 .and. ans /= 2 ) then
					write(*,*) "invalid input"
					goto 38
				endif

				write(*,*) "input lower bound"
				read(*,*) input
				OmeRL = cdsqrt(freq_conv(ans, input))

				write(*,*) "input upper bound"
				read(*,*) input
				OmeRH = cdsqrt(freq_conv(ans, input))
				
				write(*,*) "input no of steps"
				read(*,*) NR
				
				write(*,*) "input frequency shift (approx value of freq for min)"
				read(*,*) input
				OmeRShift = cdsqrt(freq_conv(ans, input))
				
				!********** Fitting **********
				call ei02_sc_fit(NR, OmeRL, OmeRH, OmeRShift, OmeC_, err_n)

				!********** Finding frequency imaginary part **********
				OmeC_sq = OmeC_**2

				write(*,*) "=========================================================="
				write(*,*) "Scattering method:"
				write(*,*) "OmeC_sq = ", OmeC_sq*c**2
				write(*,*) "Re(OmeC) = ", dreal(cdsqrt(OmeC_sq*c**2))
				write(*,*) "Re(freq) = ", dreal(cdsqrt(OmeC_sq*c**2)/2.d0/pi)
				write(*,*) "Im(OmeC) = ", dimag(cdsqrt(OmeC_sq*c**2))
				write(*,*) "Im(freq) = ", dimag(cdsqrt(OmeC_sq*c**2)/2.d0/pi)
				write(*,*) "=========================================================="
				
				if (err_n /= 0) then 
					goto 33
				else
					write(100,*) "ei02_fitting: passed"
				endif
			
		elseif (ans == 5) then
			gr_opt = 2
39				write(*,*) "Choose your input, input real values"
				write(*,*) "1. angular frequency"
				write(*,*) "2. frequency"
				write(*,*) "0. back"
				read(*,*) ans
				write(*,*) "=========================================================="

				if (ans == 0) then
					goto 33
				elseif (ans /= 1 .and. ans /= 2 ) then
					write(*,*) "invalid input"
					goto 39
				endif
				
				write(*,*) "input first value"
				read(*,*) input
				OmeRL = dreal(freq_conv(ans, input))

				write(*,*) "input second value"
				read(*,*) input
				OmeRM = dreal(freq_conv(ans, input))

				write(*,*) "input third value"
				read(*,*) input
				OmeRH = dreal(freq_conv(ans, input))

				!********** Finding minimum for real frequency **********
				OmeC_sq = brent(OmeRL,OmeRM,OmeRH,ei02_infinity_cond_R,secacc,fmin)

				!********** Finding frequency imaginary part **********
				call ei02_sc_imag(OmeC_, err_n)
				OmeC_sq = OmeC_**2

				write(*,*) "=========================================================="
				write(*,*) "Scattering method:"
				write(*,*) "OmeC_sq = ", OmeC_sq*c**2
				write(*,*) "Re(OmeC) = ", dreal(cdsqrt(OmeC_sq*c**2))
				write(*,*) "Re(freq) = ", dreal(cdsqrt(OmeC_sq*c**2)/2.d0/pi)
				write(*,*) "Im(OmeC) = ", dimag(cdsqrt(OmeC_sq*c**2))
				write(*,*) "Im(freq) = ", dimag(cdsqrt(OmeC_sq*c**2)/2.d0/pi)
				write(*,*) "=========================================================="

				write(100,*) "ei02_bracket: passed"

		elseif (ans == 6) then
				gr_opt = 2
40				write(*,*) "Choose your input"
				write(*,*) "1. angular frequency"
				write(*,*) "2. frequency"
				write(*,*) "0. back"
				read(*,*) ans
				write(*,*) "=========================================================="
				
				if (ans == 0) then
					goto 33
				elseif (ans /= 1 .and. ans /= 2 ) then
					write(*,*) "invalid input"
					goto 40
				endif
				write(*,*) "input lower bound"
				read(*,*) input
				OmeRL = cdsqrt(freq_conv(ans, input))

				write(*,*) "input upper bound"
				read(*,*) input
				OmeRH = cdsqrt(freq_conv(ans, input))
				
				write(*,*) "input no of steps"
				read(*,*) NR

				write(*,*) "=========================================================="
			!********** Scanning **********
			call ei02_scan1D(NR, OmeRL, OmeRH, 0.d0)
			write(100,*) "ei02_scan1D: passed"
			write(*,*) "Scanning complete!"
			print *, char(7)	! "BEEP"
			write(*,*) "=========================================================="
			goto 33

		elseif (ans == 7) then
				gr_opt = 3
41				write(*,*) "Choose your input"
				write(*,*) "1. angular frequency"
				write(*,*) "2. frequency"
				write(*,*) "0. back"
				read(*,*) ans
				write(*,*) "=========================================================="
				
				if (ans == 0) then
					goto 33
				elseif (ans /= 1 .and. ans /= 2 ) then
					write(*,*) "invalid input"
					goto 41
				endif
				write(*,*) "input lower bound"
				read(*,*) input
				OmeRL = cdsqrt(freq_conv(ans, input))

				write(*,*) "input upper bound"
				read(*,*) input
				OmeRH = cdsqrt(freq_conv(ans, input))
				
				write(*,*) "input no of steps"
				read(*,*) NR

				write(*,*) "=========================================================="
			!********** Scanning **********
			call ei02_scan1D(NR, OmeRL, OmeRH, 0.d0)
			write(100,*) "ei02_scan1D: passed"
			write(*,*) "Scanning complete!"
			print *, char(7)	! "BEEP"
			write(*,*) "=========================================================="
			goto 33

		elseif (ans == 8) then
				gr_opt = 2
42				write(*,*) "Choose your input"
				write(*,*) "1. angular frequency"
				write(*,*) "2. frequency"
				write(*,*) "0. back"
				read(*,*) ans
				write(*,*) "=========================================================="
				
				if (ans == 0) then
					goto 33
				elseif (ans /= 1 .and. ans /= 2 ) then
					write(*,*) "invalid input"
					goto 42
				endif

				write(*,*) "input freq"
				read(*,*) input
				OmeRL = cdsqrt(freq_conv(ans, input))
				
				OmeC_sq = OmeRL**2
		
		elseif (ans == 9) then
				write(*,*) "Finding k2:"
				zero_freq = .true.
				OmeC_sq = (0.d0, 0.d0 )

				if (eos_opt == 2) then
				    write(*,*) "Love Number Andersson NS model:"
				    write(*,*) "*********************************"
					eos_choice = 'NS_A11'
				    call lo01_iterate_fluid
				    call lo01_iterate_FE2
					call lo01_iterate_Solid2
					call lo01_iterate_EF
call lo01_iterate_FEF
call lo01_iterate_EFE2
!call lo01_iterate_fluid_W
!call lo01_iterate_FE_2Side_W
				elseif (eos_opt == 3) then
				    write(*,*) "Love Number Mannarelli QS model:"
				    write(*,*) "*********************************"
					eos_choice = 'QS'
				    call lo01_iterate_fluid
				    call lo01_iterate_FE2
					call lo01_iterate_Solid2
					call lo01_iterate_EF
!call lo01_iterate_EFE2

call lo01_iterate_fluid_W
call lo01_iterate_EF_W
!call lo01_iterate_EF_W_No_Slip
!call lo01_iterate_fluid_W
!call lo01_iterate_FE_2Side_W

			    else
				     pause 'err: ei02 opt 09 EOS mismatch'
				endif
				
				!********** Calculate Moment of Inertia **********
				call mi01_Ctrl(I_, I_bar, Eta)

			write(*,*) "=========================================================="
			goto 33
				
		elseif (ans == 0) then
			write(100,*) "=========================================================="
				
			!********** Exit Subroutine**********
			call ei02_unit(2)
			return
		else
			write(*,*) "invalid input"
			goto 33
		endif
		
		!********** Write the solution into txt files **********
		call ei02_w_soln
		write(100,*) "ei02_soln: passed"
		
		!********** Calculate Moment of Inertia **********
		call mi01_Ctrl(I_, I_bar, Eta)

		if (pe_opt == 3) then
			call ei02_o_soln_C1L4R4
			write(100,*) "ei02_o_soln: passed"
		elseif (pe_opt == 5) then
			call ei02_o_soln_C2L4R6
			write(100,*) "ei02_o_soln: passed"
		elseif (pe_opt == 6) then
			call ei02_o_soln_C3L4R6
			write(100,*) "ei02_o_soln: passed"
		elseif (pe_opt == 9) then
			call ei02_o_soln_C1L6R6
			write(100,*) "ei02_o_soln: passed"
		else
			write(*,*) "err: ei02: pe_opt"
			pause
		endif
		
		zero_freq = .false.

		write(*,*) "Solution Written!"
		print *, char(7)	! "BEEP"
		write(*,*) "=========================================================="
		goto 33

endsubroutine ei02_Ctrl





! ___     ___     ___     ___     ___     ___     ___     ___     ___     ___    
! \  \    \  \    \  \    \  \    \  \    \  \    \  \    \  \    \  \    \  \   
!  \  \    \  \    \  \    \  \    \  \    \  \    \  \    \  \    \  \    \  \  
!   )  )    )  )    )  )    )  )    )  )    )  )    )  )    )  )    )  )    )  ) 
!  /  /    /  /    /  /    /  /    /  /    /  /    /  /    /  /    /  /    /  /  
! /__/    /__/    /__/    /__/    /__/    /__/    /__/    /__/    /__/    /__/   
!__________                __    .___.___ 
!\______   \_____ ________/  |_  |   |   |
! |     ___/\__  \\_  __ \   __\ |   |   |
! |    |     / __ \|  | \/|  |   |   |   |
! |____|    (____  /__|   |__|   |___|___|

                    
!C3 A1---------------------------------------------------
subroutine ei02_unit(mode)
implicit none
integer :: mode
	if (mode == 1) then
		P = P * Grav_Const/c**4
		rho = rho * Grav_Const/c**2
		m = m * Grav_Const/c**2
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
		m = m / (Grav_Const/c**2)
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
endsubroutine ei02_unit

!C3 A2---------------------------------------------------
function freq_conv(i, input)
implicit none
complex(8) :: freq_conv, input
integer :: i
	if (i == 1) freq_conv = (input/c)**2
	if (i == 2) freq_conv = (input*(2.d0*pi/c))**2
endfunction freq_conv

!C3 A3---------------------------------------------------
subroutine ei02_shoot(x1, acc, err_n, optx2, optx3)
use Root_Finding
use Format_IO
use FWrite
implicit none
complex(8) :: x1, x2, x3
complex(8), optional :: optx2, optx3
real(8) :: acc, inix2, inix3
integer :: i, err_n
	inix2 = dreal(x1) * 1.d-4
	inix3 = cdabs(x1) * 1.d-4
	if (present(optx2)) then 
		x2 = optx2
	else 
		x2 = x1 + dcmplx(inix2, 0.d0)
	endif

	if (present(optx3)) then 
		x3 = optx3
	else 
		x3 = x1 + dcmplx(0.d0, inix3)
	endif

	OmeC_sq = rtmuller_c(ei02_muller, x1, x2, x3, acc, JMAX=muller_n_max, ERR = err_n)	! error = 0 if no problem

	write(*,*) "=========================================================="
	write(*,*) "OmeC_sq = ", OmeC_sq*c**2
	write(*,*) "Mean Density sqrt = ", (M0/R0**3)**0.5d0
	write(*,*) "=========================================================="
endsubroutine ei02_shoot

!C3 A4---------------------------------------------------
function ei02_muller(x_)
!	input OmeC_sq and return gamma
implicit none
complex(8) :: x_, beta, gamma, ei02_muller
	OmeC_sq = x_

	call puls_eqt_solve(beta, gamma, .false.)

	write(*,*) "OmeC_sq = ", OmeC_sq*c**2
	write(*,*) "Re(OmeC) = ", dreal(cdsqrt(OmeC_sq*c**2))
	write(*,*) "Re(freq) = ", dreal(cdsqrt(OmeC_sq*c**2)/2.d0/pi)
	write(*,*) "Im(OmeC) = ", dimag(cdsqrt(OmeC_sq*c**2))
	write(*,*) "Im(freq) = ", dimag(cdsqrt(OmeC_sq*c**2)/2.d0/pi)
	write(*,*) "gamma = ", gamma
	write(*,*) "=========================================================="
	ei02_muller = gamma
endfunction ei02_muller

!C3 A5---------------------------------------------------
subroutine ei02_scan1D(N, OmeRL, OmeRH, OmeI)
use Format_IO
use FWrite
implicit none
real(8) :: OmeRL, OmeRH, OmeI
real(8) :: d_OmeR, ans_R
complex(8) :: x_, ans
integer :: N, i

	open(61, file = 'data/scanning1D.txt', status='replace')
open(1234, file = 'data/scanning1D_Uz.txt', status='replace')
	d_OmeR = (OmeRH-OmeRL)/N
	do i = 1, N
		x_ = dcmplx(i*d_OmeR + OmeRL, OmeI)**2
		if (gr_opt == 1) then 
			ans = ei02_infinity_cond(x_)
			ans_R = dlog(cdabs( ans ))
		elseif (gr_opt == 2) then 
			ans = ei02_infinity_cond_R(dreal(x_))
			ans_R = dlog(cdabs( ans ))
		elseif (gr_opt == 3) then
			ans = ei02_infinity_cond_R(dreal(x_))
			ans_R = dreal( ans )
		endif
		call write3R8(61, FR8, OmeI*c/2.d0/pi, dreal(cdsqrt(x_*c**2)/2.d0/pi), ans_R)
	enddo
	close(61)
close(1234)
endsubroutine ei02_scan1D

!C3 A6---------------------------------------------------
subroutine ei02_scan2D(NR, NI, OmeRL, OmeRH, OmeIL, OmeIH)
use Format_IO
use FWrite
implicit none
real(8) :: OmeRL, OmeRH, OmeIL, OmeIH
real(8) :: d_OmeR, d_OmeI
complex(8) :: x_, ans
integer :: NR, NI, i, j

	open(61, file = 'data/scanning2D.txt', status='replace')
	d_OmeR = (OmeRH-OmeRL)/NR
	d_OmeI = (OmeIH-OmeIL)/NI

	do i = 1, NR
		do j = 1, NI
			x_ = dcmplx(i*d_OmeR + OmeRL, j*d_OmeI + OmeIL)**2
			ans = ei02_infinity_cond(x_)
			call write5R8(61, FR8, dreal(cdsqrt(x_*c**2)/2.d0/pi), dimag(cdsqrt(x_*c**2)/2.d0/pi), dreal(ans), dimag(ans), dlog(cdabs(ans)))
		enddo
	enddo
	close(61)
!	scan near (2325, 50)
endsubroutine ei02_scan2D

!C3 A7---------------------------------------------------
function ei02_infinity_cond(x_)
implicit none
complex(8) :: x_, beta, gamma, ei02_infinity_cond
	OmeC_sq = x_
	call puls_eqt_solve(beta, gamma, .false.)

	if (gr_opt == 1) then
		write(*,*) "OmeC_sq = ", OmeC_sq*c**2
		write(*,*) "Re(freq) = ", dreal(cdsqrt(OmeC_sq*c**2)/2.d0/pi)
		write(*,*) "Im(freq) = ", dimag(cdsqrt(OmeC_sq*c**2)/2.d0/pi)
		!write(*,*) "gamma/beta = ", gamma/beta
		write(*,*) "|gamma| = ", cdabs(gamma)
		write(*,*) "=========================================================="
		!ei02_infinity_cond = gamma/beta
		ei02_infinity_cond = cdabs(gamma)
	else
		write(*,*) "err: Scattering, gr_opt"
	endif
endfunction ei02_infinity_cond

!C3 A8---------------------------------------------------
function ei02_infinity_cond_R(x_)
implicit none
complex(8) :: beta, gamma
real(8) :: del, cos_del, sin_del, Amp
real(8) :: x_, ei02_infinity_cond_R
	OmeC_sq = dcmplx(x_,0.d0)
	call puls_eqt_solve(beta, gamma, .false.)

	if (gr_opt == 2) then
		write(*,*) "OmeC_sq = ", OmeC_sq*c**2
		write(*,*) "Re(freq) = ", dreal(cdsqrt(OmeC_sq*c**2)/2.d0/pi)
		write(*,*) "Im(freq) = ", dimag(cdsqrt(OmeC_sq*c**2)/2.d0/pi)
		write(*,*) "gamma^2 + beta^2 = ", gamma**2 + beta**2
		write(*,*) "=========================================================="
		ei02_infinity_cond_R = dreal(gamma**2 + beta**2)
	elseif (gr_opt == 3) then
		Amp = dsqrt( dreal(beta**2) + dreal(gamma**2))
		cos_del = dreal(beta/Amp)
		sin_del = dreal(-gamma/Amp)
		del = dabs( datan(sin_del/cos_del) )
		write(*,*) "OmeC_sq = ", OmeC_sq*c**2
		write(*,*) "Re(freq) = ", dreal(cdsqrt(OmeC_sq*c**2)/2.d0/pi)
		write(*,*) "Im(freq) = ", dimag(cdsqrt(OmeC_sq*c**2)/2.d0/pi)
		write(*,*) "atan(-beta/ gamma)/pi = ", del/pi
		write(*,*) "=========================================================="

		if (sin_del * cos_del >= 0.d0 ) then
		ei02_infinity_cond_R = del/pi
		else
		ei02_infinity_cond_R = (pi - del)/pi
		endif
	else
		write(*,*) "err: Scattering, gr_opt"
	endif
endfunction ei02_infinity_cond_R

!C3 A9---------------------------------------------------
subroutine ei02_sc_fit(N, OmeRL, OmeRH, OmeRShift, OmeC_, err_n)
use Fitting
use Format_IO
use FWrite
implicit none
integer :: N, i, err_n
real(8) :: OmeRL, OmeRH, OmeI, OmeRShift
real(8) :: d_OmeR
real(8) :: x(1:N), y(1:N), Coef_Shift(1:3), Coef(1:3), C_, Ome_0, Ome_i
real(8) :: x_, y1_, y2_
complex(8) :: OmeC_
	err_n = 0
	d_OmeR = (OmeRH-OmeRL)/N
	do i = 1, N
		x_ = (i*d_OmeR + OmeRL)**2
		y(i) = ei02_infinity_cond_R(x_)
		!¡¸Shift the curve to origin for better accuracy
		x(i) = (i*d_OmeR + OmeRL) - OmeRShift
	enddo
	call Poly_Para_Fit(n, x, y, Coef_Shift)

	Coef(1) = Coef_Shift(1) + (-Coef_Shift(2)*OmeRShift) + Coef_Shift(3)*OmeRShift**(2.d0)
	Coef(2) = Coef_Shift(2) - 2.d0*Coef_Shift(3)*OmeRShift
	Coef(3) = Coef_Shift(3)
	
	write(*,*) "y = c + bx + ax^2"
	write(*,*) "c = ", Coef(1)
	write(*,*) "b = ", Coef(2)
	write(*,*) "a = ", Coef(3)
	write(*,*) "=========================================================="
	
	open(01, file='data/sc_fit.txt', status='replace')
		do i = 1, N
			x_ = x(i)+ OmeRShift
			y1_ = Coef_Shift(1) + Coef_Shift(2)*(x(i)) + Coef_Shift(3)*(x(i))**(2.d0)
			y2_ = Coef(1) + Coef(2)*(x_) + Coef(3)*(x_)**(2.d0)			! Without shifting; Sometimes exceed the 15 digit limit of double precision 
			!¡¸compare the fitted curves with original y
			call Write4R8(01,FR8L, x_/2.d0/pi*c, y(i), y1_, y2_)		! y2 gives the approximate accuracy without shifting, ie tells if shifting is necessary or not
		enddo
	close(01)
	C_ = Coef(3)
	Ome_0 = Coef(2)/(-2.d0*C_)
	if (Coef(1)/C_ - Ome_0**2 < 0.d0) then
		pause 'err: Scattering fit; minimum below zero'
		err_n = 1
		return
	endif
	Ome_i = dsqrt(Coef(1)/C_ - Ome_0**2)
	write(*,*) "Parabola width = ", C_
	write(*,*) "Re(freq) = ", Ome_0/2.d0/pi*c
	write(*,*) "Im(freq) = ", Ome_i/2.d0/pi*c
	write(*,*) "=========================================================="
	OmeC_ = dcmplx(Ome_0,Ome_i)
endsubroutine ei02_sc_fit

!C3 A10---------------------------------------------------
subroutine ei02_sc_imag(OmeC_, err_n)
!	Fits the result to a parabola: C [(Ome - Ome_0)^2 + Ome_i^2]
!	Reference: Ferrari 1991
use FWrite
use Format_IO
implicit none
complex(8) :: beta, gamma, OmeC_
real(8) :: f(-1:1), tol_check(0:1)
real(8) :: Ome_0, Ome_i, delta
integer :: i, err_n
real(8), parameter :: factor = 5.d-5, tol = 5.d-2, iMax = 20
	
	err_n = 0
	if (dimag(OmeC_sq) /= 0) then
		write(*,*) "err: scattering imag Omega not zero"
		pause
		err_n = 1
		return
	endif
!open(01, file="data/debug_im.txt", status="replace")
	Ome_0 = dsqrt(dreal(OmeC_sq))

	delta = Ome_0 * factor

	do i = 1, iMax
		OmeC_sq = dcmplx(Ome_0**2,0.d0)
		call puls_eqt_solve(beta, gamma, .false.)
		f(0) = dreal(gamma**2 + beta**2)

		OmeC_sq = dcmplx((Ome_0-delta)**2,0.d0)
		call puls_eqt_solve(beta, gamma, .false.)
		f(-1) = dreal(gamma**2 + beta**2)

		OmeC_sq = dcmplx((Ome_0+delta)**2,0.d0)
		call puls_eqt_solve(beta, gamma, .false.)
		f(1) = dreal(gamma**2 + beta**2)

		tol_check(1) = dabs((f(1)-f(-1))/((f(1)+f(-1))/2.d0))

!call Write6R8(01, FR8, Ome_0*c/2.d0/pi, (Ome_0 - delta)*c/2.d0/pi, (Ome_0 + delta)*c/2.d0/pi, f(0), f(-1), f(1))	
!write(*,*) Ome_0
!write(*,*) delta
!write(*,*) f(1), f(0), f(-1)
!write(*,*) tol_check(0), tol_check(1)
!pause

		
		if (dabs(f(1)) < dabs(f(0)) .or. dabs(f(-1)) < dabs(f(0))) then
			write(*,*) "err: Scattering f(0) not minimum"
			pause
			err_n = 1
			exit
		endif
		if (tol_check(1) < tol .and. tol_check(1) >= tol_check(0)) then
			exit
		else
			delta = delta/2.d0
			tol_check(0) = tol_check(1)
			cycle
		endif
		
	enddo

	Ome_i = dsqrt(delta**2/(f(1)/f(0) - 1.d0))
	OmeC_ = dcmplx(Ome_0,Ome_i)
!write(*,*) "C = ", f(1)/(delta**2+Ome_i**2)
!write(*,*) "C = ", f(0)/(Ome_i**2)
!pause
!close(01)
endsubroutine ei02_sc_imag

!C3 A11---------------------------------------------------
subroutine ei02_w_soln
!	write the solution into intermediate txt files
implicit none
complex(8) :: beta, gamma, afreqC

	CoeC = 1.d0
	call puls_eqt_solve(beta, gamma, .true.)	!	Find CoeC to write in later puls_eqt_solve

	afreqC = OmeC_sq**(0.5d0)
	write(*,*) "angular freq = ", afreqC*c
	write(*,*) "freq = ", afreqC/2.d0/pi*c
	write(*,*) "=========================================================="

	call puls_eqt_solve(beta, gamma, .true.)

endsubroutine ei02_w_soln

!C3 A12---------------------------------------------------
subroutine ei02_scan1D_Love(N, OmeRL, OmeRH, OmeI)
use Format_IO
use FWrite
implicit none
real(8) :: OmeRL, OmeRH, OmeI
real(8) :: d_OmeR, ans_R
complex(8) :: x_, ans
complex(8) :: beta, gamma
integer :: N, i

	open(61, file = 'data/scanning1D_Love.txt', status='replace')
	d_OmeR = (OmeRH-OmeRL)/N
	do i = 1, N
		x_ = dcmplx(i*d_OmeR + OmeRL, OmeI)**2
		OmeC_sq = x_

		call puls_eqt_solve(beta, gamma, .false.)
		call puls_eqt_solve(beta, gamma, .true.)

		write(*,*) "OmeC_sq = ", OmeC_sq*c**2
		write(*,*) "Re(freq) = ", dreal(cdsqrt(OmeC_sq*c**2)/2.d0/pi)
		write(*,*) "Im(freq) = ", dimag(cdsqrt(OmeC_sq*c**2)/2.d0/pi)
		write(*,*) "k2 = ", k2
		write(*,*) "|gamma| = ", cdabs(gamma)
		write(*,*) "=========================================================="

		call write4R8(61, FR8, OmeI*c/2.d0/pi, dreal(cdsqrt(x_*c**2)/2.d0/pi), dreal(k2), cdabs(gamma))
	enddo
	close(61)
endsubroutine ei02_scan1D_Love

!C3 A13---------------------------------------------------
subroutine puls_eqt_solve(beta, gamma, wFile)
use puls_eqt_solve_opt03
use puls_eqt_solve_opt05
use puls_eqt_solve_opt06
use puls_eqt_solve_opt09
implicit none
complex(8) :: beta, gamma
logical :: wFile
	if (pe_opt == 3) then
		call pe03_Ctrl(beta, gamma, wFile)
	elseif (pe_opt == 5) then
		call pe05_Ctrl(beta, gamma, wFile)
	elseif (pe_opt == 6) then
		call pe06_Ctrl(beta, gamma, wFile)
	elseif (pe_opt == 9) then
		call pe09_Ctrl(beta, gamma, wFile)
	else
		write(*,*) "err: ei02: pe_opt"
		pause
	endif
endsubroutine puls_eqt_solve

endmodule eigen_freq_opt02