module EOS
!	Contains the subroutines for using EOS table or analytic EOS
use global_var
contains

	!...............................................................
	!	Part 1:
	!	Equation of states
	!...............................................................
	function sp_eos(f, x_)
	use Interpolation
	use FRead
	implicit none
	real(8), dimension(1:eos_tab_Nmax), save :: n_, nb_, rho_, P_
	real(8), save :: K, Ga
	real(8) :: sp_eos, x_
	integer :: i
	character (len=30) :: ans
	character (len=*) :: f
		if (eos_opt == 1) then
			if (sp_eos_fcall == .true.) then			! sp_eos_fcall is global variable
				call RFile4R8(eos_tab,eos_n, n_, nb_, rho_, P_)
				sp_eos_fcall = .false.
			endif
			if (f == 'p(rho)') then
				sp_eos = Int_ada(1, 1, eos_n, rho_, P_, x_, "interpolate p(rho)")
			elseif (f == 'nb(rho)') then
				sp_eos = Int_ada(1, 1, eos_n, rho_, nb_, x_, "interpolate nb(rho)")
			elseif (f == 'rho(p)') then
				sp_eos = Int_ada(1, 1, eos_n, P_, rho_, x_, "interpolate rho(p)")
			else
				write(*,*) "err: static profile eos table"
				pause
			endif
		elseif (eos_opt == 2) then
			!¡¸Polytropic EOS with ref to Ferrari 2003; uses global_var
			if (sp_eos_fcall == .true.) then
				Ga = (poly_n + 1.d0)/poly_n
				K =	Grav_Const*poly_K*1.d10/(1.d0 + drho)**Ga			! 1.d10 accounts for unit conversion from km^2 to cm^2
				if (drho /= 0.d0 ) P_t = K*(rho_t*(1.d0 + drho))**Ga
				write(*,*) "Polytropic EOS:"
				write(*,*) "K: ", K
				write(*,*) "drho: ", drho
				write(*,*) "rho_t: ", rho_t
				write(*,*) "P_t changed to: ", P_t
				write(*,*) "=========================================================="
				sp_eos_fcall = .false.
			endif
			if (f == 'p(rho)') then
				sp_eos = K*(x_)**Ga
			elseif (f == 'nb(rho)') then
				sp_eos = 0.d0
			elseif (f == 'rho(p)') then
				if (x_ >= P_t) sp_eos = (x_/K)**(1.d0/Ga)
				if (x_ < P_t) sp_eos = (x_/K)**(1.d0/Ga)/(1.d0 + drho)		
			else
				write(*,*) "err: static profile polytropic eos"
				pause
			endif
		elseif (eos_opt == 3) then
			!¡¸Quark Star EOS with ref to Pagliaroli 2014; uses global_var
			if (f == 'p(rho)') then
				sp_eos = sp_qs_P_rho(x_)
			elseif (f == 'nb(rho)') then
				sp_eos = 0.d0
			elseif (f == 'rho(p)') then
				sp_eos = sp_qs_rho_P(x_)		
			else
				write(*,*) "err: static profile quark star eos"
				pause
			endif
		elseif (eos_opt == 4) then
			!¡¸Hybrid Star: MIT Bag Model Inner Core, Neutron Star EOS For the rest
			!¡¸Quark Star EOS with ref to T Klahn 2007; uses global_var
			if (sp_eos_fcall == .true.) then			! sp_eos_fcall is global variable
				call RFile4R8(eos_tab,eos_n, n_, nb_, rho_, P_)
				do i = 1, eos_n
					if (dabs(rho_t-rho_(i))/rho_(i) <= 1.d-3) then
						P_t = P_(i)
						exit
					endif
					if (i == eos_n) then
						write(*,*) "err: eos_opt 4, can't find rho_t"
						pause
					endif
				enddo
				
				QS_a4 = 0.7d0
				QS_a2 = 0.d0
				B_eff = sp_qs_Beff(P_t, rho_t*(1.d0 + drho))
				write(*,*) "B_eff^(1/4)", B_eff**(1.d0/4.d0)
				write(*,*) "rho-", sp_qs_rho_P(P_t)
				write(*,*) "rho+", rho_t
				write(*,*) "=========================================================="

				sp_eos_fcall = .false.
			endif
			if (f == 'p(rho)') then
				if (x_ < rho_t) sp_eos = Int_ada(1, 1, eos_n, rho_, P_, x_, "interpolate p(rho)")
				if (x_ == rho_t) sp_eos = P_t
				if (x_ > rho_t) sp_eos = sp_qs_P_rho(x_)
			elseif (f == 'nb(rho)') then
				sp_eos = 0.d0
			elseif (f == 'rho(p)') then
				if (x_ < P_t) sp_eos = Int_ada(1, 1, eos_n, P_, rho_, x_, "interpolate rho(p)")
				if (x_ >= P_t) sp_eos = sp_qs_rho_P(x_)
			else
				write(*,*) "err: static profile eos table"
				pause
			endif
		elseif (eos_opt == 5) then
			!¡¸Hybrid Star with Polytropic Inner Core
			!¡¸Polytropic EOS with ref to Ferrari 2003; uses global_var
			if (sp_eos_fcall == .true.) then			! sp_eos_fcall is global variable
				call RFile4R8(eos_tab,eos_n, n_, nb_, rho_, P_)
				do i = 1, eos_n
					if (dabs(rho_t-rho_(i))/rho_(i) <= 1.d-3) then
						P_t = P_(i)
						exit
					endif
					if (i == eos_n) then
						write(*,*) "err: eos_opt 5"
						pause
					endif
				enddo
	
				Ga = (poly_n + 1.d0)/poly_n
				K =	P_t/ ( rho_t*(1.d0 + drho)	)**Ga		! 1.d10 accounts for unit conversion from km^2 to cm^2
				write(*,*) "Polytropic EOS:"
				write(*,*) "K: ", K
				write(*,*) "drho: ", drho
				write(*,*) "rho_t: ", rho_t
				write(*,*) "P_t: ", P_t
				write(*,*) "=========================================================="
				sp_eos_fcall = .false.
			endif

			if (f == 'p(rho)') then
				if (x_ < rho_t) sp_eos = Int_ada(1, 1, eos_n, rho_, P_, x_, "interpolate p(rho)")
				if (x_ == rho_t) sp_eos = P_t
				if (x_ > rho_t) sp_eos = K*(x_)**Ga
			elseif (f == 'nb(rho)') then
				sp_eos = 0.d0
			elseif (f == 'rho(p)') then
				if (x_ < P_t) sp_eos = Int_ada(1, 1, eos_n, P_, rho_, x_, "interpolate rho(p)")	
				if (x_ >= P_t) sp_eos = (x_/K)**(1.d0/Ga)
			else
				write(*,*) "err: static profile polytropic eos"
				pause
			endif
		elseif (eos_opt == 6) then
			!¡¸Hybrid Star: MIT Bag Model Inner Core, Neutron Star EOS For the rest
			!¡¸Quark Star EOS with ref to T Klahn 2007; uses global_var
			if (sp_eos_fcall == .true.) then			! sp_eos_fcall is global variable
				call RFile4R8(eos_tab,eos_n, n_, nb_, rho_, P_)
				i = 2
				do while (i <= eos_n-1)
					if ( dabs(rho_t-rho_(i)) <= dabs(rho_t-rho_(i-1)) .and. dabs(rho_t-rho_(i)) <= dabs(rho_t-rho_(i+1)))  then
						write(*,*) "Is it rho_t? (y/n)", rho_(i)
						read(*,*) ans
						if (ans == 'y' .or. ans == 'Y') then
							P_t = P_(i)
							exit
						else
							write(*,*) "What's your rho_t then?"
							read(*,*) rho_t
							i = 2
							cycle
						endif
					endif
					if (i == eos_n-1) then
						write(*,*) "err: eos_opt 6, can't find P_t"
						pause
					endif
					i = i + 1
				enddo
				QS_a4 = 0.7d0
				QS_a2 = 0.d0

				write(*,*) "B_eff^(1/4)", B_eff**(1.d0/4.d0)
				write(*,*) "rho-", sp_qs_rho_P(P_t)
				write(*,*) "rho+", rho_(i)
				write(*,*) "drho", (sp_qs_rho_P(P_t)-rho_(i))/rho_(i)
				write(*,*) "=========================================================="

				sp_eos_fcall = .false.
			endif
			if (f == 'p(rho)') then
				if (x_ < rho_t) sp_eos = Int_ada(1, 1, eos_n, rho_, P_, x_, "interpolate p(rho)")
				if (x_ == rho_t) sp_eos = P_t
				if (x_ > rho_t) sp_eos = sp_qs_P_rho(x_)
			elseif (f == 'nb(rho)') then
				sp_eos = 0.d0
			elseif (f == 'rho(p)') then
				if (x_ < P_t) sp_eos = Int_ada(1, 1, eos_n, P_, rho_, x_, "interpolate rho(p)")
				if (x_ >= P_t) sp_eos = sp_qs_rho_P(x_)
			else
				write(*,*) "err: static profile eos table"
				pause
			endif
		else
			write(*,*) "err: eos options"
			pause
		endif

	endfunction sp_eos

	function sp_qs_P_rho(x_)
	!	Using natural unit; h_bar = 1, c = 1
	implicit none
	real(8) :: sp_qs_P_rho, x_, E_
	real(8) :: C1_, A_, B_, C_, mu_2
	real(8), parameter :: e_mks = 1.60217657d-19
		
		!¡¸Quark Star EOS with ref to Pagliaroli 2014; 
		!¡¸Convert to natural units: h_bar = 1, c = 1, P in dimension [E/L^3] = [E^4/((h_bar*c)^3)]
		E_ = (h_bar*c)**3 * x_/((1.d13)*e_mks)**4 * c**2	! Energy density E_

		C1_ = 3.d0/4.d0/pi**(2.d0)
		A_ = C1_ * QS_a4
		B_ = - C1_ * QS_a2
		C_ = -B_eff
		mu_2 = -B_/6.d0/A_ + dsqrt(B_**2 + 12.d0*A_*(C_+E_))/6.d0/A_

		sp_qs_P_rho = (A_*mu_2**2 + B_*mu_2 + C_) * ((1.d13)*e_mks)**4 /(h_bar*c)**3

	endfunction sp_qs_P_rho

	function sp_qs_rho_P(x_)
	!	Using natural unit; h_bar = 1, c = 1
	implicit none
	real(8) :: sp_qs_rho_P, x_, P_
	real(8) :: C1_, A_, B_, C_, mu_2
	real(8), parameter :: e_mks = 1.60217657d-19
		
		!¡¸Quark Star EOS with ref to Pagliaroli 2014; 
		!¡¸Convert to natural units: h_bar = 1, c = 1, P in dimension [E/L^3] = [E^4/((h_bar*c)^3)]
		P_ = (h_bar*c)**3 * x_/((1.d13)*e_mks)**4

		C1_ = 3.d0/4.d0/pi**(2.d0)
		A_ = C1_ * QS_a4
		B_ = - C1_ * QS_a2
		C_ = -B_eff

!if (x_ <= P_min*2 ) then
if (B_**2 - 4.d0*A_*(C_-P_) <= 0.d0) then
write(*,*) x_, P_
pause
endif
		mu_2 = -B_/2.d0/A_ + dsqrt(B_**2 - 4.d0*A_*(C_-P_))/2.d0/A_

		sp_qs_rho_P = (3.d0*A_*mu_2**2 + B_*mu_2 - C_) * ((1.d13)*e_mks)**4/ (h_bar*c)**3/ c**2
	endfunction sp_qs_rho_P

	function sp_qs_Beff(P_, rho_)
	implicit none
	real(8) :: sp_qs_Beff, P_, rho_
	real(8) :: P_nu, E_nu
	real(8), parameter :: e_mks = 1.60217657d-19
		P_nu = (h_bar*c)**3 * P_/((1.d13)*e_mks)**4
		E_nu = (h_bar*c)**3 * rho_/((1.d13)*e_mks)**4 * c**2

		sp_qs_Beff = (E_nu - 3.d0 * P_nu)/4.d0
		
	endfunction sp_qs_Beff
	
	!...............................................................
	!	Part 2:
	!	Compositions and Shear Modulus
	!...............................................................
	subroutine sp_comp(N1, N2)
	use interpolation
	use Format_IO
	use FWrite
	use FRead
	use verify
	implicit none
	integer :: N1, N2, j
	real(8), dimension(1:comp_tab_Nmax) :: n_, nb_, rho_, nN_, Z_, A_
	real(8) :: nN_j, Z_j

		call RFile4R8(comp_tab, comp_n, n_, rho_, nN_, Z_)
		call ver_list_range(rho_(1), rho_(comp_n), rho(N1:N2), N2-N1, "range; shear_table")
		write(*,*) "range of comp_table:"
		write(*,*) "rho_i(/g cm^-3):", rho_(1)
		write(*,*) "rho_f(/g cm^-3):", rho_(comp_n)
		write(*,*) "range of crust:"
		write(*,*) "rho_i(/g cm^-3):", rho(N1)
		write(*,*) "rho_f(/g cm^-3):", rho(N2) 
		write(*,*) "=========================================================="
		open(81, file = 'data/shear.txt', status ='replace')
		do j = N1, N2
			nN_j = Int_ada(2, 1, comp_n, rho_, nN_, rho(j), "interpolate nuclei density")
			Z_j = Int_ada(2, 1, comp_n, rho_, Z_, rho(j), "interpolate charge number")
			mu(j) = sp_shear_fcn(nN_j, Z_j)
			!call Write4R8(81, FR8, sp_r(j), mu(j), dlog10(rho(j)), dlog10(mu(j)))
			call Write4R8(81, FR8, sp_r(j), mu(j), dlog10(P(j)), dlog10(mu(j)))
		enddo
		mu = mu * mu_red
		close(81)
	endsubroutine sp_comp

	function sp_shear_fcn(nN, Z)
	implicit none
	real(8) :: nN, Z, RCell, Ga, T
	real(8) :: sp_shear_fcn

		if (shear_fcn_opt == 1) then
			sp_shear_fcn = 0.3711d0*(Z*e)**2 * nN**(4.d0/3.d0)/2.d0**(1.d0/3.d0)
		elseif (shear_fcn_opt == 2) then
			T = 1.d8
			RCell = (3.d0/4.d0/pi/nN)**(1.d0/3.d0)
			Ga = (Z*e)**2/RCell/k_b/T
			sp_shear_fcn = 0.1194d0*nN * (Z*e)**2.d0/RCell / (1.d0 + 0.595d0*(173.d0/Ga)**2)
		else
			write(*,*) "err: shear fcn"
			pause
		endif
!write(*,*) 0.595d0*(173.d0/Ga)**2
!pause
	endfunction sp_shear_fcn

	subroutine sp_shear(N1, N2)
	use interpolation
	use Format_IO
	use FWrite
	use FRead
	use verify
	implicit none
	real(8), dimension(1:she_tab_Nmax) :: r_, lgrho_, lgmu_
	integer :: i, N1, N2

		if (she_format == 1) then 
			if (NCol(she_tab) /= 3) then
				write(*,*) "err: shear table column mismatch"
				pause
			endif
			call RFile3R8(she_tab, shear_n, r_, lgrho_, lgmu_)
		elseif (she_format == 2) then
			if (NCol(she_tab) /= 2) then
				write(*,*) "err: shear table column mismatch"
				pause
			endif
			call RFile2R8(she_tab, shear_n, lgrho_, lgmu_)
		endif
			write(*,*) "range of shear_table:"
			write(*,*) "rho_i(/g cm^-3):", 10.d0**lgrho_(1)
			write(*,*) "rho_f(/g cm^-3):", 10.d0**lgrho_(shear_n)
			write(*,*) "range of crust:"
			write(*,*) "rho_i(/g cm^-3):", rho(N1)
			write(*,*) "rho_f(/g cm^-3):", rho(N2) 
			write(*,*) "=========================================================="
		call ver_list_range(10.d0**lgrho_(1), 10.d0**lgrho_(shear_n), rho(N1:N2), N2-N1, "range; shear_table")
		open(81, file = 'data/shear.txt', status ='replace')
		do i = N1 , N2
			mu(i) = Int_ada(2, 1, shear_n, lgrho_, lgmu_, dlog10(rho(i)), "interpolate shear modulus")
			mu(i) = 10.d0**(mu(i))* mu_red
			!call Write4R8(81, FR8, sp_r(i), mu(i), dlog10(rho(i)), dlog10(mu(i)))
			call Write4R8(81, FR8, sp_r(i), mu(i), dlog10(P(i)), dlog10(mu(i)))
		enddo
		close(81)
	endsubroutine sp_shear

	subroutine sp_poly_shear(N1, N2)
	use Format_IO
	use FWrite
	use FRead
	implicit none
	integer :: i, N1, N2
		open(81, file = 'data/shear.txt', status ='replace')
		do i = N1 , N2
				mu(i) = rho(i)* poly_cs**2
			call Write4R8(81, FR8, sp_r(i), mu(i), dlog10(rho(i)), dlog10(mu(i)))
		enddo
		close(81)
	endsubroutine sp_poly_shear

	subroutine sp_poly_shear_Andersson(N1, N2)
	use Format_IO
	use FWrite
	use FRead
	implicit none
	integer :: i, N1, N2
		open(81, file = 'data/shear.txt', status ='replace')
		mu_poly = mu_poly
		do i = N1 , N2
				mu(i) = P(i)* Ki_poly + mu_poly* (1.d-10/(Grav_Const/c**4))
			call Write4R8(81, FR8, sp_r(i), mu(i), dlog10(rho(i)), dlog10(mu(i)))
		enddo
		close(81)
	endsubroutine sp_poly_shear_Andersson

	subroutine sp_qs_shear(N1, N2)
	use Format_IO
	use FWrite
	use FRead
	implicit none
	integer :: i, N1, N2
	real(8) :: sp_qs_P_rho, x_, E_
	real(8) :: C1_, A_, B_, C_, mu_2, mu_	! mu_ and mu_2 are chemical potential; mu(i) in global_var is shear modulus
	real(8), parameter :: e_mks = 1.60217657d-19
		
		!¡¸Quark Star EOS with ref to Pagliaroli 2014; 
		!¡¸Convert to natural units: h_bar = 1, c = 1, P in dimension [E/L^3] = [E^4/((h_bar*c)^3)]
		
		C1_ = 3.d0/4.d0/pi**(2.d0)
		A_ = C1_ * QS_a4
		B_ = - C1_ * QS_a2
		C_ = -B_eff

		open(81, file = 'data/shear.txt', status ='replace')
		do i = N1 , N2
			E_ = (h_bar*c)**3 * rho(i)/((1.d13)*e_mks)**4 * c**2	! Energy density E_
			mu_2 = -B_/6.d0/A_ + dsqrt(B_**2 + 12.d0*A_*(C_+E_))/6.d0/A_

			mu(i) = QS_mu0 *(QS_Gap/10.d0)**2 * mu_2/(400.d0)**2
			mu(i) = mu(i) * ((1.d13)*e_mks) * (1.d39)	! convert unit from MeV/fm^3 to cgs
			call Write4R8(81, FR8, sp_r(i), mu(i), dlog10(rho(i)), dlog10(mu(i)))
			!call Write4R8(81, FR8, sp_r(i), mu(i), dlog10(P(i)), dlog10(mu(i)))
		enddo
		close(81)
	endsubroutine sp_qs_shear

endmodule EOS