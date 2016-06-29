module stat_profile_opt01
use global_var
contains

!C1---------------------------------------------------
subroutine sp01_Ctrl
use EOS
implicit none
	write(100,*) "stat_profile_opt01: the most basic static profile solver "
	write(100,*) "- read eos table, solve hydrodynamic equtions, " 
	write(100,*) "- obtain P, rho, m as functions of r"
	write(100,*) "- able to cope with density jumps (pressure continuous)"
	
	sp_eos_fcall = .true.

	call sp01_sp
	write(100,*) "sp01_sp: passed"
	
	if (comp_opt == 1) then 
		call sp_comp(sp_N1, sp_N2)
		write(100,*) "sp01_comp: passed"
	elseif (comp_opt == 2) then 
		call sp_shear(sp_N1, sp_N2)
		write(100,*) "sp01_shear: passed"
	elseif (comp_opt == 3) then
		write(100,*) "0 shear"
	elseif (comp_opt == 4) then
		call sp_poly_shear(sp_N1, sp_N2)
		write(100,*) "sp01_poly_shear: passed"
	elseif (comp_opt == 5) then
		call sp_qs_shear(sp_N1, sp_N2)
		write(100,*) "sp01_qs_shear: passed"
	elseif (comp_opt == 6) then
		call sp_qs_shear(0, sp_N2)
		write(100,*) "sp02_qs_shear, solid star: passed"
	elseif (comp_opt == 7) then
		call sp_poly_shear_Andersson(0, sp_N2)
		write(100,*) "sp01_poly_shear_Andersson, solid star: passed"
	else
		write(*,*) "err: comp_opt"
		pause
	endif

	call sp01_display
	write(100,*) "=========================================================="

endsubroutine sp01_Ctrl


!C1 A1---------------------------------------------------
subroutine sp01_sp
use EOS
use sp_hyd
use Format_IO
use FWrite
implicit none
real(8) :: P_0, nb_0, x(1:3), y(1:3), dr_scale = delr_adj
integer :: i

	P_0 = sp_eos('p(rho)', rho_0)	! equation of state in "EOS" module
	nb_0 = sp_eos('nb(rho)',rho_0)
	if (P_0 <= P_t) then
		write(*,*) "err: P_0 <= P_t"
		pause
	endif

	open(201, file = 'data/Stat_Profile.txt', status = 'replace')
	call sp01_sp_ini(P_0, nb_0, x)
	i = 0
	call Write5R8(201, FR8, sp_r(i), P(i), rho(i), m(i), delr(i))
	do while (P(i) > P_min)
		sp_r(i+1) = sp_r(i) + delr(i)
		call sp_hyd_int(2, sp_r(i), sp_r(i+1), x(1:2), y(1:2))
		y(3) = sp_eos('rho(p)', y(1))
		
		call sp01_CCInt(x, y, i, dr_scale, sp_N1)

		x = y
		i = i + 1
		P(i) = x(1)
		m(i) = x(2)
		rho(i) = x(3)
		delr(i) = sp01_dr_adj(dr_scale, i)
		call Write5R8(201, FR8, sp_r(i), P(i), rho(i), m(i), delr(i))
	enddo
	sp_N2 = i
	sp_Nt = i

	R_mid = sp_r(sp_N1) 
	M0 = m(sp_Nt)
	R0 = sp_r(sp_Nt)

	close(201)

endsubroutine sp01_sp

!C1 B1---------------------------------------------------
subroutine sp01_sp_ini(P_0, nb_0, x)
implicit none
real(8) :: P_0, nb_0, x(1:3)
	P(0) = P_0
	nb(0) = nb_0
	rho(0) = rho_0
	delr(0) = delr_0
	m(0) = 0.d0
	sp_r(0) = 0.d0

	x(1) = P(0)
	x(2) = m(0)
	x(3) = rho(0)
endsubroutine sp01_sp_ini

!C1 B2---------------------------------------------------
function sp01_dr_adj(dr_, i)
use Derivatives
implicit none
real(8) :: sp01_dr_adj, dr_
integer :: i
	sp01_dr_adj = dabs(dr_/(dydx_2p(delr(i-1),m(i)-m(i-1))/m(i) - dydx_2p(delr(i-1),P(i)-P(i-1))/P(i)))
	if (delr(i-1) == 0.d0) pause 'err: sp01_dr_adj, delr = 0'
endfunction sp01_dr_adj

!C1 B3---------------------------------------------------
subroutine sp01_CCInt(x, y, i, dr_, CC_i)
!	Refine the step size near CC interface
implicit none
real(8), dimension(:) :: x, y
real(8) :: dr_
integer :: i, CC_i
integer, save :: count_ = 1

	if (P_t < x(1) .and. P_t >= y(1) .and. count_ == 1) then
		i = i - 1
		count_ = count_ + 1
		dr_ = delr_ref
		y = x
	elseif (P_t < x(1) .and. P_t >= y(1) .and. count_ == 2) then 
		CC_i = i+1
		dr_ = delr_adj
		count_ = 1
	else
		return
	endif
endsubroutine sp01_CCInt

!C1 A2---------------------------------------------------
subroutine sp01_display
implicit none
	write(*,*) "Initial Conditions:"
	write(*,*) "rho_0(/g cm^-3)", rho(0)
	write(*,*) "P_0(/dym cm^-2)", P(0)
	write(*,*) "=========================================================="
	write(*,*) "crust core transition:"
	write(*,*) "rho_-(/g cm^-3):", rho(sp_N1-1)
	write(*,*) "rho_+(/g cm^-3):", rho(sp_N1)
	write(*,*) "density jump(/rho_+):", (rho(sp_N1-1) - rho(sp_N1))/rho(sp_N1)
	write(*,*) "P_t(dyn cm^-2)", P(sp_N1)
	write(*,*) "mu_t(dyn cm^-2)", mu(sp_N1)
	write(*,*) "shear modulus jump(/P_t):", mu(sp_N1)/P(sp_N1)
	write(*,*) "=========================================================="
	write(*,*) "Summary:"
	write(*,*) "Mmid(/solar mass) = ", m(sp_N1)/M_Solar
	write(*,*) "Rmid(/km) = ", R_mid/1.d5
	write(*,*) "M(/solar mass) = ", M0/M_Solar
	write(*,*) "R(/km) = ", R0/1.d5
	write(*,*) "dMc(/solar mass) = ", (M0 - m(sp_N1))/M_Solar
	write(*,*) "dRc(/km) = ", (R0 - R_mid)/ 1.d5
	write(*,*) "Tuncate Radius(/km)", (R0 - (R0 - R_g)* pg_rs_frac)/1.d5
	write(*,*) "Compactness = ", Grav_const * M0/c**2/R0
	write(*,*) "P_cutoff(/dym cm^-2) = ", P_min
	write(*,*) "rho_cutoff(/g cm^-3) = ", rho(sp_Nt)
	write(*,*) "=========================================================="
	write(*,*) "LD formalism singularity freq(/Hz) = ", dsqrt( l_0*(l_0+1.d0)/2.d0 * Grav_Const *M0/R0**3 )/2.d0/pi
	write(*,*) "=========================================================="
endsubroutine sp01_display

endmodule stat_profile_opt01