module stat_profile_opt02
!	static profile for relativistic cowling approximation, compute nu(r)
use global_var
contains

!C1---------------------------------------------------
subroutine sp02_Ctrl
use EOS
implicit none
	write(100,*) "stat_profile_opt02: static profile solver "
	write(100,*) "- read eos table, solve hydrodynamic equtions, " 
	write(100,*) "- obtain P, rho, m as functions of r"
	write(100,*) "- able to cope with density jumps (pressure continuous)"
	write(100,*) "- sp_N1 being the position of the jump with lower density"
	write(100,*) "- solve also for the metric tensor function nu(r)"
	
	sp_eos_fcall = .true.

	call sp02_sp
	write(100,*) "sp02_sp: passed"
	
	if (comp_opt == 1) then 
		call sp_comp(sp_N1, sp_N2)
		write(100,*) "sp02_comp: passed"
	elseif (comp_opt == 2) then 
		call sp_shear(sp_N1, sp_N2)
		write(100,*) "sp02_shear: passed"
	elseif (comp_opt == 3) then
		write(100,*) "0 shear"
	elseif (comp_opt == 4) then
		call sp_poly_shear(sp_N1, sp_N2)
		write(100,*) "sp02_poly_shear: passed"
	elseif (comp_opt == 5) then
		call sp_qs_shear(sp_N1, sp_N2)
		write(100,*) "sp02_qs_shear: passed"
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

	call sp02_display
	write(100,*) "=========================================================="

endsubroutine sp02_Ctrl


!C1 A1---------------------------------------------------
subroutine sp02_sp
use EOS
use sp_hyd
use Format_IO
use FWrite
implicit none
real(8) :: P_0, nb_0, x(1:4), y(1:4), dr_scale = delr_adj
integer :: i
!integer :: position
!logical :: pass

!pass = .false.

	P_0 = sp_eos('p(rho)',rho_0)
	nb_0 = sp_eos('nb(rho)',rho_0)
	if (P_0 <= P_t) then
		write(*,*) "err: P_0 <= P_t"
		pause
	endif

	call sp02_sp_ini(P_0, nb_0, x)
	i = 0

	do while (P(i) > P_min)

		sp_r(i+1) = sp_r(i) + delr(i)
		call sp_hyd_int(3, sp_r(i), sp_r(i+1), x(1:3), y(1:3))
!if (sp_r(i) >= 7.080201094514828d5) then
!write(*,*) y(1), x(1), sp_r(i), sp_r(i+1), sp_r(i+1) - sp_r(i), delr(i), sp_r(i) + delr(i) - sp_r(i)
!pause
!endif
! Accuracy not enough: sp_r(i+1) - sp_r(i) gives 0, leads to infinity in dr[scale] -> -infinity in RK4
! Set min pressure higher
		y(4) = sp_eos('rho(p)',y(1))

		call sp02_CCInt(x, y, i, dr_scale, sp_N1)

		x = y
		i = i + 1
		P(i) = x(1)
		m(i) = x(2)
		nu(i) = x(3)
		rho(i) = x(4)
		delr(i) = sp02_dr_adj(dr_scale, i)

!if (dabs(sp_r(i) -  9.d5)  <= delr(i) .and. pass == .false.) then
!delr(i) = dabs(sp_r(i) -  9.d5)
!position = i+1
!pass = .true.
!endif
	enddo

	sp_N2 = i
	sp_Nt = i

	R_mid = sp_r(sp_N1) 
	M0 = m(sp_Nt)
	R0 = sp_r(sp_Nt)

	call sp02_nu_BC

	open(201, file = 'data/Stat_Profile.txt', status = 'replace')
open(02, file = 'data/temp_compare.txt', status = 'replace')		! write relativistic terms
	do i = 0, sp_Nt
		call Write6R8(201, FR8, sp_r(i), P(i), rho(i), m(i), nu(i), delr(i))
!if (sp_r(i) /= 0) call Write5R8(02, FR8, sp_r(i), (1.d0 - 2.d0 * Grav_Const/c**2 * m(i)/sp_r(i))**(-1), dexp(-2.d0*nu(i)), 1.d0 + 4.d0*pi*P(i)*sp_r(i)**3/m(i)/c**2, (1.d0 - 2.d0 * Grav_Const/c**2 * m(i)/sp_r(i))**(-1)/dexp(-2.d0*nu(i)) *(1.d0 + 4.d0*pi*P(i)*sp_r(i)**3/m(i)/c**2)) ! terms of c1 Omega in Yoshida
!if (sp_r(i) /= 0) call Write3R8(02, FR8, sp_r(i), rho(i)*Grav_Const*m(i)/P(i)/sp_r(i), rho(i)*Grav_Const*m(i)/P(i)/sp_r(i) * (1.d0 - 2.d0 * Grav_Const/c**2 * m(i)/sp_r(i))**(-1)*(1.d0 + 4.d0*pi*P(i)*sp_r(i)**3/m(i)/c**2)*(1.d0 + P(i)/rho(i)/c**2) )
	enddo

	close(201)
close(02)	! terms of c1 Omega in Yoshida

endsubroutine sp02_sp

!C1 B1---------------------------------------------------
subroutine sp02_sp_ini(P_0, nb_0, x)
implicit none
real(8) :: P_0, nb_0, x(1:4)
	P(0) = P_0
	nb(0) = nb_0
	rho(0) = rho_0
	delr(0) = delr_0
	m(0) = 0.d0
	nu(0) = 0.d0
	sp_r(0) = 0.d0

	x(1) = P(0)
	x(2) = m(0)
	x(3) = nu(0)
	x(4) = rho(0)
endsubroutine sp02_sp_ini

!C1 B2---------------------------------------------------
function sp02_dr_adj(dr_, i)
use Derivatives
implicit none
real(8) :: sp02_dr_adj, dr_
integer :: i
	sp02_dr_adj = dabs(dr_/(dydx_2p(delr(i-1),m(i)-m(i-1))/m(i) - dydx_2p(delr(i-1),P(i)-P(i-1))/P(i)))
	if (delr(i-1) == 0.d0) pause 'err: sp02_dr_adj, delr = 0'
endfunction sp02_dr_adj

!C1 B3---------------------------------------------------
subroutine sp02_CCInt(x, y, i, dr_, CC_i)
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
endsubroutine sp02_CCInt

!C1 B4---------------------------------------------------
subroutine sp02_nu_BC
implicit none
real(8) :: Const, nuR

	nuR = 0.5d0 *dlog(1.d0 - 2.d0*Grav_Const*M0/R0/c**2)
	Const = nuR - nu(sp_Nt)
	write(*,*) "Arbitrary constant for nu(r):"
	write(*,*) Const
	write(*,*) "=========================================================="

	nu = nu + Const
endsubroutine sp02_nu_BC


!C1 A2---------------------------------------------------
subroutine sp02_display
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
endsubroutine sp02_display

endmodule stat_profile_opt02