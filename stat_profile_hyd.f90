module sp_hyd
implicit none
real(8), parameter :: acc = 1.d-6
contains

subroutine sp_hyd_int(n, ri_, rf_, x, y )
use RK4_Set
implicit none
integer :: n
real(8) :: ri_, rf_, x(1:n), y(1:n), x_temp(1:n)
integer :: nok, nbad
real(8) :: dr_
	dr_ = rf_-ri_
	call RK4(sp_hyd_eqts, n, ri_, rf_, x(1:n), y(1:n))
!	y(1:n) = x(1:n)
!	call odeint(y(1:n),n,ri_,rf_,acc, dr_, 0.d0 ,nok,nbad,sp_hyd_eqts, rkqs)
endsubroutine sp_hyd_int

subroutine sp_hyd_eqts(n, t, xi, fcn)
use EOS
use Hydrostatic
implicit none
integer :: n
real(8), dimension(n) :: xi, fcn
real(8) :: t, x4

	x4 = sp_eos('rho(p)',xi(1))
	if (t == 0.d0) fcn(1) = 0.d0
!if (t /= 0.d0) fcn(1) = Newton_dPdr(x4, xi(2), t)
	if (t /= 0.d0) fcn(1) = TOV_dPdr(xi(1), x4, xi(2), t)

	fcn(2) = 4.d0*pi*t**2 * x4

	if (n == 3) then
		if (t == 0.d0) fcn(3) = 0.d0
		if (t /= 0.d0) fcn(3) = -TOV_dPdr(xi(1), x4, xi(2), t)/(x4*c**2 + xi(1))  !!! Metric here defined as e^(2nu)
	endif

endsubroutine sp_hyd_eqts

endmodule sp_hyd