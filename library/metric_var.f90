module metric_var
use global_var
! use global value of pi
contains
!	unit system: geometrized unit
!	metric definition: ds^2 = -e^(2nu) dt^2 + e^(2lambda) dr^2 + ...
	function met_eLamb(rho_, m_, r_)
	implicit none
	real(8) :: met_eLamb, rho_, m_, r_
		met_eLamb = (1.d0 - 2.d0*m_/r_)**(-0.5d0)
	endfunction met_eLamb

	function met_dLamb(rho_, m_, r_)
	implicit none
	real(8) :: met_dLamb, rho_, m_, r_, c1
		c1 = met_eLamb(rho_, m_, r_) **2
		met_dLamb = c1 / r_ * (4.d0*pi*r_**2*rho_ - m_/r_)
	endfunction met_dLamb

	function met_dnu(P_, rho_, m_, r_)
	implicit none
	real(8) :: met_dnu, P_, rho_, m_, r_, c1
		c1 = met_eLamb(rho_, m_, r_) **2
		met_dnu = c1 / r_ * ( 4.d0*pi*r_**(2)*P_ + m_/r_)
	endfunction met_dnu

	function met_ddnu(P_, rho_, m_, r_)
	implicit none
	real(8) :: met_ddnu, P_, rho_, m_, r_, c1
	real(8) :: dLamb, dP, dnu
	real(8) :: B(1:3)

		c1 = met_eLamb(rho_, m_, r_) **2
		dLamb = met_dLamb(rho_, m_, r_)
		dnu = met_dnu(P_, rho_, m_, r_)
		dP = -dnu *(rho_ + P_)
		
		B(1) = 2.d0*dLamb*(4.d0*pi*r_*P_ + m_/r_**2)
		B(2) = 4.d0*pi*(P_+r_*dP)
		B(3) = 4.d0*pi*rho_ - 2.d0*m_/r_**3
		met_ddnu = c1 * (B(1) +B(2) +B(3))

	endfunction met_ddnu

endmodule metric_var