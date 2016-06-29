module global_fcn
! contains temporary functions for debugging or functions for very specific uses
contains

subroutine NS_Match_r3(z, d)
!	provides a scanning scheme that contains singular point, but very useful to classify interface/ shear modes from normal fluid modes
implicit none
real(8) :: d3, a2, a3, d
real(8) :: z(1:3,1:3)
	d3 = z(1,2)*z(2,3)-z(2,2)*z(1,3)
	a2 = (z(1,1)*z(2,3)-z(2,1)*z(1,3))/d3
	a3 = (z(1,2)*z(2,1)-z(2,2)*z(1,1))/d3
	d = z(3,1) - a2*z(3,2) - a3*z(3,3)
endsubroutine NS_Match_r3

!	Debug

subroutine RC_Co_FcnBug
use global_var
use puls_eqt_set_opt02
use Format_IO
use FWrite
use Hydrostatic
implicit none
real(8) :: eLamb, dLamb, dNu,V2,gamma,U2
integer :: i
	open(01,file='data/debug.txt', status = 'replace')
	do i = 0, pg_N1*2
		eLamb = pes02_eLamb(rho_Co(i), m_Co(i), r_Co(i))
		dLamb = pes02_dLamb(rho_Co(i), m_Co(i), r_Co(i))
		dNu = pes02_dnu(P_Co(i), rho_Co(i), m_Co(i), r_Co(i))
		V2 = pes02_V2(P_Co(i), rho_Co(i), m_Co(i), r_Co(i))
		gamma = pes02_gamma(1, P_Co, rho_Co, i, 0, pg_N1*2)
		U2 = pes02_U2(rho_Co(i), m_Co(i), r_Co(i))
		call Write6R8(01,FR8,eLamb,r_Co(i)*dLamb,r_Co(i)*dNu,V2,gamma,U2)
	enddo
	close(01)
endsubroutine RC_Co_FcnBug

endmodule global_fcn
