MODULE Hydrostatic
!---------------------------------------------------------------------------
!	Contains functions for calculating dP/dr for both Newtonian & Relatistic (TOV eqt) Cases
!	Can be applied directly to any programs
!	Input:
!		- P, rho, m, r [in cgs unit]
!	Output:
!		- dPdr
!---------------------------------------------------------------------------
CONTAINS

SUBROUTINE Units(G, c)
!	Change the units here
REAL(8)	:: G, c
REAL(8), PARAMETER :: Grav_Const = 6.67428D-11, c_mks = 2.99792458D8
	G = Grav_Const * 1.D3
	c = c_mks * 1.D2
ENDSUBROUTINE Units
!(1)	Newtonian Hydrostatic Eqt (cgs unit)
FUNCTION Newton_dPdr(rho, m, r)
REAL(8)	:: rho, m, r, G, c, Newton_dPdr
REAL(8), PARAMETER :: pi = 3.14159265358979323846D0
	CALL Units(G, c)
	Newton_dPdr = - G * rho * m/ r**2
ENDFUNCTION Newton_dPdr

!(2)	Relativistic Hydrostatic Eqt (cgs unit)
!(A)	Potential
FUNCTION TOV_g(P, rho, m, r)
REAL(8)	:: P, rho, m, r, G, c, TOV_g
REAL(8), PARAMETER :: pi = 3.14159265358979323846D0
	CALL Units(G, c)
	TOV_g = G * (m + 4.D0*pi*r**3*P/c**2)/ r**2/ (1.D0 - 2.D0*G*m/c**2/r)
ENDFUNCTION TOV_g
!(B)	dPdr
FUNCTION TOV_dPdr(P, rho, m, r)
REAL(8)	:: P, rho, m, r, G, c, TOV_dPdr
REAL(8), PARAMETER :: pi = 3.14159265358979323846D0
	CALL Units(G, c)
	TOV_dPdr = - (rho + P/c**2)* TOV_g(P, rho, m, r)
ENDFUNCTION TOV_dPdr


ENDMODULE Hydrostatic
