module Fitting
contains
	subroutine Poly_Para_Fit(n, x, y, Soln)
	!	Least square soln for overdeterminate system:
	!	[A'] Soln = [B]
	!	least sq soln satisfies [A']T [A'] Soln = [A']T [B]
	!	hence, [A] Soln = Vec
	!
	!	Soln = {a,b,c}
	!	Parabola: y = a + b*x + c*x**2
	!	Large Error?
	use lin_sol_gen_int	!using IMSL; check menu of Compaq Fortran about how to activate/install
	implicit none
	integer :: n
	real(8), dimension(1:n) :: x, y
	real(8) :: A(1:3, 1:3), Vec(1:3,1), Soln(1:3,1)
	integer :: i
		A = 0.d0
		Vec = 0.d0
		do i = 1, n
			A(1,1) = A(1,1) + 1.d0
			A(1,2) = A(1,2) + x(i)
			A(1,3) = A(1,3) + x(i)**2

			A(2,1) = A(2,1) + x(i)
			A(2,2) = A(2,2) + x(i)**2
			A(2,3) = A(2,3) + x(i)**3

			A(3,1) = A(3,1) + x(i)**2
			A(3,2) = A(3,2) + x(i)**3
			A(3,3) = A(3,3) + x(i)**4

			Vec(1,1) = Vec(1,1) + y(i)
			Vec(2,1) = Vec(2,1) + x(i) * y(i)
			Vec(3,1) = Vec(3,1) + x(i)**2 * y(i)
		enddo
		call LIN_SOL_GEN(A, Vec, Soln)

	endsubroutine Poly_Para_Fit
endmodule Fitting