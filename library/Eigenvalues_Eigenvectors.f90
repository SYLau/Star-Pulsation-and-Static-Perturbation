MODULE Eigenvalues_Eigenvectors

CONTAINS

!************************************************************************************
!	Determinant of Dimension N
!************************************************************************************
SUBROUTINE Determinant_N(a, output, n)
!------------------------------------------------------------------
!	ERROR MAY BE HUGE FOR LARGE SYSTEMS 
!	USE IMSL INSTEAD!
!
!	Convert the matrix to Upper Triangular Form (UT); Then find determinant 
!	by multiplying all diagonal elements	
!
!	Input/Output:
!	a()		- input matrix
!	output	- determinant of a()
!	n		- dimension of sq matrix
!		
!	Remarks	: Uses subroutine "Tri_Matrix" to make the matrix upper tri form
!			  Should include subroutines "Tri_Matrix" & "Swap"
!------------------------------------------------------------------
REAL(8) :: a(1:n,1:n), b(1:n,1:n)
REAL(8) :: output
INTEGER :: i, n, l

	CALL Tri_Matrix(a, b, n, l)
	output = l						!	l represents the sign change due to swapping 2 rows; 
									!	It also gives 0 if the matrix cannot be UT form
	DO i = 1, n
		output = output * b(i,i)
	ENDDO

ENDSUBROUTINE Determinant_N


SUBROUTINE Determinant_N_c(a, output, n)
complex(8) :: a(1:n,1:n), b(1:n,1:n)
complex(8) :: output
INTEGER :: i, n, l

	CALL Tri_Matrix_c(a, b, n, l)
	output = l						!	l represents the sign change due to swapping 2 rows; 
									!	It also gives 0 if the matrix cannot be UT form
	DO i = 1, n
		output = output * b(i,i)
	ENDDO

ENDSUBROUTINE Determinant_N_c

!************************************************************************************
!	Converting a matrix to upper triangular form
!************************************************************************************
SUBROUTINE Tri_Matrix(a, b, n, l)
REAL(8) :: a(1:n, 1:n), b(1:n, 1:n)
REAL(8) :: Factor
INTEGER :: i, j, k, l
LOGICAL :: Det_Exist

Det_Exist = .TRUE.
!!PartI : Initialization
	l = 1
	DO i = 1,n
		DO j = 1,n
			b(i,j) = a(i,j)
		ENDDO
	ENDDO
!!PartII :
	DO i = 1,n-1
		!!PartIIA : Check existence of tri_matrix
		IF (b(i,i) == 0.d0) THEN
		Det_Exist = .FALSE.
		DO j = i+1, n
			IF (b(j,i) /= 0.d0) THEN
			Det_Exist = .TRUE.
			DO k = 1, n
				CALL SWAP(b(i,k), b(j,k))
			ENDDO
			l = (-1)*l
			EXIT
			ENDIF
		ENDDO 
		ENDIF
		IF (Det_Exist == .FALSE.) THEN
			write(*,*) "The matrix cannot be changed to UT form"
			read(*,*)
			l = 0
			RETURN
		ENDIF

		!!PartIIB : Make triangular matrix 
		DO j = i+1, n
			Factor = b(j,i)/b(i,i)
			DO k = i, n
				b(j,k) = b(j,k) - Factor * b(i,k)
			ENDDO
		ENDDO
	ENDDO

ENDSUBROUTINE Tri_Matrix

SUBROUTINE Tri_Matrix_c(a, b, n, l)
complex(8) :: a(1:n, 1:n), b(1:n, 1:n)
complex(8) :: Factor
INTEGER :: i, j, k, l
LOGICAL :: Det_Exist

Det_Exist = .TRUE.
!!PartI : Initialization
	l = 1
	DO i = 1,n
		DO j = 1,n
			b(i,j) = a(i,j)
		ENDDO
	ENDDO
!!PartII :
	DO i = 1,n-1
		!!PartIIA : Check existence of tri_matrix
		IF (b(i,i) == 0.d0) THEN
		Det_Exist = .FALSE.
		DO j = i+1, n
			IF (b(j,i) /= 0.d0) THEN
			Det_Exist = .TRUE.
			DO k = 1, n
				CALL SWAP_c(b(i,k), b(j,k))
			ENDDO
			l = (-1)*l
			EXIT
			ENDIF
		ENDDO 
		ENDIF
		IF (Det_Exist == .FALSE.) THEN
			write(*,*) "The matrix cannot be changed to UT form"
			read(*,*)
			l = 0
			RETURN
		ENDIF

		!!PartIIB : Make triangular matrix 
		DO j = i+1, n
			Factor = b(j,i)/b(i,i)
			DO k = i, n
				b(j,k) = b(j,k) - Factor * b(i,k)
			ENDDO
		ENDDO
	ENDDO

ENDSUBROUTINE Tri_Matrix_c

!************************************************************************************
!	Swapping 2 real no
!************************************************************************************
SUBROUTINE Swap(a,b)
REAL(8) :: a, b, c
		c = a
		a = b
		b = c
ENDSUBROUTINE SWAP

SUBROUTINE Swap_c(a,b)
complex(8) :: a, b, c
		c = a
		a = b
		b = c
ENDSUBROUTINE Swap_c

!************************************************************************************
!	Finding Coefficients
!************************************************************************************
SUBROUTINE Eigenvector(a, b, n)
!------------------------------------------------------------------
!	Given a set of linear combination of solution, with 1 linearly dependence relation
!	Returns the values of coefficients up to 1 normalization factor
!
!	Input/Output:
!	a()		- input matrix
!	b()		- output coefficients of the vectors; elements of eigenvector
!	n		- dimension of sq matrix
!		
!	Remarks	: Uses subroutine "Determinant_N" to output value of determinant
!			  Should include subroutines "Determinant_N", "Tri_Matrix" & "Swap"
!			  Can be used to find eigenvectors
!------------------------------------------------------------------
REAL(8)	:: a(1:n, 1:n), b(1:n), D(1:n-1,1:n-1), M(1:n-1,1:n-1)
REAL(8)	:: Determinant_M,Determinant_D
INTEGER	:: n, i,j,k

!!	Set b(1) to be normalized
b(1) = 1.D0

DO i = 1,n-1
	DO j = 1,n-1
		D(i,j) = a(i+1,j+1)
	ENDDO
ENDDO
DO k = 1, n-1
	DO i = 1,n-1
		DO j = 1,n-1
			M(i,j) = D(i,j)
		ENDDO
	ENDDO
	DO i = 1,n-1
		M(i,k) = -a(i+1,1)
	ENDDO

	CALL Determinant_N(M, Determinant_M, n-1)
	CALL Determinant_N(D, Determinant_D, n-1)
	b(k+1) = Determinant_M/Determinant_D

ENDDO

ENDSUBROUTINE Eigenvector

!************************************************************************************
!	Solving Linear System: Cramer's Rule
!************************************************************************************
!	Works quite well for a 4x4 system
!	Gives huge error for a 6s6 system (5~10%)
!
subroutine LinSys_Cramer(a, b, vec, n, err_p)
!	Solving system "[a] vec = b"
implicit none
integer :: n
real(8), dimension(1:n, 1:n) :: a
real(8), dimension(1:n) :: b, vec, row
real(8), dimension(1:n, 1:n) :: M
real(8) :: det, det_M
integer :: i, j
real(8), dimension(1:n), optional :: err_p
	CALL Determinant_N(a, det, n)
	if (det == 0.d0) then
		write(*,*) "err: LinSys_Cramer; det = 0"
		pause
	endif
	do i = 1,n
		M = a
		M(1:n, i) = b(1:n)
		CALL Determinant_N(M, det_M, n)
		vec(i) = det_M/det
	enddo
	if (present(err_p)) then
		err_p = 0.d0
		do i = 1,n
			row = 0.d0
			do j = 1,n
				row(i) = row(i) + M(i,j) * vec(j)
			enddo
			err_p(i) = dabs((b(i) - row(i)))
		enddo
	endif
endsubroutine LinSys_Cramer

subroutine LinSys_Cramer_c(a, b, vec, n, err_p)
implicit none
integer :: n
complex(8), dimension(1:n, 1:n) :: a
complex(8), dimension(1:n) :: b, vec, row
complex(8), dimension(1:n, 1:n) :: M
complex(8) :: det, det_M
integer :: i, j
real(8), dimension(1:n), optional :: err_p
	CALL Determinant_N_c(a, det, n)
	if (det == 0.d0) then
		write(*,*) "err: LinSys_Cramer_c; det = 0"
		pause
	endif
	do i = 1,n
		M = a
		M(1:n, i) = b(1:n)
		CALL Determinant_N_c(M, det_M, n)
		vec(i) = det_M/det
	enddo
	if (present(err_p)) then
		err_p = 0.d0
		do i = 1,n
			row = 0.d0
			do j = 1,n
				row(i) = row(i) + M(i,j) * vec(j)
			enddo
			err_p(i) = cdabs((b(i) - row(i)))
		enddo
	endif
endsubroutine LinSys_Cramer_c

ENDMODULE Eigenvalues_Eigenvectors