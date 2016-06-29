MODULE Root_Finding

CONTAINS

!************************************************************************************
!	Secant Method
!************************************************************************************
SUBROUTINE SECANT (fx, DL,X_In,DX_In,OUT)
! Subroutine for the root of f(x)=0 with the secant method.
! Copyright (c) Tao Pang 1997.
! 
IMPLICIT NONE
!------------------------------------------------------------------
!	Root finding with Secant Method
!	input:
!	fx		- input fcn
!	DL		- min bracket size (fractional)
!	X_In	- first trial point
!	DX_In	- first bracket size
!	output:
!	OUT		- output value
!------------------------------------------------------------------
	INTEGER :: ISTEP
	REAL(8), INTENT (IN) :: X_In,DX_In
	REAL(8) :: X0,X1,X2,D,DX, fx
	REAL(8), INTENT (IN) :: DL
	REAL(8), INTENT (OUT) :: OUT
	EXTERNAL fx
	ISTEP = 0
	X0 = X_In
	DX = DX_In
	X1 = X0+DX
!	write(*,*) X_In, DX_In
!	read(*,*)
	DO WHILE (ABS(DX)/X0 .GT. DL)
		D  = FX(X1)-FX(X0)
		X2 = X1-FX(X1)*(X1-X0)/D
		X0 = X1
		X1 = X2
		DX = X1-X0
    ISTEP = ISTEP+1
  END DO
	OUT = X2
ENDSUBROUTINE SECANT

!************************************************************************************
!	Bisection Method Numerical Recipes
!************************************************************************************
FUNCTION rtbis(func,x1,x2,xacc)
INTEGER	:: JMAX
REAL(8)	:: rtbis,x1,x2,xacc,func
EXTERNAL func
PARAMETER (JMAX=1000) !Maximum allowed number of bisections.
!Using bisection, find the root of a function func known to lie between x1 and x2. The
!root, returned as rtbis, will be refined until its accuracy is +-xacc.
INTEGER j
REAL(8)	:: dx,f,fmid,xmid

fmid=func(x2)
f=func(x1)
if(f*fmid.ge.0.D0) pause !Root must be bracketed in rtbis
if(f.lt.0.D0)then !Orient the search so that f>0 lies at x+dx.
rtbis=x1
dx=x2-x1
else
rtbis=x2
dx=x1-x2
endif
do j=1,JMAX !11 Bisection loop.
dx=dx*.5D0
xmid=rtbis+dx
fmid=func(xmid)
if(fmid.le.0.D0)rtbis=xmid
if(abs(dx).lt.xacc .or. fmid.eq.0.) return
enddo !11
pause !too many bisections in rtbis
ENDFUNCTION rtbis


!************************************************************************************
!	Muller's Method 
!************************************************************************************
function rtmuller_c(fcn, x1_, x2_, x3_, xacc, JMAX, ERR)
implicit none
complex(8) :: fcn, x1_, x2_, x3_, x1, x2, x3, f1, f2, f3, xtemp, dx1, dx2
complex(8) :: c, t, u, v, a, b, discrim, denom
real(8) :: xacc, x2_bar
complex(8) :: rtmuller_c
integer :: n
integer, optional :: JMAX, ERR
external fcn
	
	if (present(JMAX) /= .true.) JMAX=500
	if (present(ERR) == .true.) ERR = 0
	n = 0
	x1 = x1_
	x2 = x2_
	x3 = x3_

	dx1 = x2 - x1
	dx2 = x3 - x2
	x2_bar = (cdabs(x3) + cdabs(x2))/2.d0
	if (dx1 == (0.d0, 0.d0) .or. dx2 == (0.d0, 0.d0)) then
		write(*,*) "err: root-finding: rtmuller_c"
		pause
	endif
	
	do while (cdabs(dx2) > xacc* dabs(x2_bar))
		n = n + 1

		do while ((x3-x2)*(x2-x1)*(x1-x3) == 0.d0)
			if (x3-x2 == 0.d0) x3 = x3*1.00001d0
			if (x2-x1 == 0.d0) x2 = x2*1.00002d0
			if (x1-x3 == 0.d0) x1 = x1*1.00003d0
		enddo

		f1 = fcn(x1)
		f2 = fcn(x2)
		f3 = fcn(x3)

		if (f1 == 0) then
			xtemp = x1
			exit
		elseif (f2 == 0) then
			xtemp = x2
			exit
		elseif (f3 == 0) then
			xtemp = x3
			exit
		endif
				
		c = f3
		t = (f3-f2)/(x3-x2)
		u = (f2-f1)/(x2-x1)
		v = (f1-f3)/(x1-x3)
		a = (v-u)/(x3-x2)
		b = t + v - u

		discrim = sqrt(b*b - 4.d0 * a * c)
	
		IF ( abs(b+discrim) > abs(b-discrim)) THEN
			denom = b + discrim
		ELSE
			denom = b - discrim
		END IF

		xtemp = x3 - (2.d0*c) / denom

		dx1 = x2 - x1
		dx2 = x3 - x2
		x2_bar = (cdabs(x3) + cdabs(x2))/2.d0
		x1 = x2
		x2 = x3
		x3 = xtemp
		
		if (n >= JMAX) then
			write(*,*) "err: root-finding: too many steps"
			pause
			if (present(ERR)== .true.) ERR = 1
			return
		endif
	enddo
	rtmuller_c = xtemp
endfunction rtmuller_c

function vec_dot(v1, v2)
implicit none
real(8) :: vec_dot, v1(1:), v2(1:), s
integer :: i
	if (size(v1) /= size(v2)) then
		pause 'err: vector dot product, dimension not match'
	else
		s = 0.d0
		do i = 1, size(v1)
			s = s + v1(i) * v2(i)
		enddo
		vec_dot = s
	endif
endfunction vec_dot

ENDMODULE Root_Finding