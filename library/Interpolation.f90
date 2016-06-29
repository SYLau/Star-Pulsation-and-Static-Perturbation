MODULE Interpolation
!*******************************************************************
!	Interpolation:
!	#	call Int_ada to use
!	Int_ada		- control unit for interpolation
!	
!	mode		- choose option for interpolation (order/ log scale etc)
!	i1, i2		- data range, must enclose the output range
!	x, y		- array of input data
!	xi			- x value of output data
!	msg			- error msg, displace if xi is out of data range
!
!*******************************************************************
CONTAINS

!************************************************************************************
!	Interpolation control unit
!************************************************************************************
function Int_ada(mode, i1, i2, x, y, xi, msg)
implicit none
!	Instructions: 
!				- input x, y as discrete functions in j-space
!				- input i1, i2 as range of discrete functions x, y in j-space
!				- input xi as arbitrary value in i-space
!				- return Int_ada as y(xi) value in i-space
!				- return error msg if xi is not within x 
!
!	Notice: please make sure x and y are in the same j-space
!	mode:
!		1		- linear ln interpolation
!		2		- linear interpolation
!
!		4		- 4 pts polynomial interpolation
!------------------------------------------------------------------------------------
real(8) :: x(:), y(:), xi
real(8) :: Int_ada, Polint_out
integer :: mode, i1, i2, j, w
character(len=*), optional :: msg
	if (present(msg) == .false.) msg = ''
	do j = i1, i2 - 1
		if ( (xi - x(j)) * (xi - x(j+1)) <= 0.d0) then
			
			if (mode == 1) then
				Int_ada = linear_Ln_int(x(j), x(j+1), y(j), y(j+1), xi)
			elseif (mode == 2) then
				Int_ada = linear_int(x(j), x(j+1), y(j), y(j+1), xi)
			elseif (mode == 4) then
				w = j
				if (w == i1 .or. w == i1 +1) w = i1 + 2
				if (w == i2) w = i2 -1
				call polint(x(w-2:w+1),y(w-2:w+1),4,xi,Polint_out)
				Int_ada = Polint_out
			else
				write(*,*) "err: Int_ada; interpolation mode"
			endif

			exit
		endif
		if (j == i2 - 1) then
			if (present(msg)) then
				write(*,*) "err: ", msg
			else
				write(*,*) "err: interpolation pt out of range"
			endif
			pause
		endif
	enddo
endfunction Int_ada

!************************************************************************************
!	Linear Interpolation
!************************************************************************************
FUNCTION linear_int(x1,x2,y1,y2,x)
REAL(8)	:: x1, x2, y1, y2, x
REAL(8)	:: linear_int
	linear_int = (y2 - y1)/(x2 - x1) * (x - x1) + y1
ENDFUNCTION linear_int

!************************************************************************************
!	Linear Ln Interpolation
!************************************************************************************
FUNCTION linear_Ln_int(x1,x2,y1,y2,x)
REAL(8)	:: x1, x2, y1, y2, x
REAL(8)	:: linear_Ln_int
	if (x1*x2*y1*y2*x == 0.d0) then
		write(*,*) "err: linear Ln interpolation; zero arguements"
		pause
	endif
	linear_Ln_int = dexp((dlog(y2) - dlog(y1))/(dlog(x2) - dlog(x1)) * (dlog(x) - dlog(x1)) + dlog(y1))
ENDFUNCTION linear_Ln_int

!************************************************************************************
!	Polynomial Interpolation (Numerical Recipes)
!************************************************************************************
SUBROUTINE polint(xa,ya,n,x,y)
INTEGER n,NMAX
REAL(8) dy,x,y,xa(n),ya(n)
PARAMETER (NMAX=10) !	Largest anticipated value of n.
!Given arrays xa and ya, each of length n, and given a value x, this routine returns a
!value y, and an error estimate dy. If P(x) is the polynomial of degree N ? 1 such that
!P(xai ) = yai, i = 1,..., n, then the returned value y = P(x).
INTEGER i,m,ns
REAL(8) den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)

	ns=1
	dif=abs(x-xa(1))
	do i=1,n 
		dift=abs(x-xa(i))
	if (dift.lt.dif) then
		ns=i
		dif=dift
	endif
	c(i)=ya(i)
	d(i)=ya(i)
	enddo
	y=ya(ns) 
	ns=ns-1
	do m=1,n-1 
	do i=1,n-m 
		ho=xa(i)-x
		hp=xa(i+m)-x
		w=c(i+1)-d(i)
		den=ho-hp
	if(den.eq.0.)pause
		den=w/den
		d(i)=hp*den
		c(i)=ho*den
	enddo
	if (2*ns.lt.n-m)then
		dy=c(ns+1)
	else
		dy=d(ns)
		ns=ns-1
	endif
	y=y+dy
	enddo
	return
ENDSUBROUTINE polint

!!	Reference
SUBROUTINE POLINT1(XA,YA,N,X,Y,DY)
!*****************************************************
!*     Polynomial Interpolation or Extrapolation     *
!*            of a Discreet Function                 *
!* ------------------------------------------------- *
!* INPUTS:                                           *
!*   XA:    Table of abcissas  (N)                   *
!*   YA:    Table of ordinates (N)                   *
!*    N:    Number of points                         *
!*    X:    Interpolating abscissa value             *
!* OUTPUTS:                                          *
!*    Y:    Returned estimation of function for X    *
!*   DY:    Estimated error for Y                    *
!*****************************************************
PARAMETER(NMAX=25)
REAL*8 XA(N),YA(N), X,Y,DY
REAL*8 C(NMAX),D(NMAX)
REAL*8 DEN,DIF,DIFT,HO,HP,W
NS=1
DIF=DABS(X-XA(1))
DO I=1,N
  DIFT=DABS(X-XA(I))
  IF(DIFT.LT.DIF) THEN
    NS=I                 !index of closest table entry
	DIF=DIFT
  ENDIF
  C(I)=YA(I)             !Initialize the C's and D's
  D(I)=YA(I)
END DO
Y=YA(NS)                 !Initial approximation of Y
NS=NS-1
DO M=1,N-1
  DO I=1,N-M
    HO=XA(I)-X
	HP=XA(I+M)-X
    W=C(I+1)-D(I) 
    DEN=HO-HP
	IF(DEN.EQ.0.) Pause 'Error: two identical abcissas)'
	DEN=W/DEN
	D(I)=HP*DEN          !Update the C's and D's
	C(I)=HO*DEN
  END DO
  IF(2*NS.LT.N-M) THEN   !After each column in the tableau XA
    DY=C(NS+1)           !is completed, we decide which correction,
  ELSE                   !C or D, we want to add to our accumulating 
    DY=D(NS)             !value of Y, i.e. which path to take through 
	NS=NS-1              !the tableau, forking up or down. We do this
  ENDIF                  !in such a way as to take the most "straight 
  Y=Y+DY	             !line" route through the tableau to its apex,
END DO                   !updating NS accordingly to keep track of 
                         !where we are. This route keeps the partial
RETURN                   !approximations centered (insofar as possible)
ENDSUBROUTINE POLINT1    !on the target X.The last DY added is thus the
                         !error indication.
ENDMODULE Interpolation