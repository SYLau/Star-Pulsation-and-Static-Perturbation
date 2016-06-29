MODULE RK4_Set
CONTAINS
!************************************************************************************
!	Runge Kutta 4th order method
!************************************************************************************

SUBROUTINE RK4(IFcn, n, ti, tf, xi, y)
IMPLICIT NONE
!------------------------------------------------------------------
!	Solving system of ODE by RK4 method
!	input:
!	IFcn	- subroutine of the input fcn
!	n		- number of ODEs
!	ti		- initial time
!	tf		- solution time
!	xi()	- initial values
!	output:
!	y()		- output value
!------------------------------------------------------------------
INTEGER :: i, j, n
DOUBLE PRECISION, DIMENSION(1:n) :: K1, K2, K3, K4, xi, x, y, fcn
REAL(8) :: h, tf, ti, t
EXTERNAL IFcn

	h = tf - ti
	t = ti

	CALL IFcn(n, t, xi, fcn)
	DO j = 1, n
		K1(j) = fcn(j) *(h)
		x(j) = xi(j) + K1(j)/2.D0
	ENDDO

	CALL IFcn(n, t + h/2.D0, x, fcn)
	DO j = 1, n
		K2(j) = fcn(j) *(h)
		x(j) = xi(j) + K2(j)/2.D0
	ENDDO

	CALL IFcn(n, t + h/2.D0, x, fcn)
	DO j = 1, n
		K3(j) = fcn(j) *(h)
		x(j) = xi(j) + K3(j)
	ENDDO

	CALL IFcn(n, t + h, x, fcn)
	DO j = 1, n
		K4(j) = fcn(j) *(h)
		y(j) = xi(j) + (K1(j) + 2.D0*K2(j) + 2.D0*K3(j) + K4(j))/6.D0
	ENDDO

ENDSUBROUTINE RK4

!************************************************************************************
!	Runge Kutta from Numerical Recipe
!************************************************************************************

subroutine oderk(ri,re,y,n,derivs) 
INTEGER, PARAMETER :: NMAX=16
REAL(8) :: ri, re, step
REAL(8), DIMENSION(NMAX) :: y, dydx, yout
EXTERNAL derivs
           
call derivs(ri,y,dydx) 
step=re-ri 
CALL rk4_NR(y,dydx,n,ri,step,yout,derivs) 
do i=1,n 
   y(i)=yout(i) 
enddo
return 
end subroutine oderk

SUBROUTINE RK4_NR(Y,DYDX,N,X,H,YOUT,DERIVS) 
INTEGER, PARAMETER :: NMAX=16 
REAL(8) :: H,HH,XH,X,H6 
REAL(8), DIMENSION(N) :: Y, DYDX, YOUT 
REAL(8), DIMENSION(NMAX) :: YT, DYT, DYM 
EXTERNAL derivs


HH=H*0.5D0 
H6=H/6D0 
XH=X+HH
 
DO I=1,N 
   YT(I)=Y(I)+HH*DYDX(I) 
ENDDO

CALL DERIVS(XH,YT,DYT) 

DO I=1,N 
   YT(I)=Y(I)+HH*DYT(I) 
ENDDO

CALL DERIVS(XH,YT,DYM) 

DO I=1,N 
   YT(I)=Y(I)+H*DYM(I) 
   DYM(I)=DYT(I)+DYM(I) 
ENDDO

CALL DERIVS(X+H,YT,DYT) 

DO I=1,N 
   YOUT(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2*DYM(I)) 
ENDDO
            

END SUBROUTINE RK4_NR

!************************************************************************************
!	Runge Kutta Adaptive Stepsize Control from Numerical Recipe
!************************************************************************************
!
! Quality Control
!
SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
INTEGER n,NMAX
REAL(8) eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
EXTERNAL derivs
PARAMETER (NMAX=50) !   Maximum number of equations.
!   C USES derivs,rkck
!   Fifth-order Runge-Kutta step with monitoring of local truncation error to ensure accuracy
!   and adjust stepsize. Input are the dependent variable vector y(1:n) and its derivative
!   dydx(1:n) at the starting value of the independent variable x. Also input are the stepsize
!   to be attempted htry, the required accuracy eps, and the vector yscal(1:n) against
!   which the error is scaled. On output, y and x are replaced by their new values, hdid is the
!   stepsize that was actually accomplished, and hnext is the estimated next stepsize. derivs
!   is the user-supplied subroutine that computes the right-hand side derivatives.
INTEGER i
REAL(8) errmax,h,htemp,xnew,yerr(NMAX),ytemp(NMAX),SAFETY,PGROW, PSHRNK,ERRCON
PARAMETER (SAFETY=0.9d0,PGROW=-.2d0,PSHRNK=-.25d0,ERRCON=1.89d-4)
!   The value ERRCON equals (5/SAFETY)**(1/PGROW), see use below.
h=htry !   Set stepsize to the initial trial value.
1 call rkck(y,dydx,n,x,h,ytemp,yerr,derivs) !   Take a step.
errmax=0.d0 !   Evaluate accuracy.
do i=1,n
errmax=max(errmax,dabs(yerr(i)/yscal(i)))
enddo 
errmax=errmax/eps !   Scale relative to required tolerance.
if(errmax.gt.1.d0)then !   Truncation error too large, reduce stepsize.
htemp=SAFETY*h*(errmax**PSHRNK)
h=dsign(max(dabs(htemp),0.1d0*dabs(h)),h) !   No more than a factor of 10.
xnew=x+h
if(xnew.eq.x) pause 'stepsize underflow in rkqs'
goto 1 !   For another try.
else !   Step succeeded. Compute size of next step.
if(errmax.gt.ERRCON)then
hnext=SAFETY*h*(errmax**PGROW)
else !   No more than a factor of 5 increase.
hnext=5.d0*h
endif
hdid=h
x=x+h
do i=1,n
y(i)=ytemp(i)
enddo
return
endif
ENDSUBROUTINE rkqs
!
! Cash-Karp Runge-Kutta step:
! [SY Lau 25052016] Added n to derivs argument due to local definition of fcn
SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr,derivs)
INTEGER n,NMAX
REAL(8) h,x,dydx(n),y(n),yerr(n),yout(n)
EXTERNAL derivs
PARAMETER (NMAX=50) !   Set to the maximum number of functions.
!   C USES derivs
!   Given values for n variables y and their derivatives dydx known at x, use the fth-order
!   Cash-Karp Runge-Kutta method to advance the solution over an interval h and return
!   the incremented variables as yout. Also return an estimate of the local truncation error
!   in yout using the embedded fourth-order method. The user supplies the subroutine
!   derivs(x,y,dydx), which returns derivatives dydx at x.
INTEGER i
REAL(8) ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX), ytemp(NMAX),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,&
& B52,B53,B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
PARAMETER (A2=.2d0,A3=.3d0,A4=.6d0,A5=1.d0,A6=.875d0,B21=.2d0,B31=3.d0/40.d0,B32=9.d0/40.d0,B41=.3d0,B42=-.9d0,B43=1.2d0,B51=-11.d0/54.d0,B52=2.5d0,&
& B53=-70.d0/27.d0,B54=35.d0/27.d0,B61=1631.d0/55296.d0,B62=175.d0/512.d0,B63=575.d0/13824.d0,B64=44275.d0/110592.d0,B65=253.d0/4096.d0,&
& C1=37.d0/378.d0,C3=250.d0/621.d0,C4=125.d0/594.d0,C6=512.d0/1771.d0,DC1=C1-2825.d0/27648.d0,DC3=C3-18575.d0/48384.d0,&
& DC4=C4-13525.d0/55296.d0,DC5=-277.d0/14336.d0,DC6=C6-.25d0)
do i=1,n !   First step.
ytemp(i)=y(i)+B21*h*dydx(i)
enddo
call derivs(n, x+A2*h,ytemp,ak2) !   Second step.
do i=1,n
ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
enddo
call derivs(n, x+A3*h,ytemp,ak3) !   Third step.
do i=1,n
ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
enddo
call derivs(n, x+A4*h,ytemp,ak4) !   Fourth step.
do i=1,n
ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+ B54*ak4(i))
enddo
call derivs(n, x+A5*h,ytemp,ak5) !   Fifth step.
do i=1,n
ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+ B64*ak4(i)+B65*ak5(i))
enddo
call derivs(n, x+A6*h,ytemp,ak6) !   Sixth step.
do i=1,n !   Accumulate increments with proper weights.
yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+ C6*ak6(i))
enddo
do i=1,n
!   Estimate error as difference between fourth and fth order methods.
yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i) +DC6*ak6(i))
enddo
return
ENDSUBROUTINE rkck
!
! Driver ***
! [SY Lau 25052016] Added nvar to derivs argument due to local definition of fcn
SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,rkqs)
INTEGER nbad,nok,nvar,KMAXX,MAXSTP,NMAX
REAL(8) eps,h1,hmin,x1,x2,ystart(nvar),TINY
EXTERNAL derivs,rkqs
PARAMETER (MAXSTP=10000,NMAX=50,KMAXX=200,TINY=1.d-30)
!   Runge-Kutta driver with adaptive stepsize control. Integrate the starting values ystart(1:nvar)
!   from x1 to x2 with accuracy eps, storing intermediate results in the common block /path/.
!   h1 should be set as a guessed rst stepsize, hmin as the minimum allowed stepsize (can
!   be zero). On output nok and nbad are the number of good and bad (but retried and fixed) steps taken, 
!   and ystart is replaced by values at the end of the integration interval.
!   derivs is the user-supplied subroutine for calculating the right-hand side derivative, while
!   rkqs is the name of the stepper routine to be used. /path/ contains its own information
!   about how often an intermediate value is to be stored.
INTEGER i,kmax,kount,nstp
REAL(8) dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),y(NMAX), yp(NMAX,KMAXX),yscal(NMAX)
COMMON /path/ kmax,kount,dxsav,xp,yp
!   User storage for intermediate results. Preset dxsav and kmax.
x=x1
h=dsign(h1,x2-x1)
nok=0
nbad=0
kount=0
do i=1,nvar
y(i)=ystart(i)
enddo
if (kmax.gt.0) xsav=x-2.d0*dxsav !   Assures storage of rst step.
do nstp=1,MAXSTP !   Take at most MAXSTP steps.
call derivs(nvar, x,y,dydx)
do i=1,nvar
!   Scaling used to monitor accuracy. This general-purpose choice can be modified if need be.
yscal(i)=dabs(y(i))+dabs(h*dydx(i))+TINY
enddo
if(kmax.gt.0)then
if(dabs(x-xsav).gt.dabs(dxsav)) then !   Store intermediate results.
if(kount.lt.kmax-1)then
kount=kount+1
xp(kount)=x
do i=1,nvar
yp(i,kount)=y(i)
enddo
xsav=x
endif
endif
endif
if((x+h-x2)*(x+h-x1).gt.0.d0) h=x2-x !   If stepsize can overshoot, decrease.
call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
if(hdid.eq.h)then
nok=nok+1
else
nbad=nbad+1
endif
if((x-x2)*(x2-x1).ge.0.d0)then !   Are we done?
do i=1,nvar
ystart(i)=y(i)
enddo
if(kmax.ne.0)then
kount=kount+1 !   Save final step.
xp(kount)=x
do i=1,nvar
yp(i,kount)=y(i)
enddo
endif
return !   Normal exit.
endif
if(dabs(hnext).lt.hmin) pause 'stepsize smaller than minimum in odeint'
h=hnext
enddo
pause 'too many steps in odeint'
return
ENDSUBROUTINE odeint

ENDMODULE RK4_Set