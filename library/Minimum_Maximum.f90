module Minimum_Maximum
contains
	!// Brent's method:
	!// little modification from "brent" in numerical recipe P.397
	!// accuracy limited by ZEPS
FUNCTION brent(ax,bx,cx,f,tol,fmin)
INTEGER ITMAX
REAL(8) brent,ax,bx,cx,tol,fmin,f,CGOLD,ZEPS
EXTERNAL f
PARAMETER (ITMAX=100,CGOLD=.3819660,ZEPS=1.0d-40)
!Given a function f, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is
!between ax and cx, and f(bx) is less than both f(ax) and f(cx)), this routine isolates
!the minimum to a fractional precision of about tol using Brent's method. The abscissa of
!the minimum is returned as brent, and the minimum function value is returned as fmin,
!Parameters: Maximum allowed number of iterations; golden ratio; and a small number that
!protects against trying to achieve fractional accuracy for a minimum that happens to be
!exactly zero.
INTEGER iter
REAL(8) a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm

	a=min(ax,cx) !a and b must be in ascending order, though the input
	b=max(ax,cx) !abscissas need not be.
	v=bx !Initializations...
	w=v
	x=v
	e=0. !This will be the distance moved on the step before last.
	fx=f(x)
	fv=fx
	fw=fx
	do iter=1,ITMAX !Main program loop.
		xm=0.5*(a+b)
		tol1=tol*abs(x)+ZEPS
		tol2=2.*tol1
		if(abs(x-xm).le.(tol2-.5*(b-a))) then
			goto 3 !Test for done here.
		endif
		if(abs(e).gt.tol1) then !Construct a trial parabolic t.
			r=(x-w)*(fx-fv)
			q=(x-v)*(fx-fw)
			p=(x-v)*q-(x-w)*r
			q=2.*(q-r)
			if(q.gt.0.) p=-p
			q=abs(q)
			etemp=e
			if(abs(p).ge.abs(.5*q*etemp).or.p.le.q*(a-x).or.p.ge.q*(b-x)) goto 1
			!The above conditions determine the acceptability of the parabolic t. Here it is o.k.:
			d=p/q !Take the parabolic step.
			u=x+d
			if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
				goto 2 !Skip over the golden section step.
		endif
1		if(x.ge.xm) then !We arrive here for a golden section step, which we take
			e=a-x !into the  larger of the two segments.
		else
			e=b-x
		endif
		d=CGOLD*e !Take the golden section step.
2		if(abs(d).ge.tol1) then !Arrive here with d computed either from parabolic t, or else from golden section.
			u=x+d 
		else
			u=x+sign(tol1,d)
		endif
		fu=f(u) !This is the one function evaluation per iteration,
		if(fu.le.fx) then !and now we have to decide what to do with our function
			if(u.ge.x) then !evaluation. Housekeeping follows:
				a=x
			else
				b=x
			endif
			v=w
			fv=fw
			w=x
			fw=fx
			x=u
			fx=fu
		else
			if(u.lt.x) then
				a=u
			else
				b=u
			endif
			if(fu.le.fw .or. w.eq.x) then
				v=w
				fv=fw
				w=u
				fw=fu
			else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
				v=u
				fv=fu
			endif
		endif !Done with housekeeping. Back for another iteration.
	enddo
	write(*,*) "err: brent exceed maximum iterations"
	pause !'brent exceed maximum iterations'
3	fmin=fx !Arrive here ready to exit with best values.
	brent=x
	return
ENDFUNCTION brent

endmodule Minimum_Maximum