module eigen_freq_opt01_ex2
!	output options
!	normalized with crust core interface amplitude
!	C1L1R1 stands for 1 Component, 1 Left (integration forward), 1 Right (integration backwards)
use global_var

real(8), allocatable :: rho_ar(:), dis_r(:), dis_t(:), r_ar(:), grad_r(:), grad_t(:)

contains

subroutine ei01_o_soln_C1L1R1
use Format_IO
use FWrite
use FRead
implicit none
real(8), dimension(0:pg_N1+1) :: ryA_, VA_
real(8), dimension(0:pg_N2+1) :: ryFM_, VFM_
real(8), dimension(1:2, 0:pg_N1+1) :: yA_
real(8), dimension(1:2, 0:pg_N2+1) :: yFM_
real(8), dimension(1:2) :: y_norm
real(8) :: ans, norm, Q, Q_U, Q_V, Q_V_temp1, Q_V_temp2, Q_Alt
integer :: nyA, nyFM, i
	call RFile3R8(pef_yA, nyA, ryA_, yA_(1,0:pg_N1+1), yA_(2,0:pg_N1+1))
	call RFile3R8(pef_y_FM_Cr, nyFM, ryFM_, yFM_(1,0:pg_N2+1), yFM_(2,0:pg_N2+1))

	if (nyA/= pg_N1+1 .or. nyFM/= pg_N2+1) then
		write(*,*) "err: linear combination, number of puls eqt grid mismatch"
		pause
	endif

	open(72, file='data/sol_y_Co.txt', status= 'replace')
	open(73, file='data/sol_y_Cr.txt', status= 'replace')
	open(74, file='data/sol_U.txt', status= 'replace')
	open(75, file='data/sol_V.txt', status= 'replace')
	open(76, file='data/sol_stat profile.txt', status= 'replace')
	open(77, file='data/sol_divergence.txt', status= 'replace')

	do i = 0, nyA-1
		yA_(1:2, i) = yA_(1:2, i)*ryA_(i)
		if (pes_opt == 1 .and. newt_V_opt == 1) then
			VA_(i) = ei01_V_conv(1, m_Co(2*i), ryA_(i), afreq, yA_(2, i))
		elseif (pes_opt == 1 .and. newt_V_opt == 2) then
			VA_(i) = ei01_V_conv(2, m_Co(2*i), ryA_(i), afreq, yA_(2, i), P_ = P_Co(2*i), rho_ = rho_Co(2*i))
		elseif (pes_opt == 1 .and. newt_V_opt == 3) then
			VA_(i) = ei01_V_conv(1, m_Co(2*i), ryA_(i), afreq, yA_(2, i))	
		elseif (pes_opt == 2) then
			VA_(i) = ei01_V_conv(3, m_Co(2*i), ryA_(i), afreq, yA_(2, i), P_ = P_Co(2*i), rho_ = rho_Co(2*i), Nu_ = Nu_Co(2*i))
		else
			write(*,*) "err: write files, pes mismmatch"
			pause
		endif
	enddo
		y_norm(1) = dabs(yA_(1, nyA-1))
		y_norm(2) = dabs(VA_(nyA-1))

	allocate(rho_ar(0:nyA-1),dis_r(0:nyA-1), dis_t(0:nyA-1), r_ar(0:nyA-1) )
	do i = 0, nyA-1
		call Write3R8(72, FR8, ryA_(i), yA_(1,i), yA_(2,i))
		call Write2R8(74, FR8, ryA_(i), yA_(1,i)/y_norm(1))
		call Write2R8(75, FR8, ryA_(i), VA_(i)/y_norm(1))  !/y_norm(2))
		call Write2R8(76, FR8, ryA_(i), rho_Co(2*i))

		rho_ar(i) = rho_Co(2*i)
		dis_r(i) = yA_(1,i)
		dis_t(i) = VA_(i)
		r_ar(i) = ryA_(i)
		
	enddo
do i = 0, nyA-1
call Write2R8(77, FR8, r_ar(i), ei01_rho_Xi(i))
enddo
	norm = 0.d0
	call ei01_integrate(ei01_norm, 0, nyA-1, 0.d0, ans)
	norm = norm + ans

	Q = 0.d0
	Q_U = 0.d0
	Q_V = 0.d0
	Q_Alt = 0.d0
	call ei01_integrate(ei01_Q, 0, nyA-1, 0.d0, ans)
	Q = Q + ans
	call ei01_integrate(ei01_Q_U, 0, nyA-1, 0.d0, ans)
	Q_U = Q_U + ans
	call ei01_integrate(ei01_Q_V, 0, nyA-1, 0.d0, ans)
	Q_V = Q_V + ans
Q_V_temp1 = Q_V
	call ei01_integrate(ei01_Q_alt, 1, nyA-2, 0.d0, ans)
	Q_Alt = Q_Alt + ans
	deallocate(rho_ar,dis_r, dis_t, r_ar )

	do i = 0, nyFM-1
		yFM_(1:2, i) = yFM_(1:2, i)*ryFM_(i)
		if (pes_opt == 1 .and. newt_V_opt == 1) then
			VFM_(i) = ei01_V_conv(1, m_Cr(2*(nyFM- 1-i)), ryFM_(i), afreq, yFM_(2, i))
		elseif (pes_opt == 1 .and. newt_V_opt == 2) then
			VFM_(i) = ei01_V_conv(2, m_Cr(2*(nyFM- 1-i)), ryFM_(i), afreq, yFM_(2, i), P_ = P_Cr(2*(nyFM- 1-i)), rho_ = rho_Cr(2*(nyFM- 1-i)))
		elseif (pes_opt == 1 .and. newt_V_opt == 3) then
			VFM_(i) = ei01_V_conv(1, m_Cr(2*(nyFM- 1-i)), ryFM_(i), afreq, yFM_(2, i))		
		elseif (pes_opt == 2) then
			VFM_(i) = ei01_V_conv(3, m_Cr(2*(nyFM- 1-i)), ryFM_(i), afreq, yFM_(2, i), P_ = P_Cr(2*(nyFM- 1-i)), rho_ = rho_Cr(2*(nyFM- 1-i)), Nu_ = Nu_Cr(2*(nyFM- 1-i)))
		else
			write(*,*) "err: write files, pes mismmatch"
			pause
		endif

	enddo

	allocate(rho_ar(0:nyFM-1),dis_r(0:nyFM-1), dis_t(0:nyFM-1), r_ar(0:nyFM-1) )
	do i = 0, nyFM-1
		call Write3R8(73, FR8, ryFM_(i), yFM_(1,i), yFM_(2,i))
		call Write2R8(74, FR8, ryFM_(i), yFM_(1,i)/y_norm(1))
		call Write2R8(75, FR8, ryFM_(i), VFM_(i)/y_norm(1))!/y_norm(2))
		call Write2R8(76, FR8, ryFM_(i), rho_Cr(2*(nyFM- 1-i)) )
		
!VFM_(nyFM- 1-i) = VFM_(nyFM- 1-i)/1.5492d0
		rho_ar(i) = rho_Cr(2*i)
		dis_r(i) = yFM_(1,nyFM- 1-i)
		dis_t(i) = VFM_(nyFM- 1-i)
		r_ar(i) = ryFM_(nyFM- 1-i)	
	enddo
do i = 0, nyFM-1
call Write2R8(77, FR8, r_ar(i), ei01_rho_Xi(i))
enddo
	call ei01_integrate(ei01_norm, 0, nyFM-1,0.d0, ans)
	norm = norm + ans
	norm = dsqrt(M0*R0**2/norm)

	call ei01_integrate(ei01_Q, 0, nyFM-1, 0.d0, ans)
	Q = Q + ans
	call ei01_integrate(ei01_Q_U, 0, nyFM-1, 0.d0, ans)
	Q_U = Q_U + ans
	call ei01_integrate(ei01_Q_V, 0, nyFM-1, 0.d0, ans)
	Q_V = Q_V + ans
Q_V_temp2 = Q_V - Q_V_temp1
	call ei01_integrate(ei01_Q_alt, 1, nyFM-2, 0.d0, ans)
	Q_Alt = Q_Alt + ans
	deallocate(rho_ar,dis_r, dis_t, r_ar )

write(*,*) "Radial Displacement/ Tangential displacement, in = ", yA_(1,nyA -1)/VA_(nyA -1)
write(*,*) "Tangential displacement/Radial Displacement, out = ", VFM_(nyFM- 1)/yFM_(1,nyFM- 1)
write(*,*) "Tangential displacement jump = ", (VA_(nyA -1) - VFM_(nyFM-1))/ VA_(nyA -1)
write(*,*) "Fluid transverse displacement jump , 0th order =, M1 ", Grav_Const * m_Co(2* pg_N1) * yA_(1, nyA-1) /VA_(nyA -1) / afreq**2 / r_Co(2* pg_N1) ** 3 * (rho_Co(2* pg_N1) - rho_Cr(0))/ rho_Co(2* pg_N1)	! NCA rho g
write(*,*) "Fluid transverse displacement jump , 0th order, M2 = ", Grav_Const * m_Co(2* pg_N1) * yA_(1, nyA-1) /VA_(nyA -1) / afreq**2 / r_Co(2* pg_N1) ** 3 * (rho_Co(2* pg_N1) - rho_Cr(0))/ rho_Co(2* pg_N1) /(1.d0 - 2.d0 * Grav_Const * m_Co(2*pg_N1)/c**2/r_Co(2*pg_N1)) ! NCA dP/dr
write(*,*) "Fluid transverse displacement jump , 0th order, M3 = ", Grav_Const * m_Co(2* pg_N1) * yA_(1, nyA-1) /VA_(nyA -1) / afreq**2 / r_Co(2* pg_N1) ** 3 * (rho_Co(2* pg_N1) - rho_Cr(0))/ rho_Co(2* pg_N1) /(1.d0 - 2.d0 * Grav_Const * m_Co(2*pg_N1)/c**2/r_Co(2*pg_N1)) / (dexp(-2.d0* nu_Co(2*pg_N1))) ! RCA
write(*,*) "1st order correction = ", - VFM_(nyFM-1)/VA_(nyA -1) * rho_Cr(0)/rho_Co(2* pg_N1)*(rho_Co(2* pg_N1) - rho_Cr(0))/rho_Co(2* pg_N1)
write(*,*) "2nd order correction = ", - VFM_(nyFM-1)/VA_(nyA -1) * rho_Cr(0)/rho_Co(2* pg_N1)*((rho_Co(2* pg_N1) - rho_Cr(0))/rho_Co(2* pg_N1))**2
write(*,*) "U/V",yA_(1, nyA-1) /VA_(nyA -1)
write(*,*) "1/sigma^2",1.d0/ afreq**2
write(*,*) "e^(2Lambda)",1.d0 /(1.d0 - 2.d0 * Grav_Const * m_Co(2*pg_N1)/c**2/r_Co(2*pg_N1))
write(*,*) "e^(-2nu)", dexp(-2.d0* nu_Co(2*pg_N1))
write(*,*) "=========================================================="
write(*,*) "normalization constant = ", norm
write(*,*) "MR^2 ", M0*R0**2
write(*,*) "Q factor = ", (Q)*norm/(M0*R0**2)
write(*,*) "Q_U factor = ", (Q_U)*norm/(M0*R0**2)
write(*,*) "Q_V factor = ", (Q_V)*norm/(M0*R0**2)
write(*,*) "Q_V factor inner = ", (Q_V_temp1)*norm/(M0*R0**2)
write(*,*) "Q_V factor  outer = ", (Q_V_temp2)*norm/(M0*R0**2)
write(*,*) "Q_Alt factor in continuous regime = ", (Q_Alt)*norm/(M0*R0**2)
write(*,*) "Q_Alt Not Norm factor in continuous regime = ", (Q_Alt)
write(*,*) "Q_Alt_surface term", rho_Cr(2*(nyFM-1)) * yFM_(1,0) * ryFM_(0)**(l_0+2.d0) *norm/(M0*R0**2)
write(*,*) "Q_Alt jump accross interface", -ryFM_(nyFM-1) **(l_0 + 2.d0) * (rho_Cr(0) * yFM_(1,nyFM-1) - rho_Co(2*(nyA-1))* yA_(1,nyA -1)) *norm/(M0*R0**2)
write(*,*) "=========================================================="

!call ei01_poloidal_projection(nyA, yA_, VA_, ryA_, nyFM, yFM_, VFM_, ryFM_, norm)
call ei01_Greens(nyA, yA_, VA_, ryA_, nyFM, yFM_, VFM_, ryFM_, norm)
	close(72)
	close(73)
	close(74)
	close(75)
	close(76)
	close(77)
endsubroutine ei01_o_soln_C1L1R1

subroutine ei01_poloidal_projection(nyA, yA_, VA_, ryA_, nyFM, yFM_, VFM_, ryFM_, norm)
!	To find the fraction which the mode function is composed of poloidal field
use Format_IO
use FWrite
implicit none
integer :: nyA, nyFM, i, j, nt, n0
real(8), dimension(1:, 0:) :: yA_, yFM_
real(8), dimension(0:) :: VA_, ryA_, VFM_, ryFM_
real(8) :: Pol_Coe, Sph_Coe, rho_U_In, fcn10_, ans, pol_fraction, norm
real(8) :: In_A, d_A, Pol_Coe_i, Sph_Coe_i, ans_pre

	nt = nyA - 1
	n0 = nyA + nyFM - 1

	rho_U_In = 0.d0
	Pol_Coe = 0.d0
	Sph_Coe = 0.d0

	!	Core part
	allocate(rho_ar(0:nyA-1),dis_r(0:nyA-1), dis_t(0:nyA-1), r_ar(0:nyA-1) )

	do i = 0, nyA-1
		rho_ar(i) = rho_Co(2*i)
		dis_r(i) = yA_(1,i) *norm
		dis_t(i) = VA_(i) *norm
		r_ar(i) = ryA_(i)
	enddo

	call ei01_integrate(ei01_fcn_rho_U, 0, nyA-1, 0.d0, ans)
	rho_U_In = rho_U_In + ans

	Sph_Coe = Sph_Coe + r_ar(nyA-1) * rho_ar(nyA-1) * dis_t(nyA-1) - r_ar(0) * rho_ar(0) * dis_t(0)

	deallocate(rho_ar,dis_r, dis_t, r_ar )
	
	!	Crust part
	allocate(rho_ar(0:nyFM-1),dis_r(0:nyFM-1), dis_t(0:nyFM-1), r_ar(0:nyFM-1) )
	do i = 0, nyFM-1
		rho_ar(i) = rho_Cr(2*i)
		dis_r(i) = yFM_(1,nyFM- 1-i) *norm
		dis_t(i) = VFM_(nyFM- 1-i) *norm
		r_ar(i) = ryFM_(nyFM- 1-i)	
	enddo

	call ei01_integrate(ei01_fcn_rho_U, 0, nyFM-1, 0.d0, ans)
	rho_U_In = rho_U_In + ans

	Sph_Coe = Sph_Coe + r_ar(nyFM-1) * rho_ar(nyFM-1) * dis_t(nyFM-1) - r_ar(0) * rho_ar(0) * dis_t(0)

	deallocate(rho_ar,dis_r, dis_t, r_ar )

	!	Throuhout whole star
	allocate(rho_ar(0:n0),dis_r(0:n0), dis_t(0:n0), r_ar(0:n0) )
	do i = 0, nt
		rho_ar(i) = rho_Co(2*i)
		dis_r(i) = yA_(1,i) *norm
		dis_t(i) = VA_(i) *norm
		r_ar(i) = ryA_(i)
	enddo
	do i = nt + 1, n0
		j = i - (nt + 1)
		rho_ar(i) = rho_Cr(2*j)
		dis_r(i) = yFM_(1,nyFM- 1-j) *norm
		dis_t(i) = VFM_(nyFM- 1-j) *norm
		r_ar(i) = ryFM_(nyFM- 1-j)
	enddo
	
	fcn10_ = VA_(0) * rho_ar(0) * l_0 * (l_0 + 1.d0)/2.d0
	call ei01_double_integrate(ei01_fcn_Double1, ei01_fcn_Double2, 0, nt, fcn10_, 0.d0, ans)
	Pol_Coe = Pol_Coe + ans
	
	pol_fraction = (rho_U_In - Sph_Coe)/(Pol_Coe - Sph_Coe)
	write(*,*) "rho_U_In", rho_U_In
	write(*,*) "Pol_Coe", Pol_Coe
	write(*,*) "Sph_Coe", Sph_Coe

	write(*,*) "pol_fraction", pol_fraction
	write(*,*) "sph_fraction", 1.d0 - pol_fraction
	write(*,*) "=========================================================="
	
	
	open(01, file = 'data/sol_Pol_Sph_Coe.txt', status = 'replace')
	In_A = fcn10_
	d_A = VA_(0) * rho_ar(0)
	Pol_Coe_i = (  rho_ar(0) * dis_r(0) - d_A  )/(  In_A - d_A  )
	Sph_Coe_i = 1.d0 - Pol_Coe_i
	
	call Write3R8(01, FR8, r_ar(0), Pol_Coe_i, Sph_Coe_i)

	ans_pre = 0.d0
	do i = 1, nt
		call ei01_integrate(ei01_fcn_Double2, i-1, i, ans_pre, ans)
		In_A = l_0 * (l_0 + 1.d0)/r_ar(i)**2 * ans
		ans_pre = ans

		d_A = ei01_dA(i)
		Pol_Coe_i = (  rho_ar(i) * dis_r(i) - d_A  )/(  In_A - d_A  )
		Sph_Coe_i = 1.d0 - Pol_Coe_i
		call Write3R8(01, FR8, r_ar(i), Pol_Coe_i, Sph_Coe_i)
	enddo
	call ei01_integrate(ei01_fcn_Double2, nt, nt + 1, ans_pre, ans)
	ans_pre = ans
	do i = nt + 2, n0
		call ei01_integrate(ei01_fcn_Double2, i-1, i, ans_pre, ans)
		In_A = l_0 * (l_0 + 1.d0)/r_ar(i)**2 * ans
		ans_pre = ans

		d_A = ei01_dA(i)
		Pol_Coe_i = (  rho_ar(i) * dis_r(i) - d_A  )/(  In_A - d_A  )
		Sph_Coe_i = 1.d0 - Pol_Coe_i
		call Write3R8(01, FR8, r_ar(i), Pol_Coe_i, Sph_Coe_i)
	enddo
	
	close(01)
	deallocate(rho_ar,dis_r, dis_t, r_ar )

	

endsubroutine ei01_poloidal_projection

subroutine ei01_Greens(nyA, yA_, VA_, ryA_, nyFM, yFM_, VFM_, ryFM_, norm)
!	To find the fraction which the mode function is composed of poloidal field
use Format_IO
use FWrite
implicit none
integer :: nyA, nyFM, i, j, nt, n0
real(8), dimension(1:, 0:) :: yA_, yFM_
real(8), dimension(0:) :: VA_, ryA_, VFM_, ryFM_
real(8) :: ans, norm
real(8), dimension(:), allocatable :: green, dgreen

	nt = nyA - 1		! nt -> nyA -1
	n0 = nyA + nyFM - 1	!  nt + 1 -> 0; n0 -> nyFM - 1

	!	Throuhout whole star
	allocate(rho_ar(0:n0),dis_r(0:n0), dis_t(0:n0), r_ar(0:n0), green(0:n0), dgreen(0:n0) ,grad_r(0:n0), grad_t(0:n0))
	do i = 0, nt
		rho_ar(i) = rho_Co(2*i)
		dis_r(i) = yA_(1,i) *norm
		dis_t(i) = VA_(i) *norm
		r_ar(i) = ryA_(i)
	enddo
	do i = nt + 1, n0
		j = i - (nt + 1)
		rho_ar(i) = rho_Cr(2*j)
		dis_r(i) = yFM_(1,nyFM- 1-j) *norm
		dis_t(i) = VFM_(nyFM- 1-j) *norm
		r_ar(i) = ryFM_(nyFM- 1-j)
	enddo

	open(01, file = 'data/sol_Grad_Field.txt', status = 'replace')

	do i = 1, nt-1
		green(i) = 0.d0
		call ei01_integrate(ei01_Greens_Inner, 0, i, 0.d0, ans)
		green(i) = green(i) + ans /r_ar(i) ** (l_0 + 1.d0)

		call ei01_integrate(ei01_Greens_Outer, i, nt-1, 1.d0/r_ar(i)**(l_0 - 1.d0) * ei01_rho_Xi(i), ans)
		green(i) = green(i) + ans * r_ar(i) ** (l_0)

		green(i) = green(i) + r_ar(i) ** (l_0) /r_ar(nt + 1)** (l_0 - 1.d0) * dis_r(nt + 1) * (rho_ar(nt+1) - rho_ar(nt))

		call ei01_integrate(ei01_Greens_Outer, nt +2, n0, 1.d0 /r_ar(nt +2)**(l_0 - 1.d0) * ei01_rho_Xi(nt +2), ans)
		green(i) = green(i) + ans * r_ar(i) ** (l_0)

		green(i) = -green(i)/(2.d0 * l_0 + 1.d0)

		dgreen(i) = (green(i) - green(i-1) ) / (r_ar(i) - r_ar(i-1))

		grad_r(i) = dgreen(i)
		grad_t(i) = green(i)/r_ar(i)
	
		call Write6R8(01, FR8, r_ar(i), green(i), dgreen(i), green(i)/r_ar(i), rho_ar(i)*dis_r(i), rho_ar(i)*dis_t(i))
	enddo
	
	do i = nt + 2, n0
		green(i) = 0.d0
		call ei01_integrate(ei01_Greens_Inner, 0, nt-1, 0.d0, ans)
		green(i) = green(i) + ans /r_ar(i) ** (l_0 + 1.d0)

		green(i) = green(i) + r_ar(nt + 1) ** (l_0 + 2.d0) /r_ar(i)** (l_0 + 1.d0) * dis_r(nt + 1) * (rho_ar(nt+1) - rho_ar(nt))

		call ei01_integrate(ei01_Greens_Inner, nt+2, i, r_ar(nt + 2)**(l_0 + 2.d0) * ei01_rho_Xi(nt + 2), ans)
		green(i) = green(i) + ans /r_ar(i) ** (l_0 + 1.d0)

		call ei01_integrate(ei01_Greens_Outer, i, n0, 1.d0/r_ar(i)**(l_0 - 1.d0) * ei01_rho_Xi(i), ans)
		green(i) = green(i) + ans * r_ar(i) ** (l_0)

		green(i) = -green(i)/(2.d0 * l_0 + 1.d0)

		dgreen(i) = (green(i) - green(i-1) ) / (r_ar(i) - r_ar(i-1))

		grad_r(i) = dgreen(i)
		grad_t(i) = green(i)/r_ar(i)

		call Write6R8(01, FR8, r_ar(i), green(i), dgreen(i), green(i)/r_ar(i), rho_ar(i)*dis_r(i), rho_ar(i)*dis_t(i))
	enddo
	
	close(01)

	open(01, file = 'data/sol_grad_Verify.txt', status = 'replace')
	do i = 1, nt-3
		call Write3R8(01, FR8, r_ar(i), ei01_Div_grad(i), ei01_rho_Xi(i))
	enddo
	do i = nt + 3, n0
		call Write3R8(01, FR8, r_ar(i), ei01_Div_grad(i), ei01_rho_Xi(i))
	enddo
	close(01)

	deallocate(rho_ar,dis_r, dis_t, r_ar, green, dgreen,grad_r, grad_t )

	

endsubroutine ei01_Greens

subroutine ei01_o_soln_C2L1R2
use Format_IO
use FWrite
use FRead
implicit none
real(8), dimension(0:pg_N1+1) :: ryA_, VA_
real(8), dimension(0:pg_N2+1) :: rzC_, rzD_
real(8), dimension(1:2, 0:pg_N1+1) :: yA_
real(8), dimension(1:4, 0:pg_N2+1) :: zC_, zD_, z_
real(8), dimension(1:2) :: y_norm
real(8) :: norm, Q, Q_U, Q_V, ans
integer :: nyA, nzC, nzD, i
	call RFile3R8(pef_yA, nyA, ryA_, yA_(1,0:pg_N1+1), yA_(2,0:pg_N1+1))
	call RFile5R8(pef_zC, nzC, rzC_, zC_(1,0:pg_N2+1), zC_(2,0:pg_N2+1), zC_(3,0:pg_N2+1), zC_(4,0:pg_N2+1))
	call RFile5R8(pef_zD, nzD, rzD_, zD_(1,0:pg_N2+1), zD_(2,0:pg_N2+1), zD_(3,0:pg_N2+1), zD_(4,0:pg_N2+1))

	if (nyA/= pg_N1+1 .or. nzC/= pg_N2+1 .or. nzD/= pg_N2+1) then
		write(*,*) "err: linear combination, number of puls eqt grid mismatch"
		pause
	endif

	open(72, file='data/sol_y_Co.txt', status= 'replace')
	open(73, file='data/sol_zC_Cr.txt', status= 'replace')
	open(74, file='data/sol_zD_Cr.txt', status= 'replace')
	open(75, file='data/sol_z_Cr.txt', status= 'replace')
	open(76, file='data/sol_U.txt', status= 'replace')
	open(77, file='data/sol_V.txt', status= 'replace')
	
	do i = 0, nyA-1
		yA_(1:2, i) = yA_(1:2, i)*ryA_(i)
		if (pes_opt == 1 .and. newt_V_opt == 1) then
			VA_(i) = ei01_V_conv(1, m_Co(2*i), ryA_(i), afreq, yA_(2, i))
		elseif (pes_opt == 1 .and. newt_V_opt == 2) then
			VA_(i) = ei01_V_conv(2, m_Co(2*i), ryA_(i), afreq, yA_(2, i), P_ = P_Co(2*i), rho_ = rho_Co(2*i))
		elseif (pes_opt == 1 .and. newt_V_opt == 3) then
			VA_(i) = ei01_V_conv(1, m_Co(2*i), ryA_(i), afreq, yA_(2, i), P_ = P_Co(2*i), rho_ = rho_Co(2*i))
		elseif (pes_opt == 2) then
			VA_(i) = ei01_V_conv(3, m_Co(2*i), ryA_(i), afreq, yA_(2, i), P_ = P_Co(2*i), rho_ = rho_Co(2*i), Nu_ = Nu_Co(2*i))
		else
			write(*,*) "err: write files, pes mismmatch"
			pause
		endif
	enddo
		y_norm(1) = dabs(yA_(1, nyA-1))
		y_norm(2) = dabs(VA_(nyA-1))

	allocate(rho_ar(0:nyA-1),dis_r(0:nyA-1), dis_t(0:nyA-1), r_ar(0:nyA-1) )
	do i = 0, nyA-1
		call Write3R8(72, FR8, ryA_(i), yA_(1,i), yA_(2,i))
		call Write2R8(76, FR8, ryA_(i), yA_(1,i)/y_norm(1))
		call Write2R8(77, FR8, ryA_(i), VA_(i)/y_norm(1))
!call Write2R8(76, FR8, ryA_(i), yA_(1,i))
!call Write2R8(77, FR8, ryA_(i), VA_(i))
		rho_ar(i) = rho_Co(2*i)
		dis_r(i) = yA_(1,i)
		dis_t(i) = VA_(i)
		r_ar(i) = ryA_(i)
	enddo
	norm = 0.d0
	call ei01_integrate(ei01_norm, 0, nyA-1, 0.d0, ans)
	norm = norm + ans

	Q = 0.d0
	Q_U = 0.d0
	Q_V = 0.d0
	call ei01_integrate(ei01_Q, 0, nyA-1, 0.d0, ans)
	Q = Q + ans
	call ei01_integrate(ei01_Q_U, 0, nyA-1, 0.d0, ans)
	Q_U = Q_U + ans
	call ei01_integrate(ei01_Q_V, 0, nyA-1, 0.d0, ans)
	Q_V = Q_V + ans
	deallocate(rho_ar,dis_r, dis_t, r_ar )

	do i = 0, nzC-1	
		zC_(1:4, i) = zC_(1:4, i)*rzC_(i)
		zD_(1:4, i) = zD_(1:4, i)*rzD_(i)
		z_(1:4, i) = zC_(1:4, i) + zD_(1:4, i)
	enddo

	allocate(rho_ar(0:nzC-1),dis_r(0:nzC-1), dis_t(0:nzC-1), r_ar(0:nzC-1) )
	do i = 0, nzC-1	
		call Write5R8(73, FR8, rzC_(i), zC_(1,i), zC_(2,i), zC_(3,i), zC_(4,i))
		call Write5R8(74, FR8, rzD_(i), zD_(1,i), zD_(2,i), zD_(3,i), zD_(4,i))
		call Write5R8(75, FR8, rzC_(i), z_(1,i), z_(2,i), z_(3,i), z_(4,i))
		call Write2R8(76, FR8, rzC_(i), z_(1,i)/y_norm(1))
		call Write2R8(77, FR8, rzC_(i), z_(3,i)/y_norm(1))
!call Write2R8(76, FR8, rzC_(i), z_(1,i))
!call Write2R8(77, FR8, rzC_(i), z_(3,i))
		
		rho_ar(i) = rho_Cr(2*i)
		dis_r(i) = z_(1,nzC- 1-i)
		dis_t(i) = z_(3,nzC- 1-i)
		r_ar(i) = rzC_(nzC- 1-i)
	enddo
	call ei01_integrate(ei01_norm, 0, nzC-1, 0.d0, ans)
	norm = norm + ans
	norm = dsqrt(M0*R0**2/norm)
	call ei01_integrate(ei01_Q, 0, nzC-1, 0.d0, ans)
	Q = Q + ans
	call ei01_integrate(ei01_Q_U, 0, nzC-1, 0.d0, ans)
	Q_U = Q_U + ans
	call ei01_integrate(ei01_Q_V, 0, nzC-1, 0.d0, ans)
	Q_V = Q_V + ans
	deallocate(rho_ar,dis_r, dis_t, r_ar )

write(*,*) "Radial Displacement/ Tangential displacement = ", z_(3,nzC-1)/VA_(nyA -1)
write(*,*) "Tangential displacement jump = ", (VA_(nyA -1) - z_(3,nzC-1))/ VA_(nyA -1)
write(*,*) "Fluid transverse displacement jump , 0th order =, M1 ", Grav_Const * m_Co(2* pg_N1) * yA_(1, nyA-1) /VA_(nyA -1) / afreq**2 / r_Co(2* pg_N1) ** 3 * (rho_Co(2* pg_N1) - rho_Cr(0))/ rho_Co(2* pg_N1)	! NCA rho g
write(*,*) "Fluid transverse displacement jump , 0th order, M2 = ", Grav_Const * m_Co(2* pg_N1) * yA_(1, nyA-1) /VA_(nyA -1) / afreq**2 / r_Co(2* pg_N1) ** 3 * (rho_Co(2* pg_N1) - rho_Cr(0))/ rho_Co(2* pg_N1) /(1.d0 - 2.d0 * Grav_Const * m_Co(2*pg_N1)/c**2/r_Co(2*pg_N1)) ! NCA dP/dr
write(*,*) "Fluid transverse displacement jump , 0th order, M3 = ", Grav_Const * m_Co(2* pg_N1) * yA_(1, nyA-1) /VA_(nyA -1) / afreq**2 / r_Co(2* pg_N1) ** 3 * (rho_Co(2* pg_N1) - rho_Cr(0))/ rho_Co(2* pg_N1) /(1.d0 - 2.d0 * Grav_Const * m_Co(2*pg_N1)/c**2/r_Co(2*pg_N1)) / (dexp(-2.d0* nu_Co(2*pg_N1))) ! RCA
write(*,*) "1st order correction = ", - z_(3,nzC-1)/VA_(nyA -1) * rho_Cr(0)/rho_Co(2* pg_N1)*(rho_Co(2* pg_N1) - rho_Cr(0))/rho_Co(2* pg_N1)
write(*,*) "2nd order correction = ", - z_(3,nzC-1)/VA_(nyA -1) * rho_Cr(0)/rho_Co(2* pg_N1)*((rho_Co(2* pg_N1) - rho_Cr(0))/rho_Co(2* pg_N1))**2
write(*,*) "=========================================================="
write(*,*) "normalization constant = ", norm
write(*,*) "MR^2 ", M0*R0**2
write(*,*) "Q factor = ", (Q)*norm/(M0*R0**2)
write(*,*) "Q_U factor = ", (Q_U)*norm/(M0*R0**2)
write(*,*) "Q_V factor = ", (Q_V)*norm/(M0*R0**2)
write(*,*) "=========================================================="


call ei01_Greens(nyA, yA_, VA_, ryA_, nzC, z_(1:2,0:nzC-1), z_(3,0:nzC-1), rzC_, norm)
!do i = 0, 1000
!	r_ar(i) = -1.d0/1000*i
!enddo
!call ei01_integrate(trial, 0, 1000, 0.d0, ans)
!write(*,*) ans

	close(72)
	close(73)
	close(74)
	close(75)
	close(76)
	close(77)

endsubroutine ei01_o_soln_C2L1R2

subroutine ei01_o_soln_C3L1R2
use Format_IO
use FWrite
use FRead
implicit none
real(8), dimension(0:pg_N1+1) :: ryA_, VA_
real(8), dimension(0:pg_N2+1) :: rzC_, rzD_
real(8), dimension(0:pg_N3+1) :: ryB_, VB_
real(8), dimension(1:2, 0:pg_N1+1) :: yA_
real(8), dimension(1:4, 0:pg_N2+1) :: zC_, zD_, z_
real(8), dimension(1:2, 0:pg_N3+1) :: yB_
real(8), dimension(1:2) :: y_norm
integer :: nyA, nzC, nzD, nyB, i
	call RFile3R8(pef_yA, nyA, ryA_, yA_(1,0:pg_N1+1), yA_(2,0:pg_N1+1))
	call RFile5R8(pef_zC, nzC, rzC_, zC_(1,0:pg_N2+1), zC_(2,0:pg_N2+1), zC_(3,0:pg_N2+1), zC_(4,0:pg_N2+1))
	call RFile5R8(pef_zD, nzD, rzD_, zD_(1,0:pg_N2+1), zD_(2,0:pg_N2+1), zD_(3,0:pg_N2+1), zD_(4,0:pg_N2+1))
	call RFile3R8(pef_yB, nyB, ryB_, yB_(1,0:pg_N3+1), yB_(2,0:pg_N3+1))

	if (nyA/= pg_N1+1 .or. nzC/= pg_N2+1 .or. nzD/= pg_N2+1) then
		write(*,*) "err: linear combination, number of puls eqt grid mismatch"
		pause
	endif

	open(72, file='data/sol_y_Co.txt', status= 'replace')
	open(73, file='data/sol_zC_Cr.txt', status= 'replace')
	open(74, file='data/sol_zD_Cr.txt', status= 'replace')
	open(75, file='data/sol_z_Cr.txt', status= 'replace')
	open(76, file='data/sol_y_Oc.txt', status= 'replace')
	open(77, file='data/sol_U.txt', status= 'replace')
	open(78, file='data/sol_V.txt', status= 'replace')
	
	do i = 0, nyA-1
		yA_(1:2, i) = yA_(1:2, i)*ryA_(i)
		if (pes_opt == 1 .and. newt_V_opt == 1) then
			VA_(i) = ei01_V_conv(1, m_Co(2*i), ryA_(i), afreq, yA_(2, i))
		elseif (pes_opt == 1 .and. newt_V_opt == 2) then
			VA_(i) = ei01_V_conv(2, m_Co(2*i), ryA_(i), afreq, yA_(2, i))
		elseif (pes_opt == 1 .and. newt_V_opt == 3) then
			VA_(i) = ei01_V_conv(1, m_Co(2*i), ryA_(i), afreq, yA_(2, i))
		elseif (pes_opt == 2) then
			VA_(i) = ei01_V_conv(3, m_Co(2*i), ryA_(i), afreq, yA_(2, i), P_ = P_Co(2*i), rho_ = rho_Co(2*i), Nu_ = Nu_Co(2*i))
		
		else
			write(*,*) "err: write files, pes mismmatch"
			pause
		endif
	enddo
		y_norm(1) = dabs(yA_(1, nyA-1))
		y_norm(2) = dabs(VA_(nyA-1))
	do i = 0, nyA-1
		call Write3R8(72, FR8, ryA_(i), yA_(1,i), yA_(2,i))
		call Write2R8(77, FR8, ryA_(i), yA_(1,i)/y_norm(1))
		call Write2R8(78, FR8, ryA_(i), VA_(i)/y_norm(2))
	enddo

	do i = 0, nzC-1	
		zC_(1:4, i) = zC_(1:4, i)*rzC_(i)
		zD_(1:4, i) = zD_(1:4, i)*rzD_(i)
		z_(1:4, i) = zC_(1:4, i) + zD_(1:4, i)
		call Write5R8(73, FR8, rzC_(i), zC_(1,i), zC_(2,i), zC_(3,i), zC_(4,i))
		call Write5R8(74, FR8, rzD_(i), zD_(1,i), zD_(2,i), zD_(3,i), zD_(4,i))
		call Write5R8(75, FR8, rzC_(i), z_(1,i), z_(2,i), z_(3,i), z_(4,i))
		call Write2R8(77, FR8, rzC_(i), z_(1,i)/y_norm(1))
		call Write2R8(78, FR8, rzC_(i), z_(3,i)/y_norm(2))
	enddo
	do i = 0, nyB-1
		yB_(1:2, i) = yB_(1:2, i)*ryB_(i)
		if (pes_opt == 1 .and. newt_V_opt == 1) then
			VB_(i) = ei01_V_conv(1, m_Oc(2*i), ryB_(i), afreq, yB_(2, i))
		elseif (pes_opt == 1 .and. newt_V_opt == 2) then
			VB_(i) = ei01_V_conv(2, m_Oc(2*i), ryB_(i), afreq, yB_(2, i))
		elseif (pes_opt == 1 .and. newt_V_opt == 3) then
			VB_(i) = ei01_V_conv(1, m_Oc(2*i), ryB_(i), afreq, yB_(2, i))
		elseif (pes_opt == 2) then
			VB_(i) = ei01_V_conv(3, m_Oc(2*i), ryB_(i), afreq, yB_(2, i), P_ = P_Oc(2*i), rho_ = rho_Oc(2*i), Nu_ = Nu_Oc(2*i))
		else
			write(*,*) "err: write files, pes mismmatch"
			pause
		endif
		call Write3R8(76, FR8, ryB_(i), yB_(1,i), yB_(2,i))
		call Write2R8(77, FR8, ryB_(i), yB_(1,i)/y_norm(1))
		call Write2R8(78, FR8, ryB_(i), VB_(i)/y_norm(2))
	enddo

	close(72)
	close(73)
	close(74)
	close(75)
	close(76)
	close(77)
	close(78)

endsubroutine ei01_o_soln_C3L1R2


subroutine ei01_o_soln_C2L2R5
use Format_IO
use FWrite
use FRead
implicit none
real(8), dimension(0:pg_N1+1) :: ryA_
real(8), dimension(0:pg_N2+1) :: rzA_
real(8), dimension(1:2, 1:4, 0:pg_N1+1) :: yA_
real(8), dimension(1:2, 0:pg_N1+1) :: VA_
real(8), dimension(1:4, 0:pg_N1+1) :: y_
real(8), dimension(1:5, 1:6, 0:pg_N2+1) :: zA_
real(8), dimension(1:6, 0:pg_N2+1) :: z_
real(8), dimension(1:2) :: y_norm
integer :: nyA, nzA, i

	call RFile5R8(pef_FN_Co1, nyA, ryA_, yA_(1, 1,0:pg_N1+1), yA_(1, 2,0:pg_N1+1), yA_(1, 3,0:pg_N1+1), yA_(1, 4,0:pg_N1+1))
	call RFile5R8(pef_FN_Co2, nyA, ryA_, yA_(2, 1,0:pg_N1+1), yA_(2, 2,0:pg_N1+1), yA_(2, 3,0:pg_N1+1), yA_(2, 4,0:pg_N1+1))

	call RFile7R8(pef_FN_Cr1, nzA, rzA_, zA_(1, 1,0:pg_N2+1), zA_(1, 2,0:pg_N2+1), zA_(1, 3,0:pg_N2+1), zA_(1, 4,0:pg_N2+1), zA_(1, 5,0:pg_N2+1), zA_(1, 6,0:pg_N2+1))
	call RFile7R8(pef_FN_Cr2, nzA, rzA_, zA_(2, 1,0:pg_N2+1), zA_(2, 2,0:pg_N2+1), zA_(2, 3,0:pg_N2+1), zA_(2, 4,0:pg_N2+1), zA_(2, 5,0:pg_N2+1), zA_(2, 6,0:pg_N2+1))
	call RFile7R8(pef_FN_Cr3, nzA, rzA_, zA_(3, 1,0:pg_N2+1), zA_(3, 2,0:pg_N2+1), zA_(3, 3,0:pg_N2+1), zA_(3, 4,0:pg_N2+1), zA_(3, 5,0:pg_N2+1), zA_(3, 6,0:pg_N2+1))
	call RFile7R8(pef_FN_Cr4, nzA, rzA_, zA_(4, 1,0:pg_N2+1), zA_(4, 2,0:pg_N2+1), zA_(4, 3,0:pg_N2+1), zA_(4, 4,0:pg_N2+1), zA_(4, 5,0:pg_N2+1), zA_(4, 6,0:pg_N2+1))
	call RFile7R8(pef_FN_Cr5, nzA, rzA_, zA_(5, 1,0:pg_N2+1), zA_(5, 2,0:pg_N2+1), zA_(5, 3,0:pg_N2+1), zA_(5, 4,0:pg_N2+1), zA_(5, 5,0:pg_N2+1), zA_(5, 6,0:pg_N2+1))

	if (nyA/= pg_N1+1 .or. nzA/= pg_N2+1) then
		write(*,*) "err: linear combination, number of puls eqt grid mismatch"
		pause
	endif

	open(72, file='data/sol_y_Co.txt', status= 'replace')
	open(75, file='data/sol_z_Cr.txt', status= 'replace')
	open(76, file='data/sol_U.txt', status= 'replace')
	open(77, file='data/sol_V.txt', status= 'replace')
	
	do i = 0, nyA-1
		yA_(1:2, 1:4, i) = yA_(1:2, 1:4, i)*ryA_(i)
		y_(1:4, i) = yA_(1, 1:4, i) + yA_(2, 1:4, i)
		if (pes_opt == 4) then
			VA_(1, i) = ei01_V_conv(1, m_Co(2*i), ryA_(i), afreq, yA_(1, 2, i))
			VA_(2, i) = ei01_V_conv(1, m_Co(2*i), ryA_(i), afreq, yA_(2, 2, i))
		else
			write(*,*) "err: write files, pes mismmatch"
			pause
		endif
	enddo
		y_norm(1) = dabs(y_(1, nyA-1))
		y_norm(2) = dabs(VA_(1,nyA-1)+VA_(2,nyA-1))
	do i = 0, nyA-1
		call Write5R8(72, FR8, ryA_(i), y_(1,i), y_(2,i), y_(3,i), y_(4,i))
		call Write2R8(76, FR8, ryA_(i), y_(1,i)/y_norm(1))
		call Write2R8(77, FR8, ryA_(i), (VA_(1,i)+VA_(2,i))/y_norm(2))
	enddo

	do i = 0, nzA-1	
		zA_(1:5, 1:6, i) = zA_(1:5, 1:6, i)*rzA_(i)
		z_(1:6, i) = zA_(1, 1:6, i) + zA_(2, 1:6, i)+ zA_(3, 1:6, i)+ zA_(4, 1:6, i)+ zA_(5, 1:6, i)
		call Write7R8(75, FR8, rzA_(i), z_(1,i), z_(2,i), z_(3,i), z_(4,i), z_(5,i), z_(6,i))
		call Write2R8(76, FR8, rzA_(i), z_(1,i)/y_norm(1))
		call Write2R8(77, FR8, rzA_(i), z_(3,i)/y_norm(2))
	enddo

	close(72)
	close(75)
	close(76)
	close(77)

endsubroutine ei01_o_soln_C2L2R5


subroutine ei01_o_soln_C1L2R2
use Format_IO
use FWrite
use FRead
implicit none
real(8), dimension(0:pg_N1+1) :: ryA_
real(8), dimension(0:pg_N2+1) :: ryB_
real(8), dimension(1:2, 1:4, 0:pg_N1+1) :: yA_
real(8), dimension(1:2, 0:pg_N1+1) :: VA_
real(8), dimension(1:4, 0:pg_N1+1) :: y_
real(8), dimension(1:2, 1:4, 0:pg_N2+1) :: yB_
real(8), dimension(1:2, 0:pg_N2+1) :: VB_
real(8), dimension(1:4, 0:pg_N2+1) :: y2_
real(8), dimension(1:2) :: y_norm
integer :: nyA, nyB, i
real(8) :: ans, norm, Q, Q_U, Q_V, Q_V_temp1, Q_V_temp2, Q_Alt

	call RFile5R8(pef_FN_Co1, nyA, ryA_, yA_(1, 1,0:pg_N1+1), yA_(1, 2,0:pg_N1+1), yA_(1, 3,0:pg_N1+1), yA_(1, 4,0:pg_N1+1))
	call RFile5R8(pef_FN_Co2, nyA, ryA_, yA_(2, 1,0:pg_N1+1), yA_(2, 2,0:pg_N1+1), yA_(2, 3,0:pg_N1+1), yA_(2, 4,0:pg_N1+1))

	call RFile5R8(pef_FN_Cr1, nyB, ryB_, yB_(1, 1,0:pg_N2+1), yB_(1, 2,0:pg_N2+1), yB_(1, 3,0:pg_N2+1), yB_(1, 4,0:pg_N2+1))
	call RFile5R8(pef_FN_Cr2, nyB, ryB_, yB_(2, 1,0:pg_N2+1), yB_(2, 2,0:pg_N2+1), yB_(2, 3,0:pg_N2+1), yB_(2, 4,0:pg_N2+1))

	if (nyA/= pg_N1+1 .or. nyB/= pg_N2+1) then
		write(*,*) "err: linear combination, number of puls eqt grid mismatch"
		pause
	endif

	open(72, file='data/sol_y_Co.txt', status= 'replace')
	open(75, file='data/sol_z_Cr.txt', status= 'replace')
	open(76, file='data/sol_U.txt', status= 'replace')
	open(77, file='data/sol_V.txt', status= 'replace')
	
	do i = 0, nyA-1
		yA_(1:2, 1:4, i) = yA_(1:2, 1:4, i)*ryA_(i)
		y_(1:4, i) = yA_(1, 1:4, i) + yA_(2, 1:4, i)
		if (pes_opt == 4) then
			VA_(1, i) = ei01_V_conv(1, m_Co(2*i), ryA_(i), afreq, yA_(1, 2, i))
			VA_(2, i) = ei01_V_conv(1, m_Co(2*i), ryA_(i), afreq, yA_(2, 2, i))
		else
			write(*,*) "err: write files, pes mismmatch"
			pause
		endif
	enddo

		y_norm(1) =dabs( y_(1, nyA-1))
		y_norm(2) =dabs( VA_(1,nyA-1)+VA_(2,nyA-1)	)
	
	allocate(rho_ar(0:nyA-1),dis_r(0:nyA-1), dis_t(0:nyA-1), r_ar(0:nyA-1) )
	do i = 0, nyA-1
		call Write5R8(72, FR8, ryA_(i), y_(1,i), y_(2,i), y_(3,i), y_(4,i))
		call Write2R8(76, FR8, ryA_(i), y_(1,i)/y_norm(1))
		call Write2R8(77, FR8, ryA_(i), (VA_(1,i)+VA_(2,i))/y_norm(2))

		rho_ar(i) = rho_Co(2*i)
		dis_r(i) = y_(1,i)
		dis_t(i) = VA_(1,i)+VA_(2,i)
		r_ar(i) = ryA_(i)
	enddo
	norm = 0.d0
	call ei01_integrate(ei01_norm, 0, nyA-1, 0.d0, ans)
	norm = norm + ans

	Q = 0.d0
	Q_U = 0.d0
	Q_V = 0.d0
	Q_Alt = 0.d0
	call ei01_integrate(ei01_Q, 0, nyA-1, 0.d0, ans)
	Q = Q + ans
	call ei01_integrate(ei01_Q_U, 0, nyA-1, 0.d0, ans)
	Q_U = Q_U + ans
	call ei01_integrate(ei01_Q_V, 0, nyA-1, 0.d0, ans)
	Q_V = Q_V + ans
Q_V_temp1 = Q_V
	call ei01_integrate(ei01_Q_alt, 1, nyA-2, 0.d0, ans)
	Q_Alt = Q_Alt + ans
	deallocate(rho_ar,dis_r, dis_t, r_ar )


	do i = 0, nyB-1	
		yB_(1:2, 1:4, i) = yB_(1:2, 1:4, i)*ryB_(i)
		y2_(1:4, i) = yB_(1, 1:4, i) + yB_(2, 1:4, i)
		if (pes_opt == 4) then
			VB_(1, i) = ei01_V_conv(1, m_Co(2*(nyB- 1-i)), ryB_(i), afreq, yB_(1, 2, i))
			VB_(2, i) = ei01_V_conv(1, m_Co(2*(nyB- 1-i)), ryB_(i), afreq, yB_(2, 2, i))
		else
			write(*,*) "err: write files, pes mismmatch"
			pause
		endif
		
	enddo
	
	allocate(rho_ar(0:nyB-1),dis_r(0:nyB-1), dis_t(0:nyB-1), r_ar(0:nyB-1) )
	do i = 0, nyB-1
		call Write5R8(75, FR8, ryB_(i), y2_(1,i), y2_(2,i), y2_(3,i), y2_(4,i))
		call Write2R8(76, FR8, ryB_(i), y2_(1,i)/y_norm(1))
		call Write2R8(77, FR8, ryA_(i), (VB_(1,i)+VB_(2,i))/y_norm(2))

		rho_ar(i) = rho_Cr(2*i)
		dis_r(i) = y2_(1,nyB- 1-i)
		dis_t(i) = VB_(1,nyB- 1-i)+VB_(2,nyB- 1-i)
		r_ar(i) = ryB_(nyB- 1-i)	
	enddo
	call ei01_integrate(ei01_norm, 0, nyB-1,0.d0, ans)
	norm = norm + ans
	norm = dsqrt(M0*R0**2/norm)

	call ei01_integrate(ei01_Q, 0, nyB-1, 0.d0, ans)
	Q = Q + ans
	call ei01_integrate(ei01_Q_U, 0, nyB-1, 0.d0, ans)
	Q_U = Q_U + ans
	call ei01_integrate(ei01_Q_V, 0, nyB-1, 0.d0, ans)
	Q_V = Q_V + ans
Q_V_temp2 = Q_V - Q_V_temp1
	call ei01_integrate(ei01_Q_alt, 1, nyB-2, 0.d0, ans)
	Q_Alt = Q_Alt + ans
	deallocate(rho_ar,dis_r, dis_t, r_ar )


write(*,*) "normalization constant = ", norm
write(*,*) "MR^2 ", M0*R0**2
write(*,*) "Q factor = ", (Q)*norm/(M0*R0**2)
write(*,*) "=========================================================="

	close(72)
	close(75)
	close(76)
	close(77)

endsubroutine ei01_o_soln_C1L2R2

subroutine ei01_o_soln_C2L2R2
use Format_IO
use FWrite
use FRead
implicit none
real(8), dimension(0:pg_N1+1) :: ryA_
real(8), dimension(0:pg_N2+1) :: ryB_
real(8), dimension(1:2, 1:4, 0:pg_N1+1) :: yA_
real(8), dimension(1:2, 0:pg_N1+1) :: VA_
real(8), dimension(1:4, 0:pg_N1+1) :: y_
real(8), dimension(1:4, 1:4, 0:pg_N2+1) :: yB_
real(8), dimension(1:4, 0:pg_N2+1) :: y2_
real(8), dimension(1:2) :: y_norm
integer :: nyA, nyB, i

	call RFile5R8(pef_FN_Co1, nyA, ryA_, yA_(1, 1,0:pg_N1+1), yA_(1, 2,0:pg_N1+1), yA_(1, 3,0:pg_N1+1), yA_(1, 4,0:pg_N1+1))
	call RFile5R8(pef_FN_Co2, nyA, ryA_, yA_(2, 1,0:pg_N1+1), yA_(2, 2,0:pg_N1+1), yA_(2, 3,0:pg_N1+1), yA_(2, 4,0:pg_N1+1))

	call RFile5R8(pef_FN_Cr1, nyB, ryB_, yB_(1, 1,0:pg_N2+1), yB_(1, 2,0:pg_N2+1), yB_(1, 3,0:pg_N2+1), yB_(1, 4,0:pg_N2+1))
	call RFile5R8(pef_FN_Cr2, nyB, ryB_, yB_(2, 1,0:pg_N2+1), yB_(2, 2,0:pg_N2+1), yB_(2, 3,0:pg_N2+1), yB_(2, 4,0:pg_N2+1))

	if (nyA/= pg_N1+1 .or. nyB/= pg_N2+1) then
		write(*,*) "err: linear combination, number of puls eqt grid mismatch"
		pause
	endif

	open(72, file='data/sol_y_Co.txt', status= 'replace')
	open(75, file='data/sol_z_Cr.txt', status= 'replace')
	open(76, file='data/sol_U.txt', status= 'replace')
	open(77, file='data/sol_V.txt', status= 'replace')
	
	do i = 0, nyA-1
		yA_(1:2, 1:4, i) = yA_(1:2, 1:4, i)*ryA_(i)
		y_(1:4, i) = yA_(1, 1:4, i) + yA_(2, 1:4, i)
		if (pes_opt == 4) then
			VA_(1, i) = ei01_V_conv(1, m_Co(2*i), ryA_(i), afreq, yA_(1, 2, i))
			VA_(2, i) = ei01_V_conv(1, m_Co(2*i), ryA_(i), afreq, yA_(2, 2, i))
		else
			write(*,*) "err: write files, pes mismmatch"
			pause
		endif
	enddo
		y_norm(1) = dabs( y_(1, nyA-1))
		y_norm(2) = dabs(VA_(1,nyA-1)+VA_(2,nyA-1)	)
	do i = 0, nyA-1
		call Write5R8(72, FR8, ryA_(i), y_(1,i), y_(2,i), y_(3,i), y_(4,i))
		call Write2R8(76, FR8, ryA_(i), y_(1,i)/y_norm(1))
		call Write2R8(77, FR8, ryA_(i), (VA_(1,i)+VA_(2,i))/y_norm(2))
	enddo

	do i = 0, nyB-1	
		yB_(1:4, 1:4, i) = yB_(1:4, 1:4, i)*ryB_(i)
		y2_(1:4, i) = yB_(1, 1:4, i) + yB_(2, 1:4, i)+ yB_(3, 1:4, i)+ yB_(4, 1:4, i)
		call Write5R8(75, FR8, ryB_(i), y2_(1,i), y2_(2,i), y2_(3,i), y2_(4,i))
		call Write2R8(76, FR8, ryB_(i), y2_(1,i)/y_norm(1))
		call Write2R8(77, FR8, ryB_(i), y2_(3,i)/y_norm(2))
	enddo

	close(72)
	close(75)
	close(76)
	close(77)

endsubroutine ei01_o_soln_C2L2R2

function ei01_V_conv(mode, m_, r_, sigma, x_, P_, rho_, Nu_)
implicit none
real(8) :: m_, r_, sigma, x_, ei01_V_conv
real(8), optional :: P_, rho_, Nu_
integer :: mode
	if (mode == 1) then
		ei01_V_conv = x_*Grav_Const*m_/r_**3/sigma**2
	elseif (mode == 2) then
		ei01_V_conv = x_*Grav_Const*(m_ + 4.d0 * pi * r_**3 * P_ / c**2)/r_**3/(1.d0 - 2.d0 * Grav_Const*m_/r_/c**2)/sigma**2	! same as mode 3 but without redshift in frequency
	elseif (mode == 3) then
		if (present(P_) /= .true. .or. present(rho_) /= .true. .or. present(Nu_) /= .true.) then
			write(*,*) "err: ei01_V_conv"
			pause
		endif
		ei01_V_conv = x_*Grav_Const*(m_ + 4.d0 * pi * r_**3 * P_ / c**2)/r_**3/(1.d0 - 2.d0 * Grav_Const*m_/r_/c**2)/sigma**2 * dexp(2.d0*Nu_)
	endif
endfunction ei01_V_conv

subroutine ei01_integrate(fcn, i0_, N_, yi_, ans)
! fcn and r_ar in ascending order
implicit none
real(8) :: fcn, sum, yi_, ans
integer :: i0_, N_, i
external fcn
	sum = yi_

		do i = i0_, N_-1
			sum = sum + 0.5d0 * (r_ar(i+1) - r_ar(i))*(fcn(i+1) + fcn(i))
		enddo
	ans = sum
endsubroutine ei01_integrate

subroutine ei01_double_integrate(fcn1, fcn2, i0_, N_, yi1_, yi2_, ans)
! fcn and r_ar in ascending order
! Calculate the form Int[ fcn1 * Int(fcn2) ]
implicit none
real(8) :: fcn1, fcn2, sum1, yi1_, sum2, yi2_, ans
integer :: i0_, N_, i, j
external fcn1, fcn2
	sum1 = yi1_
	do i = i0_, N_-1
		sum2 = yi2_
		do j = i0_, i
			sum2 = sum2 + 0.5d0 * (r_ar(i+1) - r_ar(i))*(fcn2(i+1) + fcn2(i))
		enddo
		sum1 = sum1 + 0.5d0 * (r_ar(i+1) - r_ar(i))*(fcn1(i+1) + fcn1(i))*sum2
	enddo
	ans = sum1
endsubroutine ei01_double_integrate

function ei01_norm(i)
implicit none
real(8) :: ei01_norm
integer :: i
	ei01_norm = rho_ar(i)*r_ar(i)**2* (dis_r(i)**2 + l_0*(l_0 + 1.d0)*dis_t(i)**2 )
endfunction ei01_norm

function ei01_Q(i)
implicit none
real(8) :: ei01_Q
integer :: i
	ei01_Q = rho_ar(i)*l_0*r_ar(i)**(l_0+1.d0)* (dis_r(i) + (l_0 + 1.d0)*dis_t(i) )
endfunction ei01_Q

function ei01_Q_U(i)
implicit none
real(8) :: ei01_Q_U
integer :: i
	ei01_Q_U = rho_ar(i)*l_0*r_ar(i)**(l_0+1.d0)* (dis_r(i) )
endfunction ei01_Q_U

function ei01_Q_V(i)
implicit none
real(8) :: ei01_Q_V
integer :: i
	ei01_Q_V = rho_ar(i)*l_0*r_ar(i)**(l_0+1.d0)* ((l_0 + 1.d0)*dis_t(i) )
endfunction ei01_Q_V

function ei01_Div_dis(i)	! do not include the boundary points in the integration
implicit none
real(8) :: ei01_Div_dis, ddis_r
integer :: i

	if (i == 0) then
		ddis_r = ( dis_r(i+1) - dis_r(i) )/( r_ar(i+1) - r_ar(i) )
	elseif (i == ubound(r_ar, 1)) then
		ddis_r = ( dis_r(i) - dis_r(i-1) )/( r_ar(i) - r_ar(i-1) )
	elseif (i >0 .and. i< ubound(r_ar, 1)) then
		ddis_r = 0.5d0 * ( ( dis_r(i+1) - dis_r(i) )/( r_ar(i+1) - r_ar(i) ) + ( dis_r(i) - dis_r(i-1) )/( r_ar(i) - r_ar(i-1) ) )
	else
		write(*,*) "err: ei01_Div_dis"
		pause
	endif
	ei01_Div_dis = ddis_r + 2.d0/r_ar(i)*dis_r(i) - l_0*(l_0+1.d0)/r_ar(i)*dis_t(i)
endfunction ei01_Div_dis

function ei01_Q_alt(i)
implicit none
real(8) :: ei01_Q_alt, drho
integer :: i
	if (i == 0) then
		drho = ( rho_ar(i+1) - rho_ar(i) )/( r_ar(i+1) - r_ar(i) )
	elseif (i == ubound(r_ar, 1)) then
		drho = ( rho_ar(i) - rho_ar(i-1) )/( r_ar(i) - r_ar(i-1) )
	elseif (i >0 .and. i< ubound(r_ar, 1)) then
		drho = 0.5d0 * ( ( rho_ar(i+1) - rho_ar(i) )/( r_ar(i+1) - r_ar(i) ) + ( rho_ar(i) - rho_ar(i-1) )/( r_ar(i) - r_ar(i-1) ) )
	else
		write(*,*) "err: ei01_Q_alt"
		pause
	endif
	ei01_Q_alt = -r_ar(i) **(l_0+ 2.d0) * (rho_ar(i) * ei01_Div_dis(i) + dis_r(i) * drho)
endfunction ei01_Q_alt

function ei01_fcn_rho_U(i)
implicit none
real(8) :: ei01_fcn_rho_U
integer :: i
	ei01_fcn_rho_U = rho_ar(i) * dis_r(i)
endfunction ei01_fcn_rho_U

function ei01_fcn_Double1(i)
implicit none
real(8) :: ei01_fcn_Double1
integer :: i
	ei01_fcn_Double1 = l_0 * (l_0+1.d0)/r_ar(i)**2
endfunction ei01_fcn_Double1

function ei01_fcn_Double2(i)
implicit none
real(8) :: ei01_fcn_Double2
integer :: i
	ei01_fcn_Double2 = r_ar(i) * rho_ar(i) * dis_t(i)
endfunction ei01_fcn_Double2

function ei01_dA(i)
implicit none
real(8) :: ei01_dA
integer :: i
	if (i == 0) then
		ei01_dA = ( ei01_fcn_Double2(i+1) - ei01_fcn_Double2(i) )/( r_ar(i+1) - r_ar(i) )
	elseif (i == ubound(r_ar, 1)) then
		ei01_dA = ( ei01_fcn_Double2(i) - ei01_fcn_Double2(i-1) )/( r_ar(i) - r_ar(i-1) )
	elseif (i >0 .and. i< ubound(r_ar, 1)) then
		ei01_dA = 0.5d0 * ( ( ei01_fcn_Double2(i+1) - ei01_fcn_Double2(i) )/( r_ar(i+1) - r_ar(i) ) + ( ei01_fcn_Double2(i) - ei01_fcn_Double2(i-1) )/( r_ar(i) - r_ar(i-1) ) )
	else
		write(*,*) "err: ei01_dA"
		pause
	endif
endfunction ei01_dA

function ei01_rho_Xi(i)
implicit none
real(8) :: ei01_rho_Xi, drho
integer :: i
	if (i == 0) then
		drho = ( rho_ar(i+1) - rho_ar(i) )/( r_ar(i+1) - r_ar(i) )
	elseif (i == ubound(r_ar, 1)) then
		drho = ( rho_ar(i) - rho_ar(i-1) )/( r_ar(i) - r_ar(i-1) )
	elseif (i >0 .and. i< ubound(r_ar, 1)) then
		drho = 0.5d0 * ( ( rho_ar(i+1) - rho_ar(i) )/( r_ar(i+1) - r_ar(i) ) + ( rho_ar(i) - rho_ar(i-1) )/( r_ar(i) - r_ar(i-1) ) )
	else
		write(*,*) "err: ei01_rho_Xi"
		pause
	endif
	ei01_rho_Xi = (rho_ar(i) * ei01_Div_dis(i) + dis_r(i) * drho)
endfunction ei01_rho_Xi

function ei01_Greens_Inner(i)
implicit none
real(8) :: ei01_Greens_Inner
integer :: i
	ei01_Greens_Inner = r_ar(i)**(l_0 + 2.d0) * ei01_rho_Xi(i)
endfunction ei01_Greens_Inner

function ei01_Greens_Outer(i)
implicit none
real(8) :: ei01_Greens_Outer
integer :: i
	ei01_Greens_Outer = r_ar(i)**(-l_0 + 1.d0) * ei01_rho_Xi(i)
endfunction ei01_Greens_Outer

function ei01_Div_grad(i)	! do not include the boundary points in the integration
implicit none
real(8) :: ei01_Div_grad, dgrad_r
integer :: i

	if (i == 0) then
		dgrad_r = ( grad_r(i+1) - grad_r(i) )/( r_ar(i+1) - r_ar(i) )
	elseif (i == ubound(r_ar, 1)) then
		dgrad_r = ( grad_r(i) - grad_r(i-1) )/( r_ar(i) - r_ar(i-1) )
	elseif (i >0 .and. i< ubound(r_ar, 1)) then
		dgrad_r = 0.5d0 * ( ( grad_r(i+1) - grad_r(i) )/( r_ar(i+1) - r_ar(i) ) + ( grad_r(i) - grad_r(i-1) )/( r_ar(i) - r_ar(i-1) ) )
	else
		write(*,*) "err: ei01_Div_dis"
		pause
	endif
	ei01_Div_grad = dgrad_r + 2.d0/r_ar(i)*grad_r(i) - l_0*(l_0+1.d0)/r_ar(i)*grad_t(i)
endfunction ei01_Div_grad

!function trial(i)
!implicit none
!real(8) :: trial
!integer :: i
!	trial = r_ar(i)**2
!endfunction trial

endmodule eigen_freq_opt01_ex2
