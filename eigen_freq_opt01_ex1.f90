module eigen_freq_opt01_ex1
!	output options
!	normalized with crust core interface amplitude
!	C1L1R1 stands for 1 Component, 1 Left (integration forward), 1 Right (integration backwards)
use global_var

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

	do i = 0, nyA-1
		yA_(1:2, i) = yA_(1:2, i)*ryA_(i)
		if (pes_opt == 1 .and. newt_V_opt == 1) then
			VA_(i) = ei01_V_conv(1, m_Co(2*i), ryA_(i), afreq, yA_(2, i))
		elseif (pes_opt == 1 .and. newt_V_opt == 2) then
			VA_(i) = ei01_V_conv(2, m_Co(2*i), ryA_(i), afreq, yA_(2, i), P_ = P_Co(2*i), rho_ = rho_Co(2*i))
		elseif (pes_opt == 2) then
			VA_(i) = ei01_V_conv(3, m_Co(2*i), ryA_(i), afreq, yA_(2, i), P_ = P_Co(2*i), rho_ = rho_Co(2*i), Nu_ = Nu_Co(2*i))
		else
			write(*,*) "err: write files, pes mismmatch"
			pause
		endif
	enddo
		y_norm(1) = yA_(1, nyA-1)
		y_norm(2) = VA_(nyA-1)

	do i = 0, nyA-1
		call Write3R8(72, FR8, ryA_(i), yA_(1,i), yA_(2,i))
		call Write2R8(74, FR8, ryA_(i), yA_(1,i)/y_norm(1))
		call Write2R8(75, FR8, ryA_(i), VA_(i)/y_norm(2))
		call Write2R8(76, FR8, ryA_(i), rho_Co(2*i))
		
	enddo

	do i = 0, nyFM-1
		yFM_(1:2, i) = yFM_(1:2, i)*ryFM_(i)
		if (pes_opt == 1 .and. newt_V_opt == 1) then
			VFM_(i) = ei01_V_conv(1, m_Cr(2*(nyFM- 1-i)), ryFM_(i), afreq, yFM_(2, i))
		elseif (pes_opt == 1 .and. newt_V_opt == 2) then
			VFM_(i) = ei01_V_conv(2, m_Cr(2*(nyFM- 1-i)), ryFM_(i), afreq, yFM_(2, i), P_ = P_Cr(2*(nyFM- 1-i)), rho_ = rho_Cr(2*(nyFM- 1-i)))
		elseif (pes_opt == 2) then
			VFM_(i) = ei01_V_conv(3, m_Cr(2*(nyFM- 1-i)), ryFM_(i), afreq, yFM_(2, i), P_ = P_Cr(2*(nyFM- 1-i)), rho_ = rho_Cr(2*(nyFM- 1-i)), Nu_ = Nu_Cr(2*(nyFM- 1-i)))
		else
			write(*,*) "err: write files, pes mismmatch"
			pause
		endif

	enddo

	do i = 0, nyFM-1
		call Write3R8(73, FR8, ryFM_(i), yFM_(1,i), yFM_(2,i))
		call Write2R8(74, FR8, ryFM_(i), yFM_(1,i)/y_norm(1))
		call Write2R8(75, FR8, ryFM_(i), VFM_(i)/y_norm(2))
		call Write2R8(76, FR8, ryFM_(i), rho_Cr(2*(nyFM- 1-i)) )
	enddo


write(*,*) "Radial Displacement/ Tangential displacement = ", yA_(1,nyA -1)/VA_(nyA -1)
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


	close(72)
	close(73)
	close(74)
	close(75)
	close(76)
endsubroutine ei01_o_soln_C1L1R1

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
		if (pes_opt == 1) then
			VA_(i) = ei01_V_conv(1, m_Co(2*i), ryA_(i), afreq, yA_(2, i))
		elseif (pes_opt == 2) then
			VA_(i) = ei01_V_conv(3, m_Co(2*i), ryA_(i), afreq, yA_(2, i), P_ = P_Co(2*i), rho_ = rho_Co(2*i), Nu_ = Nu_Co(2*i))
		else
			write(*,*) "err: write files, pes mismmatch"
			pause
		endif
	enddo
		y_norm(1) = yA_(1, nyA-1)
		y_norm(2) = VA_(nyA-1)

	do i = 0, nyA-1
		call Write3R8(72, FR8, ryA_(i), yA_(1,i), yA_(2,i))
		call Write2R8(76, FR8, ryA_(i), yA_(1,i)/y_norm(1))
		call Write2R8(77, FR8, ryA_(i), VA_(i)/y_norm(2))
!call Write2R8(76, FR8, ryA_(i), yA_(1,i))
!call Write2R8(77, FR8, ryA_(i), VA_(i))
	enddo

	do i = 0, nzC-1	
		zC_(1:4, i) = zC_(1:4, i)*rzC_(i)
		zD_(1:4, i) = zD_(1:4, i)*rzD_(i)
		z_(1:4, i) = zC_(1:4, i) + zD_(1:4, i)
	enddo

	do i = 0, nzC-1	
		call Write5R8(73, FR8, rzC_(i), zC_(1,i), zC_(2,i), zC_(3,i), zC_(4,i))
		call Write5R8(74, FR8, rzD_(i), zD_(1,i), zD_(2,i), zD_(3,i), zD_(4,i))
		call Write5R8(75, FR8, rzC_(i), z_(1,i), z_(2,i), z_(3,i), z_(4,i))
		call Write2R8(76, FR8, rzC_(i), z_(1,i)/y_norm(1))
		call Write2R8(77, FR8, rzC_(i), z_(3,i)/y_norm(2))
!call Write2R8(76, FR8, rzC_(i), z_(1,i))
!call Write2R8(77, FR8, rzC_(i), z_(3,i))
	enddo

write(*,*) "Radial Displacement/ Tangential displacement = ", z_(3,nzC-1)/VA_(nyA -1)
write(*,*) "Tangential displacement jump = ", (VA_(nyA -1) - z_(3,nzC-1))/ VA_(nyA -1)
write(*,*) "Fluid transverse displacement jump , 0th order =, M1 ", Grav_Const * m_Co(2* pg_N1) * yA_(1, nyA-1) /VA_(nyA -1) / afreq**2 / r_Co(2* pg_N1) ** 3 * (rho_Co(2* pg_N1) - rho_Cr(0))/ rho_Co(2* pg_N1)	! NCA rho g
write(*,*) "Fluid transverse displacement jump , 0th order, M2 = ", Grav_Const * m_Co(2* pg_N1) * yA_(1, nyA-1) /VA_(nyA -1) / afreq**2 / r_Co(2* pg_N1) ** 3 * (rho_Co(2* pg_N1) - rho_Cr(0))/ rho_Co(2* pg_N1) /(1.d0 - 2.d0 * Grav_Const * m_Co(2*pg_N1)/c**2/r_Co(2*pg_N1)) ! NCA dP/dr
write(*,*) "Fluid transverse displacement jump , 0th order, M3 = ", Grav_Const * m_Co(2* pg_N1) * yA_(1, nyA-1) /VA_(nyA -1) / afreq**2 / r_Co(2* pg_N1) ** 3 * (rho_Co(2* pg_N1) - rho_Cr(0))/ rho_Co(2* pg_N1) /(1.d0 - 2.d0 * Grav_Const * m_Co(2*pg_N1)/c**2/r_Co(2*pg_N1)) / (dexp(-2.d0* nu_Co(2*pg_N1))) ! RCA
write(*,*) "1st order correction = ", - z_(3,nzC-1)/VA_(nyA -1) * rho_Cr(0)/rho_Co(2* pg_N1)*(rho_Co(2* pg_N1) - rho_Cr(0))/rho_Co(2* pg_N1)
write(*,*) "2nd order correction = ", - z_(3,nzC-1)/VA_(nyA -1) * rho_Cr(0)/rho_Co(2* pg_N1)*((rho_Co(2* pg_N1) - rho_Cr(0))/rho_Co(2* pg_N1))**2
write(*,*) "=========================================================="

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
		if (pes_opt == 1) then
			VA_(i) = ei01_V_conv(1, m_Co(2*i), ryA_(i), afreq, yA_(2, i))
		elseif (pes_opt == 2) then
			VA_(i) = ei01_V_conv(3, m_Co(2*i), ryA_(i), afreq, yA_(2, i), P_ = P_Co(2*i), rho_ = rho_Co(2*i), Nu_ = Nu_Co(2*i))
		else
			write(*,*) "err: write files, pes mismmatch"
			pause
		endif
	enddo
		y_norm(1) = yA_(1, nyA-1)
		y_norm(2) = VA_(nyA-1)
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
		if (pes_opt == 1) then
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
		y_norm(1) = y_(1, nyA-1)
		y_norm(2) = VA_(1,nyA-1)+VA_(2,nyA-1)
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
		y_norm(1) = y_(1, nyA-1)
		y_norm(2) = VA_(1,nyA-1)+VA_(2,nyA-1)	
	do i = 0, nyA-1
		call Write5R8(72, FR8, ryA_(i), y_(1,i), y_(2,i), y_(3,i), y_(4,i))
		call Write2R8(76, FR8, ryA_(i), y_(1,i)/y_norm(1))
		call Write2R8(77, FR8, ryA_(i), (VA_(1,i)+VA_(2,i))/y_norm(2))
	enddo

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
	do i = 0, nyB-1
		call Write5R8(75, FR8, ryB_(i), y2_(1,i), y2_(2,i), y2_(3,i), y2_(4,i))
		call Write2R8(76, FR8, ryB_(i), y2_(1,i)/y_norm(1))
		call Write2R8(77, FR8, ryA_(i), (VB_(1,i)+VB_(2,i))/y_norm(2))
	enddo

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
		y_norm(1) = y_(1, nyA-1)
		y_norm(2) = VA_(1,nyA-1)+VA_(2,nyA-1)	
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

endmodule eigen_freq_opt01_ex1
