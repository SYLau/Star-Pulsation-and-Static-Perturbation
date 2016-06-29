module eigen_freq_opt02_ex1
!	output options
!	C2L4R6 stands for 2 Component, 4 Left (integration forward), 6 Right (integration backwards)
use global_var
contains

subroutine ei02_o_soln_C2L4R6
use Format_IO
use FWrite
use FRead
implicit none
real(8), dimension(0:pg_N1+1) :: ryA_
real(8), dimension(0:pg_N2+1) :: rzA_
real(8), dimension(1:2, 1:4, 0:pg_N1+1) :: yA_R, yA_I
real(8), dimension(1:4, 0:pg_N1+1) ::  yA_full_R, yA_full_I
real(8), dimension(1:5, 1:6, 0:pg_N2+1) :: zA_R, zA_I
real(8), dimension(1:6, 0:pg_N2+1) :: zA_full_R, zA_full_I
real(8), dimension(0:pg_N1+1) :: UA_full_R, VA_full_R
real(8), dimension(0:pg_N2+1) :: UCr_full_R, VCr_full_R
real(8), dimension(1:7):: Coe_R, Coe_I
real(8) :: U_norm, V_norm
integer :: nyA, nzA, i
	!¡¸Read data from pe results
	call RFile1R82AR8(pef_LD_Co1,nyA, ryA_, yA_R(1, 1:4, 0:pg_N1+1), yA_I(1, 1:4,0:pg_N1+1), m=4)
	call RFile1R82AR8(pef_LD_Co2,nyA, ryA_, yA_R(2, 1:4, 0:pg_N1+1), yA_I(2, 1:4,0:pg_N1+1), m=4)
	call RFile1R82AR8(pef_LD_Cr1,nzA, rzA_, zA_R(1, 1:6,0:pg_N2+1), zA_I(1, 1:6,0:pg_N2+1), m=6)
	call RFile1R82AR8(pef_LD_Cr2,nzA, rzA_, zA_R(2, 1:6,0:pg_N2+1), zA_I(2, 1:6,0:pg_N2+1), m=6)
	call RFile1R82AR8(pef_LD_Cr3,nzA, rzA_, zA_R(3, 1:6,0:pg_N2+1), zA_I(3, 1:6,0:pg_N2+1), m=6)
	call RFile1R82AR8(pef_LD_Cr4,nzA, rzA_, zA_R(4, 1:6,0:pg_N2+1), zA_I(4, 1:6,0:pg_N2+1), m=6)
	call RFile1R82AR8(pef_LD_Cr5,nzA, rzA_, zA_R(5, 1:6,0:pg_N2+1), zA_I(5, 1:6,0:pg_N2+1), m=6)
	!¡¸Verify no of data points
	if (nyA/= pg_N1+1 .or. nzA/= pg_N2+1) then
		write(*,*) "err: linear combination, number of puls eqt grid mismatch"
		pause
	endif
	open(72, file='data/sol_y_Co.txt', status= 'replace')
	open(73, file='data/sol_z_Cr.txt', status= 'replace')
	open(74, file='data/sol_U.txt', status= 'replace')
	open(75, file='data/sol_V.txt', status= 'replace')
	open(76, file='data/sol_UV.txt', status= 'replace')

	Coe_R(1:7) = dreal(CoeC(1:7))
	Coe_I(1:7) = dimag(CoeC(1:7))
	
	!¡¸Write core part
	do i = 0, nyA-1
		!¡¸(CR + iCI)*(YR + iYI) = (CR*YR - CI*YI) + i(CI*YR + CR*YI)
		!cyA_R(1, 1:4, i) = yA_R(1, 1:4, i) * Coe_R(1) - yA_I(1, 1:4, i) * Coe_I(1)
		!cyA_R(2, 1:4, i) = yA_R(2, 1:4, i) * Coe_R(2) - yA_I(2, 1:4, i) * Coe_I(2)
		!cyA_I(1, 1:4, i) = yA_R(1, 1:4, i) * Coe_I(1) + yA_I(1, 1:4, i) * Coe_R(1)
		!cyA_I(2, 1:4, i) = yA_R(2, 1:4, i) * Coe_I(2) + yA_I(2, 1:4, i) * Coe_R(2)

		yA_full_R(1:4, i) = yA_R(1, 1:4, i) + yA_R(2, 1:4, i)
		yA_full_I(1:4, i) = yA_I(1, 1:4, i) + yA_I(2, 1:4, i)
UA_full_R(i) = ei02_U_conv(m_Co(2*i), ryA_(i), yA_full_R(3,i))
VA_full_R(i) = ei02_V_conv_fluid(P_Co(2*i), rho_Co(2*i), m_Co(2*i), nu_Co(2*i), ryA_(i), yA_full_R(1:4,i), yA_full_I(1:4,i))
	enddo
	U_norm = UA_full_R(nyA-1)
	V_norm = VA_full_R(nyA-1)
	do i = 0, nyA-1
		call Write1R82AR8(72, FR8,ryA_(i),yA_full_R(1:4, i), yA_full_I(1:4, i), m=4)
		call Write3R8(74, FR8, ryA_(i), UA_full_R(i)/U_norm, yA_full_R(3,i))
		call Write2R8(75, FR8, ryA_(i), VA_full_R(i)/U_norm)
		call Write3R8(76, FR8, ryA_(i), UA_full_R(i)/U_norm, VA_full_R(i)/U_norm)
	enddo

	do i = 0, nzA-1
		!czA_R(1, 1:6, i) = zA_R(1, 1:6, i) * Coe_R(3) - zA_I(1, 1:6, i) * Coe_I(3)
		!czA_R(2, 1:6, i) = zA_R(2, 1:6, i) * Coe_R(4) - zA_I(2, 1:6, i) * Coe_I(4)
		!czA_R(3, 1:6, i) = zA_R(3, 1:6, i) * Coe_R(5) - zA_I(3, 1:6, i) * Coe_I(5)
		!czA_R(4, 1:6, i) = zA_R(4, 1:6, i) * Coe_R(6) - zA_I(4, 1:6, i) * Coe_I(6)
		!czA_R(5, 1:6, i) = zA_R(5, 1:6, i) * Coe_R(7) - zA_I(5, 1:6, i) * Coe_I(7)
		!czA_I(1, 1:6, i) = zA_R(1, 1:6, i) * Coe_I(3) + zA_I(1, 1:6, i) * Coe_R(3)
		!czA_I(2, 1:6, i) = zA_R(2, 1:6, i) * Coe_I(4) + zA_I(2, 1:6, i) * Coe_R(4)
		!czA_I(3, 1:6, i) = zA_R(3, 1:6, i) * Coe_I(5) + zA_I(3, 1:6, i) * Coe_R(5)
		!czA_I(4, 1:6, i) = zA_R(4, 1:6, i) * Coe_I(6) + zA_I(4, 1:6, i) * Coe_R(6)
		!czA_I(5, 1:6, i) = zA_R(5, 1:6, i) * Coe_I(7) + zA_I(5, 1:6, i) * Coe_R(7)

		zA_full_R(1:6, i) = zA_R(1, 1:6, i) + zA_R(2, 1:6, i) + zA_R(3, 1:6, i) + zA_R(4, 1:6, i) + zA_R(5, 1:6, i)
		zA_full_I(1:6, i) = zA_I(1, 1:6, i) + zA_I(2, 1:6, i) + zA_I(3, 1:6, i) + zA_I(4, 1:6, i) + zA_I(5, 1:6, i)
UCr_full_R(i) = ei02_U_conv(m_Cr(2*(nzA-1-i)), rzA_(i), zA_full_R(4,i))
VCr_full_R(i) = -zA_full_R(5,i)*rzA_(i)/R0 ** l_0
		call Write1R82AR8(73, FR8,rzA_(i),zA_full_R(1:6, i), zA_full_I(1:6, i), m=6)
		call Write3R8(74, FR8, rzA_(i), UCr_full_R(i)/U_norm, zA_full_R(4,i))
		call Write2R8(75, FR8, rzA_(i), VCr_full_R(i)/U_norm)
		call Write3R8(76, FR8, rzA_(i), UCr_full_R(i)/U_norm, VCr_full_R(i)/U_norm)
	enddo

	close(72)
	close(73)
	close(74)
	close(75)
	close(76)
endsubroutine ei02_o_soln_C2L4R6

subroutine ei02_o_soln_C3L4R6
use Format_IO
use FWrite
use FRead
use puls_eqt_set_opt03_ex1
implicit none
real(8), dimension(0:pg_N1+1) :: ryA_
real(8), dimension(0:pg_N2+1) :: rzA_
real(8), dimension(0:pg_N3+1) :: ryB_
real(8), dimension(1:2, 1:4, 0:pg_N1+1) :: yA_R, yA_I
real(8), dimension(1:4, 0:pg_N1+1) ::  yA_full_R, yA_full_I
real(8), dimension(1:5, 1:6, 0:pg_N2+1) :: zA_R, zA_I
real(8), dimension(1:6, 0:pg_N2+1) :: zA_full_R, zA_full_I
real(8), dimension(1:3, 1:4, 0:pg_N3+1) :: yB_R, yB_I
real(8), dimension(1:4, 0:pg_N3+1) :: yB_full_R, yB_full_I
real(8), dimension(0:pg_N1+1) :: UA_full_R, VA_full_R
real(8), dimension(0:pg_N2+1) :: UCr_full_R, VCr_full_R
real(8), dimension(0:pg_N3+1) :: UB_full_R, VB_full_R
real(8), dimension(1:10):: Coe_R, Coe_I
real(8) :: U_norm, V_norm
integer :: nyA, nzA, nyB, i, j
complex(8) :: H0_R, H0_Rm, dH0_R, yLove_R, Love, E(1:2,1:4)
	!¡¸Read data from pe results
	call RFile1R82AR8(pef_LD_Co1,nyA, ryA_, yA_R(1, 1:4, 0:pg_N1+1), yA_I(1, 1:4,0:pg_N1+1), m=4)
	call RFile1R82AR8(pef_LD_Co2,nyA, ryA_, yA_R(2, 1:4, 0:pg_N1+1), yA_I(2, 1:4,0:pg_N1+1), m=4)
	call RFile1R82AR8(pef_LD_Cr1,nzA, rzA_, zA_R(1, 1:6,0:pg_N2+1), zA_I(1, 1:6,0:pg_N2+1), m=6)
	call RFile1R82AR8(pef_LD_Cr2,nzA, rzA_, zA_R(2, 1:6,0:pg_N2+1), zA_I(2, 1:6,0:pg_N2+1), m=6)
	call RFile1R82AR8(pef_LD_Cr3,nzA, rzA_, zA_R(3, 1:6,0:pg_N2+1), zA_I(3, 1:6,0:pg_N2+1), m=6)
	call RFile1R82AR8(pef_LD_Cr4,nzA, rzA_, zA_R(4, 1:6,0:pg_N2+1), zA_I(4, 1:6,0:pg_N2+1), m=6)
	call RFile1R82AR8(pef_LD_Cr5,nzA, rzA_, zA_R(5, 1:6,0:pg_N2+1), zA_I(5, 1:6,0:pg_N2+1), m=6)
	call RFile1R82AR8(pef_LD_Oc1,nyB, ryB_, yB_R(1, 1:4, 0:pg_N3+1), yB_I(1, 1:4,0:pg_N3+1), m=4)
	call RFile1R82AR8(pef_LD_Oc2,nyB, ryB_, yB_R(2, 1:4, 0:pg_N3+1), yB_I(2, 1:4,0:pg_N3+1), m=4)
	call RFile1R82AR8(pef_LD_Oc3,nyB, ryB_, yB_R(3, 1:4, 0:pg_N3+1), yB_I(3, 1:4,0:pg_N3+1), m=4)

	!¡¸Verify no of data points
	if (nyA/= pg_N1+1 .or. nzA/= pg_N2+1) then
		write(*,*) "err: linear combination, number of puls eqt grid mismatch"
		pause
	endif
	open(72, file='data/sol_y_Co.txt', status= 'replace')
	open(73, file='data/sol_z_Cr.txt', status= 'replace')
	open(74, file='data/sol_y_Oc.txt', status= 'replace')
	open(75, file='data/sol_U.txt', status= 'replace')
	open(76, file='data/sol_V.txt', status= 'replace')
	open(77, file='data/sol_UV.txt', status= 'replace')

	Coe_R(1:10) = dreal(CoeC(1:10))
	Coe_I(1:10) = dimag(CoeC(1:10))
	
	!¡¸Write core part
	do i = 0, nyA-1
		!¡¸(CR + iCI)*(YR + iYI) = (CR*YR - CI*YI) + i(CI*YR + CR*YI)
		!cyA_R(1, 1:4, i) = yA_R(1, 1:4, i) * Coe_R(1) - yA_I(1, 1:4, i) * Coe_I(1)
		!cyA_R(2, 1:4, i) = yA_R(2, 1:4, i) * Coe_R(2) - yA_I(2, 1:4, i) * Coe_I(2)
		!cyA_I(1, 1:4, i) = yA_R(1, 1:4, i) * Coe_I(1) + yA_I(1, 1:4, i) * Coe_R(1)
		!cyA_I(2, 1:4, i) = yA_R(2, 1:4, i) * Coe_I(2) + yA_I(2, 1:4, i) * Coe_R(2)

		yA_full_R(1:4, i) = yA_R(1, 1:4, i) + yA_R(2, 1:4, i)
		yA_full_I(1:4, i) = yA_I(1, 1:4, i) + yA_I(2, 1:4, i)
UA_full_R(i) = ei02_U_conv(m_Co(2*i), ryA_(i), yA_full_R(3,i))
VA_full_R(i) = ei02_V_conv_fluid(P_Co(2*i), rho_Co(2*i), m_Co(2*i), nu_Co(2*i), ryA_(i), yA_full_R(1:4,i), yA_full_I(1:4,i))
	enddo
	U_norm = UA_full_R(nyA-1)
	V_norm = VA_full_R(nyA-1)
	do i = 0, nyA-1
		call Write1R82AR8(72, FR8,ryA_(i),yA_full_R(1:4, i), yA_full_I(1:4, i), m=4)
		call Write3R8(75, FR8, ryA_(i), UA_full_R(i)/U_norm, yA_full_R(3,i))
		call Write2R8(76, FR8, ryA_(i), VA_full_R(i)/U_norm)
		call Write3R8(77, FR8, ryA_(i), UA_full_R(i)/U_norm, VA_full_R(i)/U_norm)
	enddo

	do i = 0, nzA-1
		!czA_R(1, 1:6, i) = zA_R(1, 1:6, i) * Coe_R(3) - zA_I(1, 1:6, i) * Coe_I(3)
		!czA_R(2, 1:6, i) = zA_R(2, 1:6, i) * Coe_R(4) - zA_I(2, 1:6, i) * Coe_I(4)
		!czA_R(3, 1:6, i) = zA_R(3, 1:6, i) * Coe_R(5) - zA_I(3, 1:6, i) * Coe_I(5)
		!czA_R(4, 1:6, i) = zA_R(4, 1:6, i) * Coe_R(6) - zA_I(4, 1:6, i) * Coe_I(6)
		!czA_R(5, 1:6, i) = zA_R(5, 1:6, i) * Coe_R(7) - zA_I(5, 1:6, i) * Coe_I(7)
		!czA_I(1, 1:6, i) = zA_R(1, 1:6, i) * Coe_I(3) + zA_I(1, 1:6, i) * Coe_R(3)
		!czA_I(2, 1:6, i) = zA_R(2, 1:6, i) * Coe_I(4) + zA_I(2, 1:6, i) * Coe_R(4)
		!czA_I(3, 1:6, i) = zA_R(3, 1:6, i) * Coe_I(5) + zA_I(3, 1:6, i) * Coe_R(5)
		!czA_I(4, 1:6, i) = zA_R(4, 1:6, i) * Coe_I(6) + zA_I(4, 1:6, i) * Coe_R(6)
		!czA_I(5, 1:6, i) = zA_R(5, 1:6, i) * Coe_I(7) + zA_I(5, 1:6, i) * Coe_R(7)

		zA_full_R(1:6, i) = zA_R(1, 1:6, i) + zA_R(2, 1:6, i) + zA_R(3, 1:6, i) + zA_R(4, 1:6, i) + zA_R(5, 1:6, i)
		zA_full_I(1:6, i) = zA_I(1, 1:6, i) + zA_I(2, 1:6, i) + zA_I(3, 1:6, i) + zA_I(4, 1:6, i) + zA_I(5, 1:6, i)
UCr_full_R(i) = ei02_U_conv(m_Cr(2*(nzA-1-i)), rzA_(i), zA_full_R(4,i))
VCr_full_R(i) = -zA_full_R(5,i)*rzA_(i)/R0 ** l_0
		call Write1R82AR8(73, FR8,rzA_(i),zA_full_R(1:6, i), zA_full_I(1:6, i), m=6)
		call Write3R8(75, FR8, rzA_(i), UCr_full_R(i)/U_norm, zA_full_R(4,i))
		call Write2R8(76, FR8, rzA_(i), VCr_full_R(i)/U_norm)
		call Write3R8(77, FR8, rzA_(i), UCr_full_R(i)/U_norm, VCr_full_R(i)/U_norm)
	enddo

	do i = 0, nyB -1
		yB_full_R(1:4, i) = yB_R(1, 1:4, i) + yB_R(2, 1:4, i) + yB_R(3, 1:4, i)
		yB_full_I(1:4, i) = yB_I(1, 1:4, i) + yB_I(2, 1:4, i) + yB_I(3, 1:4, i)
UB_full_R(i) = ei02_U_conv(m_Oc(2*(nyB-1-i)), ryB_(i), yB_full_R(3,i))
VB_full_R(i) = ei02_V_conv_fluid(P_Oc(2*(nyB-1-i)), rho_Oc(2*(nyB-1-i)), m_Oc(2*(nyB-1-i)), nu_Oc(2*(nyB-1-i)), ryB_(i), yB_full_R(1:4,i), yB_full_I(1:4,i))
		call Write1R82AR8(74, FR8,ryB_(i),yB_full_R(1:4, i), yB_full_I(1:4, i), m=4)
		call Write3R8(75, FR8, ryB_(i), UB_full_R(i)/U_norm, yB_full_R(4,i))
		call Write2R8(76, FR8, ryB_(i), VB_full_R(i)/U_norm)
		call Write3R8(77, FR8, ryB_(i), UB_full_R(i)/U_norm, VB_full_R(i)/U_norm)
	enddo

	!Solving for Love number
	!¡¸Solve algebraric equations from BCs to obtain H0 and K
	if (zero_freq == .true.) then

		j = nyB-1
		call pes03_E_Fluid(P_Oc(2*j), rho_Oc(2*j), m_Oc(2*j), nu_Oc(2*j), ryB_(0), E)
		H0_R = (0.d0, 0.d0)
		do i = 1,4
			H0_R = H0_R + E(1,i) * dcmplx(yB_full_R(i, j), yB_full_I(i, j)) 
		enddo
		j = nyB-1-40
		call pes03_E_Fluid(P_Oc(2*j), rho_Oc(2*j), m_Oc(2*j), nu_Oc(2*j), ryB_(40), E)
		H0_Rm = (0.d0, 0.d0)
		do i = 1,4
			H0_Rm = H0_Rm + E(1,i) * dcmplx(yB_full_R(i, j), yB_full_I(i, j))
		enddo
		
		dH0_R = (H0_R - H0_Rm)/(ryB_(0) - ryB_(40))
		yLove_R = dH0_R/H0_R*(ryB_(0))
	
		call ei02_Love(yLove_R, Love)

		write(*,*) "k2 = ", Love

	endif
	close(72)
	close(73)
	close(74)
	close(75)
	close(76)
	close(77)
endsubroutine ei02_o_soln_C3L4R6

subroutine ei02_o_soln_C1L4R4
use Format_IO
use FWrite
use FRead
use puls_eqt_set_opt03_ex1
implicit none
real(8), dimension(0:pg_N1+1) :: ryA_
real(8), dimension(0:pg_N2+1) :: ryFM_
real(8), dimension(1:2, 1:4, 0:pg_N1+1) :: yA_R, yA_I
real(8), dimension(1:4, 0:pg_N1+1) ::  yA_full_R, yA_full_I
real(8), dimension(1:3, 1:4, 0:pg_N2+1) :: yFM_R, yFM_I
real(8), dimension(1:4, 0:pg_N2+1) :: yFM_full_R, yFM_full_I
real(8), dimension(0:pg_N1+1) :: UA_full_R, VA_full_R
real(8), dimension(0:pg_N2+1) :: UFM_full_R, VFM_full_R
real(8), dimension(1:5):: Coe_R, Coe_I
real(8) :: U_norm, V_norm
integer :: nyA, nyFM, i, j
complex(8) :: H0_R, H0_Rm, dH0_R, yLove_R, Love, E(1:2,1:4)

	call RFile1R82AR8(pef_LD_Co1,nyA, ryA_, yA_R(1, 1:4, 0:pg_N1+1), yA_I(1, 1:4,0:pg_N1+1), m=4)
	call RFile1R82AR8(pef_LD_Co2,nyA, ryA_, yA_R(2, 1:4, 0:pg_N1+1), yA_I(2, 1:4,0:pg_N1+1), m=4)
	call RFile1R82AR8(pef_LD_Cr1,nyFM, ryFM_, yFM_R(1, 1:4,0:pg_N2+1), yFM_I(1, 1:4,0:pg_N2+1), m=4)
	call RFile1R82AR8(pef_LD_Cr2,nyFM, ryFM_, yFM_R(2, 1:4,0:pg_N2+1), yFM_I(2, 1:4,0:pg_N2+1), m=4)
	call RFile1R82AR8(pef_LD_Cr3,nyFM, ryFM_, yFM_R(3, 1:4,0:pg_N2+1), yFM_I(3, 1:4,0:pg_N2+1), m=4)
	!¡¸Verify no of data points
	if (nyA/= pg_N1+1 .or. nyFM/= pg_N2+1) then
		write(*,*) "err: linear combination, number of puls eqt grid mismatch"
		pause
	endif

	open(72, file='data/sol_y_Co.txt', status= 'replace')
	open(73, file='data/sol_y_Cr.txt', status= 'replace')
	open(74, file='data/sol_U.txt', status= 'replace')
	open(75, file='data/sol_V.txt', status= 'replace')
	open(76, file='data/sol_UV.txt', status= 'replace')

	Coe_R(1:5) = dreal(CoeC(1:5))
	Coe_I(1:5) = dimag(CoeC(1:5))

	!¡¸Write core part
	do i = 0, nyA-1
		!¡¸(CR + iCI)*(YR + iYI) = (CR*YR - CI*YI) + i(CI*YR + CR*YI)
		!cyA_R(1, 1:4, i) = yA_R(1, 1:4, i) * Coe_R(1) - yA_I(1, 1:4, i) * Coe_I(1)
		!cyA_R(2, 1:4, i) = yA_R(2, 1:4, i) * Coe_R(2) - yA_I(2, 1:4, i) * Coe_I(2)
		!cyA_I(1, 1:4, i) = yA_R(1, 1:4, i) * Coe_I(1) + yA_I(1, 1:4, i) * Coe_R(1)
		!cyA_I(2, 1:4, i) = yA_R(2, 1:4, i) * Coe_I(2) + yA_I(2, 1:4, i) * Coe_R(2)

		yA_full_R(1:4, i) = yA_R(1, 1:4, i) + yA_R(2, 1:4, i)
		yA_full_I(1:4, i) = yA_I(1, 1:4, i) + yA_I(2, 1:4, i)
UA_full_R(i) = ei02_U_conv(m_Co(2*i), ryA_(i), yA_full_R(3,i))
VA_full_R(i) = ei02_V_conv_fluid(P_Co(2*i), rho_Co(2*i), m_Co(2*i), nu_Co(2*i), ryA_(i), yA_full_R(1:4,i), yA_full_I(1:4,i))
	enddo
	U_norm = UA_full_R(nyA-1)
	V_norm = VA_full_R(nyA-1)
	do i = 0, nyA-1
		call Write1R82AR8(72, FR8,ryA_(i),yA_full_R(1:4, i), yA_full_I(1:4, i), m=4)
		call Write3R8(74, FR8, ryA_(i), UA_full_R(i)/U_norm, yA_full_R(3,i))
		call Write2R8(75, FR8, ryA_(i), VA_full_R(i)/U_norm)
		call Write3R8(76, FR8, ryA_(i), UA_full_R(i)/U_norm, VA_full_R(i)/U_norm)
	enddo

	!¡¸Write fluid crust part
	do i = 0, nyFM-1
		!cyFM_R(1, 1:4, i) = yFM_R(1, 1:4, i) * Coe_R(3) - yFM_I(1, 1:4, i) * Coe_I(3)
		!cyFM_R(2, 1:4, i) = yFM_R(2, 1:4, i) * Coe_R(4) - yFM_I(2, 1:4, i) * Coe_I(4)
		!cyFM_R(3, 1:4, i) = yFM_R(3, 1:4, i) * Coe_R(5) - yFM_I(3, 1:4, i) * Coe_I(5)
		!cyFM_I(1, 1:4, i) = yFM_R(1, 1:4, i) * Coe_I(3) + yFM_I(1, 1:4, i) * Coe_R(3)
		!cyFM_I(2, 1:4, i) = yFM_R(2, 1:4, i) * Coe_I(4) + yFM_I(2, 1:4, i) * Coe_R(4)
		!cyFM_I(3, 1:4, i) = yFM_R(3, 1:4, i) * Coe_I(5) + yFM_I(3, 1:4, i) * Coe_R(5)

		yFM_full_R(1:4, i) = yFM_R(1, 1:4, i) + yFM_R(2, 1:4, i) + yFM_R(3, 1:4, i)
		yFM_full_I(1:4, i) = yFM_I(1, 1:4, i) + yFM_I(2, 1:4, i) + yFM_I(3, 1:4, i)
UFM_full_R(i) = ei02_U_conv(m_Cr(2*(nyFM-1-i)), ryFM_(i), yFM_full_R(3,i))
VFM_full_R(i) = ei02_V_conv_fluid(P_Cr(2*(nyFM-1-i)), rho_Cr(2*(nyFM-1-i)), m_Cr(2*(nyFM-1-i)), nu_Cr(2*(nyFM-1-i)), ryFM_(i), yFM_full_R(1:4,i), yFM_full_I(1:4,i))

		call Write1R82AR8(73, FR8,ryFM_(i),yFM_full_R(1:4, i), yFM_full_I(1:4, i), m=4)
		call Write3R8(74, FR8, ryFM_(i), UFM_full_R(i)/U_norm, yFM_full_R(3,i))
		call Write2R8(75, FR8, ryFM_(i), VFM_full_R(i)/U_norm)
		call Write3R8(76, FR8, ryFM_(i), UFM_full_R(i)/U_norm, VFM_full_R(i)/U_norm)
	enddo

	!Solving for Love number
	!¡¸Solve algebraric equations from BCs to obtain H0 and K
	if (zero_freq == .true.) then

		j = nyFM-1
		call pes03_E_Fluid(P_Cr(2*j), rho_Cr(2*j), m_Cr(2*j), nu_Cr(2*j), ryFM_(0), E)
		H0_R = (0.d0, 0.d0)
		do i = 1,4
			H0_R = H0_R + E(1,i) * dcmplx(yFM_full_R(i, j), yFM_full_I(i, j)) 
		enddo
		j = nyFM-1-40
		call pes03_E_Fluid(P_Cr(2*j), rho_Cr(2*j), m_Cr(2*j), nu_Cr(2*j), ryFM_(40), E)
		H0_Rm = (0.d0, 0.d0)
		do i = 1,4
			H0_Rm = H0_Rm + E(1,i) * dcmplx(yFM_full_R(i, j), yFM_full_I(i, j))
		enddo
		
		dH0_R = (H0_R - H0_Rm)/(ryFM_(0) - ryFM_(40))
		yLove_R = dH0_R/H0_R*(ryFM_(0))

		call ei02_Love(yLove_R, Love)

		write(*,*) "k2 = ", Love
		
	endif

	close(72)
	close(73)
	close(74)
	close(75)
	close(76)
endsubroutine ei02_o_soln_C1L4R4

subroutine ei02_o_soln_C1L6R6
use Format_IO
use FWrite
use FRead
implicit none
real(8), dimension(0:pg_N1+1) :: ryA_
real(8), dimension(0:pg_N2+1) :: rzA_
real(8), dimension(1:3, 1:6, 0:pg_N1+1) :: yA_R, yA_I
real(8), dimension(1:6, 0:pg_N1+1) ::  yA_full_R, yA_full_I
real(8), dimension(1:5, 1:6, 0:pg_N2+1) :: zA_R, zA_I
real(8), dimension(1:6, 0:pg_N2+1) :: zA_full_R, zA_full_I
real(8), dimension(0:pg_N1+1) :: UA_full_R, VA_full_R
real(8), dimension(0:pg_N2+1) :: UCr_full_R, VCr_full_R
real(8), dimension(1:8):: Coe_R, Coe_I
real(8) :: U_norm, V_norm
integer :: nyA, nzA, i
	!¡¸Read data from pe results
	call RFile1R82AR8(pef_LD_Co1,nyA, ryA_, yA_R(1, 1:6, 0:pg_N1+1), yA_I(1, 1:6,0:pg_N1+1), m=6)
	call RFile1R82AR8(pef_LD_Co2,nyA, ryA_, yA_R(2, 1:6, 0:pg_N1+1), yA_I(2, 1:6,0:pg_N1+1), m=6)
	call RFile1R82AR8(pef_LD_Co3,nyA, ryA_, yA_R(3, 1:6, 0:pg_N1+1), yA_I(3, 1:6,0:pg_N1+1), m=6)
	call RFile1R82AR8(pef_LD_Cr1,nzA, rzA_, zA_R(1, 1:6,0:pg_N2+1), zA_I(1, 1:6,0:pg_N2+1), m=6)
	call RFile1R82AR8(pef_LD_Cr2,nzA, rzA_, zA_R(2, 1:6,0:pg_N2+1), zA_I(2, 1:6,0:pg_N2+1), m=6)
	call RFile1R82AR8(pef_LD_Cr3,nzA, rzA_, zA_R(3, 1:6,0:pg_N2+1), zA_I(3, 1:6,0:pg_N2+1), m=6)
	call RFile1R82AR8(pef_LD_Cr4,nzA, rzA_, zA_R(4, 1:6,0:pg_N2+1), zA_I(4, 1:6,0:pg_N2+1), m=6)
	call RFile1R82AR8(pef_LD_Cr5,nzA, rzA_, zA_R(5, 1:6,0:pg_N2+1), zA_I(5, 1:6,0:pg_N2+1), m=6)
	!¡¸Verify no of data points
	if (nyA/= pg_N1+1 .or. nzA/= pg_N2+1) then
		write(*,*) "err: linear combination, number of puls eqt grid mismatch"
		pause
	endif
	open(72, file='data/sol_y_Co.txt', status= 'replace')
	open(73, file='data/sol_z_Cr.txt', status= 'replace')
	open(74, file='data/sol_U.txt', status= 'replace')
	open(75, file='data/sol_V.txt', status= 'replace')
	open(76, file='data/sol_UV.txt', status= 'replace')

	Coe_R(1:8) = dreal(CoeC(1:8))
	Coe_I(1:8) = dimag(CoeC(1:8))

	!¡¸Write core part
	do i = 0, nyA-1
		!¡¸(CR + iCI)*(YR + iYI) = (CR*YR - CI*YI) + i(CI*YR + CR*YI)
		!cyA_R(1, 1:4, i) = yA_R(1, 1:4, i) * Coe_R(1) - yA_I(1, 1:4, i) * Coe_I(1)
		!cyA_R(2, 1:4, i) = yA_R(2, 1:4, i) * Coe_R(2) - yA_I(2, 1:4, i) * Coe_I(2)
		!cyA_I(1, 1:4, i) = yA_R(1, 1:4, i) * Coe_I(1) + yA_I(1, 1:4, i) * Coe_R(1)
		!cyA_I(2, 1:4, i) = yA_R(2, 1:4, i) * Coe_I(2) + yA_I(2, 1:4, i) * Coe_R(2)

		yA_full_R(1:6, i) = yA_R(1, 1:6, i) + yA_R(2, 1:6, i) + yA_R(3, 1:6, i)
		yA_full_I(1:6, i) = yA_I(1, 1:6, i) + yA_I(2, 1:6, i) + yA_I(3, 1:6, i)
UA_full_R(i) = ei02_U_conv(m_Co(2*i), ryA_(i), yA_full_R(4,i))
VA_full_R(i) = -yA_full_R(5,i)*ryA_(i)/R0 ** l_0
	enddo
	U_norm = UA_full_R(nyA-1)
	V_norm = VA_full_R(nyA-1)
	do i = 0, nyA-1
		call Write1R82AR8(72, FR8,ryA_(i),yA_full_R(1:6, i), yA_full_I(1:6, i), m=6)
		call Write3R8(74, FR8, ryA_(i), UA_full_R(i)/U_norm, yA_full_R(4,i))
		call Write2R8(75, FR8, ryA_(i), VA_full_R(i)/U_norm)
		call Write3R8(76, FR8, ryA_(i), UA_full_R(i)/U_norm, VA_full_R(i)/U_norm)
	enddo

	do i = 0, nzA-1
		!czA_R(1, 1:6, i) = zA_R(1, 1:6, i) * Coe_R(3) - zA_I(1, 1:6, i) * Coe_I(3)
		!czA_R(2, 1:6, i) = zA_R(2, 1:6, i) * Coe_R(4) - zA_I(2, 1:6, i) * Coe_I(4)
		!czA_R(3, 1:6, i) = zA_R(3, 1:6, i) * Coe_R(5) - zA_I(3, 1:6, i) * Coe_I(5)
		!czA_R(4, 1:6, i) = zA_R(4, 1:6, i) * Coe_R(6) - zA_I(4, 1:6, i) * Coe_I(6)
		!czA_R(5, 1:6, i) = zA_R(5, 1:6, i) * Coe_R(7) - zA_I(5, 1:6, i) * Coe_I(7)
		!czA_I(1, 1:6, i) = zA_R(1, 1:6, i) * Coe_I(3) + zA_I(1, 1:6, i) * Coe_R(3)
		!czA_I(2, 1:6, i) = zA_R(2, 1:6, i) * Coe_I(4) + zA_I(2, 1:6, i) * Coe_R(4)
		!czA_I(3, 1:6, i) = zA_R(3, 1:6, i) * Coe_I(5) + zA_I(3, 1:6, i) * Coe_R(5)
		!czA_I(4, 1:6, i) = zA_R(4, 1:6, i) * Coe_I(6) + zA_I(4, 1:6, i) * Coe_R(6)
		!czA_I(5, 1:6, i) = zA_R(5, 1:6, i) * Coe_I(7) + zA_I(5, 1:6, i) * Coe_R(7)

		zA_full_R(1:6, i) = zA_R(1, 1:6, i) + zA_R(2, 1:6, i) + zA_R(3, 1:6, i) + zA_R(4, 1:6, i) + zA_R(5, 1:6, i)
		zA_full_I(1:6, i) = zA_I(1, 1:6, i) + zA_I(2, 1:6, i) + zA_I(3, 1:6, i) + zA_I(4, 1:6, i) + zA_I(5, 1:6, i)
UCr_full_R(i) = ei02_U_conv(m_Cr(2*(nzA-1-i)), rzA_(i), zA_full_R(4,i))
VCr_full_R(i) = -zA_full_R(5,i)*rzA_(i)/R0 ** l_0
		call Write1R82AR8(73, FR8,rzA_(i),zA_full_R(1:6, i), zA_full_I(1:6, i), m=6)
		call Write3R8(74, FR8, rzA_(i), UCr_full_R(i)/U_norm, zA_full_R(4,i))
		call Write2R8(75, FR8, rzA_(i), VCr_full_R(i)/U_norm)
		call Write3R8(76, FR8, rzA_(i), UCr_full_R(i)/U_norm, VCr_full_R(i)/U_norm)
	enddo

	close(72)
	close(73)
	close(74)
	close(75)
	close(76)
endsubroutine ei02_o_soln_C1L6R6

subroutine ei02_o_soln_C1L6R6_2
use Format_IO
use FWrite
use FRead
implicit none
real(8), dimension(0:pg_N1+1) :: ryA_
real(8), dimension(0:pg_N2+1) :: rzA_
real(8), dimension(1:2, 1:6, 0:pg_N1+1) :: yA_R, yA_I
real(8), dimension(1:6, 0:pg_N1+1) ::  yA_full_R, yA_full_I
real(8), dimension(1:5, 1:6, 0:pg_N2+1) :: zA_R, zA_I
real(8), dimension(1:6, 0:pg_N2+1) :: zA_full_R, zA_full_I
real(8), dimension(0:pg_N1+1) :: UA_full_R, VA_full_R
real(8), dimension(0:pg_N2+1) :: UCr_full_R, VCr_full_R
real(8), dimension(1:7):: Coe_R, Coe_I
real(8) :: U_norm, V_norm
integer :: nyA, nzA, i
	!¡¸Read data from pe results
	call RFile1R82AR8(pef_LD_Co1,nyA, ryA_, yA_R(1, 1:6, 0:pg_N1+1), yA_I(1, 1:6,0:pg_N1+1), m=6)
	call RFile1R82AR8(pef_LD_Co2,nyA, ryA_, yA_R(2, 1:6, 0:pg_N1+1), yA_I(2, 1:6,0:pg_N1+1), m=6)
	call RFile1R82AR8(pef_LD_Cr1,nzA, rzA_, zA_R(1, 1:6,0:pg_N2+1), zA_I(1, 1:6,0:pg_N2+1), m=6)
	call RFile1R82AR8(pef_LD_Cr2,nzA, rzA_, zA_R(2, 1:6,0:pg_N2+1), zA_I(2, 1:6,0:pg_N2+1), m=6)
	call RFile1R82AR8(pef_LD_Cr3,nzA, rzA_, zA_R(3, 1:6,0:pg_N2+1), zA_I(3, 1:6,0:pg_N2+1), m=6)
	call RFile1R82AR8(pef_LD_Cr4,nzA, rzA_, zA_R(4, 1:6,0:pg_N2+1), zA_I(4, 1:6,0:pg_N2+1), m=6)
	call RFile1R82AR8(pef_LD_Cr5,nzA, rzA_, zA_R(5, 1:6,0:pg_N2+1), zA_I(5, 1:6,0:pg_N2+1), m=6)
	!¡¸Verify no of data points
	if (nyA/= pg_N1+1 .or. nzA/= pg_N2+1) then
		write(*,*) "err: linear combination, number of puls eqt grid mismatch"
		pause
	endif
	open(72, file='data/sol_y_Co.txt', status= 'replace')
	open(73, file='data/sol_z_Cr.txt', status= 'replace')
	open(74, file='data/sol_U.txt', status= 'replace')
	open(75, file='data/sol_V.txt', status= 'replace')
	open(76, file='data/sol_UV.txt', status= 'replace')

	Coe_R(1:7) = dreal(CoeC(1:7))
	Coe_I(1:7) = dimag(CoeC(1:7))

	!¡¸Write core part
	do i = 0, nyA-1
		!¡¸(CR + iCI)*(YR + iYI) = (CR*YR - CI*YI) + i(CI*YR + CR*YI)
		!cyA_R(1, 1:4, i) = yA_R(1, 1:4, i) * Coe_R(1) - yA_I(1, 1:4, i) * Coe_I(1)
		!cyA_R(2, 1:4, i) = yA_R(2, 1:4, i) * Coe_R(2) - yA_I(2, 1:4, i) * Coe_I(2)
		!cyA_I(1, 1:4, i) = yA_R(1, 1:4, i) * Coe_I(1) + yA_I(1, 1:4, i) * Coe_R(1)
		!cyA_I(2, 1:4, i) = yA_R(2, 1:4, i) * Coe_I(2) + yA_I(2, 1:4, i) * Coe_R(2)

		yA_full_R(1:6, i) = yA_R(1, 1:6, i) + yA_R(2, 1:6, i)
		yA_full_I(1:6, i) = yA_I(1, 1:6, i) + yA_I(2, 1:6, i)
UA_full_R(i) = ei02_U_conv(m_Co(2*i), ryA_(i), yA_full_R(4,i))
VA_full_R(i) = -yA_full_R(5,i)*ryA_(i)/R0 ** l_0
	enddo
	U_norm = UA_full_R(nyA-1)
	V_norm = VA_full_R(nyA-1)
	do i = 0, nyA-1
		call Write1R82AR8(72, FR8,ryA_(i),yA_full_R(1:6, i), yA_full_I(1:6, i), m=6)
		call Write3R8(74, FR8, ryA_(i), UA_full_R(i)/U_norm, yA_full_R(4,i))
		call Write2R8(75, FR8, ryA_(i), VA_full_R(i)/U_norm)
		call Write3R8(76, FR8, ryA_(i), UA_full_R(i)/U_norm, VA_full_R(i)/U_norm)
	enddo

	do i = 0, nzA-1
		!czA_R(1, 1:6, i) = zA_R(1, 1:6, i) * Coe_R(3) - zA_I(1, 1:6, i) * Coe_I(3)
		!czA_R(2, 1:6, i) = zA_R(2, 1:6, i) * Coe_R(4) - zA_I(2, 1:6, i) * Coe_I(4)
		!czA_R(3, 1:6, i) = zA_R(3, 1:6, i) * Coe_R(5) - zA_I(3, 1:6, i) * Coe_I(5)
		!czA_R(4, 1:6, i) = zA_R(4, 1:6, i) * Coe_R(6) - zA_I(4, 1:6, i) * Coe_I(6)
		!czA_R(5, 1:6, i) = zA_R(5, 1:6, i) * Coe_R(7) - zA_I(5, 1:6, i) * Coe_I(7)
		!czA_I(1, 1:6, i) = zA_R(1, 1:6, i) * Coe_I(3) + zA_I(1, 1:6, i) * Coe_R(3)
		!czA_I(2, 1:6, i) = zA_R(2, 1:6, i) * Coe_I(4) + zA_I(2, 1:6, i) * Coe_R(4)
		!czA_I(3, 1:6, i) = zA_R(3, 1:6, i) * Coe_I(5) + zA_I(3, 1:6, i) * Coe_R(5)
		!czA_I(4, 1:6, i) = zA_R(4, 1:6, i) * Coe_I(6) + zA_I(4, 1:6, i) * Coe_R(6)
		!czA_I(5, 1:6, i) = zA_R(5, 1:6, i) * Coe_I(7) + zA_I(5, 1:6, i) * Coe_R(7)

		zA_full_R(1:6, i) = zA_R(1, 1:6, i) + zA_R(2, 1:6, i) + zA_R(3, 1:6, i) + zA_R(4, 1:6, i) + zA_R(5, 1:6, i)
		zA_full_I(1:6, i) = zA_I(1, 1:6, i) + zA_I(2, 1:6, i) + zA_I(3, 1:6, i) + zA_I(4, 1:6, i) + zA_I(5, 1:6, i)
UCr_full_R(i) = ei02_U_conv(m_Cr(2*(nzA-1-i)), rzA_(i), zA_full_R(4,i))
VCr_full_R(i) = -zA_full_R(5,i)*rzA_(i)/R0 ** l_0
		call Write1R82AR8(73, FR8,rzA_(i),zA_full_R(1:6, i), zA_full_I(1:6, i), m=6)
		call Write3R8(74, FR8, rzA_(i), UCr_full_R(i)/U_norm, zA_full_R(4,i))
		call Write2R8(75, FR8, rzA_(i), VCr_full_R(i)/U_norm)
		call Write3R8(76, FR8, rzA_(i), UCr_full_R(i)/U_norm, VCr_full_R(i)/U_norm)
	enddo

	close(72)
	close(73)
	close(74)
	close(75)
	close(76)
endsubroutine ei02_o_soln_C1L6R6_2

function ei02_U_conv(m_, r_, W_)
implicit none
real(8) :: m_, r_, W_, ei02_U_conv
		ei02_U_conv = W_ * r_**(l_0 - 1.d0) * (1.d0 - rel * 2.d0*m_/r_)**(0.5d0)/R0 ** l_0
endfunction ei02_U_conv

function ei02_V_conv_fluid(P_, rho_, m_, nu_, r_, y_R, y_I)
use puls_eqt_set_opt03_ex1
implicit none
complex(8) :: y_(1:4), E(1:2, 1:4), V_
real(8) :: P_, rho_, m_, nu_, r_, y_R(1:4), y_I(1:4)
real(8) :: ei02_V_conv_fluid
integer :: i_
	call pes03_E_Fluid(P_, rho_, m_, nu_, r_, E)
	y_(1:4) = dcmplx(y_R(1:4), y_I(1:4))
	V_ = (0.d0, 0.d0)
	do i_ = 1,4
		V_ = V_ + E(2,i_) * y_(i_) 
	enddo
	ei02_V_conv_fluid = dreal(-r_**(l_0-2.d0) *V_) * r_ /R0 ** l_0 ! multiplying r_ since tangential displacement defined to be dimensionless in LD2 (diff from Mc Dermott)

endfunction ei02_V_conv_fluid

!Love number---------------------------------------------------
subroutine ei02_Love(yR_, k2_)
implicit none
complex(8) :: yR_, T1_, T2a_, T2b_, k2_
real(8) :: m_, r_, C_
	m_ = M0
	r_ = R0
	C_ = m_/r_

	T1_ = 8.d0/5.d0*C_**5*(1.d0 - 2.d0*C_)**2 * (2.d0 * C_ * (yR_-1.d0) + 2.d0 - yR_)
	T2a_ = 2.d0 * C_ * (4.d0 * (yR_ + 1.d0)*C_**4 + (6.d0 * yR_ - 4.d0) * C_**3 + (26.d0 - 22.d0 * yR_) * C_**2 + 3.d0 * C_ * (5.d0 * yR_ - 8.d0) - 3.d0 * yR_ + 6.d0)
	T2b_ = 3.d0 * (1.d0 - 2.d0*C_)**2 * (2.d0 * C_ * (yR_ - 1.d0) - yR_ + 2.d0) * dlog(1.d0-2.d0 * C_)

	k2_ = T1_/(T2a_+ T2b_)

endsubroutine ei02_Love



endmodule eigen_freq_opt02_ex1
