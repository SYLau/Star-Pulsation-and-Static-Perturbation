module love_opt01_P08
! Fluid Part also integrates the matter field equations
! Not Yet developed completely
use global_var

contains
     subroutine lo01_iterate_fluid_W
     use Format_IO
     use FWrite
     use RK4_Set
     use lo_match_P02
     use lo_eqt
	 implicit none
	real(8), dimension(1:8) :: x, y, z
	real(8), allocatable :: r_out(:)
	real(8), allocatable:: y_out(:,:,:)
	real(8) :: Coef(1:2)
	real(8), dimension(1:5, 1:2) :: cY
	real(8) :: r_, r0_, dr_, Pc_, rhoc_, mc_, nuc_, k2_
	integer :: i, it, mode 
	
		allocate(r_out(1:sp_N2), y_out(1:2,1:5,1:sp_N2))
		do it = 1,2
			mode = it
			call lo01_bc_Coi_W(mode, y(1:8))
			y_out(mode, 1:5, 1) = y(1:5)
			r_out(1) = sp_r(1)

			do i = 1, sp_N1-1
				r_ = sp_r(i)
				dr_ = delr(i)
				x(1:8) = y(1:8)
				call RK4(lo_fluid_W_eqt_TOV_dr, 8, r_, r_ + dr_, x(1:8), y(1:8))
				r_out(i+1) = r_ + dr_
				y_out(mode, 1:5, i+1) = y(1:5)
			enddo

Pc_ = y(6)
rhoc_ = dsign(lo_EOS(dabs(y(6)), f = 'rho(p)'), y(6))
if (Crf(2) < 0.d0) pause 'err: rho_ < 0'
r0_ = sp_r(sp_N1)

			do i = sp_N1, sp_N2-1
				r_ = sp_r(i)
				dr_ = delr(i)
				x(1:8) = y(1:8)
				call RK4(lo_fluid_W_eqt_TOV_dr, 8, r_, r_ + dr_, x(1:8), y(1:8))
				r_out(i+1) = r_ + dr_
				y_out(mode, 1:5, i+1) = y(1:5)
			enddo
			cY(1:5,it) = y(1:5)
		enddo


				! Transfer to lo_bc_Solve
				Crf(1) = y(6)
				Crf(2) = dsign(lo_EOS(dabs(y(6)), f = 'rho(p)'), y(6))
				if (Crf(2) < 0.d0) pause 'err: rho_ < 0'
				Crf(3) = y(7)
				Crf(4) = lo_mu(y(6))
				Crf(5) = y(8)
				Crf(6) = sp_r(sp_N2)
		call lo_bc_Fluid_W_Solve(cY(1:5, 1:2), k2_, Coef(1:2))

		open(01, file = pef_LD_Co2, status = 'replace')
		open(02, file = 'data/temp_lo_y_fluid_W.txt', status = 'replace')
		do i = 1, sp_N2-1
			call Write1R81AR8(01, FR8,r_out(i), Coef(1)*y_out(1,1:5, i)+Coef(2)*y_out(2,1:5, i), m=5)
			call Write1R81AR8(02, FR8,r_out(i), Coef(1)*y_out(1,1:5, i)+Coef(2)*y_out(2,1:5, i), m=5)
		enddo
		close(01)
		close(02)

		write(*,*) "Fluid W model k2", k2_
		write(*,*) "Fluid Interface W/H", (Coef(1)*y_out(1,4, sp_N1)+Coef(2)*y_out(2,4, sp_N1))/(Coef(1)*y_out(1,1, sp_N1)+Coef(2)*y_out(2,1, sp_N1))
		write(*,*) "Fluid Interface W", (Coef(1)*y_out(1,4, sp_N1)+Coef(2)*y_out(2,4, sp_N1))
			write(*,*) "Fluid Interface r^2 PF", r0_**2/2.d0*(rhoc_+Pc_) *(Coef(1)*y_out(1,1, sp_N1)+Coef(2)*y_out(2,1, sp_N1))
		write(*,*) "Lambda_bar", 2.d0/3.d0* k2_ /(M0/R0)**5
		write(*,*) "=========================================================="
		deallocate(r_out, y_out)
		Cri(1:6) = 0.d0
		Crf(1:6) = 0.d0
     endsubroutine lo01_iterate_fluid_W

	!Integrate Fluid Core Solid Envelop 2Side---------------------------------------------------
	subroutine lo01_iterate_FE_2Side_W
	! FE - Fluid Core Elastic Crust
	! This method is Flawed. I used a pure fluid star to determine the W in the fluid core, and used this to determine the W at solid interface
	use Format_IO
	use FWrite
	use RK4_Set
	use lo_match_P02
	use lo_eqt
use love_opt01_P05
	implicit none
	real(8), dimension(1:9) :: x, y, z
	real(8), allocatable :: r_out(:)
	real(8), allocatable:: y_out(:,:,:)
	real(8), allocatable:: z_out(:,:,:)
	real(8) :: Coef(1:5)
	real(8) :: H0_Co, K_Co, W_Co
	real(8), dimension(1:6, 1:7) :: cY
	real(8) :: r_, r0_, dr_, Pc_,rhoc_,  mc_, nuc_, Pf_, mf_, nuf_, k2_
	integer :: i, it, mode 
	
	allocate(r_out(1:sp_N2), y_out(1:2,1:5,1:sp_N2), z_out(1:5, 1:6, sp_N1:sp_N2))

		open(01, file = pef_LD_Co1, status = 'replace')

		cY = 0.d0
		do it = 1,2
			mode = it
			call lo01_bc_Coi_W(mode, y(1:8))
			y_out(mode, 1:5, 1) = y(1:5)
			r_out(sp_N1) = sp_r(1)

			do i = 1, sp_N1-1
				r_ = sp_r(i)
				dr_ = delr(i)
				x(1:8) = y(1:8)
				call RK4(lo_fluid_W_eqt_TOV_dr, 8, r_, r_ + dr_, x(1:8), y(1:8))
				r_out(i+1) = r_ + dr_
				y_out(mode, 1:5, i+1) = y(1:5)
			enddo
			
			if (it == 1) then
				H0_Co = y(1)
				K_Co = y(3)
				Pc_ = y(6)
				rhoc_ = dsign(lo_EOS(dabs(y(6)), f = 'rho(p)'), y(6))
				if (Crf(2) < 0.d0) pause 'err: rho_ < 0'
				mc_ = y(7) 
				nuc_ = y(8)
				
				! Transfer to lo_bc_Solve
				Cri(1) = Pc_
				Cri(2) = dsign(lo_EOS(dabs(Pc_), f = 'rho(p)'), Pc_)
				if (Cri(2) < 0.d0) pause 'err: rho_ < 0'
				Cri(3) = mc_
				Cri(4) = lo_mu(Pc_)
				Cri(5) = nuc_
				Cri(6) = sp_r(sp_N1)
			endif

			do i = sp_N1, sp_N2-1
				r_ = sp_r(i)
				dr_ = delr(i)
				x(1:8) = y(1:8)
				call RK4(lo_fluid_W_eqt_TOV_dr, 8, r_, r_ + dr_, x(1:8), y(1:8))
				r_out(i+1) = r_ + dr_
				y_out(mode, 1:5, i+1) = y(1:5)
			enddo
			cY(1:5,it) = y(1:5)
		enddo

		Pf_ = y(6)
		mf_ = y(7) 
		nuf_ = y(8)
		r0_ = sp_r(sp_N2)
		! Transfer to lo_bc_Solve
		Crf(1) = Pf_
		Crf(2) = dsign(lo_EOS(dabs(Pf_), f = 'rho(p)'), Pf_)
		if (Pf_ < 0.d0) pause 'err: rho_ < 0'
		Crf(3) = mf_
		Crf(4) = lo_mu(Pf_)
		Crf(5) = nuf_
		Crf(6) = r0_

		call lo_bc_Fluid_W_Solve(cY(1:5, 1:2), k2_, Coef(1:2))
		W_Co = Coef(1) * y_out(1, 4, sp_N1) + Coef(2) * y_out(2, 4, sp_N1)

		close(01)


		region = 'crust'
		do it = 3, 7
			mode = it - 2
			call lo01_bc_Crf(mode, z(1:6))
			z(7) = Pf_
			z(8) = mf_
			z(9) = nuf_

			z_out(mode, 1:6, sp_N2) = z(1:6)
			r_out(sp_N2) = r0_

			do i = sp_N2, sp_N1+1, -1
				r_ = sp_r(i)
				dr_ = delr(i-1)
				x(1:9) = z(1:9)
				call RK4(lo_solid_eqt_TOV_dr, 9, r_, r_ - dr_, x(1:9), z(1:9))
				r_out(i-1) = r_ - dr_
				z_out(mode, 1:6, i-1) = z(1:6)
			enddo
			cY(1:6,it) = z(1:6)
		enddo

		call lo_bc_Solve_W_2Side(H0_Co, K_Co, W_Co, cY(1:6, 3:7), k2_, Coef(1:5))

		open(01, file = pef_LD_Cr1, status = 'replace')
!		do i = sp_N1, sp_N2-1
!			call Write1R81AR8(01, FR8,r_out(i), Coef(1)*z_out(1,1:6, i)+Coef(2)*z_out(2,1:6, i)+Coef(3)*z_out(3,1:6, i)+Coef(4)*z_out(4,1:6, i)+Coef(5)*z_out(5,1:6, i), m=6)
!		enddo
		close(01)

		write(*,*) "FE W model k2", k2_
!		write(*,*) "FE W model, interface W/H", vec_dot(Coef(1:5), z_out(1:5,4, sp_N1))/vec_dot(Coef(1:5), z_out(1:5,1, sp_N1))
		write(*,*) "Lambda_bar", 2.d0/3.d0* k2_ /(M0/R0)**5
		write(*,*) "=========================================================="
		deallocate(r_out, y_out, z_out)
		Cri(1:6) = 0.d0
		Crf(1:6) = 0.d0

	endsubroutine lo01_iterate_FE_2Side_W

! Integrate Solid Core Fluid Envelop---------------------------------------------------
	subroutine lo01_iterate_EF_W
	use Format_IO
	use FWrite
	use RK4_Set
	use love_opt01_P05
	use lo_match_P02
	use lo_eqt
use Root_Finding
use puls_eqt_set_opt03_ex1
	implicit none
	real(8), dimension(1:9) :: x, y, z
	real(8), allocatable :: r_out(:)
	real(8), allocatable:: z_out(:,:,:), y_out(:,:,:)
	real(8) :: Coef(1:8)
	real(8) :: H0F, dH0F, WF, KF
	real(8), dimension(1:6, 1:3) :: cZ
	real(8), dimension(1:5, 4:8) :: cY
	real(8) :: r_, dr_, Pc_, mc_, nuc_, Pf_, mf_, nuf_, k2_, rho_
	integer :: i, it, mode 
	
	allocate(r_out(1:sp_N2), z_out(1:3, 1:6, 1:sp_N1), y_out(4:8, 1:5, 1:sp_N2))

		cY = 0.d0
		do it = 1,3
			mode = it
			call lo01_bc_Coi_Solid(mode, z(1:9))
			r_ = sp_r(1)
			dr_ = delr(1)
			r_out(1) = r_
			z_out(mode, 1:6, 1) = z(1:6)
			do i = 1, sp_N1-1
				r_ = sp_r(i)
				dr_ = delr(i)
				x(1:9) = z(1:9)
				call RK4(lo_solid_eqt_TOV_dr, 9, r_, r_ + dr_, x(1:9), z(1:9))
				r_out(i+1) = r_ + dr_
				z_out(mode, 1:6, i+1) = z(1:6)
			enddo
			cZ(1:6,it) = z(1:6)
		enddo

		! Transfer to lo_bc_Solve
		Cri(1) = z(7)
		Cri(2) = dsign(lo_EOS(dabs(z(7)), f = 'rho(p)'), z(7))
		if (Cri(2) < 0.d0) pause 'err: rho_ < 0'
		Cri(3) = z(8)
		Cri(4) = lo_mu(z(7))
		Cri(5) = z(9)
		Cri(6) = sp_r(sp_N1)
		
		do i = sp_N1, sp_N2-1
			r_ = sp_r(i)
			dr_ = delr(i)
			x(1:9) = z(1:9)
			call RK4(lo_solid_eqt_TOV_dr, 9, r_, r_ + dr_, x(1:9), z(1:9))
		enddo

		Pf_ = z(7)
		mf_ = z(8) 
		nuf_ = z(9)
		! Transfer to lo_bc_Solve
		Crf(1) = Pf_
		Crf(2) = dsign(lo_EOS(dabs(Pf_), f = 'rho(p)'), Pf_)
		if (Pf_ < 0.d0) pause 'err: rho_ < 0'
		Crf(3) = mf_
		Crf(4) = lo_mu(Pf_)
		Crf(5) = nuf_
		Crf(6) = sp_r(sp_N2)

		do it = 4, 8
			mode = it -3
			call lo01_bc_fluid_Cr_f(mode, y(1:5))
			y(6) = Pf_
			y(7) = mf_
			y(8) = nuf_
			r_ = sp_r(sp_N2)
			dr_ = delr(sp_N2-1)
			r_out(sp_N2) = r_
			y_out(it, 1:5, sp_N2) = y(1:5)
			do i = sp_N2, sp_N1+1, -1
				r_ = sp_r(i)
				dr_ = delr(i-1)
				x(1:8) = y(1:8)
				call RK4(lo_fluid_W_eqt_TOV_dr, 8, r_, r_ - dr_, x(1:8), y(1:8))
				r_out(i-1) = r_ - dr_
				y_out(it, 1:5, i-1) = y(1:5)
			enddo
			cY(1:5,it) = y(1:5)
		enddo

		call lo_bc_EF_W_Solve(cZ(1:6,1:3), cY(1:5,4:8), k2_, Coef(1:8))

		write(*,*) "EF W k2", k2_
		write(*,*) "Lambda_bar", 2.d0/3.d0* k2_ /(M0/R0)**5
		write(*,*) "=========================================================="

open(02, file = 'data/temp_lo_y_EF_W.txt', status = 'replace')
do i = 1, sp_N1
write(02,FR8,advance = 'no') r_out(i)
write(02,FR8,advance = 'no') vec_dot(Coef(1:3) ,z_out(1:3,1,i))
write(02,FR8,advance = 'no') vec_dot(Coef(1:3) ,z_out(1:3,2,i))
write(02,FR8,advance = 'no') vec_dot(Coef(1:3) ,z_out(1:3,3,i))
write(02,FR8,advance = 'no') vec_dot(Coef(1:3) ,z_out(1:3,4,i))
write(02,FR8) vec_dot(Coef(1:3) ,z_out(1:3,5,i))
enddo
do i = sp_N1, sp_N2
write(02,FR8,advance = 'no') r_out(i)
write(02,FR8,advance = 'no') vec_dot(Coef(4:8) ,y_out(4:8,1,i))
write(02,FR8,advance = 'no') vec_dot(Coef(4:8) ,y_out(4:8,2,i))
write(02,FR8,advance = 'no') vec_dot(Coef(4:8) ,y_out(4:8,3,i))
write(02,FR8,advance = 'no') vec_dot(Coef(4:8) ,y_out(4:8,4,i))
write(02,FR8) vec_dot(Coef(4:8) ,y_out(4:8,5,i))
enddo
close(02)
!write(*,*) vec_dot(Coef(1:3) ,z_out(1:3,2,sp_N1)) - vec_dot(Coef(4:8) ,y_out(4:8,2,sp_N1))
!write(*,*) 32.d0*pi* pes03_dnu(Cri(1), Cri(2), Cri(3), Cri(6)) * 2.d0* Cri(4) * (vec_dot(Coef(1:3) ,z_out(1:3,5,sp_N1)))
		deallocate(r_out, z_out, y_out)
		Cri(1:6) = 0.d0
		Crf(1:6) = 0.d0

	endsubroutine lo01_iterate_EF_W


	! Integrate Solid Core Fluid Envelop---------------------------------------------------
!#### Bug
	subroutine lo01_iterate_EF_W_No_Slip
	use Format_IO
	use FWrite
	use RK4_Set
	use love_opt01_P05
	use lo_match_P02
	use lo_eqt
use Root_Finding
use puls_eqt_set_opt03_ex1
	implicit none
	real(8), dimension(1:9) :: x, y, z
	real(8), allocatable :: r_out(:)
	real(8), allocatable:: z_out(:,:,:), y_out(:,:,:)
	real(8) :: Coef(1:8)
	real(8) :: H0F, dH0F, WF, KF
	real(8), dimension(1:6, 1:3) :: cZ
	real(8), dimension(1:5, 4:8) :: cY
	real(8) :: r_, dr_, Pc_, mc_, nuc_, Pf_, mf_, nuf_, k2_, rho_
	integer :: i, it, mode 
	
	allocate(r_out(1:sp_N2), z_out(1:3, 1:6, 1:sp_N1), y_out(4:8, 1:5, 1:sp_N2))

		cY = 0.d0
		do it = 1,3
			mode = it
			call lo01_bc_Coi_Solid(mode, z(1:9))
			r_ = sp_r(1)
			dr_ = delr(1)
			r_out(1) = r_
			z_out(mode, 1:6, 1) = z(1:6)
			do i = 1, sp_N1-1
				r_ = sp_r(i)
				dr_ = delr(i)
				x(1:9) = z(1:9)
				call RK4(lo_solid_eqt_TOV_dr, 9, r_, r_ + dr_, x(1:9), z(1:9))
				r_out(i+1) = r_ + dr_
				z_out(mode, 1:6, i+1) = z(1:6)
			enddo
			cZ(1:6,it) = z(1:6)
		enddo

		! Transfer to lo_bc_Solve
		Cri(1) = z(7)
		Cri(2) = dsign(lo_EOS(dabs(z(7)), f = 'rho(p)'), z(7))
		if (Cri(2) < 0.d0) pause 'err: rho_ < 0'
		Cri(3) = z(8)
		Cri(4) = lo_mu(z(7))
		Cri(5) = z(9)
		Cri(6) = sp_r(sp_N1)
		
		do i = sp_N1, sp_N2-1
			r_ = sp_r(i)
			dr_ = delr(i)
			x(1:9) = z(1:9)
			call RK4(lo_solid_eqt_TOV_dr, 9, r_, r_ + dr_, x(1:9), z(1:9))
		enddo

		Pf_ = z(7)
		mf_ = z(8) 
		nuf_ = z(9)
		! Transfer to lo_bc_Solve
		Crf(1) = Pf_
		Crf(2) = dsign(lo_EOS(dabs(Pf_), f = 'rho(p)'), Pf_)
		if (Pf_ < 0.d0) pause 'err: rho_ < 0'
		Crf(3) = mf_
		Crf(4) = lo_mu(Pf_)
		Crf(5) = nuf_
		Crf(6) = sp_r(sp_N2)

		do it = 4, 8
			mode = it -3
			call lo01_bc_fluid_Cr_f(mode, y(1:5))
			y(6) = Pf_
			y(7) = mf_
			y(8) = nuf_
			r_ = sp_r(sp_N2)
			dr_ = delr(sp_N2-1)
			r_out(sp_N2) = r_
			y_out(it, 1:5, sp_N2) = y(1:5)
			do i = sp_N2, sp_N1+1, -1
				r_ = sp_r(i)
				dr_ = delr(i-1)
				x(1:8) = y(1:8)
				call RK4(lo_fluid_W_eqt_TOV_dr, 8, r_, r_ - dr_, x(1:8), y(1:8))
				r_out(i-1) = r_ - dr_
				y_out(it, 1:5, i-1) = y(1:5)
			enddo
			cY(1:5,it) = y(1:5)
		enddo

		call lo_bc_EF_W_Solve_No_Slip(cZ(1:6,1:3), cY(1:5,4:8), k2_, Coef(1:8))

		write(*,*) "EF W no slip model k2", k2_
		write(*,*) "Lambda_bar", 2.d0/3.d0* k2_ /(M0/R0)**5
		write(*,*) "=========================================================="

open(02, file = 'data/temp_lo_y_EF_No_Slip.txt', status = 'replace')
do i = 1, sp_N1
write(02,FR8,advance = 'no') r_out(i)
write(02,FR8,advance = 'no') vec_dot(Coef(1:3) ,z_out(1:3,1,i))
write(02,FR8,advance = 'no') vec_dot(Coef(1:3) ,z_out(1:3,2,i))
write(02,FR8,advance = 'no') vec_dot(Coef(1:3) ,z_out(1:3,3,i))
write(02,FR8,advance = 'no') vec_dot(Coef(1:3) ,z_out(1:3,4,i))
write(02,FR8) vec_dot(Coef(1:3) ,z_out(1:3,5,i))
enddo
do i = sp_N1, sp_N2
write(02,FR8,advance = 'no') r_out(i)
write(02,FR8,advance = 'no') vec_dot(Coef(4:8) ,y_out(4:8,1,i))
write(02,FR8,advance = 'no') vec_dot(Coef(4:8) ,y_out(4:8,2,i))
write(02,FR8,advance = 'no') vec_dot(Coef(4:8) ,y_out(4:8,3,i))
write(02,FR8,advance = 'no') vec_dot(Coef(4:8) ,y_out(4:8,4,i))
write(02,FR8) vec_dot(Coef(4:8) ,y_out(4:8,5,i))
enddo
close(02)
!write(*,*) vec_dot(Coef(1:3) ,z_out(1:3,2,sp_N1)) - vec_dot(Coef(4:8) ,y_out(4:8,2,sp_N1))
!write(*,*) 32.d0*pi* pes03_dnu(Cri(1), Cri(2), Cri(3), Cri(6)) * 2.d0* Cri(4) * (vec_dot(Coef(1:3) ,z_out(1:3,5,sp_N1)))
		deallocate(r_out, z_out, y_out)
		Cri(1:6) = 0.d0
		Crf(1:6) = 0.d0

	endsubroutine lo01_iterate_EF_W_No_Slip


	 subroutine lo01_bc_Coi_W(mode, y)
	 use puls_eqt_set_opt03_ex1
	use lo_eqt
	implicit none
	real(8) :: gamma_1, P_, rho_, m_, nu_, r_, dr_, a0_, a2_, V0
	real(8) :: y(1:8)
	integer :: mode
		P_ = P(1)
		rho_ = rho(1)
		m_ = m(1)
		nu_ = nu(1)
		r_ = sp_r(1)
		gamma_1 = lo_gamma(P_)
		
		if (mode == 1) then
			a0_ = 1.d0
			V0 = 0.d0
		elseif (mode == 2) then
			a0_ = 0.d0
			V0 = 1.d0
		endif

		a2_ = -2.d0/(2.d0*l_0+3.d0) * pi * (5.d0*rho_+9.d0*P_+(rho_+P_)/(gamma_1*P_/(rho_+P_)) )
		
		
		y(1) = a0_ * r_**l_0 * (1.d0 + r_ ** 2 *  a2_ )
		y(2) = a0_ * r_**(l_0-1.d0) * (l_0 + r_ ** 2 * (l_0+2.d0) * a2_  )
		y(3) = a0_ * r_**l_0  * ( 1.d0 +  ((l_0+2.d0)*a2_ + 8.d0*pi/3.d0*(3.d0*P_+rho_))/(l_0+2.d0) )   ! need second order expansion !!
		
		y(4) = l_0 * V0 
		y(5) = 	V0

		y(6) = P_
		y(7) = m_
		y(8) = nu_
	 endsubroutine lo01_bc_Coi_W

	subroutine lo01_bc_fluid_Cr_f(mode, y)
	implicit none
	real(8) :: y(1:5)
	integer :: mode
		y = 0.d0
		if (mode == 1) y(1) = 1.d0
		if (mode == 2) y(2) = 1.d0
		if (mode == 3) y(3) = 1.d0
		if (mode == 4) y(4) = 1.d0
		if (mode == 5) y(5) = 1.d0
	endsubroutine lo01_bc_fluid_Cr_f
endmodule love_opt01_P08