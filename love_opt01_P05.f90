module love_opt01_P05
! Coupled with TOV equation
! Use the background grid to reach lowest density
! To compare with Andersson 2011, with Polytropic EOS and analytic expression of Shear modulus
! Modified the algebraic expressions of Andersson 2011 (Since it is inconsistent with Andersson 2011 and Finn 1990)
use global_var
implicit none
contains

	!Integrate Fluid Equations---------------------------------------------------
	subroutine lo01_iterate_fluid
	use Format_IO
	use FWrite
	use RK4_Set
	use lo_match
	use lo_eqt
	implicit none
	real(8), dimension(1:6) :: x, y
	real(8) :: r_, dr_, k2_, rho_
	integer :: i
write(*,*) "Multiply the shear modulus by mu_factor"
mu_factor = 1.d0
!mu_factor = 1.d-5
write(*,*) "mu_factor = "
read(*,*) mu_factor
		open(01, file = pef_LD_Co1, status = 'replace')
		open(02, file = 'data/temp_lo_y_fluid.txt', status = 'replace')
		call lo01_bc_Coi(y(1:6))

		call Write2R8(02, FR8,sp_r(1)/sp_r(sp_N2), (sp_r(1))*y(2)/y(1))

		region = 'core'
		do i = 1, sp_N1-1
				r_ = sp_r(i)
				dr_ = delr(i)
				x(1:6) = y(1:6)
				call RK4(lo_fluid_eqt_TOV_dr, 6, r_, r_ + dr_, x(1:6), y(1:6))
				call Write1R81AR8(01, FR8,r_ + dr_, (y(1:3)), m=3)
				call Write2R8(02, FR8,(r_ + dr_)/sp_r(sp_N2), (r_ + dr_)*y(2)/y(1))
		enddo
		close(01)
		
		call Write2R8(02, FR8,sp_r(sp_N1)/sp_r(sp_N2), sp_r(sp_N1)*y(2)/y(1))

		open(01, file = pef_LD_Cr1, status = 'replace')
		region = 'crust'
		call Write1R81AR8(01, FR8,r_ + dr_, (y(1:3)), m=3)
		do i = sp_N1, sp_N2-1
				r_ = sp_r(i)
				dr_ = delr(i)
				x(1:6) = y(1:6)
				call RK4(lo_fluid_eqt_TOV_dr, 6, r_, r_ + dr_, x(1:6), y(1:6))
				call Write1R81AR8(01, FR8,r_ + dr_, (y(1:3)), m=3)
				call Write2R8(02, FR8,(r_ + dr_)/sp_r(sp_N2), (r_ + dr_)*y(2)/y(1))
		enddo
		call Write2R8(02, FR8,sp_r(sp_N2)/sp_r(sp_N2), sp_r(sp_N2)*y(2)/y(1) - 4.d0*pi*lo_EOS(dabs(y(4)))*sp_r(sp_N2)**3/y(5))
		close(01)
		close(02)
rho_ = dsign(lo_EOS(dabs(y(4)), f = 'rho(p)'), y(4))
if (rho_ < 0.d0) pause 'err: rho_ < 0'
call lo_k2(sp_r(sp_N2)*y(2)/y(1), 0.d0, y(4), rho_, y(5), 0.d0, sp_r(sp_N2), k2_)

		write(*,*) "pure fluid k2", k2_
		write(*,*) "Lambda_bar", 2.d0/3.d0* k2_ /(M0/R0)**5
		write(*,*) "=========================================================="

	endsubroutine lo01_iterate_fluid
	
	!Integrate Fluid Core Solid Envelop---------------------------------------------------
	subroutine lo01_iterate_FE
	! FE - Fluid Core Elastic Crust
	use Format_IO
	use FWrite
	use RK4_Set
	use lo_match
	use lo_eqt
	implicit none
	real(8), dimension(1:9) :: x, y, z
	real(8), allocatable :: r_out(:)
	real(8), allocatable:: z_out(:,:,:)
	real(8) :: Coef(1:5)
	real(8) :: H0_Co, K_Co
	real(8), dimension(1:6, 1:5) :: cY
	real(8) :: r_, r0_, dr_, Pc_, mc_, nuc_, k2_
	integer :: i, it, mode 
	
	allocate(r_out(sp_N1:sp_N2), z_out(1:5, 1:6, sp_N1:sp_N2))

		open(01, file = pef_LD_Co1, status = 'replace')
		cY = 0.d0
		call lo01_bc_Coi(y(1:6))
		region = 'core'
		do i = 1, sp_N1-1
				r_ = sp_r(i)
				dr_ = delr(i)
			x(1:6) = y(1:6)
			call RK4(lo_fluid_eqt_TOV_dr, 6, r_, r_ + dr_, x(1:6), y(1:6))
			call Write1R81AR8(01, FR8,r_ + dr_, (y(1:3)), m=3)
		enddo
		close(01)


		H0_Co = y(1)
		K_Co = y(3)

		Pc_ = y(4)
		mc_ = y(5) 
		nuc_ = y(6)
		r0_ = sp_r(sp_N1)

! Transfer to lo_bc_Solve
Cri(1) = Pc_
!Cri(2) = dsign((dabs(Pc_)/poly_K/1.d10)**(poly_n/(poly_n + 1.d0)), Pc_)
Cri(2) = dsign(lo_EOS(dabs(Pc_), f = 'rho(p)'), Pc_)
if (Cri(2) < 0.d0) pause 'err: rho_ < 0'
Cri(3) = mc_
!Cri(4) = (Ki_poly*Pc_ + mu_poly/1.d10)/2.d0  *mu_factor
Cri(4) = lo_mu(Pc_)
Cri(5) = nuc_
Cri(6) = r0_

		region = 'crust'
		do it = 2, 6
			mode = it - 1
			call lo01_bc_Cri(mode, z(1:6))
			z(7) = Pc_
			z(8) = mc_
			z(9) = nuc_
			z_out(mode, 1:6, sp_N1) = z(1:6)
			r_out(sp_N1) = r0_

			do i = sp_N1, sp_N2-1
				r_ = sp_r(i)
				dr_ = delr(i)
				x(1:9) = z(1:9)
				call RK4(lo_solid_eqt_TOV_dr, 9, r_, r_ + dr_, x(1:9), z(1:9))
				r_out(i+1) = r_ + dr_
				z_out(mode, 1:6, i+1) = z(1:6)
			enddo
			cY(1:6,it - 1) = z(1:6)
		enddo
! Transfer to lo_bc_Solve
Crf(1) = z(7)
!Crf(2) = dsign((dabs(x(7))/poly_K/1.d10)**(poly_n/(poly_n + 1.d0)), x(7)) !!####### x or z ???
Crf(2) = dsign(lo_EOS(dabs(z(7)), f = 'rho(p)'), z(7))
if (Crf(2) < 0.d0) pause 'err: rho_ < 0'
Crf(3) = z(8)
!Crf(4) = (Ki_poly*z(7) + mu_poly/1.d10)/2.d0 *mu_factor
Crf(4) = lo_mu(z(7))
Crf(5) = z(9)
Crf(6) = sp_r(sp_N2)
		call lo_bc_Solve(H0_Co, K_Co, cY(1:6, 1:5), k2_, Coef(1:5))


		open(01, file = pef_LD_Cr1, status = 'replace')
		do i = sp_N1, sp_N2-1
			call Write1R81AR8(01, FR8,r_out(i), Coef(1)*z_out(1,1:6, i)+Coef(2)*z_out(2,1:6, i)+Coef(3)*z_out(3,1:6, i)+Coef(4)*z_out(4,1:6, i)+Coef(5)*z_out(5,1:6, i), m=6)
		enddo
		close(01)

		write(*,*) "FE model k2", k2_
		write(*,*) "Lambda_bar", 2.d0/3.d0* k2_ /(M0/R0)**5
		write(*,*) "=========================================================="
		deallocate(r_out, z_out)
		Cri(1:6) = 0.d0
		Crf(1:6) = 0.d0

	endsubroutine lo01_iterate_FE

!Integrate Fluid Core Solid Envelop 2Side---------------------------------------------------
	subroutine lo01_iterate_FE2
	! FE - Fluid Core Elastic Crust
	use Format_IO
	use FWrite
	use RK4_Set
	use lo_match
	use lo_eqt
	use Root_Finding
	implicit none
	real(8), dimension(1:9) :: x, y, z
	real(8), allocatable :: r_out(:)
	real(8), allocatable:: z_out(:,:,:)
	real(8) :: Coef(1:5)
	real(8) :: H0_Co, K_Co
	real(8), dimension(1:6, 1:5) :: cY
	real(8) :: r_, r0_, dr_, Pc_, mc_, nuc_, Pf_, mf_, nuf_, k2_
	integer :: i, it, mode 
	
	allocate(r_out(sp_N1:sp_N2), z_out(1:5, 1:6, sp_N1:sp_N2))

		open(01, file = pef_LD_Co1, status = 'replace')
		cY = 0.d0
		call lo01_bc_Coi(y(1:6))
		region = 'core'
		do i = 1, sp_N1-1
				r_ = sp_r(i)
				dr_ = delr(i)
			x(1:6) = y(1:6)
			call RK4(lo_fluid_eqt_TOV_dr, 6, r_, r_ + dr_, x(1:6), y(1:6))
			call Write1R81AR8(01, FR8,r_ + dr_, (y(1:3)), m=3)
		enddo
		close(01)

		H0_Co = y(1)
		K_Co = y(3)
		Pc_ = y(4)
		mc_ = y(5) 
		nuc_ = y(6)
		r0_ = sp_r(sp_N1)

! Transfer to lo_bc_Solve
Cri(1) = Pc_
!Cri(2) = dsign((dabs(Pc_)/poly_K/1.d10)**(poly_n/(poly_n + 1.d0)), Pc_)
Cri(2) = dsign(lo_EOS(dabs(Pc_), f = 'rho(p)'), Pc_)
if (Cri(2) < 0.d0) pause 'err: rho_ < 0'
Cri(3) = mc_
!Cri(4) = (Ki_poly*Pc_ + mu_poly/1.d10)/2.d0  *mu_factor
Cri(4) = lo_mu(Pc_)
Cri(5) = nuc_
Cri(6) = r0_

		z(7) = Pc_
		z(8) = mc_
		z(9) = nuc_
! Solve TOV equation to obtain surface P, rho, m
			do i = sp_N1, sp_N2-1
				r_ = sp_r(i)
				dr_ = delr(i)
				x(1:9) = z(1:9)
				call RK4(lo_solid_eqt_TOV_dr, 9, r_, r_ + dr_, x(1:9), z(1:9))
			enddo

		Pf_ = z(7)
		mf_ = z(8)
		nuf_ = z(9) 
		r0_ = sp_r(sp_N2)

! Transfer to lo_bc_Solve
Crf(1) = Pf_
!Crf(2) = dsign((dabs(x(7))/poly_K/1.d10)**(poly_n/(poly_n + 1.d0)), x(7)) !!####### x or z ???
Crf(2) = dsign(lo_EOS(dabs(Pf_), f = 'rho(p)'), Pf_)
if (Pf_ < 0.d0) pause 'err: rho_ < 0'
Crf(3) = mf_
!Crf(4) = (Ki_poly*z(7) + mu_poly/1.d10)/2.d0 *mu_factor
Crf(4) = lo_mu(Pf_)
Crf(5) = nuf_
Crf(6) = r0_

		region = 'crust'
		do it = 2, 6
			mode = it - 1
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
			cY(1:6,it - 1) = z(1:6)
		enddo

		call lo_bc_Solve_2Side(H0_Co, K_Co, cY(1:6, 1:5), k2_, Coef(1:5))

		open(01, file = pef_LD_Cr1, status = 'replace')
		do i = sp_N1, sp_N2-1
			call Write1R81AR8(01, FR8,r_out(i), Coef(1)*z_out(1,1:6, i)+Coef(2)*z_out(2,1:6, i)+Coef(3)*z_out(3,1:6, i)+Coef(4)*z_out(4,1:6, i)+Coef(5)*z_out(5,1:6, i), m=6)
		enddo
		close(01)

		write(*,*) "FE model k2", k2_
		write(*,*) "FE model, interface W/H", vec_dot(Coef(1:5), z_out(1:5,4, sp_N1))/vec_dot(Coef(1:5), z_out(1:5,1, sp_N1))
		write(*,*) "Lambda_bar", 2.d0/3.d0* k2_ /(M0/R0)**5
		write(*,*) "=========================================================="
		deallocate(r_out, z_out)
		Cri(1:6) = 0.d0
		Crf(1:6) = 0.d0

	endsubroutine lo01_iterate_FE2

	!Integrate Solid Eqt---------------------------------------------------
	subroutine lo01_iterate_Solid
	! Solid Star
	use Format_IO
	use FWrite
	use RK4_Set
	use lo_match
	use lo_eqt
	implicit none
	real(8), dimension(1:9) :: x, y, z
	real(8), allocatable :: r_out(:)
	real(8), allocatable:: z_out(:,:,:)
	real(8) :: Coef(1:3)
	real(8) :: H0_Co, K_Co
	real(8), dimension(1:6, 1:3) :: cY
	real(8) :: r_, r0_, dr_, Pc_, mc_, nuc_, k2_
	integer :: i, it, mode 
	
	allocate(r_out(1:sp_N2), z_out(1:3, 1:6, 1:sp_N2))

		cY = 0.d0
		do it = 1,3
			mode = it
			call lo01_bc_Coi_Solid(mode, z(1:9))
			r_ = sp_r(1)
			dr_ = delr(1)
			r_out(1) = r_
			z_out(mode, 1:6, 1) = z(1:6)
			do i = 1, sp_N2-1
				r_ = sp_r(i)
				dr_ = delr(i)
				x(1:9) = z(1:9)
				call RK4(lo_solid_eqt_TOV_dr, 9, r_, r_ + dr_, x(1:9), z(1:9))
				r_out(i+1) = r_ + dr_
				z_out(mode, 1:6, i+1) = z(1:6)
			enddo
			cY(1:6,it) = z(1:6)
		enddo
			
! Transfer to lo_bc_Solve
Crf(1) = z(7)
Crf(2) = dsign(lo_EOS(dabs(z(7)), f = 'rho(p)'), z(7))
if (Crf(2) < 0.d0) pause 'err: rho_ < 0'
Crf(3) = z(8)
Crf(4) = lo_mu(z(7))
Crf(5) = z(9)
Crf(6) = sp_r(sp_N2)

		call lo_bc_Solid_Solve(cY(1:6,1:3), k2_, Coef(1:3))

		open(01, file = pef_LD_Cr2, status = 'replace')
		do i = 1, sp_N2-1
			call Write1R81AR8(01, FR8,r_out(i), Coef(1)*z_out(1,1:6, i)+Coef(2)*z_out(2,1:6, i)+Coef(3)*z_out(3,1:6, i), m=6)
		enddo
		close(01)

		write(*,*) "pure solid outward k2", k2_
		write(*,*) "Lambda_bar", 2.d0/3.d0* k2_ /(M0/R0)**5
		write(*,*) "##solution diverges near surface due to Low Density??"
		write(*,*) "=========================================================="
		deallocate(r_out, z_out)
		Cri(1:6) = 0.d0
		Crf(1:6) = 0.d0

	endsubroutine lo01_iterate_Solid


	!Integrate Solid Eqt---------------------------------------------------
	subroutine lo01_iterate_Solid2
	! Solid Star
	use Format_IO
	use FWrite
	use RK4_Set
	use lo_match
	use lo_eqt
	use Root_Finding
	implicit none
	real(8), dimension(1:9) :: x, y, z
	real(8), allocatable :: r_out(:)
	real(8), allocatable:: z_out_Co(:,:,:), z_out_Cr(:,:,:), z_full_Co(:,:), z_full_Cr(:,:)
	real(8) :: Coef(1:8)
	real(8), dimension(1:6, 1:8) :: cY
	real(8) :: r_, r0_, dr_, Pf_, mf_, nuf_, k2_
	integer :: i, it, mode 
	
	allocate(r_out(1:sp_N2), z_out_Co(1:3, 1:6, 1:sp_N1), z_out_Cr(4:8, 1:6, sp_N1:sp_N2), z_full_Co(1:6, 1:sp_N1), z_full_Cr(1:6, sp_N1:sp_N2))

		cY = 0.d0
		do it = 1,3
			mode = it
			call lo01_bc_Coi_Solid(mode, z(1:9))
			r_ = sp_r(1)
			dr_ = delr(1)
			r_out(1) = r_
			z_out_Co(mode, 1:6, 1) = z(1:6)
			do i = 1, sp_N1-1
				r_ = sp_r(i)
				dr_ = delr(i)
				x(1:9) = z(1:9)
				call RK4(lo_solid_eqt_TOV_dr, 9, r_, r_ + dr_, x(1:9), z(1:9))
				r_out(i+1) = r_ + dr_
				z_out_Co(mode, 1:6, i+1) = z(1:6)
			enddo
			! Transfer to lo_bc_Solve
				region = 'core'
				Cri(1) = z(7)
				Cri(2) = dsign(lo_EOS(dabs(z(7)), f = 'rho(p)'), z(7))
				if (Cri(2) < 0.d0) pause 'err: rho_ < 0'
				Cri(3) = z(8)
				Cri(4) = lo_mu(z(7))
				Cri(5) = z(9)
				Cri(6) = sp_r(sp_N1)

			cY(1:6,it) = z(1:6)
		enddo
			! Solve TOV equation to obtain surface P, rho, m
			do i = sp_N1, sp_N2-1
				r_ = sp_r(i)
				dr_ = delr(i)
				x(1:9) = z(1:9)
				call RK4(lo_solid_eqt_TOV_dr, 9, r_, r_ + dr_, x(1:9), z(1:9))
			enddo

		Pf_ = z(7)
		mf_ = z(8)
		nuf_ = z(9) 

		region = 'crust'
		Cri_out(1:6) = Cri(1:6)
		Cri_out(4) = lo_mu(Cri_out(1))

! Transfer to lo_bc_Solve
Crf(1) = z(7)
Crf(2) = dsign(lo_EOS(dabs(z(7)), f = 'rho(p)'), z(7))
if (Crf(2) < 0.d0) pause 'err: rho_ < 0'
Crf(3) = z(8)
Crf(4) = lo_mu(z(7))
Crf(5) = z(9)
Crf(6) = sp_r(sp_N2)

		do it = 4, 8
			mode = it -3
			call lo01_bc_Crf(mode, z(1:6))
			z(7) = Pf_
			z(8) = mf_
			z(9) = nuf_
			r_ = sp_r(sp_N2)
			dr_ = delr(sp_N2-1)
			r_out(sp_N2) = r_
			z_out_Cr(it, 1:6, sp_N2) = z(1:6)
			do i = sp_N2, sp_N1+1, -1
				r_ = sp_r(i)
				dr_ = delr(i-1)
				x(1:9) = z(1:9)
				call RK4(lo_solid_eqt_TOV_dr, 9, r_, r_ - dr_, x(1:9), z(1:9))
				r_out(i-1) = r_ - dr_
				z_out_Cr(it, 1:6, i-1) = z(1:6)
			enddo
			cY(1:6,it) = z(1:6)
		enddo

		!¡¸ able to due with 2 Solid models, with density gap/ shear modulus gap
			call lo_bc_Solid_Solve_2Side_muGap(cY(1:6,1:8), k2_, Coef(1:8))
		!¡¸ original Solid_Solve_2Side, only for 1 component pure solid
			!call lo_bc_Solid_Solve_2Side(cY(1:6,1:8), k2_, Coef(1:8))

		z_full_Co = 0.d0
		z_full_Cr = 0.d0
		do i = 1, sp_N1
			z_full_Co(1:6, i) = z_full_Co(1:6, i) + Coef(1)*z_out_Co(1,1:6, i)+ Coef(2)*z_out_Co(2,1:6, i)+ Coef(3)*z_out_Co(3,1:6, i)
		enddo
		do i = sp_N1,sp_N2
			z_full_Cr(1:6, i) = z_full_Cr(1:6, i) + Coef(4)*z_out_Cr(4,1:6, i)+ Coef(5)*z_out_Cr(5,1:6, i)+ Coef(6)*z_out_Cr(6,1:6, i)+ Coef(7)*z_out_Cr(7,1:6, i)+ Coef(8)*z_out_Cr(8,1:6, i)
		enddo

		open(01, file = pef_LD_Cr2, status = 'replace')
		open(02, file = 'data/temp_lo_y_2solid.txt', status = 'replace')
write(*,*) "Rc/R", sp_r(sp_N1)/sp_r(sp_N2)
write(*,*) "r*dH0/H0", sp_r(sp_N1)*(z_full_Co(2,sp_N1) - z_full_Cr(2,sp_N1))/z_full_Co(1,sp_N1)
		do i = 1, sp_N1
			call Write1R81AR8(01, FR8,r_out(i), z_full_Co(1:6, i) , m=6)
!			call Write2R8(02,FR8,r_out(i)/sp_r(sp_N2),r_out(i)*z_full_Co(2, i)/z_full_Co(1, i))
write(02,FR8,advance = 'no') r_out(i)
write(02,FR8,advance = 'no') z_full_Co(1,i)
write(02,FR8,advance = 'no') z_full_Co(2,i)
write(02,FR8,advance = 'no') z_full_Co(3,i)
write(02,FR8,advance = 'no') z_full_Co(4,i)
write(02,FR8,advance = 'no') z_full_Co(5,i)
write(02,FR8) z_full_Co(6,i)
		enddo

		do i = sp_N1, sp_N2
			call Write1R81AR8(01, FR8,r_out(i), z_full_Cr(1:6, i) , m=6)
!			call Write2R8(02,FR8,r_out(i)/sp_r(sp_N2),r_out(i)*z_full_Cr(2, i)/z_full_Cr(1, i))
write(02,FR8,advance = 'no') r_out(i)
write(02,FR8,advance = 'no') z_full_Cr(1,i)
write(02,FR8,advance = 'no') z_full_Cr(2,i)
write(02,FR8,advance = 'no') z_full_Cr(3,i)
write(02,FR8,advance = 'no') z_full_Cr(4,i)
write(02,FR8,advance = 'no') z_full_Cr(5,i)
write(02,FR8) z_full_Cr(6,i)
		enddo
!			call Write2R8(02,FR8,r_out(sp_N2)/sp_r(sp_N2),solR_temp)
		close(01)
		close(02)

		write(*,*) "pure solid 2 sides k2", k2_
		write(*,*) "Solid model, interface W/H", vec_dot(Coef(1:3), z_out_Co(1:3,4, sp_N1))/vec_dot(Coef(1:3), z_out_Cr(1:3,1, sp_N1))
		write(*,*) "Lambda_bar", 2.d0/3.d0* k2_ /(M0/R0)**5
		write(*,*) "=========================================================="
		deallocate(r_out, z_out_Co, z_out_Cr, z_full_Co, z_full_Cr)
		Cri(1:6) = 0.d0
		Crf(1:6) = 0.d0

	endsubroutine lo01_iterate_Solid2

	! Integrate Solid Core Fluid Envelop---------------------------------------------------
	subroutine lo01_iterate_EF
	use Format_IO
	use FWrite
	use RK4_Set
	use lo_match
	use lo_eqt
	implicit none
	real(8), dimension(1:9) :: x, y, z
	real(8), allocatable :: r_out(:)
	real(8), allocatable:: z_out(:,:,:), z_full(:,:)
	real(8) :: Coef(1:3)
	real(8) :: H0F, dH0F, KF
	real(8), dimension(1:6, 1:3) :: cY
	real(8) :: r_, r0_, dr_, Pc_, mc_, nuc_, k2_, rho_
	integer :: i, it, mode 
	
	allocate(r_out(1:sp_N1), z_out(1:3, 1:6, 1:sp_N1), z_full(1:6,1:sp_N1))

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
			cY(1:6,it) = z(1:6)
		enddo

! Transfer to lo_bc_Solve
Cri(1) = z(7)
Cri(2) = dsign(lo_EOS(dabs(z(7)), f = 'rho(p)'), z(7))
if (Cri(2) < 0.d0) pause 'err: rho_ < 0'
Cri(3) = z(8)
Cri(4) = lo_mu(z(7))
Cri(5) = z(9)
Cri(6) = sp_r(sp_N1)
		call lo_bc_EF_Solve(cY(1:6,1:3), Coef(1:3), H0F, dH0F, KF)
		
		y(1) = H0F
		y(2) = dH0F
		y(3) = KF
		y(4) = z(7)
		y(5) = z(8)
		y(6) = z(9)
open(02, file = 'data/temp_lo_y_EF.txt', status = 'replace')
z_full = 0.d0
do i = 1, sp_N1
z_full(1:6, i) = z_full(1:6, i) + Coef(1)*z_out(1, 1:6, i) + Coef(2)*z_out(2, 1:6, i)+ Coef(3)*z_out(3, 1:6, i)
	call Write2R8(02, FR8,r_out(i)/sp_r(sp_N2), r_out(i)*z_full(2, i)/z_full(1, i))
enddo
write(*,*) "r*dH0/H0", sp_r(sp_N1)*(z_full(2,sp_N1) - dH0F)/z_full(1,sp_N1)
		open(01, file = pef_LD_Cr3, status = 'replace')
		region = 'crust'
		call Write2R8(02, FR8,sp_r(sp_N1)/sp_r(sp_N2), sp_r(sp_N1)*y(2)/y(1))
		do i = sp_N1, sp_N2-1
				r_ = sp_r(i)
				dr_ = delr(i)
			x(1:6) = y(1:6)
			call RK4(lo_fluid_eqt_TOV_dr, 6, r_, r_ + dr_, x(1:6), y(1:6))
			call Write1R81AR8(01, FR8,r_ + dr_, (y(1:3)), m=3)
			call Write2R8(02, FR8,(r_ + dr_)/sp_r(sp_N2), (r_ + dr_)*y(2)/y(1))
		enddo
		close(01)
call Write2R8(02, FR8,(r_ + dr_)/sp_r(sp_N2), (r_ + dr_)*y(2)/y(1) - 4.d0*pi*lo_EOS(dabs(y(4)))*(r_ + dr_)**3/y(5))
close(02)
		!open(01, file = pef_LD_Co2, status = 'replace')
		!do i = 1, sp_N1-1
		!	call Write1R81AR8(01, FR8,r_out(i), Coef(1)*z_out(1,1:6, i)+Coef(2)*z_out(2,1:6, i)+Coef(3)*z_out(3,1:6, i), m=6)
		!enddo
		!close(01)

		!open(01, file = 'data/temp_love_Cr1.txt', status = 'replace')
		!do i = 1, sp_N1-1
		!	call Write1R81AR8(01, FR8,r_out(i), z_out(1,1:6, i), m=6)
		!enddo
		!close(01)
		!open(01, file = 'data/temp_love_Cr2.txt', status = 'replace')
		!do i = 1, sp_N1-1
		!	call Write1R81AR8(01, FR8,r_out(i), z_out(2,1:6, i), m=6)
		!enddo
		!close(01)
		!open(01, file = 'data/temp_love_Cr3.txt', status = 'replace')
		!do i = 1, sp_N1-1
		!	call Write1R81AR8(01, FR8,r_out(i), z_out(3,1:6, i), m=6)
		!enddo
		!close(01)
rho_ = dsign(lo_EOS(dabs(y(4)), f = 'rho(p)'), y(4))
if (rho_ < 0.d0) pause 'err: rho_ < 0'
		call lo_k2(sp_r(sp_N2)*y(2)/y(1), 0.d0, y(4), rho_, y(5), 0.d0, sp_r(sp_N2), k2_)
		write(*,*) "EF model k2", k2_
		write(*,*) sp_r(sp_N2)*y(2)/y(1) - 4.d0*pi*rho_*r_**3/y(5)
		write(*,*) "Lambda_bar", 2.d0/3.d0* k2_ /(M0/R0)**5
		write(*,*) "=========================================================="

		deallocate(r_out, z_out, z_full)
		Cri(1:6) = 0.d0
		Crf(1:6) = 0.d0

	endsubroutine lo01_iterate_EF


	! Integrate Solid Core Fluid Mid Solid Envelop---------------------------------------------------
	subroutine lo01_iterate_EFE
	!## Strange Results - Cannot reduce to 2 Solid Free Slip Case
	use Format_IO
	use FWrite
	use RK4_Set
	use lo_match
	use lo_eqt
use Root_Finding
use metric_var
	implicit none
	real(8), dimension(1:9) :: x, y, z
	real(8), allocatable :: r_out(:)
	real(8), allocatable:: z_out(:,:,:), z_full_Co(:,:), z_full_Oc(:,:), y_full_Cr(:,:)
	real(8) :: Coef(1:3), CoefZ(1:5)
	real(8) :: H0F, dH0F, KF, H0_Co, K_Co
	real(8), dimension(1:6, 1:3) :: cY
	real(8), dimension(1:6, 1:5) :: cZ
	real(8) :: r_, r0_, dr_, Pc_, mc_, nuc_, k2_, rho_
	integer :: i, it, mode, j
	integer :: sp_NFE
	
	allocate(r_out(1:sp_N2), z_out(1:8, 1:6, 1:sp_N2), y_full_Cr(1:6,sp_N1:sp_N2), z_full_Co(1:6, 1:sp_N1), z_full_Oc(1:6, sp_N1:sp_N2))

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
			cY(1:6,it) = z(1:6)
		enddo

! Transfer to lo_bc_Solve
Cri(1) = z(7)
Cri(2) = dsign(lo_EOS(dabs(z(7)), f = 'rho(p)'), z(7))
if (Cri(2) < 0.d0) pause 'err: rho_ < 0'
Cri(3) = z(8)
Cri(4) = lo_mu(z(7))
Cri(5) = z(9)
Cri(6) = sp_r(sp_N1)
		call lo_bc_EF_Solve(cY(1:6,1:3), Coef(1:3), H0F, dH0F, KF)
		
		y(1) = H0F
		y(2) = dH0F
		y(3) = KF
		y(4) = z(7)
		y(5) = z(8)
		y(6) = z(9)
		
		! Finding the crust ocean transition point
		do i = sp_N1, sp_N2-1
			if ((P(i) - P_g*(Grav_Const/c**4))* ( P(i+1) - P_g*(Grav_Const/c**4)) <= 0.d0) then
				sp_NFE = i
				exit
			endif
		enddo
		if (i >= sp_N2-1) pause 'sp_NFE out of range'
			
			r_out(sp_N1) = sp_r(sp_N1)
			y_full_Cr(1:3, sp_N1) = y(1:3)
			y_full_Cr(4:6, sp_N1) = 0.d0
			do i = sp_N1, sp_NFE-1
					r_ = sp_r(i)
					dr_ = delr(i)
				x(1:6) = y(1:6)
				call RK4(lo_fluid_eqt_TOV_dr, 6, r_, r_ + dr_, x(1:6), y(1:6))
				call Write1R81AR8(01, FR8,r_ + dr_, (y(1:3)), m=3)
				r_out(i+1) = sp_r(i+1)
				y_full_Cr(1:3, i+1) = y(1:3)
				y_full_Cr(4:6, i+1) = 0.d0
			enddo
		H0_Co = y(1)
		K_Co = y(3)

		Pc_ = y(4)
		mc_ = y(5) 
		nuc_ = y(6)
		r0_ = sp_r(sp_NFE)

! Transfer to lo_bc_Solve
Cri(1) = Pc_
Cri(2) = dsign(lo_EOS(dabs(Pc_), f = 'rho(p)'), Pc_)
if (Cri(2) < 0.d0) pause 'err: rho_ < 0'
Cri(3) = mc_
Cri(4) = lo_mu(Pc_)
Cri(5) = nuc_
Cri(6) = r0_

		region = 'crust'
		do it = 4, 8
			mode = it - 3
			call lo01_bc_Cri(mode, z(1:6))
			z(7) = Pc_
			z(8) = mc_
			z(9) = nuc_
			r_out(sp_NFE) = sp_r(sp_NFE)
			z_out(it, 1:6, sp_NFE) = z(1:6)
			do i = sp_NFE, sp_N2-1
				r_ = sp_r(i)
				dr_ = delr(i)
				x(1:9) = z(1:9)
				call RK4(lo_solid_eqt_TOV_dr, 9, r_, r_ + dr_, x(1:9), z(1:9))
				r_out(i+1) = r_ + dr_
				z_out(it, 1:6, i+1) = z(1:6)
			enddo
			cZ(1:6,it - 3) = z(1:6)
		enddo
! Transfer to lo_bc_Solve
Crf(1) = z(7)
Crf(2) = dsign(lo_EOS(dabs(z(7)), f = 'rho(p)'), z(7))
if (Crf(2) < 0.d0) pause 'err: rho_ < 0'
Crf(3) = z(8)
Crf(4) = lo_mu(z(7))
Crf(5) = z(9)
Crf(6) = sp_r(sp_N2)
		call lo_bc_Solve(H0_Co, K_Co, cZ(1:6, 1:5), k2_, CoefZ(1:5))


		open(02, file = 'data/temp_lo_y_EFE.txt', status = 'replace')
		do i = 1, sp_N1
			do j = 1,6
				z_full_Co(j, i) = vec_dot(z_out(1:3, j, i), Coef(1:3))
			enddo
!			call Write2R8(02,FR8,r_out(i)/sp_r(sp_N2),r_out(i)*z_full_Co(2, i)/z_full_Co(1, i))
write(02,FR8,advance = 'no') r_out(i)
write(02,FR8,advance = 'no') z_full_Co(1,i)
write(02,FR8,advance = 'no') z_full_Co(2,i)
write(02,FR8,advance = 'no') z_full_Co(3,i)
write(02,FR8,advance = 'no') z_full_Co(4,i)
write(02,FR8,advance = 'no') z_full_Co(5,i)
write(02,FR8) z_full_Co(6,i)
		enddo
		
		do i = sp_N1, sp_NFE
!			call Write2R8(02,FR8,r_out(i)/sp_r(sp_N2),r_out(i)*y_full_Cr(2, i)/y_full_Cr(1, i))
		enddo

		do i = sp_NFE, sp_N2
			do j = 1,6
				z_full_Oc(j, i) = vec_dot(z_out(4:8, j, i), CoefZ(1:5))
			enddo
!			call Write2R8(02,FR8,r_out(i)/sp_r(sp_N2),r_out(i)*z_full_Oc(2, i)/z_full_Oc(1, i))
write(02,FR8,advance = 'no') r_out(i)
write(02,FR8,advance = 'no') z_full_Oc(1,i)
write(02,FR8,advance = 'no') z_full_Oc(2,i)
write(02,FR8,advance = 'no') z_full_Oc(3,i)
write(02,FR8,advance = 'no') z_full_Oc(4,i)
write(02,FR8,advance = 'no') z_full_Oc(5,i)
write(02,FR8) z_full_Oc(6,i)
		enddo
!			call Write2R8(02,FR8,r_out(sp_N2)/sp_r(sp_N2),solR_temp)
		close(02)
i = sp_NFE
write(*,*) "Dy from result", r_out(i)*y_full_Cr(2, i)/y_full_Cr(1, i) - r_out(i)*z_full_Oc(2, i)/z_full_Oc(1, i)
write(*,*) "Dy due to mu", 32.d0*pi*r_out(i)* (2.d0*met_dnu(Cri(1), Cri(2), Cri(3), r_out(i))) * Cri(4) * z_full_Oc(5, i)/z_full_Oc(1, i)

		deallocate(r_out, z_out, y_full_Cr, z_full_Co, z_full_Oc)
		Cri(1:6) = 0.d0
		Crf(1:6) = 0.d0

		write(*,*) "Rg", sp_r(sp_NFE)/1.d5
		write(*,*) "EFE model k2", k2_
		write(*,*) "Lambda_bar", 2.d0/3.d0* k2_ /(M0/R0)**5
		write(*,*) "=========================================================="
	endsubroutine lo01_iterate_EFE


	! Integrate Solid Core Fluid Mid Solid Envelop---------------------------------------------------
	subroutine lo01_iterate_EFE2
	use Format_IO
	use FWrite
	use RK4_Set
	use lo_match
	use lo_eqt
use Root_Finding
use metric_var
	implicit none
	real(8), dimension(1:9) :: x, y, z
	real(8), allocatable :: r_out(:)
	real(8), allocatable:: z_out(:,:,:), z_full_Co(:,:), z_full_Oc(:,:), y_full_Cr(:,:)
	real(8) :: Coef(1:3), CoefZ(1:5)
	real(8) :: H0F, dH0F, KF, H0_Co, K_Co
	real(8), dimension(1:6, 1:3) :: cY
	real(8), dimension(1:6, 1:5) :: cZ
	real(8) :: r_, r0_, dr_, Pc_, mc_, nuc_, Pf_, mf_, nuf_, k2_, rho_
	integer :: i, it, mode, j
	integer :: sp_NFE
	
	allocate(r_out(1:sp_N2), z_out(1:8, 1:6, 1:sp_N2), y_full_Cr(1:6,sp_N1:sp_N2), z_full_Co(1:6, 1:sp_N1), z_full_Oc(1:6, sp_N1:sp_N2))

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
			cY(1:6,it) = z(1:6)
		enddo

! Transfer to lo_bc_Solve
Cri(1) = z(7)
Cri(2) = dsign(lo_EOS(dabs(z(7)), f = 'rho(p)'), z(7))
if (Cri(2) < 0.d0) pause 'err: rho_ < 0'
Cri(3) = z(8)
Cri(4) = lo_mu(z(7))
Cri(5) = z(9)
Cri(6) = sp_r(sp_N1)
		call lo_bc_EF_Solve(cY(1:6,1:3), Coef(1:3), H0F, dH0F, KF)

		y(1) = H0F
		y(2) = dH0F
		y(3) = KF
		y(4) = z(7)
		y(5) = z(8)
		y(6) = z(9)
		
		! Finding the crust ocean transition point
		do i = sp_N1, sp_N2-1
			if ((P(i) - P_g*(Grav_Const/c**4))* ( P(i+1) - P_g*(Grav_Const/c**4)) <= 0.d0) then
				sp_NFE = i
				exit
			endif
		enddo
		if (i >= sp_N2-1) pause 'sp_NFE out of range'
			
			r_out(sp_N1) = sp_r(sp_N1)
			y_full_Cr(1:3, sp_N1) = y(1:3)
			y_full_Cr(4:6, sp_N1) = 0.d0
			do i = sp_N1, sp_NFE-1
					r_ = sp_r(i)
					dr_ = delr(i)
				x(1:6) = y(1:6)
				call RK4(lo_fluid_eqt_TOV_dr, 6, r_, r_ + dr_, x(1:6), y(1:6))
				call Write1R81AR8(01, FR8,r_ + dr_, (y(1:3)), m=3)
				r_out(i+1) = sp_r(i+1)
				y_full_Cr(1:3, i+1) = y(1:3)
				y_full_Cr(4:6, i+1) = 0.d0
			enddo
		H0_Co = y(1)
		K_Co = y(3)

		Pc_ = y(4)
		mc_ = y(5) 
		nuc_ = y(6)
		r0_ = sp_r(sp_NFE)
! Transfer to lo_bc_Solve2
Cri(1) = Pc_
Cri(2) = dsign(lo_EOS(dabs(Pc_), f = 'rho(p)'), Pc_)
if (Cri(2) < 0.d0) pause 'err: rho_ < 0'
Cri(3) = mc_
Cri(4) = lo_mu(Pc_)
Cri(5) = nuc_
Cri(6) = sp_r(sp_NFE)

		z(7) = Pc_
		z(8) = mc_
		z(9) = nuc_
! Solve TOV equation to obtain surface P, rho, m
			do i = sp_NFE, sp_N2-1
				r_ = sp_r(i)
				dr_ = delr(i)
				x(1:9) = z(1:9)
				call RK4(lo_solid_eqt_TOV_dr, 9, r_, r_ + dr_, x(1:9), z(1:9))
			enddo

		Pf_ = z(7)
		mf_ = z(8)
		nuf_ = z(9) 
		r0_ = sp_r(sp_N2)

! Transfer to lo_bc_Solve2
Crf(1) = Pf_
Crf(2) = dsign(lo_EOS(dabs(Pf_), f = 'rho(p)'), Pf_)
if (Pf_ < 0.d0) pause 'err: rho_ < 0'
Crf(3) = mf_
Crf(4) = lo_mu(Pf_)
Crf(5) = nuf_
Crf(6) = r0_

		region = 'crust'
		do it = 4, 8
			mode = it - 3
			call lo01_bc_Crf(mode, z(1:6))
			z(7) = Pf_
			z(8) = mf_
			z(9) = nuf_

			r_out(sp_N2) = sp_r(sp_N2)
			z_out(it, 1:6, sp_N2) = z(1:6)

			do i = sp_N2, sp_NFE+1, -1
				r_ = sp_r(i)
				dr_ = delr(i-1)
				x(1:9) = z(1:9)
				call RK4(lo_solid_eqt_TOV_dr, 9, r_, r_ - dr_, x(1:9), z(1:9))
				r_out(i-1) = r_ - dr_
				z_out(it, 1:6, i-1) = z(1:6)
			enddo
			cZ(1:6,it - 3) = z(1:6)
		enddo

		call lo_bc_Solve_2Side(H0_Co, K_Co, cZ(1:6, 1:5), k2_, CoefZ(1:5))

		open(02, file = 'data/temp_lo_y_EFE.txt', status = 'replace')
		do i = 1, sp_N1
			do j = 1,6
				z_full_Co(j, i) = vec_dot(z_out(1:3, j, i), Coef(1:3))
			enddo
!			call Write2R8(02,FR8,r_out(i)/sp_r(sp_N2),r_out(i)*z_full_Co(2, i)/z_full_Co(1, i))
write(02,FR8,advance = 'no') r_out(i)
write(02,FR8,advance = 'no') z_full_Co(1,i)
write(02,FR8,advance = 'no') z_full_Co(2,i)
write(02,FR8,advance = 'no') z_full_Co(3,i)
write(02,FR8,advance = 'no') z_full_Co(4,i)
write(02,FR8,advance = 'no') z_full_Co(5,i)
write(02,FR8) z_full_Co(6,i)
		enddo
		
		do i = sp_N1, sp_NFE
!			call Write2R8(02,FR8,r_out(i)/sp_r(sp_N2),r_out(i)*y_full_Cr(2, i)/y_full_Cr(1, i))
		enddo

		do i = sp_NFE, sp_N2
			do j = 1,6
				z_full_Oc(j, i) = vec_dot(z_out(4:8, j, i), CoefZ(1:5))
			enddo
!			call Write2R8(02,FR8,r_out(i)/sp_r(sp_N2),r_out(i)*z_full_Oc(2, i)/z_full_Oc(1, i))
write(02,FR8,advance = 'no') r_out(i)
write(02,FR8,advance = 'no') z_full_Oc(1,i)
write(02,FR8,advance = 'no') z_full_Oc(2,i)
write(02,FR8,advance = 'no') z_full_Oc(3,i)
write(02,FR8,advance = 'no') z_full_Oc(4,i)
write(02,FR8,advance = 'no') z_full_Oc(5,i)
write(02,FR8) z_full_Oc(6,i)
		enddo
!			call Write2R8(02,FR8,r_out(sp_N2)/sp_r(sp_N2),solR_temp)
		close(02)
i = sp_NFE
write(*,*) "Dy from result", r_out(i)*y_full_Cr(2, i)/y_full_Cr(1, i) - r_out(i)*z_full_Oc(2, i)/z_full_Oc(1, i)
write(*,*) "Dy due to mu", 32.d0*pi*r_out(i)* (2.d0*met_dnu(Cri(1), Cri(2), Cri(3), r_out(i))) * Cri(4) * z_full_Oc(5, i)/z_full_Oc(1, i)

		deallocate(r_out, z_out, y_full_Cr, z_full_Co, z_full_Oc)
		Cri(1:6) = 0.d0
		Crf(1:6) = 0.d0

		write(*,*) "Rg", sp_r(sp_NFE)/1.d5
		write(*,*) "EFE model k2", k2_
		write(*,*) "Lambda_bar", 2.d0/3.d0* k2_ /(M0/R0)**5
		write(*,*) "=========================================================="
	endsubroutine lo01_iterate_EFE2

!Integrate Fluid Core Solid Envelop 2Side---------------------------------------------------
	subroutine lo01_iterate_FEF
	! FE - Fluid Core Elastic Crust
	use Format_IO
	use FWrite
	use RK4_Set
!	use lo_match
	use lo_match_P02
	use lo_eqt
	use Root_Finding
	implicit none
	real(8), dimension(1:9) :: x, y, z
	real(8), allocatable :: r_out(:)
	real(8), allocatable:: z_out(:,:,:)
	real(8) :: Coef(1:5)
	real(8) :: H0_Co, K_Co, H0F, KF, dH0F
	real(8), dimension(1:6, 1:5) :: cY
	real(8) :: r_, dr_, Pc_, mc_, nuc_, Pf_, mf_, nuf_, k2_, rho_
	integer :: i, it, mode, sp_NEF
	
	allocate(r_out(sp_N1:sp_N2), z_out(1:5, 1:6, sp_N1:sp_N2))
		cY = 0.d0
		call lo01_bc_Coi(y(1:6))
		region = 'core'
		do i = 1, sp_N1-1
				r_ = sp_r(i)
				dr_ = delr(i)
			x(1:6) = y(1:6)
			call RK4(lo_fluid_eqt_TOV_dr, 6, r_, r_ + dr_, x(1:6), y(1:6))
			call Write1R81AR8(01, FR8,r_ + dr_, (y(1:3)), m=3)
		enddo

		H0_Co = y(1)
		K_Co = y(3)
		Pc_ = y(4)
		mc_ = y(5) 
		nuc_ = y(6)

		! Transfer to lo_bc_Solve
		Cri(1) = Pc_
		Cri(2) = dsign(lo_EOS(dabs(Pc_), f = 'rho(p)'), Pc_)
		if (Cri(2) < 0.d0) pause 'err: rho_ < 0'
		Cri(3) = mc_
		Cri(4) = lo_mu(Pc_)
		Cri(5) = nuc_
		Cri(6) = sp_r(sp_N1)

		! Finding the crust ocean transition point
		do i = sp_N1, sp_N2-1
			if ((P(i) - P_g*(Grav_Const/c**4))* ( P(i+1) - P_g*(Grav_Const/c**4)) <= 0.d0) then
				sp_NEF = i
				exit
			endif
		enddo
		if (i >= sp_N2-1) pause 'sp_Nmid out of range'

		region = 'crust'
		do it = 2, 6
			mode = it - 1
			call lo01_bc_Cri(mode, z(1:6))
			z(7) = Pc_
			z(8) = mc_
			z(9) = nuc_
			z_out(mode, 1:6, sp_N1) = z(1:6)
			r_out(sp_N1) = sp_r(sp_N1)

			do i = sp_N1, sp_NEF-1
				r_ = sp_r(i)
				dr_ = delr(i)
				x(1:9) = z(1:9)
				call RK4(lo_solid_eqt_TOV_dr, 9, r_, r_ + dr_, x(1:9), z(1:9))
				r_out(i+1) = r_ + dr_
				z_out(mode, 1:6, i+1) = z(1:6)
			enddo
			cY(1:6,it - 1) = z(1:6)
		enddo
		! Transfer to lo_bc_Solve
		Crf(1) = z(7)
		Crf(2) = dsign(lo_EOS(dabs(z(7)), f = 'rho(p)'), z(7))
		if (Crf(2) < 0.d0) pause 'err: rho_ < 0'
		Crf(3) = z(8)
		Crf(4) = lo_mu(z(7))
		Crf(5) = z(9)
		Crf(6) = sp_r(sp_NEF)

		call lo_bc_FEF_Solve(H0_Co, K_Co, cY(1:6,1:5), Coef(1:5), H0F, dH0F, KF)
		
		y(1) = H0F
		y(2) = dH0F
		y(3) = KF
		y(4) = z(7)
		y(5) = z(8)
		y(6) = z(9)

		open(01, file = pef_LD_Cr3, status = 'replace')
		region = 'ocean'
		do i = sp_NEF, sp_N2-1
				r_ = sp_r(i)
				dr_ = delr(i)
			x(1:6) = y(1:6)
			call RK4(lo_fluid_eqt_TOV_dr, 6, r_, r_ + dr_, x(1:6), y(1:6))
			call Write1R81AR8(01, FR8,r_ + dr_, (y(1:3)), m=3)
		enddo
		close(01)

rho_ = dsign(lo_EOS(dabs(y(4)), f = 'rho(p)'), y(4))
if (rho_ < 0.d0) pause 'err: rho_ < 0'
		call lo_k2(sp_r(sp_N2)*y(2)/y(1), 0.d0, y(4), rho_, y(5), 0.d0, sp_r(sp_N2), k2_)
		write(*,*) "FEF model k2", k2_
		write(*,*) "Lambda_bar", 2.d0/3.d0* k2_ /(M0/R0)**5
		write(*,*) "=========================================================="
		deallocate(r_out, z_out)
		Cri(1:6) = 0.d0
		Crf(1:6) = 0.d0

	endsubroutine lo01_iterate_FEF

	! Fluid Core BC---------------------------------------------------
	subroutine lo01_bc_Coi(y)
	use puls_eqt_set_opt03_ex1
	use lo_eqt
	implicit none
	real(8) :: gamma_1, P_, rho_, m_, nu_, r_, dr_, a0_, a2_
	real(8) :: y(1:6)
		P_ = P(1)
		rho_ = rho(1)
		m_ = m(1)
		nu_ = nu(1)
		r_ = sp_r(1)
gamma_1 = lo_gamma(P_)
		a0_ = 1.d0
		a2_ = -2.d0/(2.d0*l_0+3.d0) * pi * (5.d0*rho_+9.d0*P_+(rho_+P_)/(gamma_1*P_/(rho_+P_)) )

		y(1) = a0_ * r_**l_0 * (1.d0 + r_ ** 2 *  a2_ )
		y(2) = a0_ * r_**(l_0-1.d0) * (l_0 + r_ ** 2 * (l_0+2.d0) * a2_  )
		y(3) = a0_ * r_**l_0  * ( 1.d0 +  ((l_0+2.d0)*a2_ + 8.d0*pi/3.d0*(3.d0*P_+rho_))/(l_0+2.d0) )   ! need second order expansion !!

		y(4) = P_
		y(5) = m_
		y(6) = nu_
	endsubroutine lo01_bc_Coi

	!Solid Core BC---------------------------------------------------
	subroutine lo01_bc_Coi_Solid(mode, z)
	use lo_eqt
	implicit none
	real(8) :: gamma_1, P_, rho_, m_, mu_, nu_, r_, dr_
	real(8) :: H00, K0, W0, W2, V0, V2, T20, T22, T(1:3)
	real(8) :: z(1:9)
	integer :: mode
		P_ = P(1)
		rho_ = rho(1)
		m_ = m(1)
		region = 'core'
		mu_ = lo_mu(P_)
		nu_ = nu(1)
		r_ = sp_r(1)
gamma_1 = lo_gamma(P_)
		K0 = 0.d0
		V2 = 0.d0
		V0 = 0.d0
		if (mode == 1) then
			K0 = P_+rho_ !1.d0
		elseif (mode == 2) then
			K0 = -(P_+rho_) !1.d0
			V2 = 1.d0/r_**2
		elseif (mode == 3) then 
			V0 = 1.d0
		endif

		W0 = l_0 * V0
		H00 = K0 - 64.d0*pi*mu_*V0
		T20 = -4.d0*mu_*( (l_0-2.d0)*V0 + W0 )
			T(1) = 0.5d0*( (3.d0 * gamma_1 + 1.d0)*P_ + rho_ )
			T(2) = - (gamma_1 * P_ * l_0 + 2.d0 * mu_/3.d0 * (l_0 - 6.d0)) * (l_0 + 1.d0)  ! Finn (1990) missing (l_0 + 1) ???
			T(3) = (gamma_1 * P_ * (l_0 + 3.d0) + 2.d0 * mu_/3.d0 * (l_0 + 9.d0))			! Added mu*2 for Andersson (2014)
		W2 = (T(1)*(-H00) + T(2)*(-V2))/T(3)										! Negative Signs account for Notation in Andersson 2011
		T22 = -4.d0*mu_*( l_0*V2 + W2 )

		z(1) = r_**l_0 *H00
		z(2) = l_0 * r_**(l_0-1.d0) *H00
		z(3) = r_**l_0 *K0
		z(4) = r_**l_0 *(W0 + r_**2*W2)
		z(5) = r_**l_0 * (V0 + r_**2*V2)
		z(6) = r_**l_0 * (T20+ r_**2*T22)

		z(7) = P_
		z(8) = m_
		z(9) = nu_
	endsubroutine lo01_bc_Coi_Solid

	! Fluid Solid Interface---------------------------------------------------
	subroutine lo01_bc_Cri(mode, z)
	use puls_eqt_set_opt03_ex1
	!	2 components
	!	Quantity Q = {QQ(i)} [Y] {c(i)} ------ {vec} [2nd Tensor] {vec} inner product
	implicit none
	real(8) :: z(1:6)
	integer :: mode
		
		z = 0.d0
		if (mode == 1) z(1) = 1.d0 !/r_Cr(0)**2
		if (mode == 2) z(2) = 1.d0 !/r_Cr(0)**2
		if (mode == 3) z(3) = 1.d0 !/r_Cr(0)**2
		if (mode == 4) z(4) = 1.d0
		if (mode == 5) z(5) = 1.d0

	endsubroutine lo01_bc_Cri

	! Solid Crust Surface---------------------------------------------------
	subroutine lo01_bc_Crf(mode, z)
	use puls_eqt_set_opt03_ex1
	!	2 components
	!	Quantity Q = {QQ(i)} [Y] {c(i)} ------ {vec} [2nd Tensor] {vec} inner product
	implicit none
	real(8) :: z(1:6)
	integer :: mode
		
		z = 0.d0
		if (mode == 1) z(1) = 1.d0 !/r_Cr(0)**2
		if (mode == 2) z(2) = 1.d0 !/r_Cr(0)**2
		if (mode == 3) z(3) = 1.d0 !/r_Cr(0)**2
		if (mode == 4) z(4) = 1.d0
		if (mode == 5) z(5) = 1.d0

	endsubroutine lo01_bc_Crf

	! Solid Fluid Interface---------------------------------------------------
	subroutine lo01_bc_Ocf(y)
	use puls_eqt_set_opt03_ex1
	implicit none
	real(8) :: gamma_1, P_, rho_, m_, nu_, r_
	real(8) :: y(1:3)

	endsubroutine lo01_bc_Ocf

endmodule love_opt01_P05