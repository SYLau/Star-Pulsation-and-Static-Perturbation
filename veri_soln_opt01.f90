module veri_soln_opt01
use global_var
contains

!C5---------------------------------------------------
	subroutine vs01_Ctrl
	implicit none
	integer :: ans
		write(100,*) "veri_soln_opt01: verify solutions"
		write(100,*) "forward integration of core; backward integration of crust"
		write(100,*) "matching at crust core interface"

24		write(*,*) "Verify solutions:"
		write(*,*) "1. Verify Core"
		write(*,*) "2. Verify Crust"
		write(*,*) "0. Back"
		write(*,*) "=========================================================="
		read(*,*) ans

		if (ans /= 1 .and. ans /= 2 .and. ans /= 0) then
			write(*,*) "invalid input"
			goto 24
		endif
		
		if (ans == 1) then 

			call vs01_vx_Co( pef_yA ,'data/veri_Co.txt')
			write(100,*) "vs01_vx_Co: passed"
			write(*,*) "Verified core solution"
			write(*,*) "=========================================================="
		elseif (ans == 2) then 
			call vs01_vx_Cr( pef_zC ,'data/veri_Cr_C.txt')
			call vs01_vx_Cr( pef_zD ,'data/veri_Cr_D.txt')
			write(100,*) "vs01_vx_Cr: passed"
			write(*,*) "Verified crust solution"
			write(*,*) "=========================================================="
		elseif (ans == 0) then
			write(100,*) "=========================================================="
			return
		endif
		
		goto 24
	endsubroutine vs01_Ctrl
!C5 A1---------------------------------------------------
	subroutine vs01_vx_Co(FName, OName)
	use FRead
	use FWrite
	use Format_IO
	use puls_eqt_set_opt01
	use puls_eqt_set_opt02
	implicit none
	real(8), dimension (0: pg_N1 +1) :: r_, x_
	real(8), dimension (1:2, 0: pg_N1 +1) :: y_
	real(8), dimension (1:2, 0: pg_N1) :: dy, err, fcn
	integer :: ny_, i
	character (len=*) :: FName, OName
		call RFile3R8(FName,ny_, r_, y_(1,0:pg_N1 +1), y_(2,0:pg_N1 +1))

		open(1, file = OName, status = 'replace')
		do i = 0, ny_-1
			x_(i) = dlog(r_(i)/P_Co(2*i))
		enddo

		do i = 0, ny_-2
			dy(1:2,i) = (y_(1:2,i+1) - y_(1:2,i))/(x_(i+1) - x_(i))
			if (pes_opt == 1) then
				call pes01_vx_Co(2, x_(i), y_(1:2,i), fcn(1:2,i))
			elseif (pes_opt == 2) then
				call pes02_vx_Co(2, x_(i), y_(1:2,i), fcn(1:2,i))
			endif
			err(1:2,i) = (dy(1:2,i) - fcn(1:2,i))/dy(1:2,i)
			call Write3R8(1,FR8, r_(i), err(1,i), err(2,i))
		enddo
		close(1)
	endsubroutine vs01_vx_Co

!C5 A2---------------------------------------------------	
	subroutine vs01_vx_Cr(FName, OName)
	use FRead
	use FWrite
	use Format_IO
	use Interpolation
	use puls_eqt_set_opt01
	use puls_eqt_set_opt02
	implicit none
	real(8), dimension (0: pg_N2 +1) :: r_, x_
	real(8), dimension (1:4, 0: pg_N2 +1) :: z_
	real(8), dimension (1:4, 0: pg_N2) :: dz, err, fcn
	integer :: nz_, i
	character (len=*) :: FName, OName
		call RFile5R8(FName,nz_, r_, z_(1,0:pg_N2 +1), z_(2,0:pg_N2 +1), z_(3,0:pg_N2 +1), z_(4,0:pg_N2 +1))

		open(1, file = OName, status = 'replace')
		do i = 0, nz_-1
			x_(i) = dlog(r_(i)/P_Cr(2*(pg_N2-i)))
		enddo

		do i = 0, nz_-2
			dz(1:4,i) = (z_(1:4,i+1) - z_(1:4,i))/(x_(i+1) - x_(i))
			if (pes_opt == 1) then
				call pes01_vx_Cr(4, x_(i), z_(1:4,i), fcn(1:4,i))
			elseif (pes_opt == 2) then
				call pes02_vx_Cr(4, x_(i), z_(1:4,i), fcn(1:4,i))
			endif
			err(1:4,i) = (dz(1:4,i) - fcn(1:4,i))/dz(1:4,i)
			call Write5R8(1,FR8, r_(i), err(1,i), err(2,i), err(3,i), err(4,i))
		enddo
		close(1)
	endsubroutine vs01_vx_Cr

endmodule veri_soln_opt01