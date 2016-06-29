module GnuPlot
implicit none
contains

	subroutine GnuCtrl(Dir, title, xlabel, ylabel, x,y)
	real(8), intent(in), dimension(:) :: x,y
	character (len=*) :: Dir, title, xlabel, ylabel
		call GnuCF(Dir, title, xlabel, ylabel)
		call GnuDF2d(Dir, x,y)
	endsubroutine GnuCtrl
	
	subroutine GnuCF(Dir, title, xlabel, ylabel)
	character (len=*) :: Dir, title, xlabel, ylabel
		open(10,file= trim(Dir)//'GnuCF.gp',status='replace')
			write(10,*) 'set terminal postscript'
			write(10,*) 'set output '// trim(Dir)//'/' // trim(title) //'.jpg'
			write(10,*) 'set xlabel '//'"'//trim(xlabel)//'"'
			write(10,*) 'set ylabel '//'"'//trim(ylabel)//'"'
			write(10,*) 'plot '//trim(Dir)//'/GnuDF.dat'// ' using 1:2 with lines title "'//trim(title)//'"'
		close(10)

	endsubroutine GnuCF

	subroutine GnuDF2d(Dir, x,y)
	real(8), intent(in), dimension(:) :: x,y
	integer :: nx, ny, i
	character (len=*) :: Dir
		nx = size(x)
		ny = size(y)
		if (nx /= ny) then
			write(*,*) "err: 2dPlot; array size mismatch"
			pause
		endif
		open(05, file = trim(Dir)// '/GnuDF.dat', status = 'replace')
		do i = 1, nx
			write(05,*) x, y
		enddo
		close(05)
	endsubroutine GnuDF2d
	
endmodule GnuPlot