program nonradial_osc_ns
use main
implicit none
character (len=30) :: ans

	open(100, file = 'log\log.txt', status = 'replace')

01	continue
	
	write(*,*) " Main Menu:"
	write(*,*) "=========================================================="

	call main_body
	
02	write(*,*) "Redo? (y/n)"
	read(*,*) ans
	if (ans == "y") then
		write(100,*) "END. Redo."
		write(100,*) "=========================================================="
		goto 01
	elseif (ans == "n") then
		write(*,*) "Bye"
	else 
		write(*,*) "wrong output!"
		write(*,*) "=========================================================="
		goto 02
	endif
	close(100)

endprogram nonradial_osc_ns