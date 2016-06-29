module verify
contains

subroutine ver_list_range(L, H, List, N, msg)
implicit none
real(8) :: L, H, List(0:)
integer :: N
character (len=*) :: msg
integer :: i
	do i = 0, N
		if ((List(i)-L)*(List(i)-H) > 0.d0) then
			write(*,*) "err: ", msg
			pause
		endif
	enddo
endsubroutine ver_list_range


endmodule verify