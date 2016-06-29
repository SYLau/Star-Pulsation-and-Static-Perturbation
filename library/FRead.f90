module FRead
!*******************************************************************
!	Formatted input:
!	RFile(i)R8	- read i columns of real8 no.s from file
!	
!	FName		- name of the input file
!	FEntry		- count the no. of entry
!	Ei			- the ith column of Real8 entity
!
!*******************************************************************
contains

function NCol(FName)
implicit none
character (len=*) :: FName
character (len=50) :: line	!256
integer :: i, sb, n, NCol
	sb = 1
	n = 0
	open(820, file = FName, status = 'old')
	read(820,'(A)') line
	close(820)
	line = trim(line)
	do i = 1, len(line)
		if (line(i:i) /= '	' .and. sb == 1) then 
			n = n+1
			sb = 2
		elseif (line(i:i) == '	' .and. sb == 2) then
			sb = 1
		endif
	enddo
	if (n == 0) then
		write(*,*) "err: Number of columns in read file = 0"
		pause
	endif
	NCol = n
endfunction NCol

subroutine RFile2R8(FName,FEntry, E1, E2)
implicit none
real(8) :: E1(1:), E2(1:)
integer :: i, FEntry
character (len=*) :: FName
	i = 1
831	open(830, file = FName, status = 'old')
		read(830,*,end=832) E1(i), E2(i)
		i = i + 1
		goto 831
832	FEntry = i -1
	close(830)
endsubroutine RFile2R8

subroutine RFile3R8(FName,FEntry, E1, E2, E3)
implicit none
real(8) :: E1(1:), E2(1:), E3(1:)
integer :: i, FEntry
character (len=*) :: FName
	i = 1
832	open(831, file = FName, status = 'old')
		read(831,*,end=833) E1(i), E2(i), E3(i)
		i = i + 1
		goto 832
833	FEntry = i -1
	close(831)
endsubroutine RFile3R8

subroutine RFile4R8(FName,FEntry, E1, E2, E3, E4)
implicit none
real(8) :: E1(1:), E2(1:), E3(1:), E4(1:)
integer :: i, FEntry
character (len=*) :: FName
	i = 1
842	open(841, file = FName, status = 'old')
		read(841,*,end=843) E1(i), E2(i), E3(i), E4(i)
		i = i + 1
		goto 842
843	FEntry = i -1
	close(841)
endsubroutine RFile4R8

subroutine RFile5R8(FName,FEntry, E1, E2, E3, E4, E5)
implicit none
real(8) :: E1(1:), E2(1:), E3(1:), E4(1:), E5(1:)
integer :: i, FEntry
character (len=*) :: FName
	i = 1
852	open(851, file = FName, status = 'old')
		read(851,*,end=853) E1(i), E2(i), E3(i), E4(i), E5(i)
		i = i + 1
		goto 852
853	FEntry = i -1
	close(851)
endsubroutine RFile5R8


subroutine RFile6R8(FName,FEntry, E1, E2, E3, E4, E5, E6)
implicit none
real(8) :: E1(1:), E2(1:), E3(1:), E4(1:), E5(1:), E6(1:)
integer :: i, FEntry
character (len=*) :: FName
	i = 1
852	open(851, file = FName, status = 'old')
		read(851,*,end=853) E1(i), E2(i), E3(i), E4(i), E5(i), E6(i)
		i = i + 1
		goto 852
853	FEntry = i -1
	close(851)
endsubroutine RFile6R8


subroutine RFile7R8(FName,FEntry, E1, E2, E3, E4, E5, E6, E7)
implicit none
real(8) :: E1(1:), E2(1:), E3(1:), E4(1:), E5(1:), E6(1:), E7(1:)
integer :: i, FEntry
character (len=*) :: FName
	i = 1
852	open(851, file = FName, status = 'old')
		read(851,*,end=853) E1(i), E2(i), E3(i), E4(i), E5(i), E6(i), E7(i)
		i = i + 1
		goto 852
853	FEntry = i -1
	close(851)
endsubroutine RFile7R8

subroutine RFile1R82AR8(FName,FEntry, E1, A1, A2, m)
implicit none
real(8) :: E1(1:), A1(1:, 1:), A2(1:, 1:)
integer :: i, FEntry, n
integer, optional :: m
character (len=*) :: FName
	i = 1
	if (present(m) == .true.) then 
			n = m
		else
			n = size(A1,1)
	endif
852	open(851, file = FName, status = 'old')
		read(851,*,end=853) E1(i), A1(1:n,i), A2(1:n,i)
		i = i + 1
		goto 852
853	FEntry = i -1
	close(851)
endsubroutine RFile1R82AR8

endmodule FRead