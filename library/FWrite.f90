module FWrite
!*******************************************************************
!	Formatted output:
!	Codes:
!		Write(i)R8	- write i real8 no.s with FR8 format
!		Write(i)R8(j)AR8	- write i real8 no.s with FR8 format and j real8 arrays 
!	
!	Un			- unit no of output file
!	F			- format of output
!	Ei			- the ith Real8 entity
!	m			- (optional) specify the dimensions of arrays
!
!*******************************************************************
    contains
    SUBROUTINE Write2R8(Un, F,E1, E2)
    implicit none
    REAL(8) :: E1, E2
    integer :: Un
    character(len=*) :: F
        write(Un, F, advance = "no") E1
        write(Un, F) E2
    ENDSUBROUTINE Write2R8

    SUBROUTINE Write3R8(Un, F,E1, E2, E3)
    implicit none
    REAL(8) :: E1, E2, E3
    integer :: Un
    character(len=*) :: F
        write(Un, F, advance = "no") E1
        write(Un, F, advance = "no") E2
        write(Un, F) E3
    ENDSUBROUTINE Write3R8

    SUBROUTINE Write4R8(Un, F,E1, E2, E3, E4)
    implicit none
    REAL(8) :: E1, E2, E3, E4
    integer :: Un
    character(len=*) :: F
        write(Un, F, advance = "no") E1
        write(Un, F, advance = "no") E2
        write(Un, F, advance = "no") E3
        write(Un, F) E4
    ENDSUBROUTINE Write4R8

    SUBROUTINE Write5R8(Un, F,E1, E2, E3, E4, E5)
    implicit none
    REAL(8) :: E1, E2, E3, E4, E5
    integer :: Un
    character(len=*) :: F
        write(Un, F, advance = "no") E1
        write(Un, F, advance = "no") E2
        write(Un, F, advance = "no") E3
        write(Un, F, advance = "no") E4
        write(Un, F) E5
    ENDSUBROUTINE Write5R8

	SUBROUTINE Write6R8(Un, F,E1, E2, E3, E4, E5, E6)
    implicit none
    REAL(8) :: E1, E2, E3, E4, E5, E6
    integer :: Un
    character(len=*) :: F
        write(Un, F, advance = "no") E1
        write(Un, F, advance = "no") E2
        write(Un, F, advance = "no") E3
        write(Un, F, advance = "no") E4
		write(Un, F, advance = "no") E5
        write(Un, F) E6
    ENDSUBROUTINE Write6R8

	SUBROUTINE Write7R8(Un, F,E1, E2, E3, E4, E5, E6, E7)
    implicit none
    REAL(8) :: E1, E2, E3, E4, E5, E6, E7
    integer :: Un
    character(len=*) :: F
        write(Un, F, advance = "no") E1
        write(Un, F, advance = "no") E2
        write(Un, F, advance = "no") E3
        write(Un, F, advance = "no") E4
		write(Un, F, advance = "no") E5
		write(Un, F, advance = "no") E6
        write(Un, F) E7
    ENDSUBROUTINE Write7R8

	SUBROUTINE Write8R8(Un, F,E1, E2, E3, E4, E5, E6, E7, E8)
    implicit none
    REAL(8) :: E1, E2, E3, E4, E5, E6, E7, E8
    integer :: Un
    character(len=*) :: F
        write(Un, F, advance = "no") E1
        write(Un, F, advance = "no") E2
        write(Un, F, advance = "no") E3
        write(Un, F, advance = "no") E4
		write(Un, F, advance = "no") E5
		write(Un, F, advance = "no") E6
		write(Un, F, advance = "no") E7
        write(Un, F) E8
    ENDSUBROUTINE Write8R8

	SUBROUTINE Write9R8(Un, F,E1, E2, E3, E4, E5, E6, E7, E8, E9)
    implicit none
    REAL(8) :: E1, E2, E3, E4, E5, E6, E7, E8, E9
    integer :: Un
    character(len=*) :: F
        write(Un, F, advance = "no") E1
        write(Un, F, advance = "no") E2
        write(Un, F, advance = "no") E3
        write(Un, F, advance = "no") E4
		write(Un, F, advance = "no") E5
		write(Un, F, advance = "no") E6
		write(Un, F, advance = "no") E7
		write(Un, F, advance = "no") E8
        write(Un, F) E9
    ENDSUBROUTINE Write9R8

	SUBROUTINE Write11R8(Un, F,E1, E2, E3, E4, E5, E6, E7, E8, E9, E10, E11)
    implicit none
    REAL(8) :: E1, E2, E3, E4, E5, E6, E7, E8, E9, E10, E11
    integer :: Un
    character(len=*) :: F
        write(Un, F, advance = "no") E1
        write(Un, F, advance = "no") E2
        write(Un, F, advance = "no") E3
        write(Un, F, advance = "no") E4
		write(Un, F, advance = "no") E5
		write(Un, F, advance = "no") E6
		write(Un, F, advance = "no") E7
		write(Un, F, advance = "no") E8
		write(Un, F, advance = "no") E9
		write(Un, F, advance = "no") E10
        write(Un, F) E11
    ENDSUBROUTINE Write11R8

	SUBROUTINE Write12R8(Un, F,E1, E2, E3, E4, E5, E6, E7, E8, E9, E10, E11, E12)
    implicit none
    REAL(8) :: E1, E2, E3, E4, E5, E6, E7, E8, E9, E10, E11, E12
    integer :: Un
    character(len=*) :: F
        write(Un, F, advance = "no") E1
        write(Un, F, advance = "no") E2
        write(Un, F, advance = "no") E3
        write(Un, F, advance = "no") E4
		write(Un, F, advance = "no") E5
		write(Un, F, advance = "no") E6
		write(Un, F, advance = "no") E7
		write(Un, F, advance = "no") E8
		write(Un, F, advance = "no") E9
		write(Un, F, advance = "no") E10
		write(Un, F, advance = "no") E11
        write(Un, F) E12
    ENDSUBROUTINE Write12R8

	SUBROUTINE Write1R81AR8(Un, F,E1, A1, m)
    implicit none
    REAL(8) :: E1, A1(1:)
    integer :: Un, i, n
	integer, optional :: m
    character(len=*) :: F
		if (present(m) == .true.) then 
			n = m
		else
			n = size(A1)
		endif
		write(Un, F, advance = "no") E1
        do i = 1, n
			write(Un, F, advance = "no") A1(i)
		enddo
        write(Un, *)
    ENDSUBROUTINE Write1R81AR8

	SUBROUTINE Write1R82AR8(Un, F,E1, A1, A2, m)
    implicit none
    REAL(8) :: E1, A1(1:), A2(1:)
    integer :: Un, i, n
	integer, optional :: m
    character(len=*) :: F
		if (present(m) == .true.) then 
			n = m
		else
			n = size(A1)
		endif
		write(Un, F, advance = "no") E1
        do i = 1, n
			write(Un, F, advance = "no") A1(i)
		enddo
		do i = 1, n
			write(Un, F, advance = "no") A2(i)
		enddo
        write(Un, *)
    ENDSUBROUTINE Write1R82AR8

	SUBROUTINE Write4R8AN(Un, F,E1, E2, E3, E4)
    implicit none
    REAL(8) :: E1, E2, E3, E4
    integer :: Un
    character(len=*) :: F
        write(Un, F, advance = "no") E1
        write(Un, F, advance = "no") E2
        write(Un, F, advance = "no") E3
        write(Un, F, advance = "no") E4
    ENDSUBROUTINE Write4R8AN

	subroutine WriteList1R8(FName, F, i1, i2, E1)
	implicit none
	real(8) :: E1(:)
	integer :: i1, i2, i
	character (len=*) :: FName, F
		open(831, file = 'FName', status = 'replace')
		do i = i1, i2
			write(831,F) E1(i)
		enddo
		close(831)
	endsubroutine WriteList1R8

	subroutine ReverseList_R8(List)
	! Check for bug
	implicit none
	real(8) :: List(1:)
	real(8), allocatable :: Temp(:)
	integer :: i
		allocate(Temp(1: Size(List)))
		do i = 1, Size(List)
			Temp(i) = List(Size(List)-i+1)
		enddo
		List = Temp
		deallocate(Temp)
	endsubroutine ReverseList_R8
endmodule FWrite
