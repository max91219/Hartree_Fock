MODULE utils

	! Uses the atom module
	USE atom_mod
	
	IMPLICIT NONE
	
	! Private module function for checking read status values
	PRIVATE :: check_read_stat
	! Private module function for chekcing reads for comments
	PRIVATE :: check_comment

	!-----------------------------------------------------------!
	!--- Interface for a function that writes both 1D and 2D ---!
	!--- Arrays to the console in a way that is easy to read ---!
	!-----------------------------------------------------------!
	INTERFACE write_arr_console
		MODULE PROCEDURE write_arr_console_1D, write_arr_console_2D
	END INTERFACE write_arr_console

	!-----------------------------------------------------------!
	!--- Interface for a function that writes both 1D and 2D ---!
	!--- Arrays to a file in a way that is easy to read and  ---!
	!--- use as input for other programs (e.g. MATLAB)		 ---!
	!-----------------------------------------------------------!
	INTERFACE write_arr_file
		MODULE PROCEDURE write_arr_file_1D, write_arr_file_2D
	END INTERFACE write_arr_file	

CONTAINS

	!###########################################################!
	!	SUBROUTINE: read_input_file								!
	!	Arguments:												!
	!		file_name - input file (taken from CMD line)		!
	!		atoms - array that is filled with atom objects		!
	!				containing the information from the input	!
	!				file										!
	!		thresh - Variable to store the threshold value		!
	!				 from the input file						!
	!	Method:													!
	!		Starts by reading in the number of atoms then		!
	!		allocates the size of the atoms array. Rewinds		!
	!		the file stream to top of file then looks for		!
	!		the threshold value reads it in then rewinds 		!
	!		the file stream. Then searches for the atom			!
	!		positions, types and exponents to read in.			!
	!		Has commenting built into the input file (See		!
	!		SUBROUTINE check_comment below), they can be		!
	!		placed anywhere in the file as long as they are		!
	!		on their own line.									!
	!###########################################################!
	SUBROUTINE read_input_file(file_name, atoms, thresh)
		CHARACTER (len=*), INTENT(in) :: file_name
		CHARACTER (len=300) :: line
		INTEGER :: io_status = 0, N = 0, i = 0
		TYPE(atom), ALLOCATABLE :: atoms(:)
		DOUBLE PRECISION :: thresh
		
		OPEN(UNIT=1, FILE=file_name, STATUS='old', ACTION='read')		
		
		DO
			READ (1,'(A)', IOSTAT=io_status) line
			CALL check_read_stat(io_status)		
			
			IF (TRIM(line) .eq. 'NUMBER_OF_ATOMS') THEN
					CALL check_comment(1)
					READ (1,*,IOSTAT=io_status) N
					CALL check_read_stat(io_status)
					ALLOCATE (atoms(N))
					REWIND (1)
					EXIT
			END IF
		END DO
		
		DO
			READ (1,'(A)', IOSTAT=io_status) line
			CALL check_read_stat(io_status)		
			
			IF (TRIM(line) .eq. 'CUTOFF') THEN
					CALL check_comment(1)
					READ (1,*,IOSTAT=io_status) thresh
					CALL check_read_stat(io_status)
					REWIND (1)
					EXIT
			END IF
		END DO
	
		read_atoms: DO
			READ (1, '(A)',IOSTAT=io_status) line
			
			CALL check_read_stat(io_status)
			IF (TRIM(line) .eq. 'ATOM_POSITIONS') THEN
				DO i=1,N
					CALL check_comment(1)
					READ (1,*,IOSTAT=io_status) atoms(i)%typ, atoms(i)%coords(1),&
 						atoms(i)%coords(2) , atoms(i)%coords(3), atoms(i)%alpha
					CALL check_read_stat(io_status)
				END DO
				EXIT read_atoms
			END IF
		END DO read_atoms

		CLOSE(1)
	END SUBROUTINE read_input_file

	!###########################################################!
	!	SUBROUTINE: check_read_stat								!
	!	Arguments:												!
	!		rd_stat - return status of the read command			!
	!	Method:													!
	!		Checks the read commands return status, can take	!
	!		certian values.										!
	!			rd_stat < 0 - End of File						!
	!			rd_stat > 0 - Error occured						!
	!			rd_stat = 0 - Read executed fine				!
	!###########################################################!
	SUBROUTINE check_read_stat(rd_stat)
		INTEGER :: rd_stat
		IF (rd_stat .eq. 0) THEN
			RETURN
		ELSE IF (rd_stat .gt. 0) THEN
			WRITE (*,*) 'Error reading input file'
			CALL EXIT(1)
		END IF
	END SUBROUTINE check_read_stat

	!###########################################################!
	!	SUBROUTINE: write_arr_console_2D						!
	!	Arguments:												!
	!		mtx - matrix that needs to be written to console	!
	!	Method:													!
	!		Takes the 2D input array and prints it to the		!
	!		console in a way that is easy to read. uses the		!
	!		INTERFACE at the top of this module so that the		!
	!		the same function (write_arr_console) can be		!
	!		called with 1D and 2D arrays and they are dealt		!
	!		with appropriately									!
	!###########################################################!
	SUBROUTINE write_arr_console_2D(mtx)
		INTEGER :: i, j
		DOUBLE PRECISION :: mtx(:,:)
	
		DO i=1, SIZE(mtx,1)
			WRITE (*,*) (mtx(i,j), j=1, SIZE(mtx, 2))
		END DO
	
	END SUBROUTINE write_arr_console_2D

	!###########################################################!
	!	SUBROUTINE: write_arr_console_1D						!
	!	Arguments:												!
	!		mtx - matrix that needs to be written to console	!
	!	Method:													!
	!		Takes the 1D input array and prints it to the		!
	!		console in a way that is easy to read. uses the		!
	!		INTERFACE at the top of this module so that the		!
	!		the same function (write_arr_console) can be		!
	!		called with 1D and 2D arrays and they are dealt		!
	!		with appropriately									!
	!###########################################################!
	SUBROUTINE write_arr_console_1D(mtx)
		INTEGER :: i
		DOUBLE PRECISION :: mtx(:)
		WRITE (*,*) (mtx(i), i=1, SIZE(mtx, 1))
	END SUBROUTINE write_arr_console_1D

	!###########################################################!
	!	SUBROUTINE: write_arr_file_2D							!
	!	Arguments:												!
	!		mtx - matrix that needs to be written to console	!
	!		file_name - name of file to print to				!
	!	Method:													!
	!		Takes the 2D input array and prints it to the		!
	!		file in a way that is easy to read. uses the		!
	!		INTERFACE at the top of this module so that the		!
	!		the same function (write_arr_file) can be			!
	!		called with 1D and 2D arrays and they are dealt		!
	!		with appropriately									!
	!###########################################################!
	SUBROUTINE write_arr_file_2D(mtx, file_name)
		INTEGER :: i, j
		CHARACTER (len = *) :: file_name
		DOUBLE PRECISION :: mtx(:,:)
		
		OPEN (unit = 1, file = 'out/'//TRIM(file_name)//'.out', ACTION='write',&
				 STATUS='replace')
		
		DO i=1, SIZE(mtx,1)
			WRITE (1,*) (mtx(i,j), j=1, SIZE(mtx, 2))
		END DO			
		
		CLOSE (1)
	
	END SUBROUTINE write_arr_file_2D

	!###########################################################!
	!	SUBROUTINE: write_arr_file_1D							!
	!	Arguments:												!
	!		mtx - matrix that needs to be written to console	!
	!		file_name - name of file to print to				!
	!	Method:													!
	!		Takes the 1D input array and prints it to the		!
	!		file in a way that is easy to read. uses the		!
	!		INTERFACE at the top of this module so that the		!
	!		the same function (write_arr_file) can be			!
	!		called with 1D and 2D arrays and they are dealt		!
	!		with appropriately									!
	!###########################################################!
	SUBROUTINE write_arr_file_1D(mtx, file_name)
		INTEGER :: i
		CHARACTER (len = *) :: file_name
		DOUBLE PRECISION :: mtx(:)
		
		OPEN (unit = 1, file = 'out/'//TRIM(file_name)//'.out', ACTION='write',&
				 STATUS='replace')

		WRITE (1,*) (mtx(i), i=1, SIZE(mtx, 1))
		
		CLOSE (1)
	
	END SUBROUTINE write_arr_file_1D
	
	!###########################################################!
	!	SUBROUTINE: check_comment								!
	!	Arguments:												!
	!		unit_num - unit number of file to check				!
	!	Method:													!
	!		Reads next line from file to check for comment if	!
	!		comments is found then searches for end of comment	!
	!		block. Once found rewinds the file stream one step	!
	!		then returns control to calling code				!
	!###########################################################!
	SUBROUTINE check_comment(unit_num)
		CHARACTER (len = 300) line
		INTEGER :: unit_num, io_status = 0
		
		READ (1,'(A)', IOSTAT=io_status) line	
		IF (line(1:1) .eq. '#') THEN
			DO WHILE (line(1:1) .eq. '#')
				READ (1,'(A)', IOSTAT=io_status) line
				CALL check_read_stat(io_status)
			END DO
		END IF

		BACKSPACE(1)
	END SUBROUTINE check_comment

END MODULE utils
