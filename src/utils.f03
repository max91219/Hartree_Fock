MODULE utils

	USE atom_mod
	
	IMPLICIT NONE
	
	PRIVATE :: check_read_stat
	PRIVATE :: check_comment

	INTERFACE write_arr_console
		MODULE PROCEDURE write_arr_console_1D, write_arr_console_2D
	END INTERFACE write_arr_console

	INTERFACE write_arr_file
		MODULE PROCEDURE write_arr_file_1D, write_arr_file_2D
	END INTERFACE write_arr_file	

CONTAINS

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

	SUBROUTINE check_read_stat(rd_stat)
		INTEGER :: rd_stat
		IF (rd_stat .eq. 0) THEN
			RETURN
		ELSE IF (rd_stat .gt. 0) THEN
			WRITE (*,*) 'Error reading input file'
			CALL EXIT(1)
		END IF
	END SUBROUTINE check_read_stat

	SUBROUTINE write_arr_console_2D(mtx)
		INTEGER :: i, j
		DOUBLE PRECISION :: mtx(:,:)
	
		DO i=1, SIZE(mtx,1)
			WRITE (*,*) (mtx(i,j), j=1, SIZE(mtx, 2))
		END DO
	
	END SUBROUTINE write_arr_console_2D

	SUBROUTINE write_arr_console_1D(mtx)
		INTEGER :: i
		DOUBLE PRECISION :: mtx(:)
		WRITE (*,*) (mtx(i), i=1, SIZE(mtx, 1))
	END SUBROUTINE write_arr_console_1D

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

	SUBROUTINE write_arr_file_1D(mtx, file_name)
		INTEGER :: i
		CHARACTER (len = *) :: file_name
		DOUBLE PRECISION :: mtx(:)
		
		OPEN (unit = 1, file = 'out/'//TRIM(file_name)//'.out', ACTION='write',&
				 STATUS='replace')

		WRITE (1,*) (mtx(i), i=1, SIZE(mtx, 1))
		
		CLOSE (1)
	
	END SUBROUTINE write_arr_file_1D
	
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
