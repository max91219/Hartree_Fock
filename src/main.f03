PROGRAM main

!---declares the modules that are used in the program---
USE utils 
USE atom_mod
USE num_alg

IMPLICIT NONE

CHARACTER (len=300) :: file_name = '' ! holder for file name from CMD line arg
CHARACTER (len=100) :: out_file
TYPE(atom), ALLOCATABLE :: atoms(:) ! array to hold atoms of the system
DOUBLE PRECISION, ALLOCATABLE :: ovlp(:,:), m_gs(:,:), m_sl(:,:) ! matracies for overlap etc.
INTEGER :: N = 0, i = 0, j = 0, k=0, l=0 ! N number of atoms, i,j counters
DOUBLE PRECISION :: thresh =0.0D0 ! threshold value

CALL GET_COMMAND_ARGUMENT(1, file_name) ! gets input file name from CMD line

!---Checks input file was specified---
IF (TRIM(file_name) .eq. '') THEN
	WRITE (*,*) 'No input file specified'
	WRITE (*,*) 'Syntax is:'
	WRITE (*,*) './guass INPUT_FILE_NAME'
	CALL EXIT(1)
END IF


CALL read_input_file(TRIM(file_name), atoms, thresh) ! reads input file (see utils module)
N = SIZE(atoms) ! sets number of atoms
ALLOCATE(ovlp(N,N), m_gs(N,N), m_sl(N,N)) ! allocates matracies

!--------------------------------------------------!
!---Calculates the overlap matrix of the system ---!
!---only calculates the upper right half of the ---!
!---matrix. Then mirrors it to the lower left   ---!
!--------------------------------------------------!
DO i=1,N
	DO j=i,N
		ovlp(j,i) = atoms(i)%overlap(atoms(j), thresh)
		ovlp(i,j) = ovlp(j,i)
	END DO
END DO

CALL write_arr_file(ovlp, 'ovlp')
WRITE (*,*) 'The sparsity of the overlap matrix is: '
WRITE (*,*) calc_sparse(ovlp)

m_gs = ovlp
m_sl = ovlp

!--------------------------------------------------!
!---Performs gram-schmidt orthogonalisation and ---!
!---then prints the result of the orthogality 	---!
!---test to a file along with the			    ---!
!---transformation matrix 						---!
!--------------------------------------------------!

CALL gram_schmidt(m_gs)
CALL write_arr_file(m_gs,'m_gs')
CALL write_arr_file(check_orthog(m_gs, ovlp, 'G'),'gs_check')

!--------------------------------------------------!
!---Performs symetric lowdin orthogonalisation  ---!
!---then prints the result of the orthogality 	---!
!---test to a file along with the				---!
!---transformation matrix						---!
!--------------------------------------------------!
CALL sym_lowdin(m_sl)
CALL write_arr_file(m_sl,'m_sl')
CALL write_arr_file(check_orthog(m_sl, ovlp, 'S'),'sl_check')

END PROGRAM main
