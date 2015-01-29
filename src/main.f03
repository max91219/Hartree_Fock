PROGRAM main

!---declares the modules that are used in the program---
USE utils 
USE atom_mod
USE num_alg

IMPLICIT NONE

CHARACTER (len=300) :: file_name = '' 
CHARACTER (len=100) :: out_file
TYPE(atom), ALLOCATABLE :: atoms(:)
DOUBLE PRECISION, ALLOCATABLE :: ovlp(:,:), s_half(:,:), C(:,:), D(:,:)
DOUBLE PRECISION, ALLOCATABLE :: ke(:,:), nuc(:,:), H(:,:), F(:,:)
DOUBLE PRECISION, AllOCATABLE :: J_mat(:,:,:,:)
INTEGER :: N = 0, i = 0, j = 0, k=0, l=0, m=0, N_el = 0
DOUBLE PRECISION :: thresh =0.0D0, energy = 0.0D0, charge

CALL GET_COMMAND_ARGUMENT(1, file_name)

!---Checks input file was specified---
IF (TRIM(file_name) .eq. '') THEN
	WRITE (*,*) 'No input file specified'
	WRITE (*,*) 'Syntax is:'
	WRITE (*,*) './guass INPUT_FILE_NAME'
	CALL EXIT(1)
END IF


CALL read_input_file(TRIM(file_name), atoms, thresh)

N = SIZE(atoms) 
N_el = 2 
charge = 2.0D0

ALLOCATE(ovlp(N,N), ke(N,N), nuc(N,N), H(N,N), F(N,N), C(N,N))
ALLOCATE(J_mat(N,N,N,N), D(N,N))

! Caclulate the matracies S and H
DO i=1,N
	DO j=i,N
		ovlp(i,j) = atoms(i)%overlap(atoms(j), thresh)
		ovlp(j,i) = ovlp(i,j)

		ke(i,j) = atoms(i)%ke_int(atoms(j))
		ke(j,i) = ke(i,j)

		nuc(i,j) = atoms(i)%nuc_int(atoms(j), charge)
		nuc(j,i) = nuc(i,j)
	END DO
END DO

H = ke + nuc

!CAll write_arr_console(H)
!WRITE(*,*) (NEW_LINE('A'), i=1, 4)

! calculate J
DO i=1,N
	DO j=1,N
		DO k=1,N
			DO l=1,N
				J_mat(i,j,k,l) = atoms(i)%J_int(atoms(j), atoms(k), atoms(l))
			END DO
		END DO
	END DO
END DO

! make initial guess at C
C = 1.0D0 

!calculate the density matrix
CALL dense_mat(D, C, N_el)
!CAll write_arr_console(D)
!WRITE(*,*) (NEW_LINE('A'), i=1, 4)

! Calc the Fock Matrix
CALL fock_mat(F, D, J_mat, H)
!CAll write_arr_console(F)
!WRITE(*,*) (NEW_LINE('A'), i=1, 4)


!Start SFC interations
CALL scf_iter(F, ovlp, C, D, J_mat, H, N_el) 

END PROGRAM main
