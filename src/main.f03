PROGRAM main

!---declares the modules that are used in the program---
USE utils 
USE atom_mod
USE num_alg

IMPLICIT NONE

CHARACTER (len=300) :: file_name = '' 
CHARACTER (len=100) :: out_file
TYPE(atom), ALLOCATABLE :: atoms(:)
DOUBLE PRECISION, ALLOCATABLE :: ovlp(:,:), s_half(:,:), C_up(:,:), C_dwn(:,:), D_up(:,:)
DOUBLE PRECISION, ALLOCATABLE :: ke(:,:), nuc(:,:), H(:,:), F_up(:,:), F_dwn(:,:), D_dwn(:,:)
DOUBLE PRECISION, AllOCATABLE :: J_mat(:,:,:,:), D_tot(:,:)
INTEGER :: N = 0, i = 0, j = 0, k=0, l=0, m=0, N_el = 0, N_up = 0, N_dwn = 0
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
charge = DBLE(N_el)

IF ( MOD(N_el,2) .eq. 0) THEN
	N_up = MAX((charge/2),1.0)
	N_dwn = N_el - N_up
ELSE
	N_dwn = MAX((charge/2),1.0)
	N_up = N_el - N_up
END IF

ALLOCATE(ovlp(N,N), ke(N,N), nuc(N,N), H(N,N), F_up(N,N), F_dwn(N,N), C_up(N,N), C_dwn(N,N))
ALLOCATE(J_mat(N,N,N,N), D_up(N,N), D_dwn(N,N))

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
C_up = 1.0D0 
C_dwn = 1.0D0

!calculate the density matrix
CALL dense_mat(D_up, C_up, N_up)
CALL dense_mat(D_dwn, C_dwn, N_dwn)
D_tot = D_up + D_dwn

! Calc the Fock Matrix
CALL fock_mat(F_up, D_tot, D_up, J_mat, H)
CALL fock_mat(F_dwn, D_tot, D_dwn, J_mat, H)

!Start SFC interations
CALL scf_iter(F_up, F_dwn, ovlp, C_up, C_dwn, D_tot, D_up, D_dwn, J_mat, H, N_up, N_dwn) 

END PROGRAM main
