PROGRAM main

!---declares the modules that are used in the program---
USE utils 
USE atom_mod
USE num_alg

IMPLICIT NONE

CHARACTER (len=300) :: file_name = '' 
CHARACTER (len=100) :: out_file
TYPE(atom), ALLOCATABLE :: atoms(:)
DOUBLE PRECISION, ALLOCATABLE :: ovlp(:,:), m_sl(:,:), C(:), D(:,:)
DOUBLE PRECISION, ALLOCATABLE :: ke(:,:), nuc(:,:), H(:,:), F(:,:)
DOUBLE PRECISION, AllOCATABLE :: J_mat(:,:,:,:), K_mat(:,:,:,:)
INTEGER :: N = 0, i = 0, j = 0, k=0, l=0, m=0, N_el = 0
DOUBLE PRECISION :: thresh =0.0D0
DOUBLE PRECISION :: lst_E, curr_E
DOUBLE PRECISION, ALLOCATABLE :: eig_val(:), eig_vec(:,:)


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

ALLOCATE(ovlp(N,N), m_sl(N,N), ke(N,N), nuc(N,N), H(N,N), F(N,N), C(N))
ALLOCATE(J_mat(N,N,N,N), K_mat(N,N,N,N), D(N,N), eig_val(N), eig_vec(N,N))

! Caclulate the matracies S and H
DO i=1,N
	DO j=i,N
		ovlp(j,i) = atoms(i)%overlap(atoms(j), thresh)
		ovlp(i,j) = ovlp(j,i)

		ke(j,i) = atoms(i)%ke_int(atoms(j))
		ke(i,j) = ke(j,i)

		nuc(j,i) = atoms(i)%nuc_int(atoms(j))
		nuc(i,j) = nuc(j,i)
	END DO
END DO

H = ke + nuc

! calculate J and K integrals
DO i=1,N
	DO j=1,N
		DO k=1,N
			DO l=1,N
				J_mat(i,j,k,l) = atoms(i)%J_int(atoms(j), atoms(k), atoms(l))
				K_mat(i,j,k,l) = atoms(i)%K_int(atoms(j), atoms(l), atoms(j))
			END DO
		END DO
	END DO
END DO


! make initial guess at C
C = 1.0D0/DSQRT(DBLE(N))

! calculate the density matrix
D = 0.0D0
DO i=1,N_el
	DO j=1,N
		DO k=1,N
			D(j,k) = D(j,k) + (C(i) * C(j))
		END DO
	END DO
END DO

! calculate matrix to orthogonalise S and F
m_sl = ovlp
CALL sym_lowdin(m_sl)

! reduce J and K integrals into Fock matrix
DO i=1,N
	DO j=1,N
		DO k=1,(N_el/2)
			DO l=1,N
				DO m=1,N
					F(i,j) = H(i,j) + C(l) * C(m) * (2.0D0 * J_mat(i,j,l,m) - K_mat(i,m,l,j))
				END DO
			END DO			
		END DO
	END DO
END DO



END PROGRAM main
