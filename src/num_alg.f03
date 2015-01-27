MODULE num_alg

IMPLICIT NONE

CONTAINS
	
    SUBROUTINE gram_schmidt(mtx)
		DOUBLE PRECISION ::  mtx(:,:)
		INTEGER :: INFO = 0, i = 0, j = 0
			
		DO i=1,SIZE(mtx,1)
			DO j=(i+1),SIZE(mtx,1)
				mtx(i,j)= 0.0D0
			END DO
		END DO

		CALL DPOTRF('L', SIZE(mtx,1), mtx, SIZE(mtx,1), INFO)
		CALL check_lapack(INFO)

		CALL DTRTRI('L','N', SIZE(mtx,1), mtx, SIZE(mtx,1), INFO)
		CALL check_lapack(info)
	END SUBROUTINE gram_schmidt
    
	SUBROUTINE sym_lowdin(mtx)
		DOUBLE PRECISION ::  mtx(:,:)
		DOUBLE PRECISION, ALLOCATABLE :: WORK(:), e_val(:), e_vec_ct(:,:)
		INTEGER :: LWORK = -1, INFO = 0, i = 0
		
		ALLOCATE(WORK(1), e_val(SIZE(mtx,1)), e_vec_ct(SIZE(mtx, 1), SIZE(mtx, 2)))

		CALL DSYEV ('V', 'U', SIZE(mtx, 1), mtx, SIZE(e_val),&
						e_val, WORK, LWORK, INFO)
		CALL check_lapack (INFO)
		
		LWORK = INT(WORK(1))
		DEALLOCATE (WORK)
		ALLOCATE (WORK(LWORK))
		
		CALL DSYEV ('V', 'U', SIZE(mtx, 1), mtx, SIZE(e_val),&
					e_val, WORK, LWORK, INFO)	
		CALL check_lapack(INFO)
		
		e_vec_ct = TRANSPOSE(mtx)

		DO i=1, SIZE(mtx,1)
			mtx(:,i) = mtx(:,i)/ SQRT(e_val(i))
		END DO

		mtx = MATMUL(mtx, e_vec_ct)
	END SUBROUTINE sym_lowdin		

	FUNCTION calc_sparse(mtx)
		DOUBLE PRECISION :: mtx(:,:)
		INTEGER :: i = 0, j = 0, k = 0
		DOUBLE PRECISION :: calc_sparse
	
		DO i=1, SIZE(mtx, 1)
			DO j=1, SIZE(mtx, 2)
				IF (mtx(j,i) .eq. 0.0D0 ) k = k + 1
			END DO
		END DO

		calc_sparse = DBLE(k) / (DBLE(SIZE(mtx, 1)) * DBLE(SIZE(mtx, 2)))
		k = 0
	END FUNCTION calc_sparse

    SUBROUTINE check_lapack(info)
		INTEGER :: info
		IF (info .eq. 0) THEN
			RETURN
		ELSE IF (info .lt. 0) THEN
			WRITE (*,*) 'Illegal input value to LAPACK call'
			CALL EXIT(1)
		ELSE
			WRITE (*,*) 'LAPACK routine failed to converge'
			CALL EXIT(1)
		END IF
	END SUBROUTINE check_lapack	

	FUNCTION check_orthog(mtx_1, mtx_2, method)
		CHARACTER (len=1) :: method
		DOUBLE PRECISION :: mtx_1(:,:), mtx_2(:,:)
		DOUBLE PRECISION :: check_orthog(SIZE(mtx_1,1), SIZE(mtx_1,2))
	
		IF (method .eq. 'S') THEN
			check_orthog = MATMUL(MATMUL(mtx_1,mtx_1),mtx_2)
		ELSE IF (method .eq. 'G') THEN
			check_orthog = MATMUL(mtx_1,MATMUL(mtx_2,TRANSPOSE(mtx_1)))
		ELSE
			WRITE (*,*) 'Invalid argument supplied to check_orthog'
			CALL EXIT(1)
		END IF
	END FUNCTION

	SUBROUTINE trans_fock(fock, s_half, s_half_inv, C)
		DOUBLE PRECISION :: fock(:,:), s_half(:,:), s_half_inv(:,:), C(:)
		INTEGER :: INFO = 0, LWORK = -1
		INTEGER, ALLOCATABLE :: IPIV(:), WORK(:)

		ALLOCATE (IPIV(SIZE(s_half, 1)), WORK(1))

		CALL DGETRF(SIZE(s_half_inv,2), SIZE(s_half_inv,1), s_half_inv,&
				SIZE(s_half_inv,1), IPIV, INFO)
		CALL check_lapack(INFO)

		CALL DGETRI(SIZE(s_half_inv, 1), s_half_inv, SIZE(s_half_inv, 1),&
				IPIV, WORK, LWORK, INFO)
		CALL check_lapack(INFO)

		LWORK = WORK(1)
		DEALLOCATE(WORK)
		ALLOCATE (WORK(LWORK))

		CALL DGETRI(SIZE(s_half_inv,1), s_half_inv, SIZE(s_half_inv, 1),&
				IPIV, WORK, LWORK, INFO)

		fock = MATMUL(MATMUL(s_half, fock), s_half)
		C = MATMUL(s_half_inv, C)

	END SUBROUTINE trans_fock

	SUBROUTINE sfc_iter(F_mat, S_mat, C_mat, D_mat, J_mat, H_mat, N_el)
		DOUBLE PRECISION :: F_mat(:,:), S_mat(:,:), C_mat(:,:)
		DOUBLE PRECISION :: D_mat(:,:), J_mat (:,:,:,:), H_mat(:,:)
		DOUBLE PRECISION :: last_e = 1.0D10, current_e = 0.0D0
		DOUBLE PRECISION, ALLOCATABLE :: e_vals(:), C_new(:,:), WORK(:), S_work(:,:)
		INTEGER :: INFO = 0, LWORK = -1, i, j, k, l, N, N_el
		
		N = SIZE(F_mat, 1)

		ALLOCATE(e_vals(N), C_new(SIZE(C_mat, 1), SIZE(C_mat, 2)), WORK(1), S_work(N,N))
		C_new = F_mat
		S_work = S_mat

		! set up the LAPACK environment
		CALL DSYGV(1, 'V', 'L', SIZE(C_new, 1), C_new, SIZE(C_new, 1), S_work,&
			 	SIZE(S_work, 1), e_vals, WORK, LWORK, INFO)
		CALL check_lapack(INFO)

		LWORK = WORK(1)
		DEALLOCATE(WORK)
		ALLOCATE(WORK(LWORK))
		
		! start the convergence loop
		DO WHILE ( (last_e - current_e) .gt. 0.0D-5)
			
			!set the energies
			last_e = current_e

			!set up work matracies
			S_work = S_mat

			!Call lapack to get new values of c
			CALL DSYGV(1, 'V', 'L', SIZE(C_new, 1), C_new, SIZE(C_new, 1), S_work,&
				 	SIZE(S_work, 1), e_vals, WORK, LWORK, INFO)
			CALL check_lapack(INFO)

			! Recalculate the density matrix
			DO i=1, N_el/2
				DO j=1, N
					DO k=1, N
						D_mat(j,k) = 2.0D0 * C_new(j,i) * C_new(k,i)
       				END DO
			   	END DO
			END DO
			
			! Recalculate the fock Matrix
			DO i=1, N
				DO j=1, N
					DO k=1, N
						DO l=1, N
							C_new(i,j) = C_new(i,j) + D_mat(k,l) * (J_mat(i, k, j, l) + 0.5D0 * J_mat(i, k, l, j))
						END DO
					END DO
				END DO
			END DO
			C_new = C_new + H_mat
			
			! calculate the new energy
			current_e = hf_e(H_mat, D_mat, J_mat)
		END DO
		
		F_mat = C_new
	
	END SUBROUTINE sfc_iter

	FUNCTION hf_E(H_mat, D_mat, J_mat)
		DOUBLE PRECISION :: hf_E
		DOUBLE PRECISION :: H_mat(:,:), D_mat(:,:), J_mat(:,:,:,:)
		INTEGER :: i, j, k, l, N

		hf_e = 0.0D0
		N = SIZE(D_mat, 1)
		
		DO i = 1, N
			DO j = 1, N
				DO k = 1, N
					DO l = 1, N
						hf_e = hf_e + D_mat(i,k) * D_mat(j,l) * (J_mat(i,j,k,l) - 0.05D0 * J_mat(i,j,l,k))
					END DO
				END DO
			END DO
		END DO
		
		hf_e = 0.5D0 * hf_e

		DO i = 1, N
			DO j = 1, N
				hf_e = hf_e + D_mat(i,j) * H_mat(i,j)
			END DO
		END DO

	END FUNCTION hf_E

END MODULE num_alg
