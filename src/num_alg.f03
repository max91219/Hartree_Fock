MODULE num_alg

IMPLICIT NONE

CONTAINS
	
    !####################################################################
    !	SUBROUTINE: gram_schmidt										!
    !	Arguments:														!
    !		mtx - input matrix of overlaps								!
    !	Method:															!
    !		See report for method of gram schmidth orthogonalisation	!
    !####################################################################
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

    !####################################################################
    !	SUBROUTINE: sym_lowdin											!
    !	Arguments:														!
    !		mtx - input matrix of overlaps								!
    !	Method:															!
    !		See report for method of symetrix lowdin orthogonalisation	!
    !####################################################################
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

    !####################################################################
    !	FUNCTION: calc_sparse											!
    !	Arguments:														!
    !		mtx - input matrix for sparsity calculation					!
    !	Method:															!
    !		loops over all elements of matrix and keeps running			!
    !		total of how many are zero. returns the number				!
    !		of zero elements devided by the total number of				!
    !		elements													!
    !####################################################################
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

    !####################################################################
    !	SUBROUTINE: check_lapack										!
    !	Arguments:														!
    !		info - return value of lapack call to check					!
    !	Method:															!
    !		checks the value of info against known lapack error			!
    !		values to make sure that the call finished approperately	!
    !		if an error occured then it deals with it approprately		!
	!				info > 0 - failed to converge						!
	!				info < 0 - Illegal input							!
	!				info = 0 - call executed fine						!
    !####################################################################
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

    !####################################################################
    !	FUNCTION: check_orthog											!
    !	Arguments:														!
    !		mtx_1 - Orthogonal transformation Matrix					!
    !		mtx_2 - Overlap Matrix										!
	!		method - tells function which orthogonalisation method		!
	!				 to check											!
	!					S ====> Symmetric Lowdin						!
	!					G ====> Gram-Schmidt							!
    !	Method:															!
    !		Performs the orthogonality check detialed in the report		!
	!		for each of the methods used								!
    !####################################################################
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

	SUBROUTINE H_solve(H, S)
		DOUBLE PRECISION :: H(:,:), S(:,:)
		DOUBLE PRECISION, ALLOCATABLE :: WORK(:), e_val(:)
	    INTEGER :: LWORK = -1, INFO = 0, i = 0
	
		ALLOCATE(WORK(1), e_val(SIZE(H,1)))

		CALL DSYGV(1, 'V', 'L', SIZE(H,1), H, SIZE(H,1), S, SIZE(S,1), e_val,&
					WORK, LWORK, INFO)
		CALL check_lapack(INFO)

		LWORK = WORK(1)
		DEALLOCATE(WORK)
		ALLOCATE(WORK(LWORK))

		CALL DSYGV(1, 'V', 'L', SIZE(H,1), H, SIZE(H,1), S, SIZE(S,1), e_val,&
					WORK, LWORK, INFO)
		CALL check_lapack(INFO)	

		WRITE (*,*) (e_val(i), i=1, SIZE(e_val))

	END SUBROUTINE

END MODULE num_alg
