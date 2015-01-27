MODULE atom_mod

	IMPLICIT NONE
	
	DOUBLE PRECISION, PRIVATE, PARAMETER :: PI = 4.0D0 * ATAN(1.0D0) ! Pi parameter

	TYPE atom
		DOUBLE PRECISION :: coords(3)
		CHARACTER (len=3) :: typ
		DOUBLE PRECISION :: alpha
	
		CONTAINS
		
		PROCEDURE :: overlap
		PROCEDURE :: ke_int	
		PROCEDURE :: nuc_int
		PROCEDURE :: J_int
		PROCEDURE :: K_int

	END TYPE atom

	CONTAINS

		FUNCTION overlap(this, atm_2, thresh)
			CLASS(atom) :: this, atm_2
			DOUBLE PRECISION :: overlap, pre_fac, exp_fac, norm, thresh
			DOUBLE PRECISION :: diff(3)
			
			pre_fac = (PI/(this%alpha + atm_2%alpha))**(3.0D0/2.0D0)
			norm = (PI/ (this%alpha + this%alpha))**(-3.0D0/4.0D0)

			exp_fac = ((-1.0D0 * this%alpha * atm_2%alpha) / &
						(this%alpha + atm_2%alpha))
			
			diff = this%coords - atm_2%coords
			
			overlap = norm * pre_fac * DEXP(exp_fac * DOT_PRODUCT(diff,diff))
			
			IF (overlap .lt. thresh) overlap = 0.0D0
		END FUNCTION

		FUNCTION ke_int (this, atm_2)
			CLASS(atom) :: this, atm_2
			DOUBLE PRECISION :: ke_int, norm_1, norm_2

			norm_1 = (PI/(this%alpha + this%alpha))**(-3.0D0/4.0D0)
			norm_2 = (PI/(atm_2%alpha + atm_2%alpha))**(-3.0D0/4.0D0)

			ke_int = norm_1 * norm_2 * (3.0D0 * this%alpha * atm_2%alpha) /&
						 (this%alpha + atm_2%alpha)**(5.0D0/2.0D0)&
								 * PI**(3.0D0/2.0D0)

		END FUNCTION ke_int

		FUNCTION nuc_int (this, atm_2)
			class(atom) :: this, atm_2
			DOUBLE PRECISION :: nuc_int, norm_1, norm_2

			norm_1 = (PI/(this%alpha + this%alpha))**(-3.0D0/4.0D0)
			norm_2 = (PI/(atm_2%alpha + atm_2%alpha))**(-3.0D0/4.0D0)

			nuc_int = norm_1 * norm_2 * (-2.0D0 * PI * 1.0D0) / (this%alpha + atm_2%alpha)
		
		END FUNCTION nuc_int

		FUNCTION J_int(this, atm_2, atm_3, atm_4)
			CLASS(atom) :: this, atm_2, atm_3, atm_4
			DOUBLE PRECISION :: J_int, norm_1, norm_2, norm_3, norm_4

			norm_1 = (PI/(this%alpha + this%alpha))**(-3.0D0/4.0D0)
			norm_2 = (PI/(atm_2%alpha + atm_2%alpha))**(-3.0D0/4.0D0)
			norm_3 = (PI/(atm_3%alpha + atm_3%alpha))**(-3.0D0/4.0D0)
			norm_4 = (PI/(atm_4%alpha + atm_4%alpha))**(-3.0D0/4.0D0)

			J_int = norm_1 * norm_2 * norm_3 * norm_4 *&
						(2.0D0 * PI)**(5.0D0/2.0D0) / ((this%alpha + atm_2%alpha) *&
						(atm_3%alpha + atm_4%alpha) *&
						(this%alpha + atm_2%alpha + atm_3%alpha + atm_4%alpha)**(1.0D0/2.0D0))
	
		END FUNCTION J_int		

		FUNCTION K_int(this, atm_2, atm_3, atm_4)
			CLASS(atom) :: this, atm_2, atm_3, atm_4
			DOUBLE PRECISION :: K_int, norm_1, norm_2, norm_3, norm_4

			norm_1 = (PI/(this%alpha + this%alpha))**(-3.0D0/4.0D0)
			norm_2 = (PI/(atm_2%alpha + atm_2%alpha))**(-3.0D0/4.0D0)
			norm_3 = (PI/(atm_3%alpha + atm_3%alpha))**(-3.0D0/4.0D0)
			norm_4 = (PI/(atm_4%alpha + atm_4%alpha))**(-3.0D0/4.0D0)


			K_int = norm_1 * norm_2 * norm_3 * norm_4 *&
						 (2.0D0 * PI)**(5.0D0/2.0D0) / ((this%alpha + atm_2%alpha) *&
						(atm_3%alpha + atm_4%alpha) *&
						(this%alpha + atm_2%alpha + atm_3%alpha + atm_4%alpha)**(1.0D0/2.0D0))
	
		END FUNCTION K_int

END MODULE atom_mod
