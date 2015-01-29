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

	END TYPE atom

	CONTAINS

		FUNCTION overlap(this, atm_2, thresh)
			CLASS(atom) :: this, atm_2
			DOUBLE PRECISION :: overlap, pre_fac, exp_fac, norm_1, norm_2, thresh
			DOUBLE PRECISION :: diff(3)
			
			pre_fac = (PI/(this%alpha + atm_2%alpha))**(3.0D0/2.0D0)
			norm_1 = (PI/ (this%alpha + this%alpha))**(-3.0D0/4.0D0)
			norm_2 = (PI/ (atm_2%alpha + atm_2%alpha))**(-3.0D0/4.0D0)

			exp_fac = ((-1.0D0 * this%alpha * atm_2%alpha) / &
						(this%alpha + atm_2%alpha))
			
			diff = this%coords - atm_2%coords
			
			overlap = norm_1 * norm_2 * pre_fac * DEXP(exp_fac * DOT_PRODUCT(diff,diff))
			
			IF (overlap .lt. thresh) overlap = 0.0D0
		END FUNCTION

		FUNCTION ke_int (this, atm_2)
			CLASS(atom) :: this, atm_2
			DOUBLE PRECISION :: ke_int, norm_1, norm_2

			norm_1 = (PI/(this%alpha + this%alpha))**(-3.0D0/4.0D0)
			norm_2 = (PI/(atm_2%alpha + atm_2%alpha))**(-3.0D0/4.0D0)

			ke_int = norm_1 * norm_2 * ((3.0D0 * (PI)**(1.5D0) * this%alpha * atm_2%alpha))/&
				((this%alpha + atm_2%alpha)**(2.5D0))

		END FUNCTION ke_int

		FUNCTION nuc_int (this, atm_2, charge)
			class(atom) :: this, atm_2
			DOUBLE PRECISION :: nuc_int, norm_1, norm_2, charge

			norm_1 = (PI/(this%alpha + this%alpha))**(-3.0D0/4.0D0)
			norm_2 = (PI/(atm_2%alpha + atm_2%alpha))**(-3.0D0/4.0D0)

			nuc_int = norm_1 * norm_2 * (-2.0D0 * PI * charge) / (this%alpha + atm_2%alpha)
		
		END FUNCTION nuc_int

		FUNCTION J_int(this, atm_2, atm_3, atm_4)
			CLASS(atom) :: this, atm_2, atm_3, atm_4
			DOUBLE PRECISION :: J_int, norm_1, norm_2, norm_3, norm_4, norm_tot

			norm_1 = ((this%alpha + this%alpha)/PI)**(0.75D0)
			norm_2 = ((atm_2%alpha + atm_2%alpha)/PI)**(0.75D0)
			norm_3 = ((atm_3%alpha + atm_3%alpha)/PI)**(0.75D0)
			norm_4 = ((atm_4%alpha + atm_4%alpha)/PI)**(0.75D0)
			norm_tot = norm_1 * norm_2 * norm_3 * norm_4			
			
			J_int = norm_tot * 2.0D0 * (PI)**(2.5D0)/&
				((this%alpha+atm_3%alpha)*(atm_2%alpha+atm_4%alpha)*&
				(((this%alpha+atm_2%alpha+atm_3%alpha+atm_4%alpha)))**(1.0D0/2.0D0))
			 
		END FUNCTION J_int		

END MODULE atom_mod
