MODULE atom_mod

	IMPLICIT NONE
	
	DOUBLE PRECISION, PRIVATE, PARAMETER :: PI = 4.0D0 * ATAN(1.0D0) ! Pi parameter

	!####################################################################
	!	USER DECLARED TYPE: atom									    !
	!	Variables:													    !
	!		coords - array to store the coordinates of centre of atom   !
	!		typ - character array to store atomic symbol			    !
	!		alpha - variable to store the exponent of the guassian	    !
	!	Procedures:														!
	!		overlap - claculates the overlap integral of two atoms		!
	!####################################################################
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

		!####################################################################
		!	FUNCTION: overlap											    !
		!	Arguments:													    !
		!		this - implicitly passed object since type-bound procedure  !
		!		atm_2 - atom to calculate overlap with					    !
		!		thresh - threshold value, below this overlap = 0	 	    !
		!	Method:															!
		!		See report for method of calculating overlap intergral		! 
		!####################################################################
		FUNCTION overlap(this, atm_2, thresh)
			CLASS(atom) :: this, atm_2
			DOUBLE PRECISION :: overlap, pre_fac, exp_fac, norm, thresh
			DOUBLE PRECISION :: diff(3)
			
			pre_fac = (PI/(this%alpha + atm_2%alpha))**(3.0D0/2.0D0)
			norm = 1.0D0  !norm = pre_fac**(-1.0D0)
			exp_fac = ((-1.0D0 * this%alpha * atm_2%alpha) / &
						(this%alpha + atm_2%alpha))
			
			diff = this%coords - atm_2%coords
			
			overlap = norm * pre_fac * DEXP(exp_fac * DOT_PRODUCT(diff,diff))
			
			IF (overlap .lt. thresh) overlap = 0.0D0
		END FUNCTION




		FUNCTION ke_int (this, atm_2)
			CLASS(atom) :: this, atm_2
			DOUBLE PRECISION :: ke_int 

			ke_int = (3.0D0 * this%alpha * atm_2%alpha) /&
						 (this%alpha + atm_2%alpha)**(5.0D0/2.0D0)&
								 * PI**(3.0D0/2.0D0)
		
			ke_int = ke_int + 0.0D0

		END FUNCTION ke_int

		FUNCTION nuc_int (this, atm_2)
			class(atom) :: this, atm_2
			DOUBLE PRECISION :: nuc_int

			nuc_int = (-2.0D0 * PI * 1.0D0) / (this%alpha + atm_2%alpha)
		
		END FUNCTION nuc_int

		FUNCTION J_int(this, atm_2, atm_3, atm_4)
			CLASS(atom) :: this, atm_2, atm_3, atm_4
			DOUBLE PRECISION :: J_int

			J_int = (2.0D0 * PI)**(5.0D0/2.0D0) / ((this%alpha + atm_2%alpha) *&
						(atm_3%alpha + atm_4%alpha) *&
						(this%alpha + atm_2%alpha + atm_3%alpha + atm_4%alpha)**(1.0D0/2.0D0))
	
		END FUNCTION J_int		

		FUNCTION K_int(this, atm_2, atm_3, atm_4)
			CLASS(atom) :: this, atm_2, atm_3, atm_4
			DOUBLE PRECISION :: K_int

			K_int = (2.0D0 * PI)**(5.0D0/2.0D0) / ((this%alpha + atm_2%alpha) *&
						(atm_3%alpha + atm_4%alpha) *&
						(this%alpha + atm_2%alpha + atm_3%alpha + atm_4%alpha)**(1.0D0/2.0D0))
	
		END FUNCTION K_int

END MODULE atom_mod
