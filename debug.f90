MODULE MODDEBUG

  USE MODTEMPAR
  USE MODPBCLJ
  IMPLICIT NONE

  CONTAINS

  SUBROUTINE DEBUG(System)
    INTEGER :: i
    REAL*8 :: Z
    REAL*8 :: r, dr, pot, prob, expectation, restV
    REAL*8, PARAMETER :: cal2joule = 4.18400
    CLASS(LJ) :: System
  
    Z = 0.d0
    expectation = 0.d0
    dr = System%rc / 10000
    DO i=1,10000
      r = (i-1)*dr
      ! do not calculate the sum if the distance is too small
      IF (r < 0.1d0) THEN
        CONTINUE
      ELSE
        pot = System%calcpairpot(r*r)
        prob = EXP(-beta*pot)
        Z = Z + prob*(4*3.1415926536*r*r)*dr
        expectation = expectation + prob*pot*(4*3.1415926536*r*r)*dr
      END IF
    END DO
    restV = System%L**3 - 4*3.1415926536/3.d0 * System%rc**3
    Z = Z + restV
    PRINT *, "Energy expectation for 2 atoms (kJ/mol) =", expectation*(System%epsilom*cal2joule)/Z
  
  END SUBROUTINE DEBUG

END MODULE MODDEBUG
