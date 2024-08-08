MODULE MODDEBUG

  USE MODCONST
  USE MODLJ
  IMPLICIT NONE

  CONTAINS

  SUBROUTINE DEBUG(System,beta,lnZp)
    INTEGER :: i
    INTEGER :: ndim
    REAL*8 :: beta, lnZp
    REAL*8 :: Z, Eavg
    REAL*8 :: r, dr, pot, prob, integral
    CLASS(LJ) :: System
  
    ndim = 10000
    dr = System%rc / ndim
    Eavg = 0.d0
    integral = 0.d0

    DO i=1,ndim
      r = (i-1)*dr
      ! do not calculate the sum if the distance is too small
      IF (r < 0.1d0) THEN
        CONTINUE
      ELSE
        pot = System%calcpairpot(r*r)
        prob = EXP(-beta*pot)
        integral = integral + prob*(4*pi*r*r)*dr
        Eavg = Eavg + prob*pot*(4*pi*r*r)*dr
      END IF
    END DO
    Z = System%L**3 * (integral + System%L**3 - 4*pi/3.d0 * System%rc**3)
    Eavg = System%L**3 * Eavg / Z
    PRINT *, "DEBUG lnZ = ", LOG(Z)+lnZp
  
  END SUBROUTINE DEBUG

END MODULE MODDEBUG
