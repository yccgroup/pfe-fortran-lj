! Module for stand alone utilities
MODULE MODUTIL 
  IMPLICIT NONE

  CONTAINS

  ! Randomly select an atom from the system
  FUNCTION RANDOM_INTEGER(natoms) RESULT(aid)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: natoms
    INTEGER :: aid
    REAL*8 :: rand
  
    CALL RANDOM_NUMBER(rand)
    aid = FLOOR((natoms+1)*rand)
  
  END FUNCTION RANDOM_INTEGER
  
  
  ! Interpolate the data
  FUNCTION INTERPOLATE(Lenergy,LV0,nlevel,desire_energy) RESULT(V0)
    IMPLICIT NONE
    INTEGER :: i
    INTEGER :: nlevel, levelid
    REAL*8 :: desire_energy, V0, slope, val, Lenergy(nlevel), LV0(nlevel)
  
    ! find where enegy locates
    DO i = 1, nlevel-1
      val = (Lenergy(i) - desire_energy) * (Lenergy(i+1) - desire_energy)
      IF (val < 0) THEN
        levelid = i
        EXIT
      END IF
    END DO
  
    ! interpolate
    slope = (LV0(levelid)-LV0(levelid+1))/(Lenergy(levelid)-Lenergy(levelid+1))
    V0 = LV0(levelid) + slope * (desire_energy - Lenergy(levelid))
  END FUNCTION INTERPOLATE

END MODULE MODUTIL
