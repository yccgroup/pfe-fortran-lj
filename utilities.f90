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
    aid = FLOOR(natoms*rand)+1
  
  END FUNCTION RANDOM_INTEGER
  
  
  ! Interpolate the data
  FUNCTION INTERPOLATE(levels,volumes,nlevel,Etarget) RESULT(Vtarget)
    IMPLICIT NONE
    INTEGER :: i
    INTEGER :: nlevel, levelid
    REAL*8 :: Etarget, Vtarget, slope, val, levels(nlevel), volumes(nlevel)

    levelid = 0
    ! find where enegy locates
    DO i = 1, nlevel-1
      val = (levels(i) - Etarget) * (levels(i+1) - Etarget)
      IF (val < 0) THEN
        levelid = i
        EXIT
      END IF
    END DO
  
    ! Check if levelid is correctly found
    IF (levelid == 0) THEN
      PRINT *, "ERROR, levelid = 0, interpolation failed!"
      STOP 1
    END IF

    ! interpolate
    slope = (volumes(levelid)-volumes(levelid+1))/(levels(levelid)-levels(levelid+1))
    Vtarget = volumes(levelid) + slope * (Etarget - levels(levelid))
  END FUNCTION INTERPOLATE

END MODULE MODUTIL
