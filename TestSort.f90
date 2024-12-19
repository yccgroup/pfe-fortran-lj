PROGRAM TestSort

  USE MODUTIL
  USE MODSORT
  IMPLICIT NONE

  REAL*8,ALLOCATABLE :: Energy(:)
  INTEGER :: i
  REAL*8  :: prev
  LOGICAL :: check

  CALL read_array("Energy.dat", Energy)
  CALL qsort(Energy)

  check = .TRUE.
  prev = -1.d300
  DO i = 1, SIZE(Energy)
    WRITE(*,*) Energy(i)
    IF (Energy(i) < prev) THEN
      WRITE(*,*) "CHECK FAILED"
      STOP 1
    END IF
    prev = Energy(i)
  END DO

END PROGRAM TestSort
