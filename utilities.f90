! Module for stand alone utilities
MODULE MODUTIL 
  IMPLICIT NONE

  CONTAINS

  ! Set the seed for random numbers. 0 for random.
  SUBROUTINE SET_RANDOM_SEED(seed)
    INTEGER, INTENT(IN) :: seed
    INTEGER, ALLOCATABLE :: seedbuf(:)
    INTEGER :: n

    IF (SEED == 0) THEN
      CALL RANDOM_SEED()
    ELSE
      CALL RANDOM_SEED(SIZE=n)
      ALLOCATE(seedbuf(n))
      seedbuf(:) = seed
      CALL RANDOM_SEED(PUT=seedbuf)
    END IF
  END SUBROUTINE SET_RANDOM_SEED

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


  ! Read data from file into array.
  SUBROUTINE read_array(filename, data)
    CHARACTER(LEN=*),INTENT(IN) :: filename
    REAL*8, ALLOCATABLE :: data(:)
    INTEGER :: fd, io, ndata, i
    REAL*8 :: val

    ! Read file (first time to get ndata)
    OPEN(NEWUNIT=fd, FILE=TRIM(filename), STATUS='OLD', ACTION='READ')
    ndata = 0
    DO
      READ(fd, *, IOSTAT=io) val
      IF (io > 0) THEN
        WRITE(*,*) "Error reading the file: ", TRIM(filename), " !"
        STOP 1
      ELSE IF (io < 0) THEN
        EXIT
      ELSE
        ndata = ndata + 1
      END IF
    END DO
    CLOSE(fd)

    ALLOCATE(data(ndata))

    ! Read file (second time to get data)
    OPEN(NEWUNIT=fd, FILE=TRIM(filename), STATUS='OLD', ACTION='READ')
    DO i = 1, ndata
      READ(fd, *) data(i)
    END DO
    CLOSE(fd)

  END SUBROUTINE read_array

END MODULE MODUTIL
