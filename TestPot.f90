! Program to test if our potential function is correct,
! by comparing with OpenMM.

PROGRAM TestPot

  USE MODCONST
  USE MODOpenMM
  USE MODLJ
  IMPLICIT NONE

  ! hardcoded LJ parameters, OpenMM defaults
  REAL*8,PARAMETER :: mass    = 39.9d0
  REAL*8,PARAMETER :: epsilom = 0.238d0
  REAL*8,PARAMETER :: sigma   = 3.4d0
  REAL*8,PARAMETER :: rc      = 10.2d0
  REAL*8,PARAMETER :: L       = 25.0d0

  CHARACTER(LEN=80) :: dcdfilename, logfilename
  TYPE(DCDFile) :: dcd
  TYPE(LJ) :: system
  REAL*8,ALLOCATABLE :: xyz(:,:)
  REAL*8,ALLOCATABLE :: energies(:)
  INTEGER :: log_istrt, log_nsavc
  INTEGER :: i
  REAL*8 :: myenergy

  ! get command-line arguments
  CALL GET_COMMAND_ARGUMENT(1, dcdfilename)
  CALL GET_COMMAND_ARGUMENT(2, logfilename)

  ! read info from DCD file
  dcd = dcd_open(dcdfilename)
  CALL dcd%info(6)

  ! initialize LJ system
  system%natoms = dcd%natom
  system%mass = mass
  system%epsilom = epsilom
  system%sigma = sigma
  system%rc = rc
  system%L = L
  CALL system%init()

  ! space for arrays
  ALLOCATE(energies(dcd%nset))
  ALLOCATE(xyz(3,dcd%natom))

  ! read energies from the OpenMM log file
  CALL read_openmm_log(logfilename, dcd%nset, log_istrt, log_nsavc, energies)
  IF (log_istrt /= dcd%istrt) THEN
    PRINT *, 'ERROR: inconsistent istrt ', dcd%istrt, ' != ', log_istrt
    STOP 1
  END IF
  IF (log_nsavc /= dcd%nsavc) THEN
    PRINT *, 'ERROR: inconsistent nsavc ', dcd%nsavc, ' != ', log_nsavc
    STOP 1
  END IF

  ! read the DCD frame by frame, and output energies
  WRITE(6,*)
  WRITE(6,"(a8,4(a20,'     '))") '# Step', 'My_Energy', 'OpenMM_Energy', 'Delta', 'Ratio'
  DO i = 1, dcd%nset
    CALL dcd%read_frame(xyz)
    ! copy coordinates into system
    system%xyz = xyz
    system%uptodate = .FALSE.
    ! calculate energy, and fix unit
    myenergy = system%getenergy() * cal2joule
    WRITE(6,'(i8,4g25.15)') i, myenergy, energies(i), myenergy-energies(i), myenergy/energies(i)
  END DO

  ! clean up
  CALL dcd%close()
  DEALLOCATE(xyz)
  DEALLOCATE(energies)

END PROGRAM TestPot
