PROGRAM Main
  USE MODCONST
  USE MODLJ
  USE MODMC
  IMPLICIT NONE

  !REAL*8 :: Z, lnZ
  !REAL*8 :: r, dr, pot, prob, expectation, restV
  !TYPE(NSRAFEPPAR) :: NS

  TYPE(LJ) :: System
  TYPE(MC) :: MCeq, MCrun
  INTEGER :: i, j, k
  REAL*8 :: temperature, beta, kBT


  ! Read inputs from file input.dat
  OPEN(UNIT=10,FILE="input.dat",STATUS="OLD",ACTION="READ")
  
  ! temperature
  READ(10,*) temperature

  ! system parameters
  CALL System%rdinp(10)

  ! equilibrium MC parameters
  CALL MCeq%rdinp(10)

  ! production run MC parameters
  CALL MCrun%rdinp(10)

  !! Inputs of RAFEP and NS parameters
  !READ(10,*) NS%rafep_cutoff
  !READ(10,*) NS%nsamples
  !READ(10,*) NS%nsteps
  !READ(10,*) NS%stepsize
  !READ(10,*) NS%fractiom
  !READ(10,*) NS%root_energy

  CLOSE(10)

  ! Seed the random number
  CALL RANDOM_SEED()

  ! Check input parameters
  IF (System%rc > 0.5d0*System%L) THEN
    WRITE(*,*) "VDW cutoff (rc) should be smaller than half of the boxsize. Check your inputs."
    STOP
  END IF


  ! Initialize parameters
  kBT = kB*temperature
  beta = 1.d0/kBT
  MCeq%outdim = MCeq%nsteps/MCeq%outfreq
  MCrun%outdim = MCrun%nsteps/MCrun%outfreq

  ! Build LJ system and calculate the system energy
  CALL System%init()
  CALL System%calcenergy()

  ! Pre-equilibration
  CALL MCeq%Sampling(System,beta)

  ! Monte Carlo Sampling
  CALL MCrun%Sampling(System,beta)
  
  ! Output MC Trajectory
  OPEN(UNIT=20,FILE="Traj.dat",STATUS="UNKNOWN")
  DO i=1,MCrun%outdim
    WRITE(20,*) "Time = ", i*MCrun%outfreq
    DO j=1,System%natoms
      WRITE(20,*) (MCrun%Traj(k,j,i), k=1,3)
    END DO
  END DO
  CLOSE(20)

  ! Ouput MC Energy (unit: kj/mol)
  OPEN(UNIT=20,FILE="Energy.dat",STATUS="UNKNOWN")
  DO i=1,MCrun%outdim
    WRITE(20,*) MCrun%Energy(i)*cal2joule
  END DO
  CLOSE(20)
  
  ! Output MC Statistics (unit: kj/mol)
  PRINT *, "Acceptance rate =", MCrun%Accept
  PRINT *, "Average energy (kJ/mol) =", SUM(MCrun%Energy)*cal2joule/MCrun%outdim

  ! For 2 atoms debug
  !IF (System%natoms == 2) CALL DEBUG(System,beta)

  ! Remove the sampling outlier
  !NS%rafep_cutoff = NS%rafep_cutoff*kBT + MINVAL(MC%Energy)
  !CALL Truncate(MC,System%natoms,NS%rafep_cutoff)
  
  ! RAFEP
  !Z = Partition_RAFEP(System,NS,MC%NewEnergy,beta)
  !print *, "RAFEP Zest:", Z
  
STOP
END PROGRAM Main

