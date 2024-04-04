PROGRAM Main
  USE MODCONST
  USE MODLJ
  USE MODMC
  USE MODPFE
  USE MODDEBUG
  IMPLICIT NONE

  TYPE(LJ) :: System
  TYPE(MC) :: MCeq, MCrun
  TYPE(PFE) :: Parfu
  INTEGER :: i, j, k
  REAL*8 :: temperature, beta, kBT
  REAL*8 :: cutoff

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

  ! cutoff
  READ(10,*) cutoff

  ! pfe parameters
  CALL Parfu%rdinp(10)

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

  ! Remove the sampling outlier
  cutoff = cutoff*kBT + MINVAL(MCrun%Energy)
  CALL MCrun%Truncate(cutoff)
  PRINT *, "Average energy after truncation (kJ/mol) =", SUM(MCrun%Energy)*cal2joule/MCrun%outdim
  PRINT*, "cutoff = ", cutoff

  ! RAFEP
  CALL Parfu%PartFunc(System,MCrun%Energy,beta)
  print *, "RAFEP Zest:", EXP(Parfu%lnZ)
  
  ! For 2 atoms debug
  IF (System%natoms == 2) CALL DEBUG(System,beta,EXP(Parfu%lnZ))

END PROGRAM Main

