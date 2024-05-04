PROGRAM Main
  USE MODCONST
  USE MODLJ
  USE MODMC
  USE MODPFE
  USE MODDEBUG
  IMPLICIT NONE

  TYPE(LJ) :: System
  TYPE(MC) :: MCeq, MCrun
  TYPE(PFE) :: Parfu, Parfu2
  INTEGER :: i, j, k
  REAL*8 :: temperature, beta, kBT
  REAL*8 :: cutoff, lambda_th, prefactor

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
  lambda_th = (h/SQRT(2*pi*(System%mass*amu2kg)*(kBT*1000*cal2joule/NA))) * 1E10
  prefactor = -LOG(Gamma((System%natoms+1)*1d0)) - 3*System%natoms*LOG(lambda_th)
  PRINT *, "lambda_th:", lambda_th

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
  !cutoff = cutoff*kBT + MINVAL(MCrun%Energy)
  cutoff = cutoff*kBT + SUM(MCrun%Energy)/SIZE(MCrun%Energy)
  PRINT *, "Before truncation size = ", MCrun%outdim, SIZE(MCrun%Energy)
  CALL MCrun%Truncate(cutoff)
  PRINT *, "After truncation size = ", MCrun%outdim, SIZE(MCrun%Energy)
  PRINT *, "Average energy after truncation (kJ/mol) =", SUM(MCrun%Energy)*cal2joule/MCrun%outdim
  PRINT*, "cutoff = ", cutoff

  ! Partition Function for Natoms = 1
  IF (System%natoms == 1) THEN
    Parfu%lnZ = 3*LOG(System%L) + prefactor
    PRINT *, "ln(Z) = ", Parfu%lnZ
    STOP
  END IF

  ! Partition Function for Natoms /= 1
  ! RAFEP
  CALL Parfu%PartFunc(System,MCrun%Energy,beta)
  Parfu%lnZ = Parfu%lnZ + prefactor
  PRINT *, "RAFEP ln(Zest) = ", Parfu%lnZ, 3*System%natoms*LOG(lambda_th), 3*System%natoms*LOG(System%L) 

  ! NS (Emin determined by MC relaxation)
  CALL Parfu2%NSPartinit(System,MCrun,Parfu,kBT)
  CALL Parfu2%NSVolume(System)
  CALL Parfu2%NSPartition(System,beta)
  Parfu2%lnZ = Parfu2%lnZ + prefactor
  PRINT *, "NS ln(Znum) = ", Parfu2%lnZ

  ! DEBUG for Natoms == 2 only
  IF (System%natoms == 2) CALL DEBUG(System,beta,EXP(Parfu%lnZ),prefactor)

END PROGRAM Main

