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
  REAL*8 :: lambda_th, lnZ, lnQ, lnZp
  REAL*8, Parameter :: dE=0.01

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
  lnZp = -LOG(Gamma((System%natoms+1)*1d0)) - 3*System%natoms*LOG(lambda_th)
  PRINT *, "lambda_th:", lambda_th

  ! Build LJ system and calculate the system energy
  CALL System%init()
  CALL System%calcenergy()

  PRINT *, "Letitia: begin MD"
  FLUSH(6)

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

  ! Partition Function for Natoms = 1
  IF (System%natoms == 1) THEN
    lnQ = 3*LOG(System%L) 
    lnZ = lnQ + lnZp 
    PRINT *, "lnQ = ", lnQ
    PRINT *, "lnZp = ", lnZp
    PRINT *, "lnZ = ", lnZ
    STOP
  END IF

  ! Partition Function for Natoms /= 1
  IF (Parfu%Method == "Exact") THEN
    
    ! Calculate partition function by Nested Sampling (Time consuming)
    PRINT *, "Letitia: calculate partition function via NS"
    FLUSH(6)
    ! Determine Estar and Emin for Exact method
    CALL MCrun%Minimize(System) 
    CALL System%calcenergy()
    Parfu%Estar = System%Energy + kBT*0.1
    CALL Parfu%NSVolume(System,System%Energy)
    CALL Parfu%NSPartition(System,beta)
    lnQ = Parfu%lnZ
    lnZ = lnQ + lnZp
    PRINT *, "NS lnQ = ", lnQ
    PRINT *, "NS lnZp = ", lnZp
    PRINT *, "NS lnZ = ", lnZ

  ELSE IF (Parfu%Method == "PFE") THEN
    
    ! Calculate partition function by PFE (our theory)
    PRINT *, "Letitia: calculate partition function via PFE"
    FLUSH(6)
    CALL Parfu%PartFunc(System,MCrun%Energy,beta)
    lnQ = Parfu%lnZ
    lnZ = lnQ + lnZp
    PRINT *, "NS lnQ = ", lnQ
    PRINT *, "NS lnZp = ", lnZp
    PRINT *, "NS lnZ = ", lnZ
    PRINT *, "NS Err = ", SQRT(Parfu%Err2)
    PRINT *, "Estar = ", Parfu%Estar
    PRINT *, "Cutoff % = ", Parfu%percentage*100

    ! DEBUG for Natoms == 2 only
    IF (System%natoms == 2) CALL DEBUG(System,beta,lnZp)
    
  END IF


END PROGRAM Main

