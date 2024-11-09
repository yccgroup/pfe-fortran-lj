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
  INTEGER :: seed, nbin, edim, io
  REAL*8 :: temperature, beta, kBT, r2, pot, val 
  REAL*8 :: lambda_th, lnZ, lnQ, lnP
  REAL*8, Parameter :: dE=0.01
  REAL*8, ALLOCATABLE :: Energy(:)
  LOGICAL :: flag
  CHARACTER(LEN=10) :: Job
  CHARACTER(LEN=80) :: MCfilename, energyfilename

  ! Read inputs from file input.dat
  OPEN(UNIT=10,FILE="input.dat",STATUS="OLD",ACTION="READ")
  
  ! temperature
  READ(10,*) temperature

  ! system parameters
  CALL System%rdinp(10)

  ! random seed
  READ(10,*) seed

  ! equilibrium MC parameters
  CALL MCeq%rdinp(10)

  ! production run MC parameters
  CALL MCrun%rdinp(10)

  ! filename for MC data
  READ(10,*) MCfilename
  READ(10,*) energyfilename
  READ(10,*) flag
  READ(10,*) nbin

  ! pfe parameters
  CALL Parfu%rdinp(10)

  ! job type
  READ(10,*) Job

  CLOSE(10)

  ! Seed the random number
  CALL SET_RANDOM_SEED(seed)

  ! Check input parameters
  IF (System%rc > 0.5d0*System%L) THEN
    WRITE(*,*) "VDW cutoff (rc) should be smaller than half of the boxsize. Check your inputs."
    STOP
  END IF


  ! Initialize parameters
  kBT = kB*temperature
  beta = 1.d0/kBT
  lambda_th = (h/SQRT(2*pi*(System%mass*amu2kg)*(kBT*1000*cal2joule/NA))) * 1E10
  lnP = -LOG_GAMMA((System%natoms+1)*1d0) - 3*System%natoms*LOG(lambda_th)
  PRINT *, "lambda_th:", lambda_th

  ! Build LJ system and calculate the system energy
  CALL System%init()

  IF (Job == 'MC') THEN

    PRINT *, "Letitia: begin MD"
    FLUSH(6)

    ! Pre-equilibration
    CALL MCeq%Sampling(System,beta)

    ! Monte Carlo Sampling
    CALL MCrun%Sampling(System,beta)
    
    ! Save MC data
    CALL MCrun%Write(System,beta,TRIM(MCfilename))

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


  ELSE IF (Job == 'POT') THEN

    ! Print out the potential energy function
    OPEN(UNIT=20,FILE="pot.dat",STATUS="UNKNOWN")
    WRITE(20,*) "# distance (nm)  potential (kcal/mol)"
    DO i=1,300
      r2 = (0.1*i)**2
      pot = System%calcpairpot(r2)
      System%XYZ = 0
      System%XYZ(1,2) = 0.1*i
      CALL System%calcenergy()
      WRITE(20,*) 0.1*i/10, pot, System%energy
    END DO  
    CLOSE(20)

  ELSE ! Job /= MC

    Parfu%Method = Job

    !! Read MC data
    !CALL MCrun%Read(System,beta,TRIM(MCfilename))

    ! Read energy file (first time to get edim)
    OPEN(UNIT=20, FILE=energyfilename, STATUS='OLD', ACTION='READ')
    edim = 0
    DO
      READ(20, *, IOSTAT=io) val
      IF (io > 0) THEN
        WRITE(*,*) "Error in reading the energy file:", energyfilename, "!"
      ELSE IF (io < 0) THEN
        EXIT
      ELSE
        edim = edim + 1
      END IF
    END DO
    CLOSE(20)

    ALLOCATE(Energy(edim))

    ! Read energy file (second time to get energy data)
    OPEN(UNIT=20, FILE=energyfilename, STATUS='OLD', ACTION='READ')
    DO i = 1, edim
      READ(20, *) Energy(i)
    END DO
    CLOSE(20)

    ! Unit kJ/mol to kcal/mol
    Energy = Energy / cal2joule 


    ! Partition Function for Natoms = 1
    IF (System%natoms == 1) THEN
      lnQ = 3*LOG(System%L) 
      lnZ = lnQ + lnP 
      PRINT *, "lnQ = ", lnQ
      PRINT *, "lnP = ", lnP
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
      CALL Parfu%NSVolume(System,System%Energy,.FALSE.)
      CALL Parfu%NSPartition(System,beta)
      lnQ = Parfu%lnQ
      lnZ = lnQ + lnP
      PRINT *, "NS lnQ = ", lnQ
      PRINT *, "NS lnP = ", lnP
      PRINT *, "NS lnZ = ", lnZ
  
    ELSE IF (Parfu%Method == "PFE") THEN
      
      ! Calculate partition function by PFE (our theory)
      PRINT *, "Letitia: calculate partition function via PFE"
      FLUSH(6)
      IF (flag) THEN
        PRINT *, "call PartFunc3"
        CALL Parfu%PartFunc3(System,Energy,beta,nbin)
      ELSE
        CALL Parfu%PartFunc(System,Energy,beta)
      END IF
      lnQ = Parfu%lnQ
      lnZ = lnQ + lnP
      PRINT *, "PFE lnQ = ", lnQ
      PRINT *, "PFE lnP = ", lnP
      PRINT *, "PFE lnZ = ", lnZ
      PRINT *, "PFE Err VErr Err+VErr = ", SQRT(Parfu%Err2), SQRT(Parfu%VErr2), SQRT(Parfu%Err2+Parfu%VErr2)
      PRINT *, "Estar = ", Parfu%Estar*cal2joule, "kJ/mol = ", Parfu%Estar, "kcal/mol"
      PRINT *, "Cutoff % = ", Parfu%percentage*100
  
      ! DEBUG for Natoms == 2 only
      IF (System%natoms == 2) CALL DEBUG(System,beta,lnP)
      
    END IF ! Method

  END IF ! Job

END PROGRAM Main

