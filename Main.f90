PROGRAM Main
  USE MODTEMPAR
  USE MODMC
  USE MODNSRAFEP
  USE MODPBCLJ
  IMPLICIT NONE

  INTEGER :: i, j, k
  REAL*8 :: Z, lnZ
  REAL*8 :: r, dr, pot, prob, expectation, restV
  REAL*8, PARAMETER :: cal2joule = 4.18400
  TYPE(LJ) :: System
  TYPE(MCPAR) :: MCeq, MC
  TYPE(NSRAFEPPAR) :: NS


  ! Read inputs from file input.dat
  OPEN(UNIT=10,FILE="input.dat",STATUS="OLD",ACTION="READ")

    ! Inputs of System parameters
    READ(10,*) System%natoms
    READ(10,*) System%mass
    READ(10,*) System%epsilom
    READ(10,*) System%sigma
    READ(10,*) System%L
    READ(10,*) System%rc

    ! Inputs of TEMP parameters
    READ(10,*) temperature

    ! Inputs of MCeq parameters
    READ(10,*) MCeq%nsteps
    READ(10,*) MCeq%stepsize
    READ(10,*) MCeq%outfreq

    ! Inputs of MC parameters
    READ(10,*) MC%nsteps
    READ(10,*) MC%stepsize
    READ(10,*) MC%outfreq

    ! Inputs of RAFEP and NS parameters
    READ(10,*) NS%rafep_cutoff
    READ(10,*) NS%nsamples
    READ(10,*) NS%nsteps
    READ(10,*) NS%stepsize
    READ(10,*) NS%fractiom
    READ(10,*) NS%root_energy

  CLOSE(10)

  ! Seed the random number
  CALL RANDOM_SEED()

  ! Check input parameters
  IF (System%rc > 0.5d0*System%L) THEN
    WRITE(*,*) "VDW cutoff (rc) should be smaller than half of the boxsize. Check your inputs."
    STOP
  END IF

  ! Unit transform (to reduce units)
  System%L = System%L/System%sigma
  System%rc = System%rc/System%sigma 
  temperature = kB*temperature/System%epsilom
  MCeq%stepsize = MCeq%stepsize/System%sigma
  MC%stepsize = MC%stepsize/System%sigma
  NS%stepsize = NS%stepsize/System%sigma
  NS%root_energy = NS%root_energy/System%epsilom

  ! Initialize parameters
  System%rc2 = System%rc**2
  System%Ec = 4*(1.d0/System%rc**12 - 1.d0/System%rc**6)
  beta = 1.d0/temperature
  kBT = temperature
  MCeq%outdim = MCeq%nsteps/MCeq%outfreq
  MC%outdim = MC%nsteps/MC%outfreq

  ! Build LJ system and calculate the system energy
  ALLOCATE(System%XYZ(3,System%natoms))
  CALL System%genXYZ()
  CALL System%calcenergy()

  ! Pre-equilibration
  CALL MCSampling(System, MCeq)

  ! Monte Carlo Sampling
  CALL MCSampling(System, MC)
  
  ! Output MC Trajectory
  OPEN(UNIT=10,FILE="Traj.dat",STATUS="UNKNOWN")
  DO i=1,MC%outdim
    WRITE(10,*) "Time = ", i*MC%outfreq
    DO j=1,System%natoms
      WRITE(10,*) (MC%Traj(k,j,i)*System%sigma, k=1,3)
    END DO
  END DO
  CLOSE(10)

  ! Ouput MC Energy
  OPEN(UNIT=10,FILE="Energy.dat",STATUS="UNKNOWN")
  DO i=1,MC%outdim
    WRITE(10,*) MC%Energy(i)*System%epsilom
  END DO
  CLOSE(10)
  
  ! Output MC Statistics
  PRINT *, "Acceptance rate =", MC%Accept
  PRINT *, "Average energy (kJ/mol) =", SUM(MC%Energy)*(System%epsilom*cal2joule)/MC%outdim

  ! For 2 atoms debug
  Z = 0.d0
  expectation = 0.d0
  dr = System%rc / 10000
  DO i=1,10000
    r = (i-1)*dr
    ! do not calculate the sum if the distance is too small
    IF (r < 0.1d0) THEN
      CONTINUE
    ELSE
      pot = System%calcpairpot(r*r)
      prob = EXP(-beta*pot)
      Z = Z + prob*(4*3.1415926536*r*r)*dr
      expectation = expectation + prob*pot*(4*3.1415926536*r*r)*dr
    END IF
  END DO
  restV = System%L**3 - 4*3.1415926536/3.d0 * System%rc**3
  Z = Z + restV
  PRINT *, "Energy expectation for 2 atoms (kJ/mol) =", expectation*(System%epsilom*cal2joule)/Z



  STOP
  ! Remove the sampling outlier
  NS%rafep_cutoff = NS%rafep_cutoff*kBT + MINVAL(MC%Energy)
  CALL Truncate(MC,System%natoms,NS%rafep_cutoff)
  
  ! RAFEP
  Z = Partition_RAFEP(System,NS,MC%NewEnergy)
  print *, "RAFEP Zest:", Z
  
STOP
END PROGRAM Main
