! Module for Partition Function Estimator (PFE)
MODULE MODPFE
    USE MODCONST
    USE MODLJ
    USE MODMC 
    USE MODUTIL
  IMPLICIT NONE

  TYPE PFE
    INTEGER :: nsamples, nsteps, nlevel
    REAL*8  :: stepsize, fract, percentage
    REAL*8  :: Estar, Eroot, Err2, VErr2
    REAL*8  :: lnZ, volume
    REAL*8, ALLOCATABLE  :: levels(:), volumes(:)
    CHARACTER(LEN=10) :: Method
    CONTAINS
      PROCEDURE :: rdinp => pfe_rdinp
      PROCEDURE :: PartFunc
      PROCEDURE :: NSVolume 
      PROCEDURE :: NSPartition 
  END TYPE PFE

  CONTAINS

  ! read input parameters for PFE
  SUBROUTINE pfe_rdinp(self,fd)
    CLASS(PFE) :: self
    INTEGER :: fd

    READ(fd,*) self%nsamples
    READ(fd,*) self%nsteps
    READ(fd,*) self%stepsize
    READ(fd,*) self%fract
    READ(fd,*) self%Eroot
    READ(fd,*) self%Method
  END SUBROUTINE pfe_rdinp

  ! Calculate the partition function via Partition Function Estimator under PBC
  SUBROUTINE PartFunc(self,System,Energy,beta)
    CLASS(PFE) :: self
    TYPE(LJ) :: System
    INTEGER :: i
    INTEGER :: ndim
    REAL*8 :: beta, Emax, Emin, Estar_prev, Estar, difference 
    REAL*8 :: Avg, Avg2, Err2, LogOmega
    REAL*8 :: Energy(:)
    REAL*8, ALLOCATABLE  :: Work(:), Heaviside(:), Func(:), Func2(:)

    ! Initialization
    ndim = SIZE(Energy)
    ALLOCATE(Work(ndim),Heaviside(ndim),Func(ndim),Func2(ndim))
    Emax = MAXVAL(Energy)
    Emin = MINVAL(Energy)
    Func(:) = EXP(beta*(Energy-Emax))
    Func2(:) = EXP(2*beta*(Energy-Emax))

    Estar = Emax
    Estar_prev = Estar
    difference = 1

    ! Solve Estar iteratively
    DO WHILE (difference > 1E-6)
      ! Generate the Heaviside function
      Heaviside(:) = 1.d0
      DO i = 1, ndim
        IF (Energy(i) > Estar) THEN 
          Heaviside(i) = 0.d0
        END IF
      END DO

      ! Calculate Avg and Avg2 (Energy shifted)
      Avg = SUM(Func*Heaviside)/ndim
      Avg2= SUM(Func2*Heaviside)/ndim

      ! Update E*, Err2, difference
      Estar = (LOG(2.d0) + LOG(Avg2) - LOG(Avg))/beta + Emax
      Err2 = (Avg2/Avg**2-1)/ndim
      difference = ABS(Estar-Estar_prev)
      

      ! For next iteration
      Estar_prev = Estar
    END DO

    ! For Estar = 0, choose Estar=Emin+4*kBT
    IF (Estar == 0.d0) Estar = MINVAL(Energy) + 4.d0/beta

    ! Assign Estar
    self%Estar = Estar

    ! Recalculate Avg, Avg2, Err based on the final Estar, also calculate the cutoff percentage for statistics
    Heaviside(:) = 1.d0
    self%percentage = 0
    DO i = 1, ndim
      IF (Energy(i) > Estar) THEN
        Heaviside(i) = 0.d0
        self%percentage = self%percentage + 1
      END IF
    END DO 
    Avg = SUM(Func*Heaviside)/ndim
    Avg2= SUM(Func2*Heaviside)/ndim
    Err2= (Avg2/Avg**2 -1)/ndim
    self%percentage = self%percentage/ndim
    self%Err2 = Err2

    ! calculate the volume by the nested sampling
    CALL self%NSVolume(System,Emin)
  
    ! calculate the partition function (Omega = L**3N * volume)
    LogOmega = 3*System%natoms*LOG(System%L) + LOG(self%volume) 
    self%lnZ = LogOmega - LOG(Avg) - beta*Emax

    !! DEALLOCATE 
    !DEALLOCATE(Work(ndim),Heaviside(ndim),Func(ndim),Func2(ndim))
  END SUBROUTINE PartFunc
  
  
  ! Calculate the volume ratio via the nested sampling
  SUBROUTINE NSVolume(self,System,Emin)
    CLASS(PFE) :: self
    TYPE(LJ) :: System
    TYPE(LJ), ALLOCATABLE :: Samples(:)
    INTEGER :: i
    INTEGER :: niter, nsamples, nsteps 
    INTEGER :: iterid, nrelaxsteps, tnrelaxsteps, noutliers, ninliers
    REAL*8 :: Emin, Estar, Eroot, fract, Elevel
    REAL*8 :: stepsize, volume
    
    ! initialization
    Estar = self%Estar
    Eroot = self%Eroot
    fract = self%fract
    nsamples = self%nsamples
    nsteps = self%nsteps
    stepsize = self%stepsize
    nrelaxsteps = 0
    self%VErr2 = 0

    ! determine niter (number of iterations required to reach Estar = self%nlevel)
    niter = 0
    Elevel = Eroot
    DO WHILE (Elevel > Estar)
      ! introduce Emin to ensure the threshold is correctly scaled
      Elevel = (Elevel - Emin) * fract + Emin
      niter = niter + 1
    END DO
    self%nlevel = niter
    
    ! allocate arrays
    ALLOCATE(self%levels(niter),self%volumes(niter))
    ALLOCATE(Samples(nsamples))
  
    ! initialization
    DO i = 1, nsamples
      Samples(i) = System
      CALL Samples(i)%genXYZ()
    END DO
  
    ! to start, make sure all points are within the root energy
    DO i = 1, nsamples
      CALL Samples(i)%calcenergy()
      IF (Samples(i)%energy > Eroot) THEN
        CALL relax(Samples(i), Eroot, stepsize, nrelaxsteps)
        CALL propagate(Samples(i), Eroot, nsteps, stepsize)
      END IF
    END DO
  
    ! ouput for NS performance statistics
    OPEN(UNIT=10,FILE="Statistics.dat",STATUS="UNKNOWN") 

    iterid = 0
    Elevel = Eroot
    volume = 1.d0

    ! calculate the relative volume iteratively
    DO WHILE (Elevel > Estar)
      Elevel = (Elevel - Emin) * fract + Emin
  
      ! data for performance statistics
      iterid = iterid + 1
      noutliers = 0
      ninliers = 0
      tnrelaxsteps = 0

      ! push the samples outside to the area under Elevel, count ninsider
      DO i = 1, nsamples
        CALL Samples(i)%calcenergy()
        IF (Samples(i)%energy <= Elevel) THEN
           ninliers = ninliers + 1
        ELSE
           CALL relax(Samples(i), Elevel, stepsize, nrelaxsteps)
           CALL propagate(Samples(i), Elevel, nsteps, stepsize)
           ! data for performance statistics
           noutliers = noutliers + 1
           tnrelaxsteps = tnrelaxsteps + nrelaxsteps
        END IF
      END DO
      

      ! data for performance statistics
      WRITE(10,*) iterid, noutliers, tnrelaxsteps, nsteps*noutliers

      ! calculate the current relative volume and save it
      volume = volume * (1.d0*ninliers/nsamples)
      self%levels(iterid) = Elevel
      self%volumes(iterid) = volume

      ! calculate the cummulated error VErr2
      IF (ninliers /= 0) THEN
        self%VErr2 = self%VErr2 + (1.d0/ninliers - 1.d0/nsamples)
      END IF

    END DO

    CLOSE(10)

    IF (Estar == 0.d0) THEN
      ! Estar = 0.d0 is a discontinuity in accessible volume space. It cannot be treated using linear interpolation.
      self%volume = self%volumes(niter-1)
    ELSE
      ! interpolate to find the relative volume under Estar
      self%volume = INTERPOLATE(self%levels,self%volumes,niter,Estar)
    END IF
    PRINT *, "volume = ", self%volume, "err = ", SQRT(self%VErr2)

    ! Ouput levels (unit: kj/mol)
    OPEN(UNIT=20,FILE="Levels.dat",STATUS="UNKNOWN")
    DO i = 1, niter
      WRITE(20,*) self%levels(i)*cal2joule, self%volumes(i)
    END DO
    CLOSE(20)

  END SUBROUTINE NSVolume
  
  
  ! Calculate the partition function via nested sampling (use NSVolume routine)
  SUBROUTINE NSPartition(self,System,beta)
    CLASS(PFE) :: self
    TYPE(LJ) :: System
    INTEGER :: i
    REAL*8 :: beta
    REAL*8 :: energy, Emin
    REAL*8 :: summation
 
    ! integration
    summation = 0.d0
    Emin = MINVAL(self%levels)
    self%levels = self%levels - Emin
    DO i = 1, self%nlevel-1
      energy = 0.5d0*(self%levels(i)+self%levels(i+1))
      summation = summation + EXP(-beta*energy)*(self%volumes(i)-self%volumes(i+1))
    END DO  

    ! lnZ = 3N*ln(L) - beta * Emin + ln(sum(EXP(-beta*(E-Emin))*dV)) 
    self%lnZ = 3*System%natoms*LOG(System%L) - beta*Emin + LOG(summation)

  END SUBROUTINE NSPartition
  
  
  ! Relax the system to energy lower than threshold
  SUBROUTINE relax(System, threshold, stepsize, nrelaxsteps)
    TYPE(LJ) System
    INTEGER :: aid, nrelaxsteps
    REAL*8 :: threshold, stepsize, E1, E2
    REAL*8 :: coords(3)
  
    nrelaxsteps = 0
    coords = 0.d0
    CALL System%calcenergy()
    E1 = System%energy
    DO WHILE (E1 > threshold)
      aid = RANDOM_INTEGER(System%natoms)
      coords(:) = System%XYZ(:,aid) 
      CALL System%move(System%XYZ(:,aid),stepsize)
      CALL System%calcenergy()
      E2 = System%energy
      IF (E2 > E1) THEN
        ! reject the move if energy gets higher
        System%XYZ(:,aid) = coords(:)
        System%energy = E1
      ELSE
        ! accept the move if energy gets lower
        E1 = E2
      END IF
      nrelaxsteps = nrelaxsteps + 1
    END DO
  
  END SUBROUTINE relax
  
  ! Propagate the system while maintaining the energy under the threshold
  SUBROUTINE propagate(System, threshold, nsteps, stepsize)
    TYPE(LJ) System
    INTEGER :: i
    INTEGER :: aid, nsteps
    REAL*8 :: threshold, stepsize, E1, E2
    REAL*8 :: coords(3)
  
    coords = 0.d0
    CALL System%calcenergy()
    E1 =System%energy
    DO i = 1, nsteps
      aid = RANDOM_INTEGER(System%natoms)
      coords(:) = System%XYZ(:,aid)
      CALL System%move(System%XYZ(:,aid),stepsize)
      CALL System%calcenergy()
      E2 = System%energy
      IF (E2 > threshold) THEN
        ! reject the move, restore the coordinates
        System%XYZ(:,aid) = coords(:)
        System%energy = E1
      ELSE
        ! accept the move, update the energy
        E1 = E2
      END IF
    END DO
  
  END SUBROUTINE propagate

END MODULE MODPFE
