! Module for Partition Function Estimator (PFE)
MODULE MODPFE
    USE MODCONST
    USE MODLJ
    USE MODMC
    USE MODUTIL
  IMPLICIT NONE

  TYPE PFE
    INTEGER :: nsamples, nsteps, nlevel, nextrasteps
    REAL*8  :: stepsize, fract, percentage
    REAL*8  :: Edagg, Estar, Err2, Avg, Avg2
    REAL*8  :: Eroot, VErr2
    REAL*8  :: lnQ, logVolume, logVstar, logVdagg
    REAL*8, ALLOCATABLE  :: levels(:), logVolumes(:)
    CHARACTER(LEN=10) :: Method
    CONTAINS
      PROCEDURE :: rdinp => pfe_rdinp
      PROCEDURE :: PartFunc
      PROCEDURE :: PartFunc2
      PROCEDURE :: NSVolume
      PROCEDURE :: NSPartition
  END TYPE PFE

  CONTAINS

  ! read input parameters for PFE
  SUBROUTINE pfe_rdinp(self,fd)
    CLASS(PFE) :: self
    INTEGER :: fd

    READ(fd,*) self%nsamples
    READ(fd,*) self%nextrasteps
    READ(fd,*) self%nsteps
    READ(fd,*) self%stepsize
    READ(fd,*) self%fract
    READ(fd,*) self%Eroot
  END SUBROUTINE pfe_rdinp

  ! Calculate the partition function via Partition Function Estimator under PBC
  SUBROUTINE PartFunc(self,System,Energy,beta)
    CLASS(PFE) :: self
    TYPE(LJ),INTENT(IN) :: System
    REAL*8,INTENT(IN) :: Energy(:)
    REAL*8,INTENT(IN) :: beta
    INTEGER :: i, edim
    REAL*8 :: Emax, Emin, Estar_prev, Estar, difference
    REAL*8 :: Avg, Avg2, Err2, LogOmega
    REAL*8, ALLOCATABLE  :: Work(:), Heaviside(:), Func(:), Func2(:)

    ! Initialization
    edim = SIZE(Energy)
    ALLOCATE(Work(edim),Heaviside(edim),Func(edim),Func2(edim))
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
      DO i = 1, edim
        IF (Energy(i) > Estar) THEN
          Heaviside(i) = 0.d0
        END IF
      END DO

      ! Calculate Avg and Avg2 (Energy shifted)
      Avg = SUM(Func*Heaviside)/edim
      Avg2= SUM(Func2*Heaviside)/edim

      ! Update E*, Err2, difference
      Estar = (LOG(2.d0) + LOG(Avg2) - LOG(Avg))/beta + Emax
      Err2 = (Avg2/Avg**2-1)/edim
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
    DO i = 1, edim
      IF (Energy(i) > Estar) THEN
        Heaviside(i) = 0.d0
        self%percentage = self%percentage + 1
      END IF
    END DO
    Avg = SUM(Func*Heaviside)/edim
    Avg2= SUM(Func2*Heaviside)/edim
    Err2= (Avg2/Avg**2 -1)/edim
    self%percentage = self%percentage/edim
    self%Err2 = Err2

    ! calculate the volume by the nested sampling
    CALL self%NSVolume(System,Emin,.FALSE.)

    ! calculate the partition function (Omega = L**3N * volume)
    LogOmega = 3*System%natoms*LOG(System%L) + self%logVolume
    self%lnQ = LogOmega - LOG(Avg) - beta*Emax

  END SUBROUTINE PartFunc


  ! Calculate the partition function via Partition Function Estimator (both side cutoff)
  SUBROUTINE PartFunc2(self,System,Energy,beta,nbin)
    CLASS(PFE) :: self
    TYPE(LJ),INTENT(IN) :: System
    REAL*8,INTENT(IN) :: Energy(:)
    REAL*8,INTENT(IN) :: beta
    INTEGER,INTENT(IN) :: nbin
    INTEGER :: i, j, k, edim
    REAL*8 :: Emax, Emin, dE
    REAL*8 :: Estar, Edagg
    REAL*8 :: Avg, Avg2, Err2, Err2min, LogOmega
    REAL*8, ALLOCATABLE  :: Work(:), Heaviside(:), Func(:), Func2(:), Egrid(:)

    ! Initialization
    edim = SIZE(Energy)
    ALLOCATE(Work(edim),Heaviside(edim),Func(edim),Func2(edim),Egrid(nbin))
    Emax = MAXVAL(Energy)
    Emin = MINVAL(Energy)
    PRINT *, 'DEBUG PartFunc2: Emax  = ', Emax, 'Emin  = ', Emin
    Func(:) = EXP(beta*(Energy-Emax))
    Func2(:) = EXP(2*beta*(Energy-Emax))

    ! Calculate the energy grid for the distribution histogram
    dE = (Emax-Emin)/(nbin-1)
    DO i = 1, nbin
       Egrid(i) = Emin + (i-1)*dE
    END DO
    Egrid(1) = Emin      ! for numerical stability
    Egrid(nbin) = Emax   ! for numerical stability

    ! Calculate the initial guess for Err2min, use the absolute value for numerical stability
    Avg = SUM(Func)/edim
    Avg2= SUM(Func2)/edim
    Err2min = ABS((Avg2/Avg**2-1)/edim)

    ! Find Estar and Edagg based on Egrid
    OPEN(UNIT=15, FILE='errmap.dat', STATUS='unknown', ACTION='write')
    DO i = 1, nbin

      Estar = Egrid(nbin-i+1)

      DO j = 1, nbin

        Edagg = Egrid(j)
        IF (Edagg >= Estar) THEN
          WRITE(15,*) Edagg, Estar, 'NaN'
          CYCLE
        END IF

        ! Calculate the averages and error
        Heaviside(:) = 1.d0
        DO k = 1, edim
          IF ((Energy(k) < Edagg) .OR. (Energy(k) > Estar)) THEN
            Heaviside(k) = 0.d0
          END IF
        END DO

        Avg = SUM(Func*Heaviside)/edim
        Avg2= SUM(Func2*Heaviside)/edim
        Err2 = ABS((Avg2/Avg**2-1)/edim)

        ! Update
        IF (Err2 < Err2min) THEN
          self%Estar = Estar
          self%Edagg = Edagg
          self%Err2 = Err2
          self%Avg  = Avg
          self%Avg2 = Avg2
          Err2min = Err2
        END IF
        WRITE(15,*) Edagg, Estar, SQRT(Err2)

      END DO
      WRITE(15,*)
    END DO
    CLOSE(UNIT=15)

    ! Reassign Estar and Edagg
    Estar = self%Estar
    Edagg = self%Edagg
    PRINT *, 'DEBUG PartFunc2: Estar = ', Estar, 'Edagg = ', Edagg, ' Err2 = ', Err2min

    ! Calculate the cutoff percentage for statistics
    Heaviside(:) = 1.d0
    self%percentage = 0
    DO i = 1, edim
      IF ((Energy(i) < Edagg) .OR. (Energy(i) > Estar)) THEN
        self%percentage = self%percentage + 1
      END IF
    END DO
    self%percentage = self%percentage/edim

    ! calculate the volume by the nested sampling
    IF (Edagg == Emin) THEN
      CALL self%NSVolume(System,Emin,.FALSE.)
    ELSE
      CALL self%NSVolume(System,Emin,.TRUE.)
    END IF

    ! calculate the partition function (Omega = L**3N * volume)
    LogOmega = 3*System%natoms*LOG(System%L) + self%logVolume
    self%lnQ = LogOmega - LOG(self%Avg) - beta*Emax

  END SUBROUTINE PartFunc2


  ! Calculate the volume ratio via the nested sampling
  SUBROUTINE NSVolume(self,System,Emin,flag)
    CLASS(PFE) :: self
    TYPE(LJ),INTENT(IN) :: System
    REAL*8,INTENT(IN) :: Emin
    LOGICAL,INTENT(IN) :: flag ! .True. for both end cutting
    TYPE(LJ), ALLOCATABLE :: Samples(:)
    INTEGER :: i
    INTEGER :: ndagg, nstar, niter, nsamples, nsteps, nextrasteps
    INTEGER :: iterid, nrelaxsteps, tnrelaxsteps, noutliers, ninliers
    INTEGER :: naccept, tnaccept_relax, tnaccept_prop
    REAL*8 :: Edagg, Estar, Eroot, fract, Elevel, Elevel_further
    REAL*8 :: stepsize, logVolume, VErr2, Einitmin

    ! initialization
    Edagg = self%Edagg
    Estar = self%Estar
    Eroot = self%Eroot
    fract = self%fract
    nsamples = self%nsamples
    nextrasteps = self%nextrasteps
    nsteps = self%nsteps
    stepsize = self%stepsize
    nrelaxsteps = 0
    VErr2 = 0

    ! determine nstar (number of iterations required to reach Estar)
    nstar = 0
    Elevel = Eroot
    DO WHILE (Elevel > Estar)
      ! introduce Emin to ensure the threshold is correctly scaled
      Elevel = (Elevel - Emin) * fract + Emin
      IF (Elevel < Estar)  Elevel = Estar
      nstar = nstar + 1
    END DO

    ! determine ndagg (number of extra iterations required to reach Edagg)
    ndagg = 0
    IF (flag) THEN
      DO WHILE (Elevel > Edagg)
        ! introduce Emin to ensure the threshold is correctly scaled
        Elevel = (Elevel - Emin) * fract + Emin
        IF (Elevel < Edagg)  Elevel = Edagg
        ndagg = ndagg + 1
      END DO
    END IF

    ! assign the number of iteration to self%nlevel
    niter = nstar + ndagg
    self%nlevel = niter

    ! allocate arrays
    ALLOCATE(self%levels(niter),self%logVolumes(niter))
    ALLOCATE(Samples(nsamples))

    ! to start, generate samples all within the root energy,
    ! via rejection sampling
    Einitmin = Eroot
    i = 1
    DO WHILE (i <= nsamples)
      Samples(i) = System
      CALL Samples(i)%genXYZ()
      IF (Samples(i)%getenergy() <= Eroot) THEN
        Einitmin = Samples(i)%getenergy()
        i = i + 1 
      END IF
    END DO

    ! ouput for NS performance statistics
    OPEN(UNIT=10,FILE="Statistics.dat",STATUS="UNKNOWN")

    Elevel = Eroot
    logVolume = 0.d0

    ! calculate the log of the relative volume iteratively
    DO iterid = 1, niter

      Elevel = (Elevel - Emin) * fract + Emin
      Elevel_further = (Elevel - Emin) * fract + Emin
      ! fix the level to Edagg or Estar to avoid interpolation
      IF (iterid == niter)  Elevel = Edagg
      IF (iterid == nstar)  Elevel = Estar

      ! data for performance statistics
      noutliers = 0
      ninliers = 0
      tnrelaxsteps = 0
      tnaccept_relax = 0
      tnaccept_prop = 0

      ! push the samples outside to the area under Elevel, count inliers
      DO i = 1, nsamples
        IF (Samples(i)%getenergy() <= Elevel) THEN
           ninliers = ninliers + 1
        ELSE
           CALL relax(Samples(i), Elevel, stepsize, nrelaxsteps, nextrasteps, naccept)
           !CALL relax(Samples(i), Elevel_further, stepsize, nrelaxsteps, nextrasteps, naccept)
           ! data for performance statistics
           tnrelaxsteps = tnrelaxsteps + nrelaxsteps
           tnaccept_relax = tnaccept_relax + naccept

           CALL propagate(Samples(i), Elevel, nsteps, stepsize, naccept)
           ! data for performance statistics
           tnaccept_prop = tnaccept_prop + naccept
           noutliers = noutliers + 1
        END IF
      END DO

      ! data for performance statistics
      WRITE(10,'(i8,i8,2i12,2f8.3)') &
        iterid, noutliers, tnrelaxsteps, nsteps*noutliers, &
        1.d0*tnaccept_relax/tnrelaxsteps, 1.d0*tnaccept_prop/(nsteps*noutliers)
      FLUSH(10)

      PRINT *, 'DEBUG: step ', iterid, '/', niter, 'inliers', ninliers, '/', nsamples, 'Elevel / Einitmin = ', Elevel/Einitmin

      ! abort if no inliers (volume becomes zero)
      IF (ninliers == 0) THEN
        PRINT *, 'NSVolume: volume=0 after ',iterid,' iterations -- increase frac?'
        STOP 1
      END IF

      ! calculate the current relative volume and save it
      logVolume = logVolume + LOG(1.d0*ninliers/nsamples)
      self%levels(iterid) = Elevel
      self%logVolumes(iterid) = logVolume

      ! calculate the cumulated error VErr2
      VErr2 = VErr2 + (1.d0/ninliers - 1.d0/nsamples)

      ! update when iterid == nstar
      IF (iterid == nstar) THEN
        self%logVstar = logVolume
        self%VErr2 = VErr2
      END IF

    END DO

    ! update when iterid == niter
    IF (flag) THEN
      self%logVdagg = logVolume
      !self%VErr2 = self%VErr2 + VErr2 ! not sure about this...
    END IF

    CLOSE(10)


    IF (flag) THEN
      ! Volume = Vstar - Vdagg
      ! but Vdagg <<< Vstar, so the correction hardly matters...
      self%logVolume = LOG(EXP(self%logVstar) - EXP(self%logVdagg))
    ELSE
      self%logVolume = self%logVstar
    END IF

    PRINT *, "logVolume = ", self%logVolume, "volume = ", EXP(self%logVolume), "err = ", SQRT(self%VErr2)

    ! Output levels (unit: kj/mol)
    OPEN(UNIT=20,FILE="Levels.dat",STATUS="UNKNOWN")
    DO i = 1, niter
      WRITE(20,*) self%levels(i)*cal2joule, self%logVolumes(i), EXP(self%logVolumes(i))
    END DO
    CLOSE(20)

  END SUBROUTINE NSVolume


  ! Calculate the partition function via nested sampling (use NSVolume routine)
  SUBROUTINE NSPartition(self,System,beta)
    CLASS(PFE) :: self
    TYPE(LJ),INTENT(IN) :: System
    REAL*8,INTENT(IN) :: beta
    INTEGER :: i
    REAL*8 :: energy, Emin
    REAL*8 :: summation

    ! integration
    summation = 0.d0
    Emin = MINVAL(self%levels)
    self%levels = self%levels - Emin
    DO i = 1, self%nlevel-1
      energy = 0.5d0*(self%levels(i)+self%levels(i+1))
      summation = summation + EXP(-beta*energy)*(EXP(self%logVolumes(i))-EXP(self%logVolumes(i+1)))
    END DO

    ! lnQ = 3N*ln(L) - beta * Emin + ln(sum(EXP(-beta*(E-Emin))*dV))
    self%lnQ = 3*System%natoms*LOG(System%L) - beta*Emin + LOG(summation)

  END SUBROUTINE NSPartition


  ! Relax the system to energy lower than threshold
  SUBROUTINE relax(System, threshold, stepsize, nrelaxsteps, nextrasteps, naccept)
    TYPE(LJ) System
    REAL*8,INTENT(IN) :: threshold, stepsize
    INTEGER,INTENT(IN) :: nextrasteps
    INTEGER,INTENT(OUT) :: nrelaxsteps, naccept
    INTEGER :: i, aid
    REAL*8 :: E1, E2
    REAL*8 :: coords(3)

    naccept = 0
    nrelaxsteps = 0
    coords = 0.d0
    E1 = System%getenergy()

    i = 0 ! count extra steps
    DO WHILE ((E1 > threshold) .OR. (i < nextrasteps))
      IF (E1 <= threshold)  i = i + 1
      aid = RANDOM_INTEGER(System%natoms)
      coords(:) = System%XYZ(:,aid)
      CALL System%move(aid,stepsize)
      E2 = System%getenergy()
      IF (E2 > E1) THEN
        ! reject the move if energy gets higher
        System%XYZ(:,aid) = coords(:)
        System%energy = E1
      ELSE
        ! accept the move if energy gets lower
        E1 = E2
        naccept = naccept + 1
      END IF
      nrelaxsteps = nrelaxsteps + 1
    END DO

  END SUBROUTINE relax


  ! Propagate the system while maintaining the energy under the threshold
  SUBROUTINE propagate(System, threshold, nsteps, stepsize, naccept)
    TYPE(LJ) System
    REAL*8,INTENT(IN) :: threshold, stepsize
    INTEGER,INTENT(IN) :: nsteps
    INTEGER,INTENT(OUT) :: naccept
    INTEGER :: i, aid
    REAL*8 :: E1, E2
    REAL*8 :: coords(3)

    naccept = 0
    coords = 0.d0
    E1 = System%getenergy()
    DO i = 1, nsteps
      aid = RANDOM_INTEGER(System%natoms)
      coords(:) = System%XYZ(:,aid)
      CALL System%move(aid,stepsize)
      E2 = System%getenergy()
      IF (E2 > threshold) THEN
        ! reject the move, restore the coordinates
        System%XYZ(:,aid) = coords(:)
        System%energy = E1
      ELSE
        ! accept the move, update the energy
        E1 = E2
        naccept = naccept + 1
      END IF
    END DO

  END SUBROUTINE propagate

END MODULE MODPFE
