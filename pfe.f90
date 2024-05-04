! Module for Partition Function Estimator (PFE)
MODULE MODPFE
    USE MODCONST
    USE MODLJ
    USE MODMC 
    USE MODUTIL
  IMPLICIT NONE

  TYPE PFE
    INTEGER :: nsamples, nsteps, nlevel
    REAL*8  :: stepsize, fractiom
    REAL*8  :: Emax, Emin, Eroot
    REAL*8  :: lnZ, volume
    REAL*8, ALLOCATABLE  :: levels(:), volumes(:)
    CONTAINS
      PROCEDURE :: rdinp => pfe_rdinp
      PROCEDURE :: PartFunc
      PROCEDURE :: NSVolume 
      PROCEDURE :: NSPartinit
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
    READ(fd,*) self%fractiom
    READ(fd,*) self%Eroot
  END SUBROUTINE pfe_rdinp

  ! Calculate the partition function via Partition Function Estimator under PBC
  SUBROUTINE PartFunc(self,System,Energy,beta)
    CLASS(PFE) :: self
    TYPE(LJ) :: System
    INTEGER :: outdim
    REAL*8 :: beta 
    REAL*8 :: Energy(:)
    REAL*8, ALLOCATABLE  :: Work(:)

    ! initialization
    outdim = SIZE(Energy)
    self%Emax = MAXVAL(Energy)
    self%Emin = MINVAL(Energy)
    print *, "Emin =", self%Emin
    print *, "Emax =", self%Emax

    ! calculate the relative volume by the nested sampling 
    CALL self%NSVolume(System)
  
    ! calculate the partition function
    ALLOCATE(Work(outdim))
    Work(:) = EXP(beta*(Energy(:) - self%Emax))
    self%lnZ = (LOG(self%volume) + 3*System%natoms * LOG(System%L)) - (LOG(SUM(Work)/outdim) + beta*self%Emax)
  
  END SUBROUTINE PartFunc
  
  
  ! Calculate the volume ratio via the nested sampling
  SUBROUTINE NSVolume(self,System)
    CLASS(PFE) :: self
    TYPE(LJ) :: System
    TYPE(LJ), ALLOCATABLE :: Samples(:)
    INTEGER :: i
    INTEGER :: nlevel, mycount
    REAL*8 :: level, volume, ratio
    
    ! determine nlevel (number of energy levels required until Emax is reached)
    nlevel = 0
    level = self%Eroot
    DO WHILE (level > self%Emax)
      nlevel = nlevel + 1
      level = (level - self%Emin) * self%fractiom + self%Emin
    END DO
    self%nlevel = nlevel
    
    ! allocate arrays
    ALLOCATE(self%levels(nlevel),self%volumes(nlevel))
    ALLOCATE(Samples(self%nsamples))
  
    ! initialization
    DO i = 1, self%nsamples
      Samples(i) = System
      CALL Samples(i)%genXYZ()
    END DO
  
    ! to start, make sure all points are within the root energy
    DO i = 1, self%nsamples
      CALL Samples(i)%calcenergy()
      IF (Samples(i)%energy > self%Eroot) THEN
        CALL relax(Samples(i), self%Eroot, self%stepsize)
        CALL propagate(Samples(i), self%Eroot, self%nsteps, self%stepsize)
      END IF
    END DO
  
    ! calculate the relative volume iteratively
    nlevel = 0
    level = self%Eroot
    volume = 1.d0
    DO WHILE (level > self%Emax)
      ! introduce the Emin to ensure the threshold is correctly scaled
      level = (level - self%Emin) * self%fractiom + self%Emin
  
      ! push down (push the samples outside to inside the threshold, while counting how many samples are already inside)
      mycount = 0
      DO i = 1, self%nsamples
        CALL Samples(i)%calcenergy()
        IF (Samples(i)%energy <= level) THEN
           mycount = mycount +1
        ELSE
           CALL relax(Samples(i), level, self%stepsize)
           CALL propagate(Samples(i), level, self%nsteps, self%stepsize)
        END IF
      END DO
      
      ! calculate the relative volume under the current energy level
      ratio = mycount*1.d0/self%nsamples
      volume = volume * ratio
      nlevel = nlevel + 1
      self%levels(nlevel) = level
      self%volumes(nlevel) = volume
  
    END DO
   

    IF (self%Emax == 0.d0) THEN
      ! Emax = 0.d0 is a discontinuity in accessible volume space. It cannot be treated using linear interpolation.
      self%volume = self%volumes(self%nlevel-1)
    ELSE
      ! interpolate to find target relative volume at Emax
      self%volume = INTERPOLATE(self%levels,self%volumes,self%nlevel,self%Emax)
    END IF
    PRINT *, "volume = ", self%volume

    ! Ouput levels (unit: kj/mol)
    OPEN(UNIT=20,FILE="Levels.dat",STATUS="UNKNOWN")
    DO i=1,self%nlevel
      WRITE(20,*) self%levels(i)*cal2joule, self%volumes(i)
    END DO
    CLOSE(20)

  END SUBROUTINE NSVolume
  
  
  ! Initialize parameters for the NSPartition (Emax and Emin differ from PFE method)
  SUBROUTINE NSPartinit(self,System,MCrun,Parfu,kBT)
    CLASS(PFE) :: self, Parfu
    TYPE(LJ) :: System
    TYPE(MC) :: MCrun
    REAL*8 :: kBT

    self%nsamples = Parfu%nsamples 
    self%nsteps   = Parfu%nsteps
    self%stepsize = Parfu%stepsize
    self%fractiom = Parfu%fractiom
    self%Eroot    = Parfu%Eroot

    CALL MCrun%Minimize(System) 
    CALL System%calcenergy()
    self%Emin = System%Energy
    self%Emax = self%Emin + kBT*0.01

  END SUBROUTINE        


  ! Calculate the partition function via the nested sampling (need NSVolume)
  SUBROUTINE NSPartition(self,System,beta)
    CLASS(PFE) :: self
    TYPE(LJ) :: System
    INTEGER :: i
    REAL*8 :: beta
    REAL*8 :: energy, Eshift
    REAL*8 :: summation
 
    ! integration
    summation = 0.d0
    Eshift = MINVAL(self%levels)
    self%levels = self%levels - Eshift
    DO i = 1, self%nlevel-1
      energy = 0.5d0*(self%levels(i)+self%levels(i+1))
      summation = summation + EXP(-beta*energy)*(self%volumes(i)-self%volumes(i+1))
    END DO  

    ! lnZ = 3N*ln(L) - beta * Emin + ln(sum(EXP(-beta*(E-Emin))*dV)) 
    self%lnZ = 3*System%natoms*LOG(System%L) - beta*Eshift + LOG(summation)

    ! output
    OPEN(UNIT=20,FILE="Levels2.dat",STATUS="UNKNOWN")
    DO i=1,self%nlevel
      WRITE(20,*) self%levels(i)*cal2joule, self%volumes(i)
    END DO
    CLOSE(20)
    
  END SUBROUTINE NSPartition
  
  
  ! Relax the system to energy lower than threshold
  SUBROUTINE relax(System, threshold, stepsize)
    TYPE(LJ) System
    INTEGER :: aid
    REAL*8 :: threshold, stepsize, E1, E2
    REAL*8 :: coords(3)
  
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
