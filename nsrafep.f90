! Module for NS/RAFEP methods
MODULE MODNSRAFEP
    USE MODTEMPAR
    USE MODPBCLJ
    USE MODMC 
    USE MODUTIL
  IMPLICIT NONE

  TYPE NSRAFEPPAR
    INTEGER :: nsamples, nsteps, nlevel
    REAL*8  :: stepsize, fractiom
    REAL*8  :: rafep_cutoff
    REAL*8  :: desire_energy, zero_energy, root_energy
    REAL*8, ALLOCATABLE  :: Lenergy(:), LV0(:)
  END TYPE NSRAFEPPAR

  CONTAINS

  !---------------------------------------------------------------------------------
  !
  ! Subroutines related to NS/RAFEP Algorithms
  !
  !---------------------------------------------------------------------------------
  
  SUBROUTINE Truncate(MC, natoms, cutoff)
    TYPE(MCPAR) :: MC
    INTEGER :: natoms
    INTEGER :: i, mycount, outdim, newoutdim
    REAL*8 :: cutoff
  
  
    ! reassign parameters
    outdim = MC%outdim
  
    ! calculate newoutdim (number of data with energy < threshold)
    newoutdim = 0
    DO i = 1, outdim
      IF (MC%Energy(i) < cutoff) THEN
        newoutdim = newoutdim + 1
      END IF
    END DO
    MC%newoutdim = newoutdim
  
    ! allocate the arrays
    ALLOCATE(MC%NewTraj(3,natoms,newoutdim),MC%NewEnergy(newoutdim))
  
    ! truncate the unwanted energy and associated trajectory
    mycount = 0
    DO i = 1, outdim
      IF (MC%Energy(i) < cutoff) THEN
        mycount = mycount + 1
        MC%NewEnergy(mycount) = MC%Energy(i)
        MC%NewTraj(:,:,mycount) = MC%Traj(:,:,i)
      END IF
    END DO
  
  END SUBROUTINE Truncate
  
  
  ! Calculate the partition function via RAFEP under PBC
  FUNCTION Partition_RAFEP(System,NS,Energy) RESULT(Z)
    TYPE(LJ) :: System
    TYPE(NSRAFEPPAR) :: NS
    INTEGER :: nsamples, nsteps, nlevel, outdim
    REAL*8  :: stepsize, fractiom
    REAL*8  :: V0, volume, Avg, Z
    REAL*8  :: Energy(:)
    REAL*8, ALLOCATABLE  :: Work(:)

    ! initialization
    outdim = SIZE(Energy)
    NS%desire_energy = MAXVAL(Energy)
    NS%zero_energy = MINVAL(Energy)
    volume = System%L ** (3*System%natoms)
    ALLOCATE(Work(outdim))

    ! calculate the volume by the nested sampling 
    CALL NSVolume(System, NS)
    
    ! interpolate to find desired V0 and turn it to volume
    V0 = INTERPOLATE(NS%Lenergy,NS%LV0,NS%nlevel,NS%desire_energy)
    volume = volume * V0
  
    ! calculate the partition function
    Work(:) = EXP(beta*(Energy(:) - NS%desire_energy))
    Avg=(SUM(Work)/outdim)*EXP(beta*NS%desire_energy)
    Z = volume / Avg
  
  END FUNCTION Partition_RAFEP
  
  
  ! Calculate the volume ratio via the nested sampling
  SUBROUTINE NSVolume(System, NS)
    TYPE(LJ) :: System
    TYPE(LJ), ALLOCATABLE :: Samples(:)
    TYPE(NSRAFEPPAR) :: NS
    INTEGER :: i
    INTEGER :: nsamples, nsteps, nlevel, mycount
    REAL*8 :: zero_energy, desire_energy, stepsize, fractiom, root_energy, V0
    REAL*8 :: threshold, energy, ratio
    
    ! assign parameters for easy typing
    nsamples = NS%nsamples
    nsteps = NS%nsteps
    stepsize = NS%stepsize
    fractiom = NS%fractiom
    zero_energy = NS%zero_energy
    desire_energy = NS%desire_energy
    root_energy = NS%root_energy

    ! determine nlevel
    nlevel = 0
    threshold = root_energy
    DO WHILE (threshold > desire_energy)
      nlevel = nlevel + 1
      threshold = (threshold - zero_energy) * fractiom + zero_energy
    END DO
    NS%nlevel = nlevel
    
    ! allocate arrayis
    ALLOCATE(NS%Lenergy(nlevel),NS%LV0(nlevel))
    ALLOCATE(Samples(nsamples))
  
    ! initialization
    DO i = 1, nsamples
      Samples(i) = System
      CALL Samples(i)%genXYZ()
    END DO
  
    ! to start, make sure all points are within the root energy
    DO i = 1, nsamples
      CALL Samples(i)%calcenergy()
      energy = Samples(i)%energy
      IF (energy > root_energy) THEN
        CALL relax(Samples(i), root_energy, stepsize)
      END IF
    END DO
  
    ! calculate the volume iteratively
    V0 = 1.d0
    nlevel = 0
    threshold = root_energy
    DO WHILE (threshold > desire_energy)
      ! introduce the zero_energy to ensure the threshold is correctly scaled
      threshold = (threshold - zero_energy) * fractiom + zero_energy
  
      ! push down (push the samples outside to inside the threshold, while counting how many samples are already inside)
      mycount = 0
      DO i = 1, nsamples
        CALL Samples(i)%calcenergy()
        energy = Samples(i)%energy
        IF (energy <= threshold) THEN
           mycount = mycount +1
        ELSE
           CALL relax(Samples(i), threshold, stepsize)
           CALL propagate(Samples(i), threshold, nsteps, stepsize)
        END IF
      END DO
      
      ! calculate the volume of the current energy level
      ratio = mycount*1.d0/nsamples
      V0 = V0 * ratio
      nlevel = nlevel + 1
      NS%Lenergy(nlevel) = threshold
      NS%LV0(nlevel) = V0
  
    END DO
  
  
  END SUBROUTINE NSVolume
  
  
  ! Relax the system to energy lower than threshold
  SUBROUTINE relax(System, threshold, stepsize)
    TYPE(LJ) System
    INTEGER :: aid
    REAL*8 :: threshold, stepsize, E1, E2
    REAL*8 :: coords(3)
  
    coords = 0.d0
    DO WHILE (E1 > threshold)
      CALL System%calcenergy()
      E1 = System%energy
      aid = RANDOM_INTEGER(System%natoms)
      coords(:) = System%XYZ(:,aid) 
      CALL System%move(System%XYZ(:,aid),stepsize)
      CALL System%calcenergy()
      E2 = System%energy
      IF (E2 > E1) THEN
        ! reject the move if energy gets higher
        System%XYZ(:,aid) = coords(:)
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
      ELSE
        ! accept the move, update the energy
        E1 = E2
      END IF
    END DO
  
  END SUBROUTINE propagate

END MODULE MODNSRAFEP
