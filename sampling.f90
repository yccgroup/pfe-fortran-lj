! Module for Monte Carlo Algorithm
MODULE MODMC
  USE MODCONST
  USE MODLJ
  USE MODUTIL
  IMPLICIT NONE

  TYPE MC
    INTEGER :: nsteps, outfreq, outdim
    REAL*8  :: stepsize
    REAL*8  :: Accept
    REAL*8, ALLOCATABLE :: Traj(:,:,:), Energy(:)
    CONTAINS
      PROCEDURE :: rdinp => mc_rdinp
      PROCEDURE :: Sampling
      PROCEDURE :: Truncate
      PROCEDURE :: Minimize
      PROCEDURE :: Write
      PROCEDURE :: Read
  END TYPE MC 

  CONTAINS
  !-----------------------------------------------------------------------------
  !
  ! Subroutines for Monte Carlo Sampling 
  !
  !-----------------------------------------------------------------------------
 
  ! read input parameters 
  SUBROUTINE mc_rdinp(self,fd)
    CLASS(MC) :: self
    INTEGER :: fd

    READ(fd,*) self%nsteps
    READ(fd,*) self%stepsize
    READ(fd,*) self%outfreq
    ! calculate derived parameters
    self%outdim = self%nsteps / self%outfreq

  END SUBROUTINE mc_rdinp

  ! Monte Carlo Sampling
  SUBROUTINE Sampling(self,System,beta)
    CLASS(MC) :: self
    TYPE(LJ) System
    INTEGER :: i
    INTEGER :: acount, aid, outindex
    REAL*8 :: beta
    REAL*8 :: Eold, Enew
    REAL*8 :: Xold(3), Xnew(3)
  
    ! initialize and allcoate arrays
    self%Accept = 0
    outindex = 0
    ALLOCATE(self%Traj(3,System%natoms,self%outdim),self%Energy(self%outdim))
    
    ! propagation
    DO i = 1, self%nsteps
  
      ! pick one atom to move
      aid = RANDOM_INTEGER(System%natoms)

      ! back up energy and coordinates
      Xold = System%XYZ(:,aid)
      Eold = System%getenergy()
    
      ! move the selected atom and calculate the new energy
      CALL System%move(aid,self%stepsize)
      Xnew = System%XYZ(:,aid)
      Enew = System%getenergy()
    
      ! Metropolis (Xnew, Enew updated)
      CALL Metropolis(Xold,Eold,Xnew,Enew,acount,beta)
      System%XYZ(:,aid) = Xnew
      System%energy = Enew
      self%Accept = self%Accept + acount
    
      ! output
      IF (MOD(i,self%outfreq) == 0) THEN
         outindex = outindex +1 
         self%Energy(outindex) = Enew
         self%Traj(:,:,outindex) = System%XYZ
      END IF
  
    END DO
  
    ! for statisitical output
    self%Accept = self%Accept*1.d0/self%nsteps
  
  END SUBROUTINE Sampling
  
  
  ! Determine Xnew and Enew by Metropolis algorithm
  SUBROUTINE Metropolis(Xold,Eold,Xnew,Enew,acount,beta)
    INTEGER :: acount
    REAL*8 :: Eold, Enew, acceptance, rand, beta
    REAL*8 :: Xold(3), Xnew(3)
  
    ! accept if Enew <= Eold
    IF (Enew <= Eold) THEN
      acount = 1
      RETURN
    END IF
  
    ! calculate the acceptance ratio
    acceptance = EXP(-beta*(Enew-Eold))
    CALL RANDOM_NUMBER(rand)
    IF (rand <= acceptance) THEN
      acount = 1
    ELSE
      acount = 0
      Enew = Eold
      Xnew = Xold
    END IF
    
  END SUBROUTINE Metropolis

  ! Truncate MC points with high energy
  SUBROUTINE Truncate(self, cutoff)
    CLASS(MC) :: self
    INTEGER :: i, mycount
    INTEGER :: natoms, newoutdim
    REAL*8 :: cutoff
    REAL*8, ALLOCATABLE :: Energy(:), Traj(:,:,:) 
  
    ! get value of natoms
    natoms = SIZE(self%Traj,2)

    ! calculate newoutdim (number of data with energy < threshold)
    newoutdim = 0
    DO i = 1, self%outdim
      IF (self%Energy(i) < cutoff) THEN
        newoutdim = newoutdim + 1
      END IF
    END DO
  
    ! allocate the arrays
    ALLOCATE(Energy(newoutdim),Traj(3,natoms,newoutdim))
  
    ! truncate the unwanted energy and associated trajectory
    mycount = 0
    DO i = 1, self%outdim
      IF (self%Energy(i) < cutoff) THEN
        mycount = mycount + 1
        Energy(mycount) = self%Energy(i)
        Traj(:,:,mycount) = self%Traj(:,:,i)
      END IF
    END DO

    ! assign the truncated data to self
    DEALLOCATE(self%Energy,self%Traj)
    CALL MOVE_ALLOC(Energy,self%Energy)
    CALL MOVE_ALLOC(Traj,self%Traj)
    self%outdim = newoutdim
  
  END SUBROUTINE Truncate


  ! Energy minimization for nsteps
  SUBROUTINE Minimize(self,System)
    CLASS(MC) :: self
    TYPE(LJ) System
    INTEGER :: i, aid
    REAL*8 :: E1, E2
    REAL*8 :: coords(3)
  
    coords = 0.d0
    E1 = System%getenergy()
    DO i = 1, self%nsteps
      aid = RANDOM_INTEGER(System%natoms)
      coords(:) = System%XYZ(:,aid) 
      CALL System%move(aid,self%stepsize)
      E2 = System%getenergy()
      IF (E2 > E1) THEN
        ! reject the move if energy gets higher
        System%XYZ(:,aid) = coords(:)
        System%energy = E1
      ELSE
        ! accept the move if energy gets lower
        E1 = E2
      END IF
    END DO
  
  END SUBROUTINE Minimize
  

  ! Write out the Sampling data
  SUBROUTINE Write(self,System,beta,filename)
    CLASS(MC) :: self
    TYPE(LJ),INTENT(IN) :: System
    REAL*8,INTENT(IN) :: beta
    CHARACTER(LEN=*),INTENT(IN) :: filename
    INTEGER :: fd

    OPEN(FILE=filename, NEWUNIT=fd, ACTION='write', FORM='unformatted', STATUS='unknown', ERR=500)
    WRITE(UNIT=fd, ERR=510) System%natoms, System%mass, System%epsilom, System%sigma, System%L, System%rc, beta
    WRITE(UNIT=fd, ERR=510) self%nsteps, self%outfreq, self%stepsize, self%Accept
    WRITE(UNIT=fd, ERR=510) self%Traj
    WRITE(UNIT=fd, ERR=510) self%Energy
    CLOSE(UNIT=fd)
    RETURN

    500 CONTINUE
    PRINT *, 'ERROR: could not open file ', filename
    STOP 1

    510 CONTINUE
    PRINT *, 'ERROR: could not write to file ', filename
    STOP 1

  END SUBROUTINE Write


  ! Read in the Sampling data
  SUBROUTINE Read(self,System,beta,filename)
    CLASS(MC) :: self
    TYPE(LJ),INTENT(IN) :: System
    REAL*8,INTENT(IN) :: beta
    CHARACTER(LEN=*),INTENT(IN) :: filename
    INTEGER :: fd
    INTEGER :: natoms, nsteps, outfreq, outdim
    REAL*8  :: mass, epsilom, sigma, L, rc, mybeta, stepsize

    OPEN(FILE=filename, NEWUNIT=fd, ACTION='read', FORM='unformatted', STATUS='old', ERR=600)
    ! read and check System parameters
    READ(UNIT=fd, ERR=610) natoms, mass, epsilom, sigma, L, rc, mybeta
    IF (natoms /= System%natoms) CALL COMPLAIN('System natoms')
    IF (mass /= System%mass) CALL COMPLAIN('System mass')
    IF (epsilom /= System%epsilom) CALL COMPLAIN('System epsilon')
    IF (sigma /= System%sigma) CALL COMPLAIN('System sigma')
    IF (L /= System%L) CALL COMPLAIN('System boxsize')
    IF (rc /= System%rc) CALL COMPLAIN('System LJ cutoff')
    IF (mybeta /= beta) CALL COMPLAIN('temperature')
    ! read and check MC parameters
    READ(UNIT=fd, ERR=610) nsteps, outfreq, stepsize, self%Accept
    IF (outfreq /= self%outfreq) CALL WARN('MC output frequency')
    IF (stepsize /= self%stepsize) CALL WARN('MC stepsize')
    outdim = nsteps / outfreq
    IF (outdim < self%outdim) THEN
      PRINT *, 'ERROR: trajectory in ', filename, ' only has ', outdim, ' samples, but you want ', self%outdim
      STOP 1
    ELSE IF (outdim > self%outdim) THEN
      PRINT *, 'INFO: reading only ', self%outdim, ' samples from the ', outdim, ' in ', filename
    END IF
    ! read in the trajectory -- note, partial reads of records is ok in Fortran
    ALLOCATE(self%Traj(3,natoms,self%outdim))
    ALLOCATE(self%Energy(self%outdim))
    READ(UNIT=fd, ERR=610) self%Traj
    READ(UNIT=fd, ERR=610) self%Energy
    CLOSE(UNIT=fd)
    RETURN

    600 CONTINUE
    PRINT *, 'ERROR: could not open file ', filename
    STOP 1

    610 CONTINUE
    PRINT *, 'ERROR: could not read from file ', filename
    STOP 1

  CONTAINS

    SUBROUTINE COMPLAIN(what)
      CHARACTER(LEN=*) :: what
      PRINT *, 'ERROR: parameter "', what, '" differs between input and ', filename
      STOP 1
    END SUBROUTINE COMPLAIN

    SUBROUTINE WARN(what)
      CHARACTER(LEN=*) :: what
      PRINT *, 'WARNING: parameter "', what, '" differs between input and ', filename
    END SUBROUTINE WARN

  END SUBROUTINE Read

END MODULE MODMC
