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
      Eold = System%energy
    
      ! move the selected atom and calculate the new energy
      CALL System%move(System%XYZ(:,aid),self%stepsize)
      CALL System%calcenergy()
      Xnew = System%XYZ(:,aid)
      Enew = System%energy
    
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
  

END MODULE MODMC
