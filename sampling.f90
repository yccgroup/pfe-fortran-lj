! Module for Monte Carlo Algorithm
MODULE MODMC
  USE MODTEMPAR
  USE MODPBCLJ
  USE MODUTIL
  IMPLICIT NONE

  TYPE MCPAR
    INTEGER :: nsteps, outfreq, outdim, newoutdim
    REAL*8  :: stepsize
    REAL*8  :: Accept
    REAL*8, ALLOCATABLE :: Traj(:,:,:), Energy(:), NewTraj(:,:,:), NewEnergy(:)
  END TYPE MCPAR    

  CONTAINS
  !-----------------------------------------------------------------------------
  !
  ! Subroutines for Monte Carlo Sampling 
  !
  !-----------------------------------------------------------------------------
  
  ! Monte Carlo Sampling
  SUBROUTINE MCSampling(System,MC)
    TYPE(LJ) System
    TYPE(MCPAR) MC
    INTEGER :: i
    INTEGER :: natoms, nsteps, outdim, outfreq
    INTEGER :: acount, aid, outindex
    REAL*8 :: stepsize, Accept
    REAL*8 :: Eold, Enew
    REAL*8 :: Xold(3), Xnew(3)
  
    ! assign parameters for simplicity
    natoms = System%natoms
    nsteps = MC%nsteps
    outdim = MC%outdim
    outfreq = MC%outfreq
    stepsize = MC%stepsize
  
    ! initialize and allcoate arrays
    Accept = 0
    outindex = 0
    ALLOCATE(MC%Traj(3,natoms,outdim),MC%Energy(outdim))
    
    ! propagation
    DO i = 1, nsteps
  
      ! pick one atom to move
      aid = RANDOM_INTEGER(natoms)

      ! back up energy and coordinates
      !CALL System%calcenergy()  ! DEBUG: see if turn this on/off, the result change
      Xold = System%XYZ(:,aid)
      Eold = System%energy
    
      ! move the selected atom and calculate the new energy
      CALL System%move(System%XYZ(:,aid),stepsize)
      CALL System%calcenergy()
      Xnew = System%XYZ(:,aid)
      Enew = System%energy
    
      ! Metropolis (Xnew, Enew updated)
      CALL Metropolis(Xold,Eold,Xnew,Enew,acount)
      System%XYZ(:,aid) = Xnew
      System%energy = Enew
      Accept = Accept + acount
    
      ! output
      IF (MOD(i,outfreq) == 0) THEN
         outindex = outindex +1 
         MC%Energy(outindex) = Enew
         MC%Traj(:,:,outindex) = System%XYZ
      END IF
  
    END DO
  
    ! for statisitical output
    MC%Accept = Accept*1.d0/nsteps
  
  END SUBROUTINE MCSampling
  
  
  ! Determine Xnew and Enew by Metropolis algorithm
  SUBROUTINE Metropolis(Xold,Eold,Xnew,Enew,acount)
    INTEGER :: acount
    REAL*8 :: Eold, Enew, acceptance, rand
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


END MODULE MODMC
