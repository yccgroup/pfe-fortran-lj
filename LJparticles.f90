MODULE MODLJ
  USE MODCONST
  IMPLICIT NONE
  
  TYPE :: LJ
      INTEGER :: natoms
      REAL*8  :: mass, epsilom, sigma
      REAL*8  :: L, rc, rc2, Ec, energy 
      REAL*8, ALLOCATABLE :: XYZ(:,:)
      LOGICAL :: uptodate ! whether energy needs updating after move
    CONTAINS
      PROCEDURE :: rdinp => lj_rdinp
      PROCEDURE :: init => lj_init
      PROCEDURE :: genXYZ
      PROCEDURE :: getenergy
      PROCEDURE :: calcr2
      PROCEDURE :: calcpairpot
      PROCEDURE :: calcenergy
      PROCEDURE :: partenergy
      PROCEDURE :: move
      !PROCEDURE :: partition_function_integral
      !PROCEDURE :: volume_integral 
  END TYPE LJ
  
  
  CONTAINS
  
    ! read input parameters
    SUBROUTINE lj_rdinp(self,fd)
      CLASS(LJ) :: self
      INTEGER :: fd

      READ(fd,*) self%natoms
      READ(fd,*) self%mass
      READ(fd,*) self%epsilom
      READ(fd,*) self%sigma
      READ(fd,*) self%L
      READ(fd,*) self%rc   
    END SUBROUTINE lj_rdinp

    ! initialize the system parameters and coordinates
    SUBROUTINE lj_init(self)
      CLASS(LJ) :: self
      REAL*8 :: sig_rc

      self%rc2 = self%rc**2
      sig_rc = self%sigma/self%rc
      self%Ec = 4*self%epsilom*(sig_rc**12 - sig_rc**6)

      ALLOCATE(self%XYZ(3,self%natoms))
      CALL self%genXYZ()
    END SUBROUTINE lj_init

    ! initialize the LJ particles' coordinates
    SUBROUTINE genXYZ(self)
      CLASS(LJ) :: self

      CALL RANDOM_NUMBER(self%XYZ)
      self%XYZ = self%XYZ * self%L
      self%uptodate = .FALSE.
    END SUBROUTINE genXYZ

    ! return the energy
    FUNCTION getenergy(self) RESULT(energy)
      CLASS(LJ) :: self
      REAL*8 :: energy
      IF (.NOT. self%uptodate) CALL self%calcenergy()
      energy = self%energy
    END FUNCTION getenergy

    ! calculate r2 between particle i and j
    FUNCTION calcr2(self,i,j) RESULT(r2)
      CLASS(LJ) :: self
      INTEGER :: i, j, k
      REAL*8 :: r2
      REAL*8 :: dist(3)
    
      dist = self%XYZ(:,i) - self%XYZ(:,j)

      DO k = 1, 3
        IF (dist(k) > 0.5d0*self%L) THEN
          dist(k) = dist(k) - self%L
        ELSE IF (dist(k) < -0.5d0*self%L) THEN
          dist(k) = dist(k) + self%L
        END IF
      END DO
     
      r2 = DOT_PRODUCT(dist,dist)

    END FUNCTION calcr2


    ! calculate a single pair potential
    FUNCTION calcpairpot(self,r2) RESULT(pot)
      CLASS(LJ) :: self
      REAL*8 :: r2, r_2
      REAL*8 :: r_6, r_12, pot
      
      ! potential cutoff
      IF (r2 >= self%rc2) THEN
        pot = 0.d0
        RETURN
      END IF

      r_2 = self%sigma**2 / r2
      r_6  = r_2 * r_2 * r_2
      r_12 = r_6 * r_6
      pot  = 4 * self%epsilom * (r_12 - r_6) - self%Ec
    END FUNCTION calcpairpot


    ! calculate the total potential energy of the system
    SUBROUTINE calcenergy(self)
      CLASS(LJ) :: self
      INTEGER :: i, j
      REAL*8 :: r2

      self%energy =0.d0
      DO i = 1, self%natoms
        DO j = 1, i-1
          r2 = self%calcr2(i,j)
          self%energy = self%energy + self%calcpairpot(r2)
        END DO
      END DO
      self%uptodate = .TRUE.
    END SUBROUTINE calcenergy

    ! calculate the potential energy contribution from atom aid
    FUNCTION partenergy(self,aid) RESULT(myenergy)
      CLASS(LJ) :: self
      INTEGER,INTENT(IN) :: aid
      REAL*8 :: myenergy, r2
      INTEGER :: i

      myenergy = 0.d0
      DO i = 1, self%natoms
        IF (i == aid) CYCLE
        r2 = self%calcr2(i,aid)
        myenergy = myenergy + self%calcpairpot(r2)
      END DO
    END FUNCTION partenergy

    ! move particles (3D)    
    SUBROUTINE move(self,aid,stepsize)
      CLASS(LJ) :: self
      INTEGER,INTENT(IN) :: aid
      REAL*8,INTENT(IN) :: stepsize
      INTEGER :: i
      REAL*8 :: Xold(3), Xnew(3), rand(3), step(3)
      REAL*8 :: partEold, partEnew

      Xold = self%XYZ(:,aid)
      partEold = self%partenergy(aid)

      CALL RANDOM_NUMBER(rand)
      step = 2 * (rand - 0.5d0) * stepsize
      Xnew = Xold + step

      ! place the particle back to the box if it runs out
      DO i = 1, 3
        IF (Xnew(i) > self%L) THEN
          Xnew(i) = Xnew(i) - self%L
        ELSE IF (Xnew(i) < 0) THEN
          Xnew(i) = Xnew(i) + self%L
        END IF
      END DO

      ! assign the final coordinate
      self%XYZ(:,aid) = Xnew

      ! update the system energy
      partEnew = self%partenergy(aid)
      self%energy = self%energy + (partEnew - partEold)
      self%uptodate = .TRUE.

    END SUBROUTINE move
    
    ! Debug: Up to here 

!   ! calculate the partition function via direct numerical integration    
!   FUNCTION partition_function_integral(self,upper,lower,dx) RESULT(Z)
!     CLASS(LJ) :: self
!     INTEGER :: i
!     INTEGER :: ndim, id1, id2, ndof
!     INTEGER, ALLOCATABLE :: indexs(:) 
!     REAL*8 :: upper, lower, dx
!     REAL*8 :: Z, dv, r2, V
!     REAL*8 :: x1(3), x2(3), dist(3)
!     REAL*8, ALLOCATABLE :: coords(:)
!     
!     ! initialization
!     ndim = NINT((upper-lower)/dx)
!     ndof = 3*self%natoms
!     Z = 0.d0
!     dv = dx**ndof

!     ! allocate the arrays
!     ALLOCATE(indexs(ndof),coords(ndof))
!     indexs(:) = 0
!     coords(:) = lower

!     ! integration
!     DO WHILE (.TRUE.)
!       V = 0.d0
!       DO id1 = 0, self%natoms-1
!         x1 = coords(3*id1+1 : 3*id1+3)
!         DO id2 = 0, id1-1
!           x2 = coords(3*id2+1 : 3*id2+3)
!           dist = x2 - x1
!           r2 = self%calcr2(dist)
!           V = V + self%calcpairpot(r2)
!         END DO
!       END DO
!       Z  = Z + EXP(-beta*V)
!       DO i = 1, ndof
!          indexs(i) = indexs(i) + 1
!          coords(i) = lower + indexs(i)*dx
!          IF (indexs(i) >= ndim) THEN
!            indexs(i) = 0
!            coords(i) = lower
!          ELSE
!            EXIT
!          END IF
!       END DO
!       IF (ALL(indexs==0)) THEN
!         EXIT
!       END IF

!     END DO

!     Z = Z * dv

!   END FUNCTION partition_function_integral


!   ! calculate the volume via direct numerical integration    
!   FUNCTION volume_integral(self,upper,lower,dx,threshold) RESULT(Volume)
!     CLASS(LJ) :: self
!     INTEGER :: i
!     INTEGER :: ndim, id1, id2, ndof
!     INTEGER, ALLOCATABLE :: indexs(:) 
!     REAL*8 :: upper, lower, dx, threshold
!     REAL*8 :: dv, r2, V, Volume
!     REAL*8 :: x1(3), x2(3), dist(3)
!     REAL*8, ALLOCATABLE :: coords(:)
!     
!     ! initialization
!     ndim = NINT((upper-lower)/dx)
!     ndof = 3*self%natoms
!     Volume = 0.d0
!     dv = dx**ndof

!     ! allocate the arrays
!     ALLOCATE(indexs(ndof),coords(ndof))
!     indexs(:) = 0
!     coords(:) = lower

!     ! integration
!     DO WHILE (.TRUE.)
!       V = 0.d0
!       DO id1 = 0, self%natoms-1
!         x1 = coords(3*id1+1 : 3*id1+3)
!         DO id2 = 0, id1-1
!           x2 = coords(3*id2+1 : 3*id2+3)
!           dist = x2 - x1
!           r2 = self%calcr2(dist)
!           V = V + self%calcpairpot(r2)
!         END DO
!       END DO
!       ! count the volume if the energy < threshold
!       IF (V < threshold) THEN
!         Volume = Volume + 1
!       END IF
!       DO i = 1, ndof
!          indexs(i) = indexs(i) + 1
!          coords(i) = lower + indexs(i)*dx
!          IF (indexs(i) >= ndim) THEN
!            indexs(i) = 0
!            coords(i) = lower
!          ELSE
!            EXIT
!          END IF
!       END DO
!       IF (ALL(indexs==0)) THEN
!         EXIT
!       END IF

!     END DO

!     Volume = Volume * dv

!   END FUNCTION volume_integral
         



END MODULE MODLJ
