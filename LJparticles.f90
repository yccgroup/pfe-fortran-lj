MODULE MODLJ
  USE MODTEMPAR 
  IMPLICIT NONE
  
  TYPE :: LJ
      INTEGER :: natoms
      REAL*8  :: mass, epsilom, sigma
      REAL*8  :: L, rc, rc2, Ec, energy 
      REAL*8, ALLOCATABLE :: XYZ(:,:)
    CONTAINS
      PROCEDURE :: genXYZ
      PROCEDURE :: calcpairpot
      PROCEDURE :: calcenergy
      PROCEDURE :: partition_function_integral
      PROCEDURE :: volume_integral 
      PROCEDURE :: move
      PROCEDURE :: calcr2
  END TYPE LJ
  
  
  CONTAINS
  

    ! initialize the LJ particles' coordinates
    SUBROUTINE genXYZ(self)
      CLASS(LJ) :: self

      CALL RANDOM_NUMBER(self%XYZ)
      self%XYZ = self%XYZ * self%L
    END SUBROUTINE genXYZ


    ! calculate r2
    FUNCTION calcr2(self,dist) RESULT(r2)
      CLASS(LJ) :: self
      INTEGER :: k
      REAL*8 :: r2
      REAL*8 :: dist(3)
    
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
      REAL*8 :: r2
      REAL*8 :: r_6, r_12, pot
      
      ! potential cutoff
      IF (r2 >= self%rc2) THEN
        pot = 0.d0
        RETURN
      END IF

      r_6  = 1.d0 / (r2*r2*r2)
      r_12 = r_6 * r_6
      pot  = 4 * (r_12 - r_6) - self%Ec
    END FUNCTION calcpairpot


    ! calculate the total potential energy of the system
    SUBROUTINE calcenergy(self)
      CLASS(LJ) :: self
      INTEGER :: i, j
      REAL*8 :: dist(3), r2

      self%energy =0.d0
      DO i = 1, self%natoms
        DO j = 1, i-1
          dist = self%XYZ(:,i) - self%XYZ(:,j)
          r2 = self%calcr2(dist)
          self%energy = self%energy + self%calcpairpot(r2)
        END DO
      END DO
    END SUBROUTINE calcenergy


    ! move particles (3D)    
    SUBROUTINE move(self,x,stepsize) 
      CLASS(LJ) :: self
      INTEGER :: i
      REAL*8 :: stepsize
      REAL*8 :: x(3), pos(3), rand(3), step(3)

      CALL RANDOM_NUMBER(rand)
      step = 2 * (rand - 0.5d0) * stepsize
      pos = x + step

      ! place the particle back to the box if it runs out
      DO i = 1, 3
        IF (pos(i) > self%L) THEN
          pos(i) = pos(i) - self%L
        ELSE IF (pos(i) < 0) THEN
          pos(i) = pos(i) + self%L
        END IF
      END DO

      ! assign the final coordinate for output
      x = pos

    END SUBROUTINE move
    
    ! Debug: Up to here 

    ! calculate the partition function via direct numerical integration    
    FUNCTION partition_function_integral(self,upper,lower,dx) RESULT(Z)
      CLASS(LJ) :: self
      INTEGER :: i
      INTEGER :: ndim, id1, id2, ndof
      INTEGER, ALLOCATABLE :: indexs(:) 
      REAL*8 :: upper, lower, dx
      REAL*8 :: Z, dv, r2, V
      REAL*8 :: x1(3), x2(3), dist(3)
      REAL*8, ALLOCATABLE :: coords(:)
      
      ! initialization
      ndim = NINT((upper-lower)/dx)
      ndof = 3*self%natoms
      Z = 0.d0
      dv = dx**ndof

      ! allocate the arrays
      ALLOCATE(indexs(ndof),coords(ndof))
      indexs(:) = 0
      coords(:) = lower

      ! integration
      DO WHILE (.TRUE.)
        V = 0.d0
        DO id1 = 0, self%natoms-1
          x1 = coords(3*id1+1 : 3*id1+3)
          DO id2 = 0, id1-1
            x2 = coords(3*id2+1 : 3*id2+3)
            dist = x2 - x1
            r2 = self%calcr2(dist)
            V = V + self%calcpairpot(r2)
          END DO
        END DO
        Z  = Z + EXP(-beta*V)
        DO i = 1, ndof
           indexs(i) = indexs(i) + 1
           coords(i) = lower + indexs(i)*dx
           IF (indexs(i) >= ndim) THEN
             indexs(i) = 0
             coords(i) = lower
           ELSE
             EXIT
           END IF
        END DO
        IF (ALL(indexs==0)) THEN
          EXIT
        END IF

      END DO

      Z = Z * dv

    END FUNCTION partition_function_integral


    ! calculate the volume via direct numerical integration    
    FUNCTION volume_integral(self,upper,lower,dx,threshold) RESULT(Volume)
      CLASS(LJ) :: self
      INTEGER :: i
      INTEGER :: ndim, id1, id2, ndof
      INTEGER, ALLOCATABLE :: indexs(:) 
      REAL*8 :: upper, lower, dx, threshold
      REAL*8 :: dv, r2, V, Volume
      REAL*8 :: x1(3), x2(3), dist(3)
      REAL*8, ALLOCATABLE :: coords(:)
      
      ! initialization
      ndim = NINT((upper-lower)/dx)
      ndof = 3*self%natoms
      Volume = 0.d0
      dv = dx**ndof

      ! allocate the arrays
      ALLOCATE(indexs(ndof),coords(ndof))
      indexs(:) = 0
      coords(:) = lower

      ! integration
      DO WHILE (.TRUE.)
        V = 0.d0
        DO id1 = 0, self%natoms-1
          x1 = coords(3*id1+1 : 3*id1+3)
          DO id2 = 0, id1-1
            x2 = coords(3*id2+1 : 3*id2+3)
            dist = x2 - x1
            r2 = self%calcr2(dist)
            V = V + self%calcpairpot(r2)
          END DO
        END DO
        ! count the volume if the energy < threshold
        IF (V < threshold) THEN
          Volume = Volume + 1
        END IF
        DO i = 1, ndof
           indexs(i) = indexs(i) + 1
           coords(i) = lower + indexs(i)*dx
           IF (indexs(i) >= ndim) THEN
             indexs(i) = 0
             coords(i) = lower
           ELSE
             EXIT
           END IF
        END DO
        IF (ALL(indexs==0)) THEN
          EXIT
        END IF

      END DO

      Volume = Volume * dv

    END FUNCTION volume_integral
         



END MODULE MODLJ
