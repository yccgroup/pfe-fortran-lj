MODULE MODSORT

  IMPLICIT NONE

CONTAINS

  ! Sort the array in ascending order.
  RECURSIVE SUBROUTINE qsort(data)
    REAL*8,INTENT(INOUT) :: data(:)
    REAL*8 :: pivot
    INTEGER :: n, mid, i, j

    n = SIZE(data)
    IF (n <= 1) RETURN
    IF (n == 2) THEN
      IF (data(2) < data(1))  CALL swap(1,2)
      RETURN
    END IF

    mid = (n+1)/2  ! middle index
    ! use Sedgewick's "median of three" pivot
    pivot = median3(data(1), data(mid), data(n))
    ! partition the array
    i = 1
    j = n
    DO
      DO WHILE (data(i) < pivot)
        i = i + 1
      END DO
      DO WHILE (data(j) > pivot)
        j = j - 1
      END DO
      IF (i >= j) EXIT
      CALL swap(i,j)
    END DO
    ! recurse
    CALL qsort(data(1:j))
    CALL qsort(data(j+1:n))

  CONTAINS
  
    SUBROUTINE swap(i,j)
      INTEGER,INTENT(in) :: i,j
      REAL*8 :: tmp
      tmp = data(i)
      data(i) = data(j)
      data(j) = tmp
    END SUBROUTINE swap

    FUNCTION median3(a,b,c) RESULT (med)
      REAL*8,INTENT(IN) :: a,b,c
      REAL*8 :: med
      ! Is a the median? [b,a,c] or [c,a,b]
      IF ( ((b<=a).AND.(a<=c)) .OR. ((c<=a).AND.(a<=b)) ) THEN
        med = a
      ! Is c the median? [b,c,a] or [a,c,b]
      ELSEIF ( ((b<=c).AND.(c<=a)) .OR. ((a<=c).AND.(c<=b)) ) THEN
        med = c
      ELSE
        med = b
      END IF
    END FUNCTION median3

  END SUBROUTINE qsort


  ! Compute the desired quantile (p between 0 and 1) of the data.
  FUNCTION quantile(data, p) RESULT (q)
    REAL*8,INTENT(IN) :: data(:)
    REAL*8,INTENT(IN) :: p
    REAL*8 :: q
    INTEGER :: ndata, idx
    REAL*8,ALLOCATABLE :: sorted(:)

    ndata = SIZE(data)
    ! copy the data and sort it
    ALLOCATE(sorted, SOURCE=data)
    CALL qsort(sorted)
    ! calculate the quantile index
    idx = NINT(1 + p*(ndata-1))
    ! return the result
    q = sorted(idx)
  END FUNCTION quantile

END MODULE MODSORT
