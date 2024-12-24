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

    ! use Sedgewick's "median of three" pivot
    mid = (n+1)/2  ! middle index
    IF (data(mid) < data(n))  CALL swap(n,mid)
    IF (data(1)   < data(n))  CALL swap(1,n)
    IF (data(mid) < data(1))  CALL swap(mid,1)
    pivot = data(1)
    ! partition the array
    i = 0
    j = n+1
    DO
      DO
        i = i + 1
        IF (data(i) >= pivot) EXIT
      END DO
      DO
        j = j - 1
        IF (data(j) <= pivot) EXIT
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
