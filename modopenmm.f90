! Module for reading OpenMM DCD & log files.
!
! Some other variant of DCD is documented at
! https://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dcdplugin.html
! but OpenMM's is somewhat different (e.g. REAL*4 instead of REAL*8)

MODULE MODOpenMM

  IMPLICIT NONE

  CHARACTER(LEN=4), PARAMETER :: MAGIC = 'CORD'
  INTEGER, PARAMETER :: MAXTITL = 80

  TYPE :: DCDFile
    INTEGER :: fd     ! Fortran unit number
    INTEGER :: nset   ! number of frame
    INTEGER :: istrt  ! step of first frame
    INTEGER :: nsavc  ! step interval between frames
    REAL*4  :: delta  ! time step delta (???)
    INTEGER :: ntitle ! number of title lines
    CHARACTER(LEN=MAXTITL),ALLOCATABLE :: title(:)
    INTEGER :: natom  ! number of atoms
  CONTAINS
    PROCEDURE :: close => dcd_close
    PROCEDURE :: info => dcd_info
    PROCEDURE :: read_frame => dcd_read_frame
  END TYPE DCDFile

CONTAINS

  ! Open DCD file and read header information.
  FUNCTION dcd_open(filename) RESULT(dcd)
    CHARACTER(LEN=*),INTENT(IN) :: filename
    TYPE(DCDFile) :: dcd
    INTEGER :: fd, i, c
    CHARACTER(LEN=4) :: hdr
    INTEGER :: junk(6)

    OPEN(FILE=filename, NEWUNIT=fd, FORM='unformatted', STATUS='old', ACTION='read', ERR=500)
    dcd%fd = fd

    ! read header
    READ(UNIT=fd, ERR=510) hdr, dcd%nset, dcd%istrt, dcd%nsavc, junk, dcd%delta
    IF (hdr /= MAGIC) THEN
      PRINT *, 'ERROR: wrong DCD header in ', filename
      STOP 1
    END IF

    ! read title
    READ(UNIT=fd, ERR=510) dcd%ntitle
    BACKSPACE(UNIT=fd)
    ALLOCATE(dcd%title(dcd%ntitle))
    READ(UNIT=fd, ERR=510) dcd%ntitle, dcd%title
    ! fix up null bytes
    DO i = 1, dcd%ntitle
      DO c = 1, MAXTITL
        IF (dcd%title(i)(c:c) == ACHAR(0)) dcd%title(i)(c:c) = ' '
      END DO
    END DO

    ! read natom
    READ(UNIT=fd, ERR=510) dcd%natom

    ! file is now positioned at the first frame
    RETURN

    500 CONTINUE
    PRINT *, 'ERROR: could not open ', filename
    STOP 1

    510 CONTINUE
    PRINT *, 'ERROR: could not read from ', filename
    STOP 1

  END FUNCTION dcd_open


  ! Close the DCD file.
  SUBROUTINE dcd_close(self)
    CLASS(DCDFile) :: self

    DEALLOCATE(self%title)
    CLOSE(UNIT=self%fd)
    self%fd = -1
  END SUBROUTINE dcd_close


  ! Print information about DCD file to given unit.
  SUBROUTINE dcd_info(self, fd)
    CLASS(DCDFile) :: self
    INTEGER,INTENT(IN) :: fd
    INTEGER :: i

    DO i = 1, self%ntitle
      WRITE(UNIT=fd, FMT='(a,i1,2a)') '# Title ', i, ': ', TRIM(self%title(i))
    END DO
    WRITE(UNIT=fd, FMT='(a,i0)') '# Number of frames   : ', self%nset
    WRITE(UNIT=fd, FMT='(a,i0)') '# Step of 1st frame  : ', self%istrt
    WRITE(UNIT=fd, FMT='(a,i0)') '# Step between frames: ', self%nsavc
    !WRITE(UNIT=fd, FMT='(a,f15.6)') 'dt : ', self%delta
    WRITE(UNIT=fd, FMT='(a,i0)') '# Number of atoms    : ', self%natom
   END SUBROUTINE dcd_info


  ! Read a frame from the DCD file into the provided array.
  SUBROUTINE dcd_read_frame(self, xyz)
    CLASS(DCDFile) :: self
    REAL*8,INTENT(INOUT) :: xyz(:,:)
    REAL*4,ALLOCATABLE :: coords(:)
    INTEGER :: i, k

    ALLOCATE(coords(self%natom))
    DO k = 1, 3
      READ(UNIT=self%fd) coords
      DO i = 1, self%natom
        xyz(k,i) = DBLE(coords(i))
      END DO
    END DO
    DEALLOCATE(coords)
  END SUBROUTINE dcd_read_frame


  ! Read OpenMM log file (StateDataReporter).
  ! Assuming the entries are:
  ! #"Step","Potential Energy (kJ/mole)",...
  SUBROUTINE read_openmm_log(filename, nset, istrt, nsavc, energies)
    CHARACTER(LEN=*),INTENT(IN) :: filename
    INTEGER,INTENT(IN) :: nset
    INTEGER,INTENT(OUT) :: istrt, nsavc
    REAL*8,INTENT(OUT) :: energies(nset)
    INTEGER :: fd, i, step
    REAL*8 :: energy
    CHARACTER(LEN=240) :: line

    OPEN(FILE=filename, NEWUNIT=fd, STATUS='old', ACTION='read', ERR=600)
    istrt = -1
    nsavc = -1
    i = 1
    DO WHILE (i <= nset)
      READ(UNIT=fd, FMT='(a)', ERR=610) line
      IF (line(1:1) == '#') CYCLE
      READ(line, *) step, energy
      IF (istrt == -1) THEN
        istrt = step
      ELSEIF (nsavc == -1) THEN
        nsavc = step - istrt
      ENDIF
      energies(i) = energy
      i = i + 1
    END DO
    CLOSE(UNIT=fd)
    RETURN

    600 CONTINUE
    PRINT *, 'ERROR: could not open ', filename
    STOP 1

    610 CONTINUE
    PRINT *, 'ERROR: could not read from ', filename
    STOP 1

  END SUBROUTINE read_openmm_log

END MODULE MODOpenMM
