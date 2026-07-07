 PROGRAM acf
  use iso_c_binding, only: sp=>c_float, dp=>c_double
  IMPLICIT NONE
  INTEGER :: stdout
  REAL(DP), ALLOCATABLE :: p(:)
  REAL(DP), ALLOCATABLE :: c(:,:)
  REAL(DP), ALLOCATABLE :: d(:)
  REAL(DP), ALLOCATABLE :: a(:)
  REAL(DP), ALLOCATABLE :: Cm(:)
  REAL(DP) :: sum,tau
  INTEGER :: i, j, k, n, n_in, mm
  integer :: nargs, ierr
  CHARACTER(LEN=80) :: arg
  !
  character(len=80) :: filename
  integer :: DCinput1,DCinput2,Pinput,ACFoutput
  !
  !   default parameter 
  !
  n_in = 100
  !
  nargs = command_argument_count()
  if(nargs ==10) then
  do i = 1, nargs - 1
    CALL get_command_argument(i, arg)
    IF( TRIM( arg ) == '-n' ) THEN
       CALL get_command_argument(i+1, arg)
       READ( arg, * ) n_in
    END IF
    IF( TRIM( arg ) == '-i1' ) THEN
       CALL get_command_argument(i+1, arg)
       READ( arg, * ) filename
       DCinput1=86
       OPEN (UNIT=DCinput1, FILE=trim(filename), STATUS='old', &
        FORM='formatted', iostat=ierr)
       if(ierr /= 0) then
        write(6,*) ' Error opening filename '//trim(filename), ierr
        stop 1
       endif
    END IF
    IF( TRIM( arg ) == '-i2' ) THEN
       CALL get_command_argument(i+1, arg)
       READ( arg, * ) filename
       DCinput2=87
       OPEN (UNIT=DCinput2, FILE=trim(filename), STATUS='old', &
        FORM='formatted', iostat=ierr)
       if(ierr /= 0) then
        write(6,*) ' Error opening filename '//trim(filename), ierr
        stop 1
       endif
    END IF
    IF( TRIM( arg ) == '-i3' ) THEN
       CALL get_command_argument(i+1, arg)
       READ( arg, * ) filename
       Pinput=88
       OPEN (UNIT=Pinput, FILE=trim(filename), STATUS='old', &
        FORM='unformatted', iostat=ierr)
       if(ierr /= 0) then
        write(6,*) ' Error opening filename '//trim(filename), ierr
        stop 1
       endif
    END IF
    IF( TRIM( arg ) == '-o' ) THEN
       CALL get_command_argument(i+1, arg)
       READ( arg, * ) filename
       ACFoutput=89
       OPEN (UNIT=ACFoutput, FILE=trim(filename), STATUS='new', &
        FORM='formatted', iostat=ierr)
       if(ierr /= 0) then
        write(6,*) ' Error opening filename '//trim(filename), ierr
        stop 1
       endif
    END IF
  end do
  else
   write(6,*) 'Required 5 parameters:'
   write(6,*) 'n (-n), eigenvalues (-i1), eigenvectors (-i2), projections (-i3), '
   write(6,*) 'output table A,rate (-o)'
   stop 1
  end if

  n = n_in

  ALLOCATE( d( n ) )
  ALLOCATE( c( n, n ) )
  ALLOCATE( p( n ) )
  ALLOCATE( a( n ) )
  ALLOCATE( Cm( n ) )

  do i=1,n
   read(DCinput1,*) d(i)
  enddo  
  do i=1,n
   do j=1,n
    read(DCinput2,*) c(i,j)
!    write(6,'(f10.5)') c(i,j)
   enddo
  enddo
  close(DCinput1)
  close(DCinput2)
!  do i=1,n
!   read(Pinput,*) p(i)
!  enddo  
  read(Pinput) p
!  do i=1,n
!   write(6,*) 'p= ',p(i)
!  enddo
  close(Pinput)

  tau=0.d0
  a=0.d0
  mm=5.
  do i=1,N
   Cm = c(:,i)
   sum=DOT_PRODUCT(Cm,p)
   a(i)=mm*sum*sum
   write(89,*) 'A/lambda= ',a(i),d(i)
   tau = tau+a(i)/d(i)
!   write(6,*) a(i)/d(i)
  enddo
!  write(6,*) 'tau= ',tau

  close(ACFoutput)
  deallocate( p, d, c, Cm )

end program acf
