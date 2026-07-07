program many_files

implicit none
integer :: ierr
integer :: ifile
integer :: n1,n2,dimI,dim12,dimII,n    ! matrix dimensions
integer :: nf,nfiles,start,freq        ! number of files to be read
real*8, allocatable :: f1(:),s1(:),f12(:),s12(:),f2(:),s2(:)
real*8, allocatable :: p(:)
real*8, allocatable :: f1g(:),s1g(:),f12g(:),s12g(:),f2g(:),s2g(:)
real*8, allocatable :: pg(:)
real*8 :: fact
character(len=80) :: fname,ffname,buffer,pfname
integer :: i_arg, argc
character(len=32) :: argv
integer :: i,iproj

!-------------------------------------------------------------------------------
!  Reading dimensions from standard input

n1 = 1
n2 = 1
nf = 1
start=1
freq=10
ffname = 'fsI.dump'
pfname = 'proj.dump'
argc = command_argument_count()
if (argc == 7) then
 call get_command_argument(1, argv)
 read(argv, *) n1
 call get_command_argument(2, argv)
  read(argv, *) n2
 call get_command_argument(3, argv)
  read(argv, *) ffname
 call get_command_argument(4, argv)
  read(argv, *) pfname
 call get_command_argument(5, argv)
  read(argv, *) nf
 call get_command_argument(6, argv)
  read(argv, *) start
 call get_command_argument(7, argv)
  read(argv, *) freq
else
 write(6,*) 'Required 7 parameters:'
 write(6,*) 'n1, n2, F/S file name, proj file name, number of directories, starting directory, and directories to skip.'
 stop 1
end if

dimI = n1*(n1+1)/2
dim12 = n1*n2
dimII = n2*(n2+1)/2
n=n1+n2

! Allocate matrices 

allocate( f1(dimI) )
allocate( s1(dimI) )
allocate( f1g(dimI) )
allocate( s1g(dimI) )
if(n2 /= 0) then
 allocate( f12(dim12) )
 allocate( s12(dim12) )
 allocate( f2(dimII) )
 allocate( s2(dimII) )
 allocate( f12g(dim12) )
 allocate( s12g(dim12) )
 allocate( f2g(dimII) )
 allocate( s2g(dimII) )
endif
allocate( p(n) )
allocate( pg(n) )

f1g = 0.d0
s1g = 0.d0
pg= 0.d0
if(n2 /= 0) then
 f12g = 0.d0
 s12g = 0.d0
 f2g = 0.d0
 s2g = 0.d0
endif
nfiles=0
do i=start,nf,freq
 nfiles=nfiles+1
 write(buffer,"(i4)") i
 fname = trim(buffer)//"/"//trim(ffname)
 fname = trim(adjustl(fname))
 write(6,"(a)") fname
 OPEN(1,file=fname,status="old",form="unformatted",iostat=ierr)
 ifile=1
 read(ifile) f1
 f1g = f1g + f1
 read(ifile) s1
 s1g = s1g + s1
 if(n2 /= 0) then
  read(ifile) f12
  f12g = f12g + f12
  read(ifile) s12
  s12g = s12g + s12
  read(ifile) f2
  f2g = f2g + f2
  read(ifile) s2
  s2g = s2g + s2
 endif
 CLOSE(ifile)

 iproj=0
 fname = trim(buffer)//"/"//trim(pfname)
 fname = trim(adjustl(fname))
 write(6,"(a)") fname
 OPEN(1,file=fname,status="old",form="unformatted",iostat=ierr)
 if(ierr==0) then
  iproj=1
  ifile=1
  read(ifile) p
  pg = pg + p
  CLOSE(ifile)
 else
  write(6,"(a)") "continuing without projections..."
 endif

enddo

!nfiles=4
fact = 1./DBLE(nfiles)
f1g = fact*f1g
s1g = fact*s1g
if (iproj==1) then
 pg = fact*pg
endif
if(n2 /= 0) then
 f12g = fact*f12g
 s12g = fact*s12g
 f2g = fact*f2g
 s2g = fact*s2g
endif

OPEN(2,file="fs.dump",status="new",form="unformatted",iostat=ierr)
write(2) f1g
write(2) s1g
if(n2 /= 0) then
 write(2) f12g
 write(2) s12g
 write(2) f2g
 write(2) s2g
endif
close(2)
if(iproj==1) then
 OPEN(2,file="proj.dump",status="new",form="unformatted",iostat=ierr)
 write(2) pg
 close(2)
endif

end program
