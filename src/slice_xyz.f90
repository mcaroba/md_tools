! Copyright 2020-2021 by Miguel A. Caro (mcaroba@gmail.com)

program slice_xyz

! This code takes in XYZ data and makes a slice along a desired axis,
! between a minimum and maximum value of the positions along that axis.
! It does not perform any checks related to periodic boundary conditions,
! i.e., a wrapped trajectory must be provided. Since the axes are the
! good old Cartesian ones (i.e., they're not the lattice vectors), this
! works best with orthorhombic unit cells or unit cells where the slice
! axis is parallel to one of the lattice vectors (e.g., the c axis in
! hexagonal systems).

  implicit none

  integer :: natoms0, natoms, iostatus, i_frame, i, axis
  character*8192 :: comment_line
  character*1024 :: file_in, file_out
  real*8, allocatable :: pos(:,:), vel(:,:)
  real*8 :: z_min, z_max
  character*256 :: argument
  character*8, allocatable :: element(:)
  character*1 :: axis_c = "z"
  logical, allocatable :: keep_atom(:)

!  read(*,*) file_in, file_out, z_min, z_max

  do i = 1, 8
    call getarg(i, argument)
    if( argument(1:8) == "file_in=" )then
      read(argument(9:256), *) file_in
    else if( argument(1:9) == "file_out=" )then
      read(argument(10:256), *) file_out
    else if( argument(1:5) == "axis=" )then
      read(argument(6:7), *) axis_c
    else if( argument(1:4) == "min=" )then
      read(argument(5:256), *) z_min
    else if( argument(1:4) == "max=" )then
      read(argument(5:256), *) z_max
    else if( argument == "" )then
      continue
    else
      write(0,*) "ERROR: I don't understand command line argument ", argument
      stop
    end if
  end do

! Check if axis is fine
  if( axis_c == "X" .or. axis_c == "x" )then
    axis = 1
  else if( axis_c == "Y" .or. axis_c == "y" )then
    axis = 2
  else if( axis_c == "Z" .or. axis_c == "z" )then
    axis = 3
  else
    write(0,*) "ERROR: the value of 'axis' must be 'x', 'y' or 'z' (default is 'z')"
    stop
  end if

  open(unit=10, file=file_in, status="old")
  open(unit=20, file=file_out, status="new")

  iostatus = 0
  i_frame = 0
  do
    i_frame = i_frame + 1
    read(10,*, iostat=iostatus) natoms0
    if( iostatus /= 0 )exit
    write(*,*) "Doing frame no. ", i_frame
    read(10, '(A)') comment_line
    allocate( pos(1:3, 1:natoms0) )
    allocate( vel(1:3, 1:natoms0) )
    allocate( element(1:natoms0) )
    allocate( keep_atom(1:natoms0) )
    keep_atom = .false.
    do i = 1, natoms0
      read(10, *) element(i), pos(1:3, i), vel(1:3, i)
    end do
    natoms = 0
    do i = 1, natoms0
      if( pos(axis, i) >= z_min .and. pos(axis,i) <= z_max )then
        keep_atom(i) = .true.
        natoms = natoms + 1
      end if
    end do
    write(20, *) natoms
    write(20, *) trim(adjustl(comment_line))
    do i = 1, natoms0
      if( keep_atom(i) )then
        write(20,*) element(i), pos(1:3, i), vel(1:3, i)
      end if
    end do
    deallocate( pos, vel, element, keep_atom )
  end do

  close(10)
  close(20)

end program
