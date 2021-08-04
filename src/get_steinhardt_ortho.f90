! Copyright (c) 2020-2021 by Miguel A. Caro

program get_steinhardt_ortho

  use neighbors
  use potentials
  use analyze

  implicit none

  real*8 :: Lv(1:3), rmax = 0.d0
  integer :: Np, iostatus, i, image, lmax = 10, l, every = 1, counter
  character*1024 :: filename = ""
  character*64 :: el
! Modes can be:
!
! atoms_mode = all|each
!   * all: it returns the average Ql for all atoms
!   * each: it returns the Ql for each atom
!
! trajectory_mode = all|each
!   * all: it returns the Ql averaged over all trajectory snapshots
!   * each: it returns the Ql for each trajectory snapshots
  character*64 :: atoms_mode = "all", trajectory_mode = "all"
  real*8, allocatable :: pos(:,:), Ql(:,:), Qlsum(:,:)
  integer, allocatable :: nn(:), nn_list(:,:)
  logical :: PBC(1:3)
  real*8 :: time1, time2
  integer :: n_args, j
  character*1024 :: cmd_string, junk
  logical :: Lv_from_file = .true.

  PBC = .true.
  Lv = 0.d0

!  This was a shitty interface; I'm obsoleting it
!  read(*,*) filename, Lv(1:3), rmax, lmax, every, atoms_mode, trajectory_mode

!  This is more elegant :-D
  n_args = command_argument_count()
  do i = 1, n_args
    call getarg(i, cmd_string)
    write(cmd_string, '(A)') adjustl(cmd_string)
    if( cmd_string(1:5) == "rmax=" )then
      read(cmd_string(6:1024), *) rmax
    else if( cmd_string(1:5) == "lmax=" )then
      read(cmd_string(6:1024), *) lmax
    else if( cmd_string(1:6) == "every=" )then
      read(cmd_string(7:1024), *) every
    else if( cmd_string(1:15) == "atoms_filename=" )then
      read(cmd_string(16:1024), '(A)') filename
      filename = adjustl(filename)
    else if( cmd_string(1:11) == "atoms_mode=" )then
      read(cmd_string(12:1024), '(A)') atoms_mode
    else if( cmd_string(1:16) == "trajectory_mode=" )then
      read(cmd_string(17:1024), '(A)') trajectory_mode
    else if( cmd_string(1:4) == "box=" )then
      read(cmd_string(5:1024), '(A)') cmd_string
      read(cmd_string, *) Lv(1:3)
      Lv_from_file = .false.
    end if
  end do

! Checks
  if( filename == "" )then
    write(*,*) "ERROR: you need to specify an atoms_filename!"
    stop
  else if( rmax <= 0.d0 )then
    write(*,*) "ERROR: you need to specify an rmax > 0!"
    stop
  else if( .not. atoms_mode == "all" .and. .not. atoms_mode == "each" )then
    write(*,*) "ERROR: atoms_mode must be ""all"" or ""each"""
    stop
  else if( .not. trajectory_mode == "all" .and. .not. trajectory_mode == "each" )then
    write(*,*) "ERROR: trajectory_mode must be ""all"" or ""each"""
    stop
  end if

  open(unit=10, file=filename, status="old", iostat=iostatus)
  image = 0
  counter = 0
  do while( iostatus == 0 )
    image = image + 1
    read(10, *) Np
!   Allocate arrays
    if( counter == 0 )then
      if( atoms_mode == "all" )then
        allocate( Ql(0:lmax, 1) )
        allocate( Qlsum(0:lmax, 1) )
      else if( atoms_mode == "each" )then
        allocate( Ql(0:lmax, Np) )
        allocate( Qlsum(0:lmax, Np) )
      end if
      Qlsum = 0.d0
    end if
    read(10, '(A)') cmd_string
    if( Lv_from_file )then
      do i = 1, 1024-8
        if( cmd_string(i:i+4) == "box=""" )then
          do j = i+5, 1024
            if( cmd_string(j:j) == """" )exit
          end do
          read(cmd_string(i+5:j-1), *) Lv(1:3)
          exit
        end if
      end do
    end if
    if( any( Lv == 0.d0 ) )then
      write(*,*) "ERROR: box size either zero or unspecified!"
      stop
    end if
    if( image == 1 )then
      allocate( pos(1:3, 1:Np) )
      allocate( nn(1:Np) )
      allocate( nn_list(1:Np, 1:Np) )
    end if
!call cpu_time(time1)
    do i = 1, Np
      read(10, *, iostat=iostatus) el, pos(1:3, i)
      read(10, *, iostat=iostatus)
      backspace(10)
    end do
!call cpu_time(time2)
!write(*,*) "Read: ", time2-time1
    if( mod(image, every) == 0 )then
!call cpu_time(time1)
      call build_neighbors(pos, Np, Lv, PBC, rmax, nn, nn_list) 
!call cpu_time(time2)
!write(*,*) "Neighbors: ", time2-time1
      if( atoms_mode == "all" )then
!$omp parallel do
        do l = 0, lmax
          call get_Ql(pos, 0, Lv, PBC, nn, nn_list, rmax, l, Ql(l,1))
        end do
!$omp end parallel do
      else if( atoms_mode == "each" )then
!$omp parallel do private(l)
        do i = 1, Np
          do l = 0, lmax
            call get_Ql(pos, i, Lv, PBC, nn, nn_list, rmax, l, Ql(l,i))
          end do
        end do
!$omp end parallel do
      end if
!     Done for each trajectory frame
      if( trajectory_mode == "each" )then
        do i = 1, size(Ql, 2)
          write(*,*) Ql(0:lmax, i)
        end do
        write(*,*)
      else
        Qlsum = Qlsum + Ql
      end if
      counter = counter + 1
    end if
  end do
  close(10)

! Only done at the end of the trajectory
  if( trajectory_mode == "all" )then
    Qlsum = Qlsum / dfloat(counter)
    if( atoms_mode == "all" )then
      write(*,*) Qlsum(0:lmax, 1)
    else if( atoms_mode == "each" )then
      do i = 1, Np
        write(*,*) Qlsum(0:lmax, i)
      end do
    end if
  end if

end program
