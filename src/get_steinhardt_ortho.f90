! Copyright (c) 2020-2021 by Miguel A. Caro

program get_steinhardt_ortho

  use neighbors
  use potentials
  use analyze

  implicit none

  real*8 :: Lv(1:3), rmax
  integer :: Np, iostatus, i, image, lmax, l, every, counter
  character*1024 :: filename
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

  PBC = .true.

  read(*,*) filename, Lv(1:3), rmax, lmax, every, atoms_mode, trajectory_mode

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
    read(10, *)
    if( image == 1 )then
      allocate( pos(1:3, 1:Np) )
      allocate( nn(1:Np) )
      allocate( nn_list(1:Np, 1:Np) )
    end if
    do i = 1, Np
      read(10, *, iostat=iostatus) el, pos(1:3, i)
      read(10, *, iostat=iostatus)
      backspace(10)
    end do
    if( mod(image, every) == 0 )then
      call build_neighbors(pos, Np, Lv, PBC, rmax, nn, nn_list) 
      if( atoms_mode == "all" )then
        do l = 0, lmax
          call get_Ql(pos, 0, Lv, PBC, nn, nn_list, rmax, l, Ql(l,1))
        end do
      else if( atoms_mode == "each" )then
        do i = 1, Np
          do l = 0, lmax
            call get_Ql(pos, i, Lv, PBC, nn, nn_list, rmax, l, Ql(l,i))
          end do
        end do
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
