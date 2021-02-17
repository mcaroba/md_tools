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
  real*8, allocatable :: pos(:,:), Ql(:), Qlsum(:)
  integer, allocatable :: nn(:), nn_list(:,:)
  logical :: PBC(1:3)

  PBC = .true.

  read(*,*) filename, Lv(1:3), rmax, lmax, every

  allocate( Ql(0:lmax) )
  allocate( Qlsum(0:lmax) )
  Qlsum = 0.d0

  open(unit=10, file=filename, status="old", iostat=iostatus)
  image = 0
  counter = 0
  do while( iostatus == 0 )
    image = image + 1
    read(10, *) Np
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
      do l = 0, lmax
        call get_Ql(pos, Lv, PBC, nn, nn_list, rmax, l, Ql(l))
      end do
      Qlsum = Qlsum + Ql
      counter = counter + 1
    end if
  end do
  close(10)

  Qlsum = Qlsum / dfloat(counter)
  write(*,*) Qlsum(0:lmax)

end program
