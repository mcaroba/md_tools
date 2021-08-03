! Copyright (c) 2018 by Miguel A. Caro

module neighbors

use potentials

contains

!**********************************************************************************************
subroutine build_neighbors(pos, Np, L, PBC, R, nn, list)

! This routine 

  implicit none

  real*8, intent(in) :: pos(:,:), R, L(1:3)
  integer, intent(in) :: Np
  logical, intent(in) :: PBC(1:3)
  integer, intent(inout) :: list(:,:)
  integer, intent(inout) :: nn(:)
  integer :: i, j, k, mx, my, mz, i2, j2, k2, i3, j3, k3
  real*8 :: d, dist(1:3), tol = 1.d-10
  integer, allocatable :: head(:), this_list(:)

  if( any( L < R ) )then
    write(*,*) "ERROR: cutoff for neighbors build is longer than shortest lattice vector!"
    stop
  end if

  mx = int( L(1) / R )
  my = int( L(2) / R )
  mz = int( L(3) / R )

  allocate( head(1:mx*my*mz) )
  head = 0
  allocate( this_list(1:Np) )

  nn = 0
  do i = 1, Np
    call get_distance(L/2.d0, pos(1:3, i), L, PBC, dist, d)
    dist = dist + L/2.d0
    j = 1 + modulo(int( dist(1) / (L(1)+tol) * mx ), mx) &
          + modulo(int( dist(2) / (L(2)+tol) * my ), my) * mx &
          + modulo(int( dist(3) / (L(3)+tol) * mz ), mz) * my * mx
    this_list(i) = head(j)
    head(j) = i
  end do
  do i = 1, Np
    call get_distance(L/2.d0, pos(1:3, i), L, PBC, dist, d)
    dist = dist + L/2.d0
    i2 = 1 + int( dist(1) / (L(1)+tol) * mx )
    j2 = 1 + int( dist(2) / (L(2)+tol) * my )
    k2 = 1 + int( dist(3) / (L(3)+tol) * mz )
    do k3 = k2-1, k2+1
      if( mz == 1 .and. k3 /= 1 )cycle
      if( mz == 2 .and. k2 == 1 .and. k3 == 0 )cycle
      if( mz == 2 .and. k2 == 2 .and. k3 == 3 )cycle
      do j3 = j2-1, j2+1
        if( my == 1 .and. j3 /= 1 )cycle
        if( my == 2 .and. j2 == 1 .and. j3 == 0 )cycle
        if( my == 2 .and. j2 == 2 .and. j3 == 3 )cycle
        do i3 = i2-1, i2+1
          if( mx == 1 .and. i3 /= 1 )cycle
          if( mx == 2 .and. i2 == 1 .and. i3 == 0 )cycle
          if( mx == 2 .and. i2 == 2 .and. i3 == 3 )cycle
          j = 1 + modulo(i3-1,mx) + modulo(j3-1,my)*mx + modulo(k3-1,mz)*mx*my
          k = head(j)
          do while( k /= 0 )
            if( k /= i )then
              call get_distance(pos(1:3, i), pos(1:3, k), L, PBC, dist, d)
              if( d < R )then
                nn(i) = nn(i) + 1
                list(nn(i), i) = k
              end if
            end if
            k = this_list(k)
          end do
        end do
      end do
    end do
  end do
  deallocate(head, this_list)

end subroutine
!**********************************************************************************************

end module

