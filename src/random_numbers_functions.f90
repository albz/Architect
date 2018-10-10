!*****************************************************************************************************!
!             Copyright 2014-2016 Alberto Marocchino, Francesco Massimo                               !
!*****************************************************************************************************!

!*****************************************************************************************************!
!  This file is part of architect.                                                                    !
!                                                                                                     !
!  Architect is free software: you can redistribute it and/or modify                                  !
!  it under the terms of the GNU General Public License as published by                               !
!  the Free Software Foundation, either version 3 of the License, or                                  !
!  (at your option) any later version.                                                                !
!                                                                                                     !
!  Architect is distributed in the hope that it will be useful,                                       !
!  but WITHOUT ANY WARRANTY; without even the implied warranty of                                     !
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                      !
!  GNU General Public License for more details.                                                       !
!                                                                                                     !
!  You should have received a copy of the GNU General Public License                                  !
!  along with architect.  If not, see <http://www.gnu.org/licenses/>.                                 !
!*****************************************************************************************************!

 module random_numbers_functions
 implicit none

 !--- --- ---!
 contains
 !--- --- ---!

 !--- Box-Muller for a Norm(0,1)
 !--- random distributed variable
 subroutine boxmuller(randnormal)
 real(8), intent(inout) :: randnormal
 real(8) :: x,y,s,r

 s=10.
 do while( s >= 1. )
  call random_number(x)
  call random_number(y)
  x = 2.* x -1.
  y = 2.* y -1.
  s = x**2+y**2
 end do
 r=sqrt(-2.*log(s)/s)
 x = x*r
 randnormal=x
 end subroutine boxmuller


 subroutine boxmuller_vector(randnormal,len)
 integer, intent(in) :: len
 real(8), intent(inout) :: randnormal(len)
 real(8) :: x,y,s,r
 real(8) :: mu,std
 integer :: i

 do i=1,len
  s=10.
  do while( s >= 1. )
   call random_number(x)
   call random_number(y)
   x = 2.* x -1.
   y = 2.* y -1.
   s = x**2+y**2
  end do
  r=sqrt(-2.*log(s)/s)
  randnormal(i)=x*r
 enddo
 !--- convergence to N(0,1) ---!
 mu=sum(randnormal)/(1.*len)
 std=sqrt( sum((randnormal-mu)**2) / (1.*len-1.) )
 randnormal = (randnormal-mu)/std
 end subroutine boxmuller_vector

 subroutine boxmuller_vector_cut(randnormal,len,cut)
 integer, intent(in) :: len
 real(8), intent(inout) :: randnormal(len)
 real(8), intent(in) :: cut
 real(8) :: x,y,s,r
 integer :: i

 do i=1,len
112 continue !if the particle is cut :: recalculate particle position
  s=10.
  do while( s >= 1.)
   call random_number(x)
   call random_number(y)
   x = 2.* x -1.
   y = 2.* y -1.
   s = x**2+y**2
  end do
  r=sqrt(-2.*log(s)/s)
  randnormal(i)=x*r
  if(abs(randnormal(i))>cut) goto 112
 enddo
 randnormal = randnormal-sum(randnormal)/(1.*max(1,size(randnormal)))
 end subroutine boxmuller_vector_cut


 !---------------------------------------------------!
 !---*** BoxMuller with cut for multidimension ***---!
 !---------------------------------------------------!
 subroutine boxmuller_vector_cutDimensional(randnormal,len,s_cut,dimension)
 integer, intent(in) :: len,dimension
 real(8), intent(inout) :: randnormal(3,len)
 real(8), intent(in) :: s_cut
 real(8) :: x,y,z,r
 integer :: i

 select case(dimension)
 case(3) !case 3D
    do i=1,len
      call boxmuller(x)
      call boxmuller(y)
      call boxmuller(z)
      Do While( x**2 + y**2 + z**2 > s_cut**2 )
        call boxmuller(x)
        call boxmuller(y)
        call boxmuller(z)
      enddo
      randnormal(1,i)=x
      randnormal(2,i)=y
      randnormal(3,i)=z
    enddo
  case(2) !case 2D
    do i=1,len
      call boxmuller(x)
      call boxmuller(y)
      Do While( x**2 + y**2 > s_cut**2 )
        call boxmuller(x)
        call boxmuller(y)
      enddo
      randnormal(1,i)=x
      randnormal(2,i)=y
    enddo
  case(1) !case 1D
    do i=1,len
      115 continue
      call boxmuller(x)
      if( (x/s_cut)**2 > 1) goto 115
      randnormal(1,i)=x
    enddo
  end select
 end subroutine boxmuller_vector_cutDimensional



 subroutine random_uniform_vector(randuniform,len)
 integer, intent(in) :: len
 real(8), intent(inout) :: randuniform(len)
 real(8) :: x
 integer :: i

 do i=1,len
    call random_number(x)
    randuniform(i)=x
 enddo

 end subroutine random_uniform_vector


 subroutine random_INTeger_uniform_vector(randuniform,len,min,max)
 integer, intent(in) :: len,min,max
 integer, intent(inout) :: randuniform(len)
 real(8) :: x
 integer :: i

 do i=1,len
    call random_number(x)
    randuniform(i)=min + Nint(x*(max-min))
 enddo

 end subroutine random_INTeger_uniform_vector


 !--- function: uniform distribution between 'min' and 'max'
  real(8) function random_number_range(minimum,maximum)
  real(8), intent(in) :: minimum,maximum
  real(8) :: x
    call random_number(x)
    random_number_range = (maximum-minimum)*x+minimum
  end function random_number_range



 !--- --- ---!
 end module random_numbers_functions
 !--- --- ---!
