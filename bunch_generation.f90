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

 module bunch_generation
 implicit none

 !--- --- ---!
 contains
 !--- --- ---!

 subroutine generate_bunch(x_cm,y_cm,z_cm,s_x,s_y,s_z,gamma_m,eps_x,eps_y,dgamma,nparticles,generated_bunch)

 integer,intent(in)   :: nparticles
 real(8),intent(in)      :: x_cm,y_cm,z_cm,s_x,s_y,s_z,gamma_m,eps_x,eps_y,dgamma
 real(8),intent(inout)   :: generated_bunch(:,:)
 real(8) :: rnumber(nparticles)

	call boxmuller_vector(rnumber,nparticles)
	generated_bunch(1,:)=rnumber*s_x + x_cm
	call boxmuller_vector(rnumber,nparticles)
	generated_bunch(2,:)=rnumber*s_y + y_cm
	call boxmuller_vector(rnumber,nparticles)
	generated_bunch(3,:)=rnumber*s_z + z_cm
	call boxmuller_vector(rnumber,nparticles)
	generated_bunch(4,:)=rnumber*eps_x/s_x
	call boxmuller_vector(rnumber,nparticles)
	generated_bunch(5,:)=rnumber*eps_y/s_y
	call boxmuller_vector(rnumber,nparticles)
	generated_bunch(6,:)=rnumber * sqrt(3.0)*0.01*dgamma*gamma_m + gamma_m
!	bunch(7,n1:n2)=weight

 end subroutine generate_bunch

!--- *** ---!
 subroutine generate_bunch_twiss(x_cm,y_cm,z_cm,s_x,s_y,s_z,gamma_m,eps_x,eps_y,dgamma,nparticles,generated_bunch,alphaTwiss,betaTwiss)

 integer,intent(in)   :: nparticles
 real(8),intent(in)      :: x_cm,y_cm,z_cm,s_x,s_y,s_z,gamma_m,eps_x,eps_y,dgamma
 real(8),intent(in)      :: alphaTwiss,betaTwiss
 real(8),intent(inout)   :: generated_bunch(:,:)
 real(8) :: rnumber(nparticles)
 real(8)  :: ax11,ax12,ay11,ay12
 integer :: i

 !twiss-matrix
 ax11=sqrt( eps_x*betaTwiss/(s_x**2+s_x**2*alphaTwiss**2) )
 ay11=sqrt( eps_y*betaTwiss/(s_y**2+s_y**2*alphaTwiss**2) )
 ax12=-ax11*alphaTwiss*s_x**2/eps_x
 ay12=-ay11*alphaTwiss*s_y**2/eps_y

	call boxmuller_vector(rnumber,nparticles)
	generated_bunch(1,:)=rnumber*s_x
	call boxmuller_vector(rnumber,nparticles)
	generated_bunch(2,:)=rnumber*s_y
	call boxmuller_vector(rnumber,nparticles)
	generated_bunch(3,:)=rnumber*s_z + z_cm
	call boxmuller_vector(rnumber,nparticles)
	generated_bunch(4,:)=rnumber*eps_x/s_x
	call boxmuller_vector(rnumber,nparticles)
	generated_bunch(5,:)=rnumber*eps_y/s_y
	call boxmuller_vector(rnumber,nparticles)
	generated_bunch(6,:)=rnumber * sqrt(3.0)*0.01*dgamma*gamma_m + gamma_m
!	bunch(7,n1:n2)=weight

  !twiss-shifting
  DO i=1,nparticles
   generated_bunch(1,i)=ax11*generated_bunch(1,i)+ax12*generated_bunch(4,i) + x_cm
   generated_bunch(2,i)=ay11*generated_bunch(2,i)+ay12*generated_bunch(5,i) + y_cm
   generated_bunch(3,i)=generated_bunch(3,i)
   generated_bunch(4,i)=generated_bunch(4,i)/ax11
   generated_bunch(5,i)=generated_bunch(5,i)/ay11
   generated_bunch(6,i)=generated_bunch(6,i)
  ENDDO
 end subroutine generate_bunch_twiss


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
 randnormal = randnormal-sum(randnormal)/(1.*max(1,size(randnormal)))

 end subroutine boxmuller_vector


 subroutine read_from_external_file(name_file,nparticles,generated_bunch)
 character,intent(in) :: name_file*255
 integer,intent(in)   :: nparticles
 real(8),intent(inout)   :: generated_bunch(:,:)
 integer		      :: i

 open(9,file=trim(name_file),status='old')

 write(*,*) 'loading: ',trim(name_file)
 DO i=1,nparticles
 	read(9,*) generated_bunch(1,i), generated_bunch(2,i), generated_bunch(3,i), generated_bunch(4,i), generated_bunch(5,i), generated_bunch(6,i)
 ENDDO

 close(9)

 write(*,*) trim(name_file),'   loaded'
 write(*,*) ''

 end subroutine read_from_external_file





 !--- --- ---!
 end module bunch_generation
 !--- --- ---!
