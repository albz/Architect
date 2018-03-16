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

 USE random_numbers_functions
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
	generated_bunch(1,:)=rnumber*s_x !+ x_cm
	call boxmuller_vector(rnumber,nparticles)
	generated_bunch(2,:)=rnumber*s_y !+ y_cm
	call boxmuller_vector(rnumber,nparticles)
	generated_bunch(3,:)=rnumber*s_z !+ z_cm
	call boxmuller_vector(rnumber,nparticles)
	generated_bunch(4,:)=rnumber*eps_x/s_x
	call boxmuller_vector(rnumber,nparticles)
	generated_bunch(5,:)=rnumber*eps_y/s_y
	call boxmuller_vector(rnumber,nparticles)
	generated_bunch(6,:)=rnumber *0.01*dgamma*gamma_m + gamma_m
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
	generated_bunch(3,:)=rnumber*s_z !+ z_cm
	call boxmuller_vector(rnumber,nparticles)
	generated_bunch(4,:)=rnumber*eps_x/s_x
	call boxmuller_vector(rnumber,nparticles)
	generated_bunch(5,:)=rnumber*eps_y/s_y
	call boxmuller_vector(rnumber,nparticles)
	generated_bunch(6,:)=rnumber * 0.01*dgamma*gamma_m + gamma_m
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

!--- *** ---!
 subroutine generate_hollow_bunch(x_cm,y_cm,z_cm,s_x,s_y,s_z,gamma_m,eps_x,eps_y,dgamma,nparticles,generated_bunch)

 integer,intent(in)   :: nparticles
 real(8),intent(in)      :: x_cm,y_cm,z_cm,s_x,s_y,s_z,gamma_m,eps_x,eps_y,dgamma
 real(8),intent(inout)   :: generated_bunch(:,:)
 real(8) :: rnumber(nparticles),particle_radius(nparticles),particle_theta(nparticles)

  call random_uniform_vector(particle_radius,nparticles)
  particle_radius=s_x+particle_radius*(s_y-s_x)
  call random_uniform_vector(particle_theta,nparticles)
  particle_theta=2.*3.1415*particle_theta
	! call boxmuller_vector(rnumber,nparticles)
	generated_bunch(1,:)=particle_radius*cos(particle_theta)
	! call boxmuller_vector(rnumber,nparticles)
	generated_bunch(2,:)=particle_radius*sin(particle_theta)
	call random_uniform_vector(rnumber,nparticles)
	generated_bunch(3,:)=(rnumber*s_z)-s_z/2. !+ z_cm
	call boxmuller_vector(rnumber,nparticles)
	generated_bunch(4,:)=rnumber*eps_x/s_x
	call boxmuller_vector(rnumber,nparticles)
	generated_bunch(5,:)=rnumber*eps_y/s_y
	call boxmuller_vector(rnumber,nparticles)
	generated_bunch(6,:)=rnumber * 0.01*dgamma*gamma_m + gamma_m


 end subroutine generate_hollow_bunch


!--- *** ---!
 subroutine generate_triangularZ_uniformR_bunch(x_cm,y_cm,z_cm,s_x,s_y,s_z,gamma_m,eps_x,eps_y,dgamma,nparticles,generated_bunch,Charge_right,Charge_left)
 integer,intent(in)   :: nparticles
 real(8),intent(in)      :: x_cm,y_cm,z_cm,s_x,s_y,s_z,gamma_m,eps_x,eps_y,dgamma,Charge_right,Charge_left
 real(8),intent(inout)   :: generated_bunch(:,:)
 real(8) :: rnumber(nparticles),particle_radius(nparticles),particle_theta(nparticles)
 integer :: i
 real(8) :: z,y,x,a

  do i=1,nparticles
    call random_number(z)
    call random_number(a)

    if (Charge_right>=Charge_left) then
      Do while(a>Charge_left+(Charge_right-Charge_left)*z)
        call random_number(z)
        call random_number(a)
      enddo
    else
      Do while(a>Charge_right+(Charge_left-Charge_right)*z)
        call random_number(z)
        call random_number(a)
      enddo
    endif



    x=random_number_range(-1.d0,1.d0)
    y=random_number_range(-1.d0,1.d0)
    Do while(sqrt(x**2+y**2)>1.d0)
      x=random_number_range(-1.d0,1.d0)
      y=random_number_range(-1.d0,1.d0)
    enddo
    generated_bunch(1,i)=x*s_x!+x_cm
    generated_bunch(2,i)=y*s_y!+y_cm
    generated_bunch(3,i)=z*s_z!+z_cm
  enddo

	call boxmuller_vector(rnumber,nparticles)
	generated_bunch(4,:)=rnumber*eps_x/s_x
	call boxmuller_vector(rnumber,nparticles)
	generated_bunch(5,:)=rnumber*eps_y/s_y
	call boxmuller_vector(rnumber,nparticles)
	generated_bunch(6,:)=rnumber * 0.01*dgamma*gamma_m + gamma_m

 end subroutine generate_triangularZ_uniformR_bunch

!--- *** Cylindrical Bunch ***----!
subroutine generate_cylindrical_bunch(x_cm,y_cm,z_cm,s_x,s_y,s_z,gamma_m,eps_x,eps_y,dgamma,nparticles,generated_bunch)
integer,intent(in)   :: nparticles
real(8),intent(in)      :: x_cm,y_cm,z_cm,s_x,s_y,s_z,gamma_m,eps_x,eps_y,dgamma
real(8),intent(inout)   :: generated_bunch(:,:)
real(8) :: rnumber(nparticles),particle_radius(nparticles),particle_theta(nparticles)
integer :: i
real(8) :: z,y,x,a

  do i=1,nparticles
     x=random_number_range(-1.d0,1.d0)
     y=random_number_range(-1.d0,1.d0)
     z=random_number_range(-1.d0,1.d0)
     Do while(sqrt(x**2+y**2)>1.d0)
       x=random_number_range(-1.d0,1.d0)
       y=random_number_range(-1.d0,1.d0)
     enddo
     generated_bunch(1,i)=x*s_x!+x_cm
     generated_bunch(2,i)=y*s_y!+y_cm
     generated_bunch(3,i)=z*s_z!+z_cm
   enddo

 call boxmuller_vector(rnumber,nparticles)
 generated_bunch(4,:)=rnumber*eps_x/s_x
 call boxmuller_vector(rnumber,nparticles)
 generated_bunch(5,:)=rnumber*eps_y/s_y
 call boxmuller_vector(rnumber,nparticles)
 generated_bunch(6,:)=rnumber * 0.01*dgamma*gamma_m + gamma_m
end subroutine generate_cylindrical_bunch


 subroutine generate_triangularZ_normalR_bunch(x_cm,y_cm,z_cm,s_x,s_y,s_z,gamma_m,eps_x,eps_y,dgamma,nparticles,generated_bunch,Charge_right,Charge_left)
 integer,intent(in)   :: nparticles
 real(8),intent(in)      :: x_cm,y_cm,z_cm,s_x,s_y,s_z,gamma_m,eps_x,eps_y,dgamma,Charge_right,Charge_left
 real(8),intent(inout)   :: generated_bunch(:,:)
 real(8) :: rnumber(nparticles),particle_radius(nparticles),particle_theta(nparticles)
 integer :: i
 real(8) :: z,y,x,a

  do i=1,nparticles
    call random_number(z)
    call random_number(a)
    if(Charge_left >= Charge_right) then !triangular shape left-right
      Do while(a*Charge_left>Charge_left+(Charge_right-Charge_left)*z)
        call random_number(z)
        call random_number(a)
      enddo
    elseif(Charge_left < Charge_right) then !triangular shape right-left
      Do while(a*Charge_right>Charge_left+(Charge_right-Charge_left)*z)
        call random_number(z)
        call random_number(a)
      enddo
    endif
    generated_bunch(3,i)=z*s_z!+z_cm
  enddo

  call boxmuller_vector(rnumber,nparticles)
	generated_bunch(1,:)=rnumber*s_x !+ x_cm
	call boxmuller_vector(rnumber,nparticles)
	generated_bunch(2,:)=rnumber*s_y !+ y_cm

	call boxmuller_vector(rnumber,nparticles)
	generated_bunch(4,:)=rnumber*eps_x/s_x
	call boxmuller_vector(rnumber,nparticles)
	generated_bunch(5,:)=rnumber*eps_y/s_y
	call boxmuller_vector(rnumber,nparticles)
	generated_bunch(6,:)=rnumber * 0.01*dgamma*gamma_m + gamma_m

end subroutine generate_triangularZ_normalR_bunch


!--- *** ---!
 subroutine read_from_external_file(name_file,number_of_particles,generated_bunch)
 character,intent(in) :: name_file*255
 integer,intent(in)   :: number_of_particles
 real(8),intent(out)   :: generated_bunch(:,:)
 real(8) :: dummy_charge
 integer :: i,dummy_nparticles
 character(199) :: string

 open(9,file=trim(name_file),status='old')

 write(*,'(A,A)') 'loading particle from external file: ',trim(name_file)

!--- read: header ---!
read(9,*) string
read(9,*) dummy_charge
read(9,*) string
read(9,*) dummy_nparticles
read(9,*) string

 DO i=1,number_of_particles
 	read(9,*) generated_bunch(1,i), generated_bunch(2,i), generated_bunch(3,i), generated_bunch(4,i), generated_bunch(5,i), generated_bunch(6,i)
 ENDDO

 close(9)

 write(*,'(A,A)') trim(name_file),'   loaded'

 end subroutine read_from_external_file


 !--- *** ---!
  subroutine read_header_external_file(name_file,charge,number_of_particles)
    character,intent(in) :: name_file*255
    real(8),intent(out)   :: charge
    integer,intent(out)   :: number_of_particles
    character(199) :: string

    open(9,file=trim(name_file),status='unknown')
    !--- read: header ---!
    read(9,*) string
    read(9,*) charge
    read(9,*) string
    read(9,*) number_of_particles
    write(*,'(A,f6.3,A)') 'reading header :: Bunch Charge::  ',charge,' (nC)'
    write(*,'(A,1I8)') 'Reading header :: Bunch number of particles:: ',number_of_particles
    close(9)
  end subroutine read_header_external_file




 !--- --- ---!
 end module bunch_generation
 !--- --- ---!
