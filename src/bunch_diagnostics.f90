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

MODULE Diagnostics_on_Bunches

  ! USE class_species
  USE my_types
  USE use_my_types
  USE moments
  USE random_numbers_functions
  USE shapiro_wilks

  IMPLICIT NONE

CONTAINS

  ! --- BUNCH DIAGNOSTIC MANAGER ---!
  SUBROUTINE bunch_diagnostics_Integrated_AllBunches
  integer :: bunch_number, n,n1,n2,ier

  do bunch_number=1,bunchip%n_total_bunches
    call bunch_integrated_diagnostics(bunch_number)
    if(sim_parameters%L_SW_test) call bunch_shapiro_wilks(bunch_number)
  enddo

  if (sim_parameters%diagnostics_with_dcut.eq.1) then
    call apply_Sigma_cut(sim_parameters,plasma%k_p)
    do bunch_number=1,bunchip%n_total_bunches
      call bunch_integrated_diagnostics_for_dcutparticles(bunch_number)
    enddo
  endif

END SUBROUTINE bunch_diagnostics_Integrated_AllBunches
! --- --- --- --- --- --- --- --- !



!--- --- --- --- --- --- ---!
SUBROUTINE apply_Sigma_cut(sim_par,plasma_k_0)

USE my_types
USE class_species
USE class_particle

  TYPE(simul_param), INTENT(IN) :: sim_par
  INTEGER :: ss, npc, ni, nf, bunch_number, iparticle
  REAL(8) avgz,plasma_k_0,radius,drelZ,max_ellipsoidal
  REAL(8) :: s_x(1),s_y(1),s_z(1)
  REAL(8) :: mu_x(1),mu_y(1),mu_z(1)

  do bunch_number=1,sim_par%Nbunches

    mu_x(1)  = calculate_nth_moment(bunch_number,1,1,'nocentral')
  	mu_y(1)  = calculate_nth_moment(bunch_number,1,2,'nocentral')
  	mu_z(1)  = calculate_nth_moment(bunch_number,1,3,'nocentral')
    s_x(1)  = sqrt(calculate_nth_moment(bunch_number,2,1,'central'))
  	s_y(1)  = sqrt(calculate_nth_moment(bunch_number,2,2,'central'))
  	s_z(1)  = sqrt(calculate_nth_moment(bunch_number,2,3,'central'))

    do iparticle=1,size(bunch(bunch_number)%part(:))
 			if( bunch(bunch_number)%part(iparticle)%cmp(8) .eq. 1.) then
        if( abs(bunch(bunch_number)%part(iparticle)%cmp(1)-mu_x(1)) > 4.*bunchip%sx_um(bunch_number) ) bunch(bunch_number)%part(iparticle)%cmp(8)=0.
        if( abs(bunch(bunch_number)%part(iparticle)%cmp(2)-mu_y(1)) > 4.*bunchip%sy_um(bunch_number) ) bunch(bunch_number)%part(iparticle)%cmp(8)=0.
        if( abs(bunch(bunch_number)%part(iparticle)%cmp(3)-mu_z(1)) > 4.*bunchip%sz_um(bunch_number) ) bunch(bunch_number)%part(iparticle)%cmp(8)=0.
 			endif
 		enddo

 	enddo

END SUBROUTINE



!--- --- --- --- --- --- ---!
SUBROUTINE bunch_shapiro_wilks(bunch_number)
integer, intent(in) :: bunch_number
integer :: n,n1,n2,ier,sub_iter,component
real(8) :: w, pw, w_tmp,pw_tmp, SW(6), PSW(6)
real(8), allocatable :: x(:),a(:),selected(:),bunch_mask(:)
integer, allocatable :: order(:),randuniform(:)
logical :: init
character(1) :: num2str
character*90 :: filename


n=sim_parameters%SW_sample_dimension
n1 = n
n2 = n/2
ALLOCATE( x(n1), a(n2), order(n1), selected(n), randuniform(n))

do component=1,6
  w_tmp=0.d0
  pw_tmp=0.d0
  do sub_iter=1,sim_parameters%SW_sub_iter
    call random_INTeger_uniform_vector(randuniform,n,1,size(bunch(bunch_number)%part(:)))
    selected=bunch(bunch_number)%part(randuniform)%cmp(component)
    !--- Sort ascending order ---!
    CALL quick_sort(selected, order)
    !--- Shapiro-Wilks test ---!
    init = .FALSE.
    CALL swilk(init, selected, n, n1, n2, a, w, pw, ier)
    w_tmp=w_tmp+w
    pw_tmp=pw_tmp+pw
  enddo
  SW(component)=w_tmp/real(sim_parameters%SW_sub_iter)
  PSW(component)=pw_tmp/real(sim_parameters%SW_sub_iter)
enddo

write(num2str,'(I1)') bunch_number
filename=TRIM(sim_parameters%path_integrated_diagnostics)//'bunch_ShapiroWilks_'//num2str//'.dat'
call open_file(OSys%macwin,filename)
  !--- ShapiroWilks W-statistics for: X,Y,Z,Px,Py,Pz --- P-value for the same quantities ---!
  write(11,'(1p100e14.5)') sim_parameters%sim_time*c,SW(:),PSW(:)
close(11)


DEALLOCATE( x, a, order, selected, randuniform )
end subroutine bunch_shapiro_wilks


END MODULE Diagnostics_on_Bunches
