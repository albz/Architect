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

  ! USE pstruct_data
  USE my_types
  USE use_my_types
  USE moments
  USE random_numbers_functions
  USE shapiro_wilks

  IMPLICIT NONE

CONTAINS

  ! --- BUNCH DIAGNOSTIC MANAGER ---!
  SUBROUTINE bunch_diagnostics_Integrated_AllBunches
  integer :: i, n,n1,n2,ier
  real(8) :: mu_z,sigma_z

  do i=1,bunch_initialization%n_total_bunches
    call bunch_sliced_diagnostics(i)
    call bunch_integrated_diagnostics(i)
    if(sim_parameters%L_SW_test) call bunch_shapiro_wilks(i)
  enddo

  if (sim_parameters%diagnostics_with_dcut.eq.1) then
    call apply_Sigma_cut(sim_parameters,plasma%k_p)
    do i=1,bunch_initialization%n_total_bunches
      mu_z    = calculate_nth_moment_bunch_dcut(i,1,3)
      sigma_z = sqrt( calculate_nth_central_moment_bunch_dcut(i,2,3) )

      call bunch_sliced_diagnostics_dcut(i,mu_z,sigma_z)
      call bunch_integrated_diagnostics_dcut(i)
    enddo
  endif

END SUBROUTINE bunch_diagnostics_Integrated_AllBunches
! --- --- --- --- --- --- --- --- !






!--- --- --- --- --- --- ---!
SUBROUTINE apply_Sigma_cut(sim_par,plasma_k_0)

USE my_types
USE pstruct_data
USE architect_class_structure

  TYPE(simul_param), INTENT(IN) :: sim_par
  INTEGER :: ss, npc, ni, nf, ibunch, iparticle
  REAL(8) avgz,plasma_k_0,radius,drelZ,max_ellipsoidal
  REAL(8) :: s_x(1),s_y(1),s_z(1)
  REAL(8) :: mu_x(1),mu_y(1),mu_z(1)

  ! do ibunch=1,sim_par%Nbunches
  !   mu_x(1)  = calculate_nth_moment(ibunch,1,1,'nocentral')
  ! 	mu_y(1)  = calculate_nth_moment(ibunch,1,2,'nocentral')
  ! 	mu_z(1)  = calculate_nth_moment(ibunch,1,3,'nocentral')
  !   s_x(1)  = sqrt(calculate_nth_moment(ibunch,2,1,'central'))
  ! 	s_y(1)  = sqrt(calculate_nth_moment(ibunch,2,2,'central'))
  ! 	s_z(1)  = sqrt(calculate_nth_moment(ibunch,2,3,'central'))
  !
  !   do iparticle=1,size(bunch(ibunch)%part(:))
 	! 		if( bunch(ibunch)%part(iparticle)%cmp(8) .eq. 1.) then
  !       if( abs(bunch(ibunch)%part(iparticle)%cmp(1)-mu_x(1)) > 4.*s_x(1)  ) bunch(ibunch)%part(iparticle)%cmp(8)=0.
  !       if( abs(bunch(ibunch)%part(iparticle)%cmp(2)-mu_y(1)) > 4.*s_y(1)  ) bunch(ibunch)%part(iparticle)%cmp(8)=0.
  !       if( abs(bunch(ibunch)%part(iparticle)%cmp(3)-mu_z(1)) > 4.*s_z(1)  ) bunch(ibunch)%part(iparticle)%cmp(8)=0.
 	! 		endif
 	! 	enddo
  !
 	! enddo


  do ibunch=1,sim_par%Nbunches

    mu_x(1)  = calculate_nth_moment(ibunch,1,1,'nocentral')
  	mu_y(1)  = calculate_nth_moment(ibunch,1,2,'nocentral')
  	mu_z(1)  = calculate_nth_moment(ibunch,1,3,'nocentral')
    s_x(1)  = sqrt(calculate_nth_moment(ibunch,2,1,'central'))
  	s_y(1)  = sqrt(calculate_nth_moment(ibunch,2,2,'central'))
  	s_z(1)  = sqrt(calculate_nth_moment(ibunch,2,3,'central'))

    do iparticle=1,size(bunch(ibunch)%part(:))
 			if( bunch(ibunch)%part(iparticle)%cmp(8) .eq. 1.) then
        if( abs(bunch(ibunch)%part(iparticle)%cmp(1)-mu_x(1)) > 4.*bunch_initialization%bunch_s_x(ibunch) ) bunch(ibunch)%part(iparticle)%cmp(8)=0.
        if( abs(bunch(ibunch)%part(iparticle)%cmp(2)-mu_y(1)) > 4.*bunch_initialization%bunch_s_y(ibunch) ) bunch(ibunch)%part(iparticle)%cmp(8)=0.
        if( abs(bunch(ibunch)%part(iparticle)%cmp(3)-mu_z(1)) > 4.*bunch_initialization%bunch_s_z(ibunch) ) bunch(ibunch)%part(iparticle)%cmp(8)=0.
 			endif
 		enddo

 	enddo

END SUBROUTINE




!--- --- --- --- --- --- ---!
SUBROUTINE write_integrated_diagnostics(sim_par,plasma_k_0)


USE my_types
USE pstruct_data
USE architect_class_structure


   TYPE(simul_param), INTENT(IN) :: sim_par
   INTEGER :: ss, npc, npl, nplost, ni, nf, ii
   REAL(8) :: avgz,plasma_k_0,radius,drelZ
   REAL(8) :: mu_x,mu_y,mu_z,mu_px,mu_py,mu_pz,sigma_x,sigma_y,sigma_z
   REAL(8) :: sigma_px,sigma_py,sigma_pz,Corr_x_px,Corr_y_py
   REAL(8) :: eps_x,eps_y,gamma_mean,sigma_g,dgsug
   REAL(8), DIMENSION(:), ALLOCATABLE :: sigma__x,sigma__y,sigma__z,sigma__pz,sigma__py,sigma__px,gamma_sp
   CHARACTER(1) :: num2str

	! --- --- --- !
	! bunch(es) diagnostic
	! --- --- --- !


	do ii=1,sim_par%Nbunches

!---
		npc=sum(bunch(ii)%part(:)%cmp(8))
		npl=size(bunch(ii)%part(:))
		nplost=npl-npc
		allocate( sigma__x(1:npl),sigma__y(1:npl),sigma__z(1:npl),sigma__pz(1:npl),sigma__py(1:npl),sigma__px(1:npl),gamma_sp(1:npl))
!---
		mu_x  = sum(bunch(ii)%part(:)%cmp(1)*bunch(ii)%part(:)%cmp(8))/npc
		mu_y  = sum(bunch(ii)%part(:)%cmp(2)*bunch(ii)%part(:)%cmp(8))/npc
		mu_z  = sum(bunch(ii)%part(:)%cmp(3)*bunch(ii)%part(:)%cmp(8))/npc
		mu_px = sum(bunch(ii)%part(:)%cmp(4)*bunch(ii)%part(:)%cmp(8))/npc
		mu_py = sum(bunch(ii)%part(:)%cmp(5)*bunch(ii)%part(:)%cmp(8))/npc
		mu_pz = sum(bunch(ii)%part(:)%cmp(6)*bunch(ii)%part(:)%cmp(8))/npc

		sigma__x  = (bunch(ii)%part(:)%cmp(1)-mu_x)*bunch(ii)%part(:)%cmp(8)
		sigma_x   = sqrt( dot_product(dble(sigma__x),dble(sigma__x))/dble(npc) )
		sigma__y  = (bunch(ii)%part(:)%cmp(2)-mu_y)*bunch(ii)%part(:)%cmp(8)
		sigma_y   = sqrt( dot_product(dble(sigma__y),dble(sigma__y))/dble(npc) )
		sigma__z  = (bunch(ii)%part(:)%cmp(3)-mu_z)*bunch(ii)%part(:)%cmp(8)
		sigma_z   = sqrt( dot_product(dble(sigma__z),dble(sigma__z))/dble(npc) )
		sigma__px = (bunch(ii)%part(:)%cmp(4)-mu_px)*bunch(ii)%part(:)%cmp(8)
		sigma_px  = sqrt( dot_product(dble(sigma__px),dble(sigma__px))/dble(npc) )
		sigma__py = (bunch(ii)%part(:)%cmp(5)-mu_py)*bunch(ii)%part(:)%cmp(8)
		sigma_py  = sqrt( dot_product(dble(sigma__py),dble(sigma__py))/dble(npc) )
		sigma__pz = (bunch(ii)%part(:)%cmp(6)-mu_pz)*bunch(ii)%part(:)%cmp(8)
		sigma_pz  = sqrt( dot_product(dble(sigma__pz),dble(sigma__pz))/dble(npc) )

		Corr_x_px = dot_product(dble(sigma__x),dble(sigma__px))/dble(npc)
		Corr_y_py = dot_product(dble(sigma__y),dble(sigma__py))/dble(npc)

		eps_x = sqrt(sigma_x**2 * sigma_px**2 -Corr_x_px**2)
		eps_y = sqrt(sigma_y**2 * sigma_py**2 -Corr_y_py**2)

		gamma_sp=real(bunch(ii)%part(:)%cmp(8))* &
		sqrt(1.+(bunch(ii)%part(:)%cmp(4))**2+(bunch(ii)%part(:)%cmp(5))**2+(bunch(ii)%part(:)%cmp(6))**2)

		gamma_mean=sum(gamma_sp)/real(npc)
		gamma_sp=(gamma_sp-gamma_mean)*real(bunch(ii)%part(:)%cmp(8))
		sigma_g=sqrt(dot_product(gamma_sp,gamma_sp)/real(npc))
		dgsug=sigma_g/gamma_mean

		write(num2str,'(I1)') ii
		open(15,file='bunch_integrated_diagnostic_'//num2str//'.dat',status='unknown',position='append')
			write(15,'(1p100e14.5)') sim_par%sim_time,-mu_z,mu_x,mu_y,mu_z,mu_px,mu_py,mu_pz,sigma_x,sigma_y,sigma_z,Corr_x_px,Corr_y_py,eps_x,eps_y,sigma_g,gamma_mean,dgsug,real(npc),real(nplost)
		close(15)

		deallocate( sigma__x,sigma__y,sigma__z,sigma__pz,sigma__py,sigma__px,gamma_sp)

	ENDDO

END SUBROUTINE


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
