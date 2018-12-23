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

 module moments

 USE class_species
 USE my_types
 USE use_my_types
 USE utilities
 USE random_numbers_functions
 USE shapiro_wilks
 USE Compute_plasma_current

 implicit none

 !--- --- ---!
 contains
 !--- --- ---!

!--- weighted moments ---!
 real(8) FUNCTION calculate_nth_moment(number_bunch,nth,component,central)
 integer, intent(in) :: nth, component, number_bunch
 character(len=*), intent(in) :: central
 integer :: np,i
 real(8) :: mu_mean(1),moment(1),weight(1)
 logical, allocatable, dimension(:) :: maskbunch

 !--- create logical mask---!
 np = size(bunch(number_bunch)%part(:))
 allocate(maskbunch(np))
 maskbunch=.true.
 DO i=1,np
   if(bunch(number_bunch)%part(i)%cmp(14)<1.) maskbunch(i)=.false.
 enddo
 !--- ---!

 !--- mean calculation
 if(component<7) then
    if(trim(bunchip%PWeights(number_bunch))=='weighted' .and. (component==1 .or. component==4) ) then
            mu_mean=0.D0
            weight   = sum( bunch(number_bunch)%part(:)%cmp(13),                                            mask=maskbunch )
    else
            mu_mean  = sum( bunch(number_bunch)%part(:)%cmp(component)*bunch(number_bunch)%part(:)%cmp(13), mask=maskbunch )
            weight   = sum( bunch(number_bunch)%part(:)%cmp(13),                                            mask=maskbunch )
            mu_mean  = mu_mean / weight
    endif
 else if(component==7) then !---gamma
   mu_mean  = sum( &
   sqrt(1. + bunch(number_bunch)%part(:)%cmp(4)**2  &
           + bunch(number_bunch)%part(:)%cmp(5)**2  &
           + bunch(number_bunch)%part(:)%cmp(6)**2) &
           * bunch(number_bunch)%part(:)%cmp(13), mask=maskbunch )
   weight   = sum( bunch(number_bunch)%part(:)%cmp(13),           mask=maskbunch )
   mu_mean  = mu_mean / weight
 endif

 !--- moment calculation
 if     ( trim(central)=='central' .and. component<7) then
   moment = sum( (bunch(number_bunch)%part(:)%cmp(component)-mu_mean(1))**nth *bunch(number_bunch)%part(:)%cmp(13), mask=maskbunch )
   moment = moment / weight
 elseif ( trim(central)=='central' .and. component==7) then
   moment  = sum( &
   (sqrt(1. + bunch(number_bunch)%part(:)%cmp(4)**2  &
          + bunch(number_bunch)%part(:)%cmp(5)**2   &
          + bunch(number_bunch)%part(:)%cmp(6)**2)  &
          - mu_mean(1))**nth                        &
          * bunch(number_bunch)%part(:)%cmp(13), mask=maskbunch )
   moment  = moment / weight
 elseif ( trim(central)=='nocentral' .and. component<7) then
        if(trim(bunchip%PWeights(number_bunch))=='weighted' .and. (component==1 .or. component==4) ) then
            moment=0.D0
        else
            moment = sum( bunch(number_bunch)%part(:)%cmp(component)**nth *bunch(number_bunch)%part(:)%cmp(13), mask=maskbunch )
            moment = moment / weight
        endif
 elseif ( trim(central)=='nocentral' .and. component==7) then
   moment  = sum( &
   (sqrt(1. + bunch(number_bunch)%part(:)%cmp(4)**2   &
            + bunch(number_bunch)%part(:)%cmp(5)**2   &
            + bunch(number_bunch)%part(:)%cmp(6)**2)  &
            )**nth                                    &
            * bunch(number_bunch)%part(:)%cmp(13), mask=maskbunch )
   moment  = moment / weight
 endif

 !---
 deallocate(maskbunch)
 if(trim(bunchip%PWeights(number_bunch))=='weighted' .and. (component==1) ) moment=moment/2.D0
 calculate_nth_moment = moment(1)
 END FUNCTION calculate_nth_moment


 !--- weighted correlation ---!
 real(8) FUNCTION calculate_correlation(number_bunch,component1,component2)
 integer, intent(in) :: number_bunch, component1, component2
 real(8) :: mu_component1(1),mu_component2(1)
 real(8) :: corr,weight
 integer :: np,i
 logical, allocatable, dimension(:) :: maskbunch

 !--- create logical mask---!
 np = size(bunch(number_bunch)%part(:))
 allocate(maskbunch(np))
 maskbunch=.true.
 DO i=1,np
   if(bunch(number_bunch)%part(i)%cmp(14)<1.) maskbunch(i)=.false.
 enddo
 !--- ---!

 mu_component1(1) = calculate_nth_moment(number_bunch,1,component1,'nocentral')
 mu_component2(1) = calculate_nth_moment(number_bunch,1,component2,'nocentral')
 weight           = sum( bunch(number_bunch)%part(:)%cmp(13), mask=maskbunch )

 corr = sum(                                                     &
 (bunch(number_bunch)%part(:)%cmp(component1)-mu_component1(1))* &
 (bunch(number_bunch)%part(:)%cmp(component2)-mu_component2(1))* &
  bunch(number_bunch)%part(:)%cmp(13), mask=maskbunch)
 corr = corr / weight

 !---
 deallocate(maskbunch)
 calculate_correlation = corr
 END FUNCTION calculate_correlation


 !--- --- ---!
 SUBROUTINE bunch_integrated_diagnostics(number_bunch)
 integer, intent(in) :: number_bunch
 real(8) :: mu_x(1),mu_y(1),mu_z(1)          !spatial meam
 real(8) :: mu_px(1),mu_py(1),mu_pz(1)       !momenta mean
 real(8) :: s_x(1),s_y(1),s_z(1)                         !spatial variance
 real(8) :: s_px(1),s_py(1), s_pz(1)                      !momenta variance
 real(8) :: m4_x(1),m4_y(1),m4_z(1),m4_px(1),m4_py(1),m4_pz(1) !4th central moments for kurtosis
 real(8) :: mu_gamma(1),s_gamma(1),dgamma_su_gamma(1)     !gamma mean-variance
 real(8) :: corr_y_py(1),corr_z_pz(1),corr_x_px(1) !correlation transverse plane
 real(8) :: emittance_x(1),emittance_y(1) !emittance variables
 real(8) :: u1(1),u2(1),u3(1),u4(1),u5(1)
 character(1) :: b2str
 character*90 :: filename
 TYPE(simul_param) :: sim_par

  !---mask :: particle selection :: all particles---!
  bunch(number_bunch)%part(:)%cmp(14)=1.
  !--- *** ---!

	mu_x(1)  = calculate_nth_moment(number_bunch,1,1,'nocentral')
	mu_y(1)  = calculate_nth_moment(number_bunch,1,2,'nocentral')
	mu_z(1)  = calculate_nth_moment(number_bunch,1,3,'nocentral')
	mu_px(1) = calculate_nth_moment(number_bunch,1,4,'nocentral')
	mu_py(1) = calculate_nth_moment(number_bunch,1,5,'nocentral')
	mu_pz(1) = calculate_nth_moment(number_bunch,1,6,'nocentral')

	s_x(1)  = sqrt(calculate_nth_moment(number_bunch,2,1,'central'))
	s_y(1)  = sqrt(calculate_nth_moment(number_bunch,2,2,'central'))
	s_z(1)  = sqrt(calculate_nth_moment(number_bunch,2,3,'central'))
	s_px(1) = sqrt(calculate_nth_moment(number_bunch,2,4,'central'))
	s_py(1) = sqrt(calculate_nth_moment(number_bunch,2,5,'central'))
	s_pz(1) = sqrt(calculate_nth_moment(number_bunch,2,6,'central'))

  m4_x(1)  = calculate_nth_moment(number_bunch,4,1,'central')
	m4_y(1)  = calculate_nth_moment(number_bunch,4,1,'central')
	m4_z(1)  = calculate_nth_moment(number_bunch,4,1,'central')
	m4_px(1) = calculate_nth_moment(number_bunch,4,1,'central')
	m4_py(1) = calculate_nth_moment(number_bunch,4,1,'central')
	m4_pz(1) = calculate_nth_moment(number_bunch,4,1,'central')

	corr_x_px(1) = calculate_correlation(number_bunch,1,4)
	corr_y_py(1) = calculate_correlation(number_bunch,2,5)
	corr_z_pz(1) = calculate_correlation(number_bunch,3,6)

  !---!
  mu_gamma(1) = calculate_nth_moment(number_bunch,1,7,'nocentral')
	s_gamma(1)  = sqrt(calculate_nth_moment(number_bunch,2,7,'central'))
	dgamma_su_gamma = s_gamma(1)/mu_gamma(1)

	!---!
	emittance_x = sqrt( s_x(1)**2 *s_px(1)**2 - corr_x_px(1)**2 )
	emittance_y = sqrt( s_y(1)**2 *s_py(1)**2 - corr_y_py(1)**2 )

  write(b2str,'(I1.1)') number_bunch
  filename=TRIM(sim_parameters%path_integrated_diagnostics)//'bunch_integrated_quantity_'//b2str//'.dat'
  call open_file(OSys%macwin,filename)
  !1  2   3   4   5    6    7     8      9      10     11      12      13     14    15    16     17       18       19       20        21
  !t,<X>,<Y>,<Z>,<Px>,<Py>,<Pz>,<rmsX>,<rmsY>,<rmsZ>,<rmsPx>,<rmsPy>,<rmsPz>,<Emx>,<Emy>,<Gam>,DGam/Gam,cov<xPx>,cov<yPy>,cov<zPz>,n_over_ne
  write(11,'(100e14.5)') sim_parameters%sim_time*c,mu_x,mu_y,mu_z+mesh_par%z_max_moving_um,mu_px,mu_py,mu_pz,s_x,s_y,s_z,s_px,s_py, &
   s_pz,emittance_x,emittance_y,mu_gamma,dgamma_su_gamma,corr_x_px,corr_y_py,corr_z_pz, &
   background_density_value(1+int((mu_z(1)*plasma%k_p-mesh_par%z_min_moving)*one_over_dz),2)
  close(11)

  filename=TRIM(sim_parameters%path_integrated_diagnostics)//'kurtosis_'//b2str//'.dat'
  call open_file(OSys%macwin,filename)
  !1,  2,  3,  4,   5,   6,   7
  !t,K_x,K_y,K_z,K_px,K_py,K_pz
  write(11,'(100e14.5)') sim_parameters%sim_time*c, &
                         m4_x/s_x**4-3.D0,   m4_y/s_y**4-3.D0,   m4_z/s_z**4-3.D0, &
                         m4_px/s_px**4-3.D0, m4_py/s_py**4-3.D0, m4_pz/s_pz**4-3.D0
  close(11)

 END SUBROUTINE


  !--- --- --- for dcut particles --- --- ---!
 SUBROUTINE bunch_integrated_diagnostics_for_dcutparticles(number_bunch)
 integer, intent(in) :: number_bunch
 real(8) :: mu_x(1),mu_y(1),mu_z(1)          !spatial meam
 real(8) :: mu_px(1),mu_py(1),mu_pz(1)       !momenta mean
 real(8) :: s_x(1),s_y(1),s_z(1)                         !spatial variance
 real(8) :: s_px(1),s_py(1), s_pz(1)                      !momenta variance
 real(8) :: m4_x(1),m4_y(1),m4_z(1),m4_px(1),m4_py(1),m4_pz(1) !4th central moments for kurtosis
 real(8) :: mu_gamma(1),s_gamma(1),dgamma_su_gamma(1)     !gamma mean-variance
 real(8) :: corr_y_py(1),corr_z_pz(1),corr_x_px(1) !correlation transverse plane
 real(8) :: emittance_x(1),emittance_y(1) !emittance variables
 real(8) :: u1(1),u2(1),u3(1),u4(1),u5(1)
 character(1) :: b2str
 character*90 :: filename
 TYPE(simul_param) :: sim_par

  !---mask :: particle selection :: need to diagnose only not-cut particles, so only keeping dcut==1---!
  bunch(number_bunch)%part(:)%cmp(14)=bunch(number_bunch)%part(:)%cmp(8)
  !--- *** ---!

	mu_x(1)  = calculate_nth_moment(number_bunch,1,1,'nocentral')
	mu_y(1)  = calculate_nth_moment(number_bunch,1,2,'nocentral')
	mu_z(1)  = calculate_nth_moment(number_bunch,1,3,'nocentral')
	mu_px(1) = calculate_nth_moment(number_bunch,1,4,'nocentral')
	mu_py(1) = calculate_nth_moment(number_bunch,1,5,'nocentral')
	mu_pz(1) = calculate_nth_moment(number_bunch,1,6,'nocentral')

	s_x(1)   = sqrt(calculate_nth_moment(number_bunch,2,1,'central'))
	s_y(1)   = sqrt(calculate_nth_moment(number_bunch,2,2,'central'))
	s_z(1)   = sqrt(calculate_nth_moment(number_bunch,2,3,'central'))
	s_px(1)  = sqrt(calculate_nth_moment(number_bunch,2,4,'central'))
	s_py(1)  = sqrt(calculate_nth_moment(number_bunch,2,5,'central'))
	s_pz(1)  = sqrt(calculate_nth_moment(number_bunch,2,6,'central'))

  m4_x(1)  = calculate_nth_moment(number_bunch,4,1,'central')
	m4_y(1)  = calculate_nth_moment(number_bunch,4,1,'central')
	m4_z(1)  = calculate_nth_moment(number_bunch,4,1,'central')
	m4_px(1) = calculate_nth_moment(number_bunch,4,1,'central')
	m4_py(1) = calculate_nth_moment(number_bunch,4,1,'central')
	m4_pz(1) = calculate_nth_moment(number_bunch,4,1,'central')

	corr_x_px(1) = calculate_correlation(number_bunch,1,4)
	corr_y_py(1) = calculate_correlation(number_bunch,2,5)
	corr_z_pz(1) = calculate_correlation(number_bunch,3,6)

  !---!
  mu_gamma(1)     = calculate_nth_moment(number_bunch,1,7,'nocentral')
	s_gamma(1)      = sqrt(calculate_nth_moment(number_bunch,2,7,'central'))
	dgamma_su_gamma = s_gamma(1)/mu_gamma(1)

	!---!
	emittance_x = sqrt( s_x(1)**2 *s_px(1)**2 - corr_x_px(1)**2 )
	emittance_y = sqrt( s_y(1)**2 *s_py(1)**2 - corr_y_py(1)**2 )

  write(b2str,'(I1.1)') number_bunch
  filename=TRIM(sim_parameters%path_integrated_diagnostics)//'bunch_integrated_quantity_dcut_'//b2str//'.dat'
  call open_file(OSys%macwin,filename)
  !1  2   3   4   5    6    7     8      9      10     11      12      13     14    15    16     17       18       19       20        21
  !t,<X>,<Y>,<Z>,<Px>,<Py>,<Pz>,<rmsX>,<rmsY>,<rmsZ>,<rmsPx>,<rmsPy>,<rmsPz>,<Emx>,<Emy>,<Gam>,DGam/Gam,cov<xPx>,cov<yPy>,cov<zPz>,n_over_ne
  write(11,'(100e14.5)') sim_parameters%sim_time*c,mu_x,mu_y,mu_z+mesh_par%z_max_moving_um,mu_px,mu_py,mu_pz,s_x,s_y,s_z,s_px,s_py, &
   s_pz,emittance_x,emittance_y,mu_gamma,dgamma_su_gamma,corr_x_px,corr_y_py,corr_z_pz, &
   background_density_value(1+int((mu_z(1)*plasma%k_p-mesh_par%z_min_moving)*one_over_dz),2)
  close(11)

  filename=TRIM(sim_parameters%path_integrated_diagnostics)//'kurtosis_'//b2str//'.dat'
  call open_file(OSys%macwin,filename)
  !1,  2,  3,  4,   5,   6,   7
  !t,K_x,K_y,K_z,K_px,K_py,K_pz
  write(11,'(100e14.5)') sim_parameters%sim_time*c, &
                         m4_x/s_x**4-3.D0,   m4_y/s_y**4-3.D0,   m4_z/s_z**4-3.D0, &
                         m4_px/s_px**4-3.D0, m4_py/s_py**4-3.D0, m4_pz/s_pz**4-3.D0
  close(11)

 END SUBROUTINE


 !----------------------------------------------------------!
 !--- functions to calculate emittance and energy spread ---!
 !----------------------------------------------------------!
 real(8) FUNCTION calculate_energy_spread(number_bunch)
 integer, intent(in) :: number_bunch
 real(8) :: mu_gamma(1),s_gamma(1),dgamma_su_gamma(1)

  !---mask :: particle selection---!
  bunch(number_bunch)%part(:)%cmp(14)=1.
  !--- *** ---!
 mu_gamma(1) = calculate_nth_moment(number_bunch,1,7,'nocentral')
 s_gamma(1)  = sqrt(calculate_nth_moment(number_bunch,2,7,'central'))
 dgamma_su_gamma = s_gamma(1)/mu_gamma(1)
 !--- *** ---!
 calculate_energy_spread = dgamma_su_gamma(1)
END FUNCTION calculate_energy_spread


real(8) FUNCTION calculate_emittance_x(number_bunch)
integer, intent(in) :: number_bunch
real(8) :: s_x(1),s_px(1),corr_x_px(1),emittance_x(1)

 !---mask :: particle selection---!
 bunch(number_bunch)%part(:)%cmp(14)=1.
 !--- *** ---!

s_x(1)  = sqrt(calculate_nth_moment(number_bunch,2,1,'central'))
s_px(1) = sqrt(calculate_nth_moment(number_bunch,2,4,'central'))
corr_x_px(1) = calculate_correlation(number_bunch,1,4)
!---!
emittance_x = sqrt( s_x(1)**2 *s_px(1)**2 - corr_x_px(1)**2 )
calculate_emittance_x = emittance_x(1)
END FUNCTION calculate_emittance_x


real(8) FUNCTION calculate_emittance_y(number_bunch)
integer, intent(in) :: number_bunch
real(8) :: s_y(1),s_py(1),corr_y_py(1),emittance_y(1)

 !---mask :: particle selection---!
 bunch(number_bunch)%part(:)%cmp(14)=1.
 !--- *** ---!
s_y(1)  = sqrt(calculate_nth_moment(number_bunch,2,2,'central'))
s_py(1) = sqrt(calculate_nth_moment(number_bunch,2,5,'central'))
corr_y_py(1) = calculate_correlation(number_bunch,2,5)
!---!
emittance_y = sqrt( s_y(1)**2 *s_py(1)**2 - corr_y_py(1)**2 )
calculate_emittance_y = emittance_y(1)
END FUNCTION calculate_emittance_y


 !--- --- ---!
 end module moments
 !--- --- ---!
