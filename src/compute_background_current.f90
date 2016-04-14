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


MODULE Compute_plasma_current

USE my_types
USE use_my_types

IMPLICIT NONE

CONTAINS

! Compute_plasma_3current_FDTD
! Computes the backgound plasma electron current in the
! longitudinal and transverse  direction.

    SUBROUTINE Compute_plasma_electron_current

   	USE my_types

   	IMPLICIT NONE

   	INTEGER i,j,iter
   	REAL(8), DIMENSION(mesh_par%Nzm,mesh_par%Nxm) :: Jpe_r, Jpe_z,Jpe_r_for_FDTD,Jpe_z_for_FDTD
    REAL(8), DIMENSION(mesh_par%Nzm,mesh_par%Nxm) :: beta_x,beta_z
    REAL(8), DIMENSION(mesh_par%Nxm) :: r_mesh_EM
    REAL(8) :: threshold_factor=1e-6
    REAL(8) :: test_ne
!---------------------------------------------------------------------------------!
!                           Z axis                                                !
!                                                                                 !
!    ghost cell        physical domain         ghost cell                         !
!  |______________|__________________________|______________|                     !
!  |              |                          |              |                     !
!  1         Node_min_z                 Node_max_z       Node_end_z               !
!                                                                                 !
!---------------------------------------------------------------------------------!
!              Yee Lattice (See notes for centering details)                      !
!---------------------------------------------------------------------------------!
!                           R axis                                                !
!                                                                                 !
!      Axis    physical domain                   ghost cell                       !
!  |________|________________________________|______________|                     !
!  |        |                                |              |                     !
!  1     Node_min_r=2                    Node_max_r       Node_end_r              !
!                                                                                 !
!---------------------------------------------------------------------------------!




!---------------------------------------------------------------------------------!
!                                     Z axis                                      !
!         								                                                        !
!     ghost cell                  physical domain                ghost cell       !
!  |______________|___________________________________________|____________|      !
!  |              |            |                |             |            |      !
!       1           Node_min_z                     Node_max_z   Node_end_z        !
!---------------------------------------------------------------------------------!
!                                    Fluid Mesh				                            !
!---------------------------------------------------------------------------------!
!                                      R axis                                     !
!									                                                                !
!     ghost cell                  physical domain                ghost cell       !
!  |______________|___________________________________________|____________|      !
!  |              |            |                |             |            |      !
!       1         ^ Node_min_z                     Node_max_z   Node_end_z        !
!                 |                                                               !
!               Z Axis                                                            !
!---------------------------------------------------------------------------------!



  ! Initialize auxiliary vectors and variables
  Jpe_r              = 0.
  Jpe_z              = 0.
  beta_x             = 0.
  beta_z             = 0.
  Jpe_r_for_FDTD		 = 0.
  Jpe_z_for_FDTD		 = 0.

	! ux, uz, beta_x, beta_z and ne are centered at the center of the cell

	!---------------------!
  !  Betas computation  !
  !---------------------!

  do i= 1,Node_end_z
  do j= 1,Node_end_r

    if ( mesh(i,j)%n_plasma_e.le.threshold_factor) then
      Jpe_r(i,j) = 0.
      Jpe_z(i,j) = 0.
    else
      beta_x(i,j) = mesh(i,j)%ux/sqrt( 1. + mesh(i,j)%ux**2+ mesh(i,j)%uz**2 + 1e-10 )
      beta_z(i,j) = mesh(i,j)%uz/sqrt( 1. + mesh(i,j)%ux**2+ mesh(i,j)%uz**2 + 1e-10 )

      Jpe_r(i,j) = mesh(i,j)%n_plasma_e*beta_x(i,j) !face centered
      Jpe_z(i,j) = mesh(i,j)%n_plasma_e*beta_z(i,j) !face centered
    endif

  enddo
  enddo

  !--------------------------------------------------------!
	!  Centers the currents as needed by FDTD field solver   !
  !--------------------------------------------------------!

	! Jr, Jz for FDTD are centered as in the Yee cell
    do i= 2,Node_max_z
    do j= 2,Node_max_r
      Jpe_r_for_FDTD(i,j) = &
                 ( Jpe_r(i  ,j)*(j-2+0.5)+Jpe_r(i  ,j+1)*(j+1-2+0.5) ) / (4.*(j+1-2)) &
            +    ( Jpe_r(i+1,j)*(j-2+0.5)+Jpe_r(i+1,j+1)*(j+1-2+0.5) ) / (4.*(j+1-2))
            !Surfaces are at: inner-surface (DeltaR*(j-2)+DeltaR/2.)
            !                 outer-surface (DeltaR*(j+1-2)+DeltaR/2.)
      Jpe_z_for_FDTD(i,j) =  Jpe_z(i,j)
    enddo
    enddo

	!--------------------------!
  !    Boundary conditions   !
  !      for currents and    !
	!      electron density    !
  !--------------------------!

	! lower boundary
	do i= Node_min_z,Node_max_z
		Jpe_r_for_FDTD(i,1)  = 0.
		Jpe_z_for_FDTD(i,1)  = Jpe_z_for_FDTD(i,Node_min_r)
		mesh(i,1)%n_plasma_e = mesh(i,Node_min_r)%n_plasma_e
  enddo
  ! upper boundary
  Jpe_r_for_FDTD(:,Node_end_r)  = Jpe_r_for_FDTD(:,Node_max_r)
  Jpe_z_for_FDTD(:,Node_end_r)  = Jpe_z_for_FDTD(:,Node_max_r)
  mesh(:,Node_end_r)%n_plasma_e = mesh(:,Node_max_r)%n_plasma_e

  !---backward substitution---!
  mesh(:,:)%Jpe_r      = Jpe_r_for_FDTD(:,:)
  mesh(:,:)%Jpe_z      = Jpe_z_for_FDTD(:,:)

   	return

   	END SUBROUTINE




    SUBROUTINE set_initial_plasma_density
    ! define capillary density at first iteration
    REAL :: test_ne
		INTEGER j

		mesh(:,:)%ux 			    = 0.
		mesh(:,:)%uz 			    = 0.
    mesh(:,:)%n_plasma_e 	= 0.

		if (sim_parameters%ramps_order.eq.1) then ! linear ramp at capillary entrance

			if (sim_parameters%I_start_ramp_peak_position.gt.2) then
				! ramp initially inside the window

				if (sim_parameters%order_capillary_density_z.eq.0) then ! uniform density profile along z
          mesh(1:sim_parameters%I_start_ramp_peak_position,1:mesh_par%NRmax_plasma)%n_plasma_e = 1.
				else if (sim_parameters%order_capillary_density_z.eq.2) then ! parabolic density profile along z

				endif

				do j=sim_parameters%I_start_ramp_peak_position,(sim_parameters%I_start_ramp_peak_position+sim_parameters%I_start_ramp_length_in_window)
					test_ne				 = 1.-1./sim_parameters%I_start_ramp_length_in_window*(1.*(j-sim_parameters%I_start_ramp_peak_position) )
					if (test_ne .le. 0.) test_ne = 0.
          mesh(j,1:1:mesh_par%NRmax_plasma)%n_plasma_e  = test_ne
				enddo


			else if (sim_parameters%I_start_ramp_peak_position.eq.2)	then
				! ramp initially in part inside the window, in part outside the window
				do j=2,(1+sim_parameters%I_start_ramp_length_in_window)
					test_ne				 = sim_parameters%start_n_peak_in_window*(1.-1./sim_parameters%I_start_ramp_length_in_window*(1.*j) )
					if (test_ne .le. 0.) test_ne = 0.
          mesh(j,1:mesh_par%NRmax_plasma)%n_plasma_e  = test_ne
				enddo
        mesh(1,1:mesh_par%NRmax_plasma)%n_plasma_e = mesh(2,1:mesh_par%NRmax_plasma)%n_plasma_e
			endif

		else if (sim_parameters%ramps_order.eq.2) then ! longitudinal parabolic ramp
			!DO NOT DELETE UNTIL PARABOLIC Z PROFILE IS IMPLEMENTED!!!

			!~do j=2,(sim_parameters%I_parabola_length_in_window+1)
			!~
			!~	ne (j,1:Node_end_r)  = - sim_parameters%I_parabola_normalization_factor* &
			!~	 (j-sim_parameters%I_origin_z_axis - 1 + sim_parameters%I_parabola_start) &
			!~	*(j-sim_parameters%I_origin_z_axis - 1 + sim_parameters%I_parabola_end)
			!~
			!~enddo
			!~
			!~ne (1,1:Node_end_r) = ne (2,1:Node_end_r)

		endif

		if (sim_parameters%order_capillary_density_r.eq.2) then
			! multiplication by radial factor
			do j= 2,Node_max_r
				mesh(:,j)%n_plasma_e 	  = mesh(:,j)%n_plasma_e*radial_factor(j)
			enddo
			! BC
			mesh(:,Node_end_r)%n_plasma_e = mesh(:,Node_max_r)%n_plasma_e ! upper boundary
			mesh(:,1         )%n_plasma_e = mesh(:,2         )%n_plasma_e ! lower boundary
		endif


	END SUBROUTINE set_initial_plasma_density


END MODULE
