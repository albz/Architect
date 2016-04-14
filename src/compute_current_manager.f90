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

MODULE ComputeCurrentFDTD

USE my_types
USE use_my_types
USE Compute_beam_current_FDTD
USE Compute_plasma_current
USE pstruct_data
USE architect_class_structure
USE utilities



IMPLICIT NONE

CONTAINS




   SUBROUTINE Kernel_ComputeCurrent_FDTD
   IMPLICIT NONE
   INTEGER inc,tt
   REAL(8) :: k_0,n_0,conv

   k_0 = plasma%k_p
   n_0 = plasma%n0

   call Compute_beam_3current_FDTD

   ! Conversion factor
   conv     = (1./(2.*pi*mesh_par%dzm*mesh_par%dxm))*(k_0*1.e4)**3   	!Divide by cell volume (1/r factor included in subroutine Compute_beam_3current_FDTD)
   conv     = conv/n_0     						!dimensionless

   mesh%rho =   mesh%rho*conv
   mesh%Jz  =   mesh%Jz*conv
   mesh%Jr  =   mesh%Jr*conv

   if ( sim_parameters%Fluid_Scheme==0 ) then
	     write(*,*) 'error: old upwind centering is not supported in this version'
	     stop
   !call Compute_plasma_electron_current_FDTD
   else if (sim_parameters%Fluid_Scheme==2   .or.   sim_parameters%Fluid_Scheme==1) then
       call Compute_plasma_electron_current
   endif

   return
   END SUBROUTINE


   SUBROUTINE define_capillary_ramp
		!INTEGER s,e
		write(*,*) 'defining plasma capillary parameters'
		call define_capillary_ramps
		call define_capillary_density_z
		call define_capillary_density_r
		write(*,*) 'plasma capillary parameters defined'
		write(*,*)
    END SUBROUTINE define_capillary_ramp


	SUBROUTINE define_capillary_ramps
		write(*,*) 'Defining capillary ramps'
		if(sim_parameters%ramps_order.eq.0) then
			write(*,*) 'Error: value of ramps_order not valid'
			stop
		else if (sim_parameters%ramps_order.eq.1) then ! linear ramps
			write(*,*) 'Linear ramps at entrance and exit of capillary'
			call define_capillary_linear_ramps
		else if (sim_parameters%ramps_order.eq.2) then ! parabolic ramps
			write(*,*) 'Error: value of ramps_order not yet supported'
			stop
		endif

		write(*,*) 'Entrance ramp length (lp) = ',sim_parameters%start_ramp_length
		write(*,*) 'Exit ramp length (lp)     = ',sim_parameters%end_ramp_length
		write(*,*)
		write(*,*) 'Capillary lengths defined'
		write(*,*)

	END SUBROUTINE define_capillary_ramps


	SUBROUTINE define_capillary_linear_ramps
		INTEGER :: j

		! whole initial ramp outside window
			if ( sim_parameters%distance_capillary.ge.sim_parameters%Left_Domain_boundary) then
				write(*,*) sim_parameters%distance_capillary,sim_parameters%Left_Domain_boundary
				write(*,*) 'Error: capillary ramp outside window'
				stop
			endif

			! whole ramp inside window
			if ( (sim_parameters%distance_capillary+sim_parameters%start_ramp_length).le.sim_parameters%Left_Domain_boundary) then
				sim_parameters%start_ramp_length_in_window = sim_parameters%start_ramp_length
			else ! ramp partially outside window
				sim_parameters%start_ramp_length_in_window = sim_parameters%Left_Domain_boundary - sim_parameters%distance_capillary
			endif

			sim_parameters%start_ramp_peak_position = sim_parameters%Left_Domain_boundary- &
					(sim_parameters%distance_capillary+sim_parameters%start_ramp_length)

			if (sim_parameters%start_ramp_peak_position.ge.0.) then
				! whole ramp inside window
				sim_parameters%I_start_ramp_peak_position		= Nint(sim_parameters%start_ramp_peak_position/mesh_par%dzm*2.*pi)
			else
				! ramp partially outside window
				sim_parameters%I_start_ramp_peak_position		= 2
			endif

			sim_parameters%I_start_ramp_length            		= max(1,Nint(sim_parameters%start_ramp_length/mesh_par%dzm*2.*pi))
			sim_parameters%I_start_ramp_length_in_window  		= max(1,Nint(sim_parameters%start_ramp_length_in_window/mesh_par%dzm*2.*pi))-2
			sim_parameters%I_start_ramp_length_out_window 		= sim_parameters%I_start_ramp_length-sim_parameters%I_start_ramp_length_in_window

			if ((sim_parameters%distance_capillary+sim_parameters%start_ramp_length).gt.sim_parameters%Left_Domain_boundary) then
			sim_parameters%start_n_peak_in_window 	= (sim_parameters%start_ramp_length_in_window)/sim_parameters%start_ramp_length
			endif

			sim_parameters%start_exit_ramp	= plasma%l_acc+( sim_parameters%distance_capillary + &
															 sim_parameters%start_ramp_length    &
															-sim_parameters%Left_Domain_boundary &
															)*plasma%lambda_p ! zg value for which the window boundary reaches the exit ramp

			sim_parameters%end_exit_ramp	= sim_parameters%start_exit_ramp + sim_parameters%end_ramp_length*plasma%lambda_p

			if(sim_parameters%order_capillary_density_r.eq.2) then
				sim_parameters%R_capillary_adim_squared = (sim_parameters%R_capillary*plasma%k_p)**2
				allocate( radial_factor(mesh_par%Nxm) )
				do j=1, mesh_par%Nxm
					radial_factor(j) = 1.- ( x_mesh_shifted(j) )**2/sim_parameters%R_capillary_adim_squared + &
									   sim_parameters%fraction_n0_capillary_wall*( x_mesh_shifted(j) )**2/sim_parameters%R_capillary_adim_squared
				enddo
			endif

	END SUBROUTINE define_capillary_linear_ramps



	SUBROUTINE define_capillary_density_z

		if 		(sim_parameters%order_capillary_density_z.eq.0) then ! uniform density along z
			write(*,*) 'Uniform capillary density along z, n0 (cm^-3) = ', plasma%n0
			write(*,*)
		else if (sim_parameters%order_capillary_density_z.eq.1) then ! not valid
			write(*,*) 'Error: value of order_capillary_density_z not valid'
			stop
		else if (sim_parameters%order_capillary_density_z.eq.2) then ! parabolic density along z
			write(*,*) 'Error: value of order_capillary_density_z not yet supported'
			stop
		endif

	END SUBROUTINE define_capillary_density_z


	SUBROUTINE define_capillary_density_r
		REAL(8) :: rad_factor_window_boundary
		if 		(sim_parameters%order_capillary_density_r.eq.0) then ! uniform density along r
			write(*,*) 'Uniform capillary density along r'
			write(*,*)
		else if (sim_parameters%order_capillary_density_r.eq.1) then ! not valid
			write(*,*) 'Error: value of order_capillary_density_r not valid'
			stop
		else if (sim_parameters%order_capillary_density_r.eq.2) then ! parabolic density along r
			write(*,*) 'Parabolic capillary density along r'
			write(*,*) 'Capillary radius                              = ',sim_parameters%R_capillary,' um'
			write(*,*) 'n_e on axis            (if uniform z profile) = ',plasma%n0,' cm^-3'
			rad_factor_window_boundary = 1.- ( x_mesh_shifted(Node_max_r) )**2/sim_parameters%R_capillary_adim_squared + &
									   sim_parameters%fraction_n0_capillary_wall*( x_mesh_shifted(Node_max_r) )**2/sim_parameters%R_capillary_adim_squared
			write(*,*) 'n_e at window boundary (if uniform z profile) = ',plasma%n0*rad_factor_window_boundary,' cm^-3'
			write(*,*) 'n_e at capillary wall  (if uniform z profile) = ',sim_parameters%fraction_n0_capillary_wall*plasma%n0,' cm^-3'
			write(*,*)
		endif

	END SUBROUTINE define_capillary_density_r


END MODULE


!~		PREVIOUS VERSION CAPILLARY DEFINITION. DO NOT DELETE UNTIL PARABOLIC Z PROFILE IS IMPLEMENTED!!!
!~		write(*,*) 'defining plasma capillary parameters'
!~		if (sim_parameters%longitudinal_density_profile.eq.0) then ! linear ramp at capillary entrance
!~			write(*,*) 'plasma initial density shape in the capillary: linear ramp'
!~			! Initial density ramp parameters
!~
!~			! whole initial ramp outside window
!~			if ( sim_parameters%distance_capillary.ge.sim_parameters%Left_Domain_boundary) then
!~				write(*,*) 'Error: capillary ramp outside window'
!~				stop
!~			endif
!~
!~			! whole ramp inside window
!~			if ( (sim_parameters%distance_capillary+sim_parameters%ramp_length).le.sim_parameters%Left_Domain_boundary) then
!~				sim_parameters%ramp_length_in_window = sim_parameters%ramp_length
!~			else ! ramp partially outside window
!~				sim_parameters%ramp_length_in_window = sim_parameters%Left_Domain_boundary - sim_parameters%distance_capillary
!~
!~			endif
!~
!~			sim_parameters%ramp_peak_position = sim_parameters%Left_Domain_boundary- &
!~					(sim_parameters%distance_capillary+sim_parameters%ramp_length)
!~
!~			if (sim_parameters%ramp_peak_position.ge.0.) then
!~				! whole ramp inside window
!~				sim_parameters%I_ramp_peak_position		= Nint(sim_parameters%ramp_peak_position/mesh_par%dzm*2.*pi)
!~			else
!~				! ramp partially outside window
!~				sim_parameters%I_ramp_peak_position		= 2
!~			endif
!~
!~			sim_parameters%I_ramp_length            	= max(1,Nint(sim_parameters%ramp_length/mesh_par%dzm*2.*pi))
!~			sim_parameters%I_ramp_length_in_window  	= max(1,Nint(sim_parameters%ramp_length_in_window/mesh_par%dzm*2.*pi))
!~			sim_parameters%I_ramp_length_out_window 	= sim_parameters%I_ramp_length-sim_parameters%I_ramp_length_in_window
!~
!~			if ((sim_parameters%distance_capillary+sim_parameters%ramp_length).gt.sim_parameters%Left_Domain_boundary) then
!~			sim_parameters%n_peak_in_window = sim_parameters%ramp_length_in_window/sim_parameters%ramp_length
!~			endif
!~
!~			sim_parameters%start_exit_ramp	= plasma%l_acc+(sim_parameters%distance_capillary + &
!~															sim_parameters%ramp_length        &
!~															)*plasma%lambda_p
!~
!~			sim_parameters%end_exit_ramp	= sim_parameters%start_exit_ramp + sim_parameters%L_exit_ramp*plasma%lambda_p
!~
!~		else if (sim_parameters%longitudinal_density_profile.eq.1) then ! longitudinal parabola
!~			write(*,*) 'plasma initial density shape in the capillary: longitudinal parabola'
!~			sim_parameters%I_origin_z_axis					= max(1,Nint(sim_parameters%Left_Domain_boundary/mesh_par%dzm*2.*pi)   ) !computed excluding ghost cell
!~			sim_parameters%I_parabola_length_in_window		= &
!~			max(1,Nint((sim_parameters%Left_Domain_boundary - sim_parameters%distance_capillary)/mesh_par%dzm*2.*pi) )
!~			sim_parameters%I_parabola_start					= max(1,Nint(sim_parameters%distance_capillary/mesh_par%dzm*2.*pi) )
!~			sim_parameters%I_parabola_end					= &
!~			max(1,Nint((sim_parameters%distance_capillary)/mesh_par%dzm*2.*pi + plasma%l_acc/plasma%lambda_p/mesh_par%dzm*2.*pi  )     )
!~
!~			s = (sim_parameters%I_origin_z_axis + 1 - sim_parameters%I_parabola_start)
!~			e = (sim_parameters%I_origin_z_axis + 1 - sim_parameters%I_parabola_end  )
!~			sim_parameters%I_parabola_normalization_factor	= 1./(0.25*(s**2+e**2 )-0.5*e*s)
!~
!~		endif
