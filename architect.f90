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

PROGRAM Architect

USE my_types
USE use_my_types
USE read_input_module
USE Make_a_Mesh
USE Diagnostics_on_Bunches
USE Fields_FDTD
USE Move_Window_FDTD
USE ComputeCurrentFDTD
USE Compute_beam_current_FDTD
USE MoveParticle_FDTD
USE FluidAdvance_FDTD
USE bunch_generation
USE Data_dumping
USE pstruct_data
USE architect_class_structure
USE moments
USE initialize_bunch
USE utilities
USE linear_algebra
USE init_fields



IMPLICIT NONE

LOGICAL getFileFlag
INTEGER i,j,ii,ss
REAL(8) ::  mu_z,sigma_z,total_run_distance
INTEGER Lapl_dim






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            LOAD INPUT DATA, SET PARAMETERS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   sim_parameters%iter=1

   call read_input

   sim_parameters%dim_Laplacian=2*bunch_initialization%init_width_z*mesh_par%Nsample_z* &
   (bunch_initialization%init_width_r*mesh_par%Nsample_r/2)

   total_run_distance = plasma%l_acc + & !in um
   						plasma%lambda_p * (sim_parameters%distance_capillary+sim_parameters%start_ramp_length &
                                          +sim_parameters%end_ramp_length+sim_parameters%distance_after_end_ramp )

   call generate_output_tree

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            BUNCH INIT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   call init_bunch

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            WINDOW INIT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Initializes simulation window parameters
   call init_window

   ! Creates the mesh and capillary ramp
   call Kernel_Make_a_mesh
   call define_capillary_ramp


	! Initialize simulation time
	sim_parameters%sim_time	= 0.

	! First global print
	call first_print_at_screen

	! Write Started Flag
	!File Flags: begin
	call SetFileFlag  ('==started==')


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            TIME EVOLUTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Sets plasma background density,
   ! computes charge and fields at first iteration

   !Self-consistent-bunch field initialization
  if( bunch_initialization%self_consistent_field_bunch==0 ) then 		!no field initialization, only for comparison purposes
		call init_null_EM_fields
		call init_external_Bfields
		call set_initial_plasma_density
		call Kernel_ComputeCurrent_FDTD
	elseif( bunch_initialization%self_consistent_field_bunch==1 ) then 	!coax shells
		call set_initial_plasma_density
		call Kernel_ComputeCurrent_FDTD
		call init_EM_fields_coax_shells
		call init_external_Bfields
	elseif(bunch_initialization%self_consistent_field_bunch>1) then 	!LU or SOR: LU option (2), SOR option (3)
		call init_EM_fields
		call set_initial_plasma_density
		call Kernel_ComputeCurrent_FDTD
		!call analytic_rho_only_first_bunch
		call init_external_Bfields
	endif


   ! New bunch integrated diagnostic without the cut + sliced diagnostics
   do i=1,bunch_initialization%n_total_bunches
		mu_z    = calculate_nth_moment_bunch(i,1,3)
		sigma_z = calculate_nth_central_moment_bunch(i,2,3)

		call bunch_sliced_diagnostics(i)
		call bunch_integrated_diagnostics(i)
   enddo

	if (sim_parameters%diagnostics_with_dcut.eq.1) then
		call apply_Sigma_cut(sim_parameters,plasma%k_p)
		do i=1,bunch_initialization%n_total_bunches
			mu_z    = calculate_nth_moment_bunch_dcut(i,1,3)
			sigma_z = calculate_nth_central_moment_bunch_dcut(i,2,3)

			call bunch_sliced_diagnostics_dcut(i,mu_z,sigma_z)
			call bunch_integrated_diagnostics_dcut(i)
		enddo
	endif

   ! Plasma and bunch diagnostics
   call final_data_dump
   call dump_input_file

! ----------- MAIN TIME LOOP -----------------------------
! --------------------------------------------------------

   main_loop: do

     !write(*,'(A,100e14.5)') 'time >',sim_parameters%dt,Dt,sim_parameters%sim_time

		!---check suspension flag---!
		if(getFileFlag('==suspend==')) then
			call SetFileFlag('==suspended==')
			call UnsetFileFlag('==suspend==')
			call UnsetFileFlag('==started==')
			write(*,*) 'Simulation suspended, at:',sim_parameters%zg
			stop
		endif


    ! ---- BULK
		! Computes plasma and beam current
		call Kernel_ComputeCurrent_FDTD
		! Advances electromagnetic fields
		call Kernel_Fields_FDTD_COMB
		! Fluid advance
		call Kernel_FluidAdvance_FDTD
		! Particle pusher
		call Kernel_MoveParticle_FDTD

	  ! New bunch integrated diagnostic without the cut + sliced diagnostics
	  sim_parameters%IntDeltaOutput=abs(sim_parameters%sim_time*c-sim_parameters%IntLastOutput)
      if(mod(sim_parameters%iter,sim_parameters%output_Integrated_params_nstep).eq.0 &
      .or. sim_parameters%IntDeltaOutput>sim_parameters%output_Integrated_params_dist) then

		do i=1,bunch_initialization%n_total_bunches
			mu_z    = calculate_nth_moment_bunch(i,1,3)
			sigma_z = calculate_nth_central_moment_bunch(i,2,3)

			call bunch_sliced_diagnostics(i)
			call bunch_integrated_diagnostics(i)
		enddo

		if (sim_parameters%diagnostics_with_dcut.eq.1) then
			call apply_Sigma_cut(sim_parameters,plasma%k_p)
			do i=1,bunch_initialization%n_total_bunches
				mu_z    = calculate_nth_moment_bunch_dcut(i,1,3)
				sigma_z = calculate_nth_central_moment_bunch_dcut(i,2,3)

				call bunch_sliced_diagnostics_dcut(i,mu_z,sigma_z)
				call bunch_integrated_diagnostics_dcut(i)
			enddo
		endif
		sim_parameters%IntLastOutput=sim_parameters%sim_time*c
	  endif

      ! Plasma and bunch data dump
	  call data_dump

      ! Advances simulation iteration
      sim_parameters%iter 	= sim_parameters%iter+1


      !--------------------------------------------------------!
      !---> Re-initialise bunch density whethere necessary <---!
      !--------------------------------------------------------!
      call ReInitialise_EB_Bunch_fields
      !---> <---!

	  ! Advances simulation window
	  call control_window

      ! continous diagnostic
      if(mod(sim_parameters%iter,10)==0) then
			write(*,*) 'Iteration =',sim_parameters%iter,' - Moving window position =',sim_parameters%zg
      endif

      ! Print at screen, when data are dumped
      if (mod(sim_parameters%iter,sim_parameters%output_grid_nstep).eq.0 ) then
			call print_at_screen
      endif

	  ! Advance simulation time
      !sim_parameters%sim_time=sim_parameters%sim_time+sim_parameters%dt
      sim_parameters%sim_time=sim_parameters%iter*sim_parameters%dt

      if(abs(calculate_nth_moment_bunch(1,1,3)).gt.total_run_distance) then
        write(*,*) 'End of run'
        write(*,'(A,1e14.5)') 'total distance run:',abs(calculate_nth_moment_bunch(1,1,3))
        goto 123
    	endif

   enddo main_loop

123 continue

   ! Final dump
   call Compute_bunch_density
   call final_data_dump


   ! Deallocate variables
   DEALLOCATE(x_mesh,z_mesh,mesh)
   if(sim_parameters%order_capillary_density_r.eq.2) then
		 DEALLOCATE(radial_factor)
   endif
	!end of run flags
	call setFileFlag('==completed==')
	call UnsetFileFlag('==started==')

   stop
   END
