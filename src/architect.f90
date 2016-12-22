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
USE dump_status
USE ion_background
USE ionisation_module



IMPLICIT NONE

LOGICAL getFileFlag
INTEGER i,j,ii,ss
REAL(8) ::  mu_z,sigma_z,total_run_distance
INTEGER Lapl_dim






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            LOAD INPUT DATA, SET PARAMETERS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   sim_parameters%iter=0

   call SetFileFlag('==started==')
   call read_input
   call write_read_nml

  !  total_run_distance = plasma%l_acc + & !in um
  !  						plasma%lambda_p * (sim_parameters%distance_capillary+sim_parameters%start_ramp_length &
  !                               +sim_parameters%end_ramp_length+sim_parameters%distance_after_end_ramp )
     total_run_distance = plasma%l_acc !in um

   call generate_output_tree

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    INITIALISE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call dimension_first_bunch
   call init_window
   call Kernel_Make_a_mesh
   call init_bunch
   call dt_calculation

	! Initialize simulation time
	sim_parameters%sim_time	= 0.

	! First global print
	call first_print_at_screen



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            TIME EVOLUTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Sets plasma background density,
   ! computes charge and fields at first iteration

   !Self-consistent-bunch field initialization
  if( bunch_initialization%self_consistent_field_bunch==0 ) then 		!no field initialization, only for comparison purposes
		call init_null_EM_fields
		call init_external_Bfields
		call set_initial_background_condition
    call initialise_ion_background
		call Kernel_ComputeCurrent_FDTD
	elseif( bunch_initialization%self_consistent_field_bunch==1 ) then 	!coax shells
		call set_initial_background_condition
    call initialise_ion_background
		call Kernel_ComputeCurrent_FDTD
		call init_EM_fields_coax_shells
		call init_external_Bfields
	elseif(bunch_initialization%self_consistent_field_bunch>1) then 	!LU or SOR: LU option (2), SOR option (3)
		call init_EM_fields
		call set_initial_background_condition
    call initialise_ion_background
		call Kernel_ComputeCurrent_FDTD
		!call analytic_rho_only_first_bunch
		call init_external_Bfields
	endif


  !--- Bunch Diagnostics ---!
  call bunch_diagnostics_Integrated_AllBunches
  ! Plasma and bunch diagnostics
  call final_data_dump
  call write_read_nml

! ----------- MAIN TIME LOOP -----------------------------
! --------------------------------------------------------

   main_loop: do

		!---check suspension flag---!
		if(getFileFlag('==suspend==')) then
			call SetFileFlag('==suspended==')
			call UnsetFileFlag('==suspend==')
			call UnsetFileFlag('==started==')
			write(*,*) 'Simulation suspended, at:',sim_parameters%zg
			stop
		endif

    call ionise

    ! --- --- --- BULK --- --- --- !
    if(sim_parameters%L_plasma_evolution) then
      ! *** full PLASMA evolution ***!
  		call Kernel_ComputeCurrent_FDTD ! Computes plasma and beam current
  		call Kernel_Fields_FDTD_bck    ! Advances electromagnetic fields background
      if(sim_parameters%L_Bunch_evolve) call Kernel_Fields_FDTD_bunch    ! Advances electromagnetic field bunch(es)
  		call Kernel_FluidAdvance_FDTD   ! Fluid advance
  		call Kernel_MoveParticle_FDTD   ! Particle pusher (bunch electrons)
      if(ionisation%L_ionisation) call ion_advection  ! Particle pusher: background ions
    else
      ! *** only Particle Tracking ***!
  		call Kernel_ComputeCurrent_FDTD ! Computes plasma and beam current
  		call Kernel_MoveParticle_FDTD   ! Particle pusher
    endif
    ! --- --- --- BULK --- --- --- !



    ! --- --- --- OUTPUT --- --- --- !
    call print_integrated_diagnostics
	  call data_dump
    call print_at_screen
    call print_cpu_file
    ! --- --- --- output --- --- --- !

    !--- un poco sporchetta da mettere a posto---!
    if(mod(sim_parameters%iter,sim_parameters%output_Integrated_params_nstep).eq.0 &
    .or. sim_parameters%IntDeltaOutput>sim_parameters%output_Integrated_params_dist.or.&
    sim_parameters%iter==1) then

    open(unit=19, file='file.out', form='formatted')
    do i=1,ionisation%tot_ions
      if(static_ion(i)%cmp(7)>0.) write(19,'(3e14.5)') static_ion(i)%cmp(1)/plasma%k_p,static_ion(i)%cmp(2)/plasma%k_p
    enddo
    close(19)
  endif
  ! --- --- --- output --- --- --- !


    !--------------------------------------------------------!
    !---> Re-initialise bunch density whethere necessary <---!
    !--------------------------------------------------------!
    call ReInitialise_EB_Bunch_fields
    !---> <---!

	  ! Advances simulation window
	  call control_window


	  ! Advance simulation time
      !sim_parameters%sim_time=sim_parameters%sim_time+sim_parameters%dt
      sim_parameters%sim_time=sim_parameters%iter*sim_parameters%dt

      if(abs(calculate_nth_moment_bunch(1,1,3)).gt.total_run_distance) then
        write(*,'(A,1e14.5)') ' > total run distance:',abs(calculate_nth_moment_bunch(1,1,3))
        write(*,'(A)') '--- --- --- End of run --- --- ---'
        goto 123
    	endif

      !--- Advances simulation iteration ---!
      sim_parameters%iter = sim_parameters%iter+1

   enddo main_loop

123 continue


   !--- Final dump ---!
   call Compute_bunch_density
   call final_data_dump


    ! Deallocate variables
    DEALLOCATE(x_mesh,z_mesh,mesh)
    if( allocated(mesh_util%Bphi_BC_Left) ) DEALLOCATE(mesh_util%Bphi_BC_Left)
    !end of run flags
    call setFileFlag('==completed==')
    call UnsetFileFlag('==started==')

   stop
   END
