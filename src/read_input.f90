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

module read_input_module

USE my_types
USE use_my_types
USE pstruct_data
USE architect_class_structure



 implicit none

 !--- --- ---!
 contains
 !--- --- ---!

subroutine read_input
    integer		:: i
    real :: Delta

    NAMELIST / in_phys_pars / plasma, sim_parameters, mesh_par

    !pre-definitions
    sim_parameters%time_step 	= .5
		sim_parameters%Output_format= 1 ! Binary PS and grid output
		sim_parameters%jump_grid	= 1 ! no jump of grid points in output
		sim_parameters%reduced_PS 	= 1
		sim_parameters%jump_PS		= 1 ! no jump of particles in PSoutput

		!--- ---!
      write(*,*) "preset and read twiss nml"
      call preset_twiss_nml
      call read_twiss_nml
		!--- ---!

    !--- ---!
      write(*,*) "preset and read poloidal magnetic field values"
      call preset_Bpoloidal
      call read_Bpoloidal_nml
		!--- ---!

    !--- ---!
      write(*,*) "preset capillary and read input file"
      call preset_capillary
      open(9,file='architect.nml',status='old')
      READ(9,NML=in_phys_pars,ERR=11,END=11)
  11 continue
      close(9)
    !--- ---!
      if( mesh_par%Rmax_plasma .eq. -1.D0 ) mesh_par%Rmax_plasma=mesh_par%Rmax
    !--- ---!

    !--- ---!
		write(*,*) "Initialise bunch(es) and self consistent fields"
		call preset_bunch_initialization
		!~write(*,*) bunch_initialization%ChargeB(1)
		call read_nml_bunch_initialization
		!~write(*,*) bunch_initialization%ChargeB(1)
		!~stop
    !--- ---!

    !--- ---!
		write(*,*) "preset and read operating system nml"
		call preset_OS_nml
		call read_OS_nml
    !--- ---!

        if(sim_parameters%jump_grid.le.0) then
            write(*,*) 'jump_grid parameter lesser than 0.'
            write(*,*) 'Setting it to 1'
            sim_parameters%jump_grid=1
        endif

        if(sim_parameters%jump_PS.le.0) then
            write(*,*) 'jump_PS parameter lesser than 0.'
            write(*,*) 'Setting it to 1'
            sim_parameters%jump_PS=1
        endif



		!---parameters control for Laplacian matrix---!
  !~      bunch_initialization%init_width_r	=	int(max(mesh_par%Rmax,1.*bunch_initialization%init_width_r))
  !~      bunch_initialization%init_width_z	=	min(6,bunch_initialization%init_width_z)


		bunch_initialization%init_width_r	=	int(min(mesh_par%Rmax,1.*bunch_initialization%init_width_r))
    bunch_initialization%init_width_z	=	min(17,bunch_initialization%init_width_z)

		!--- ---!



        ! simulation and plasma parameters
        sim_parameters%Nbunches	= bunch_initialization%n_total_bunches 	! number of bunches
        plasma%lambda_p			= sqrt(1.11491e21/(plasma%n0))			! Plasma wavelength (um)
        plasma%k_p				= 2.*pi/plasma%lambda_p					! Plasma wavenumber (rad/um)
        plasma%omega_p 			= plasma%k_p*c   						! Plasma pulsation (rad/fs)

        ! number of particles
        sim_parameters%Np=sum(bunch_initialization%n_particles(1:bunch_initialization%n_total_bunches))

        !-from CFL to Dt in (fs)-!
        Delta = min( 2.*bunch_initialization%bunch_s_x(1)*1e-6/real(mesh_par%Nsample_r) , &
                        bunch_initialization%bunch_s_z(1)*1e-6/real(mesh_par%Nsample_z)   )
        sim_parameters%dt = sim_parameters%time_step * Delta / 3e8 / 1e-15

        write(*,*) "input file read"
        write(*,*)
end subroutine read_input





subroutine generate_output_tree
        character	:: command*200


		select case (OSys%macwin)

			case (0)
				sim_parameters%out_dir   = 'out/'
	            sim_parameters%out_root  = ''

    	        sim_parameters%path_PS 	 = 'out/PS/'
        	    sim_parameters%path_grid = 'out/2D/'
            	call system('mkdir -p out/PS')
            	call system('mkdir -p out/2D')

			case (1)
				sim_parameters%out_dir   = TRIM(ADJUSTL(Osys%output_dir))//'\'
            	sim_parameters%out_root  = TRIM(ADJUSTL(Osys%output_dir))//'\'
            	sim_parameters%path_PS 	 = TRIM(ADJUSTL(Osys%output_dir))//'PS\'
            	sim_parameters%path_grid = TRIM(ADJUSTL(Osys%output_dir))//'2D\'
            	sim_parameters%path_1D   = TRIM(ADJUSTL(Osys%output_dir))//'1D\'

            command='mkdir '//TRIM(ADJUSTL(sim_parameters%path_PS))
            write(*,*) command
            call system(command)

            command='del /S /Q '//TRIM(ADJUSTL(sim_parameters%path_PS))//'*.arch'
            write(*,*) command
            call system(command)

            command='mkdir '//TRIM(ADJUSTL(sim_parameters%path_grid))
            write(*,*) command
            call system(command)

            command='del /S /Q '//TRIM(ADJUSTL(sim_parameters%path_grid))//'*.arch'
            write(*,*) command
            call system(command)

            command='mkdir '//TRIM(ADJUSTL(sim_parameters%path_1D))
            write(*,*) command
            call system(command)

            command='del /S /Q '//TRIM(ADJUSTL(sim_parameters%path_1D))//'*.dat'
            write(*,*) command
            call system(command)

            command='del /Q '//TRIM(ADJUSTL(Osys%output_dir))//'*.dat'
            write(*,*) command
            call system(command)

        end select


end subroutine generate_output_tree


!CAPILLARY
subroutine preset_capillary

	!Capillary channel density profile
	sim_parameters%order_capillary_density_z 	= 0		! uniform density along z
	sim_parameters%order_capillary_density_r 	= 0		! uniform density along r
	sim_parameters%distance_capillary			= .1	! distance between the center of z axis and the initial capillary ramp, in lambda_p

	!Ramps parameters
	sim_parameters%ramps_order					= 1 	! linear ramps
	sim_parameters%start_ramp_length			= 0. 	! total entrance ramp length
	sim_parameters%end_ramp_length 				= 0.	! ramp length at the end of the capillary
	sim_parameters%distance_after_end_ramp 		= 1.	! distance traveled by driver after exit ramp

  !Plasma extension
  mesh_par%Rmax_plasma = -1.D0

end subroutine preset_capillary

!TWISS
subroutine preset_twiss_nml
	twiss%L_TWISS=.false.
    twiss%alpha_new_factor(:)=1.
    twiss%beta_new_factor(:) =1.
end subroutine preset_twiss_nml

subroutine read_twiss_nml
		NAMELIST / twiss_par / twiss
        open(9,file='architect.nml',status='old')
            READ(9,NML=twiss_par,ERR=12,END=12)
        12 continue
        close(9)
end subroutine read_twiss_nml



!Poloidal magnetic field namelist control!TWISS
subroutine preset_Bpoloidal
	Bpoloidal%L_Bpoloidal=.false.
  Bpoloidal%B_ex_poloidal=0.0
  Bpoloidal%radius_poloidal=1.0
end subroutine preset_Bpoloidal

subroutine read_Bpoloidal_nml
  NAMELIST / Bpoloidal_par / Bpoloidal
  open(11,file='architect.nml',status='old')
  READ(11,NML=Bpoloidal_par,ERR=13,END=13)
  13 continue
  close(11)
  !convert
  !Bpoloidal%B_ex_poloidal=Bpoloidal%B_ex_poloidal / (electron_mass*c*c*plasma%omega_p/electron_charge)
  !Bpoloidal%radius_poloidal=Bpoloidal%radius_poloidal*plasma%plasma%k_p
end subroutine read_Bpoloidal_nml



!OSystem
subroutine preset_OS_nml
    OSys%macwin=0 !0=mac
    OSys%output_dir='.'
end subroutine preset_OS_nml


subroutine read_OS_nml
		NAMELIST / OS / OSys

        open(9,file='architect.nml',status='old')
            READ(9,NML=OS,ERR=13,END=13)
        13 continue
        close(9)

end subroutine read_OS_nml


subroutine preset_bunch_initialization
    bunch_initialization%l_bunch_internal_init=.true.
    !1=coax-shells,2=LU,3=SOR
    bunch_initialization%self_consistent_field_bunch=1
    !width(in sigma_r)of the initialization domain around bunch
    bunch_initialization%init_width_r=1000
    !width(in sigma_z)of the initialization domain around bunch
    bunch_initialization%init_width_z=3
		bunch_initialization%iter_max=2000
		bunch_initialization%maxnorm=1e-4
		bunch_initialization%wsor=1.45

		!just one driver
		bunch_initialization%n_total_bunches=1
		bunch_initialization%chargeB(:)= 0.200
		bunch_initialization%n_particles(:)= 50000
		bunch_initialization%bunch_s_x(:)=8.0
		bunch_initialization%bunch_s_y(:)=8.0
		bunch_initialization%bunch_s_z(:)=50.0
		bunch_initialization%bunch_gamma_m(:)=200.
		bunch_initialization%bunch_eps_x(:)=1.0
		bunch_initialization%bunch_eps_y(:)=1.0
		bunch_initialization%bunch_dgamma(:)=0.1
		bunch_initialization%db(:)=0.0
end subroutine preset_bunch_initialization


subroutine read_nml_bunch_initialization
		NAMELIST / BUNCHINIT / bunch_initialization
        open(9,file='architect.nml',status='old')
            READ(9,NML=BUNCHINIT,ERR=13,END=13)
		13 continue
        close(9)
end subroutine read_nml_bunch_initialization






 !--- --- ---!
 end module read_input_module
 !--- --- ---!
