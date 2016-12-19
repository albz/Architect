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

integer, public :: iounit,ierr
character(100), public :: error_message
data iounit,ierr,error_message /9,0,' '/

 !--- --- ---!
 contains
 !--- --- ---!

subroutine read_input
    integer		:: i
    real :: Delta

    ! NAMELIST / in_phys_pars / plasma, sim_parameters

    !--- ---!
      write(*,*) 'preset :: read :: plasma parameters'
      call preset_plasma_parameters
      call read_plasma_parameters
		!--- ---!

    !--- ---!
      write(*,*) 'preset :: read :: sim parameters'
      call preset_sim_parameters
      call read_sim_parameters
    !--- ---!

    !--- ---!
      write(*,*) 'preset :: read :: mesh_parameters'
      call preset_mesh_parameters
      call read_mesh_parameters
		!--- ---!

    !--- ---!
      write(*,*) 'preset :: read :: background plasma profile'
      call preset_background_plasma_profile
      call read_background_plasma_profile
		!--- ---!

		!--- ---!
      write(*,*) 'preset :: read :: twiss nml'
      call preset_twiss_nml
      call read_twiss_nml
		!--- ---!

    !--- ---!
      write(*,*) 'preset :: read :: poloidal magnetic field values'
      call preset_Bpoloidal
      call read_Bpoloidal_nml
		!--- ---!

    !--- ---!
      write(*,*) 'preset :: read :: ionisation module'
      call preset_ionisation
      call read_ionisation_nml
		!--- ---!

    !--- ---!
		write(*,*) 'preset :: read :: bunch initilisation and self-consistent fields initialisation technique'
		call preset_bunch_initialization
		call read_nml_bunch_initialization
    !--- ---!

    !--- ---!
		write(*,*) 'preset :: read :: operating system'
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



		! !---parameters control for Laplacian matrix---!
    ! !--- not clear procedure :: to be changed
		! bunch_initialization%init_width_r	=	int(min(mesh_par%Rmax,1.*bunch_initialization%init_width_r))
    ! bunch_initialization%init_width_z	=	min(17,bunch_initialization%init_width_z)

		!--- ---!



        !---simulation and plasma parameters---!
        sim_parameters%Nbunches	= bunch_initialization%n_total_bunches 	! number of bunches
        plasma%lambda_p			= sqrt(1.11491e21/(plasma%n0))			! Plasma wavelength (um)
        plasma%k_p				= 2.*pi/plasma%lambda_p					! Plasma wavenumber (rad/um)
        plasma%omega_p 			= plasma%k_p*c   						! Plasma pulsation (rad/fs)

        ! number of particles
        sim_parameters%Np=sum(bunch_initialization%n_particles(1:bunch_initialization%n_total_bunches))

        ! !-from CFL to Dt in (fs)-!
        ! Delta = min( 2.*bunch_initialization%bunch_s_x(1)*1e-6/real(mesh_par%Nsample_r) , &
        !                 bunch_initialization%bunch_s_z(1)*1e-6/real(mesh_par%Nsample_z)   )
        ! sim_parameters%dt = sim_parameters%CFL * Delta / 3e8 / 1e-15
        !
        write(*,'(A)') 'input file read'
        ! write(*,*)
        ! write(*,'(A,f6.3)') 'CFL     >',sim_parameters%CFL
        ! write(*,'(A,f6.3)') 'CFL (fs)>',sim_parameters%dt
        ! write(*,'(A,f6.3)') 'CFL (um)>',sim_parameters%dt*c
        ! write(*,*)
end subroutine read_input





subroutine generate_output_tree
  character	:: command*200

  select case (OSys%macwin)

    case (0)
      sim_parameters%out_dir   = 'out/'
	    sim_parameters%out_root  = ''
      sim_parameters%path_PS 	 = 'out/PS/'
	    sim_parameters%path_grid = 'out/2D/'
      sim_parameters%path_integrated_diagnostics = 'out/integrated_diagnostics/'
      sim_parameters%path_dumprestart = 'dumprestart'
    	call system('mkdir -p out/PS')
    	call system('mkdir -p out/2D')
      call system('mkdir -p out/integrated_diagnostics')
      call system('mkdir -p dumprestart')

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
!--- --- --- --- --- --- ---!


!--- --- --- --- --- --- ---!
!--- PLASMA PARAMETERS
!--- --- --- --- --- --- ---!
subroutine preset_plasma_parameters
  plasma%l_acc=3d4
  plasma%n0=1d16
end subroutine preset_plasma_parameters


subroutine read_plasma_parameters
  NAMELIST/plasma_parameters/plasma
  open(iounit,file='architect.nml',status='old')
  READ(iounit,NML=plasma_parameters,iostat=ierr)
  error_message='plasma parameters'
  close(iounit)
  if(ierr/=0) call print_at_screen_nml_error
end subroutine read_plasma_parameters
!--- --- --- --- --- --- ---!


!--- --- --- --- --- --- ---!
!--- SIM PARAMETERS
!--- --- --- --- --- --- ---!
subroutine preset_sim_parameters
  sim_parameters%CFL 	= .9
  sim_parameters%Output_format= 1 ! Binary PS and grid output
  sim_parameters%jump_grid	= 1 ! no jump of grid points in output
  sim_parameters%reduced_PS 	= 1
  sim_parameters%jump_PS		= 1 ! no jump of particles in PSoutput

  !--- LineOut ---!
  sim_parameters%L_lineout = .FALSE.

  !--- Shapiro Wilks Test---!
  sim_parameters%L_SW_test=.FALSE.
  sim_parameters%SW_sample_dimension=5000
  sim_parameters%SW_sub_iter=10

  !--- bunch field treatment ---!
  sim_parameters%L_BunchREinit=.false.
  sim_parameters%L_plasma_evolution=.true.
  sim_parameters%L_Bunch_evolve=.false.
end subroutine preset_sim_parameters


subroutine read_sim_parameters
    NAMELIST/simulation_parameters/sim_parameters
    open(iounit,file='architect.nml',status='old')
    READ(iounit,NML=simulation_parameters,iostat=ierr)
    error_message='sim_parameters'
    close(iounit)
    if(ierr/=0) call print_at_screen_nml_error
end subroutine read_sim_parameters
!--- --- --- --- --- --- ---!



!--- --- --- --- --- --- ---!
!--- MESH
!--- --- --- --- --- --- ---!
subroutine preset_mesh_parameters
  mesh_par%dxm=-1.D0           !transverse mesh size        in [um] -> converted in kp
  mesh_par%dzm=-1.D0           !longitudinal mesh size      in [um] -> converted in kp
  mesh_par%Nsample_z=-1        !longitudinal number of points per sigma_z
  mesh_par%Nsample_r=-1        !transverse number of points per 2*sigma_r
  mesh_par%R_mesh=-1.D0        !transverse domain dimension in [um] -> converted in kp
  mesh_par%R_mesh_plasma=-1.D0 !transverse plasma extension in [um] -> converted in kp
  mesh_par%Left_mesh=-1.D0     !longitudinal Left  plasma extension in [um] -> converted in kp
  mesh_par%Right_mesh=-1.D0    !longitudinal Right plasma extension in [um] -> converted in kp
end subroutine preset_mesh_parameters

subroutine read_mesh_parameters
    NAMELIST / mesh_parameters / mesh_par
    open(iounit,file='architect.nml',status='old')
    READ(iounit,NML=mesh_parameters,iostat=ierr)
    error_message='mesh parameters'
    if(ierr/=0) call print_at_screen_nml_error
    close(iounit)
end subroutine read_mesh_parameters

!--- --- --- --- --- --- ---!
!--- Ramp
!--- --- --- --- --- --- ---!
subroutine preset_background_plasma_profile
  bck_plasma%order_logitudinal=0 !0: constant
                                 !1: linear ramp
  bck_plasma%order_radial=0      !0: constant
  bck_plasma%n_over_n0=1.d0      !n/n0
  bck_plasma%z_coordinate_um(1)=0
  bck_plasma%z_coordinate_um(2)=plasma%l_acc
  bck_plasma%perturbation_amplitude=0.d0 !amplitude perturbation background density
  !bck_plasma%z_coordinate(8)
end subroutine preset_background_plasma_profile

subroutine read_background_plasma_profile
    NAMELIST / BACKGROUND_PLASMA_PROFILE / bck_plasma
    open(iounit,file='architect.nml',status='old')
    READ(iounit,NML=BACKGROUND_PLASMA_PROFILE,iostat=ierr)
    error_message='background plasma profile'
    close(iounit)
    if(ierr/=0) call print_at_screen_nml_error
end subroutine read_background_plasma_profile


!TWISS
subroutine preset_twiss_nml
    twiss%L_TWISS=.false.
    twiss%alpha_new_factor(:)=1.
    twiss%beta_new_factor(:) =1.
end subroutine preset_twiss_nml

subroutine read_twiss_nml
    NAMELIST / twiss_par / twiss
    open(iounit,file='architect.nml',status='old')
    READ(iounit,NML=twiss_par,iostat=ierr)
    error_message='twiss parameters'
    close(iounit)
    if(ierr/=0) call print_at_screen_nml_error
end subroutine read_twiss_nml



!Poloidal magnetic field namelist control!TWISS
subroutine preset_Bpoloidal
	Bpoloidal%L_Bpoloidal          =.false. !Externally imposed magnetic field
	Bpoloidal%L_BfieldfromV        =.false. !self-generated B-field
	!---> obsolete:: Bpoloidal%B_ex_poloidal=0.0D0
  Bpoloidal%Bprofile             =1
                                    ! 1: Linear+Cubic
                                    ! 2: Linear+FlatSaturation
                                    ! 3: Exponential Profile
                                    ! 4: x^a
                                    ! 5: a^-1 + (1-a^-1) * 1/x^a
	Bpoloidal%capillary_radius_um  =1.0D3 !Capillary Radius in [um]
  Bpoloidal%capillary_radius     =1.0D3 !Capillary Radius in dimensionless units > Self-Calculated
	Bpoloidal%background_current_A =0.D0 !Current in Ampere [A]
  Bpoloidal%a_shape=0.0D0
                          ! 1: weightening Linear and Cubic (1)
                          ! 2: Saturation radius (um)
                          ! 3: Saturation radius Exponential Profile (um)
                          ! 4: Exponent (1)
                          ! 5: Exponent and coefficient
  Bpoloidal%z_coordinate_um      =-1.0D10 !longitudinal coordinate for Bpoloidal fields
end subroutine preset_Bpoloidal

subroutine read_Bpoloidal_nml
    NAMELIST / Bpoloidal_par / Bpoloidal
    open(iounit,file='architect.nml',status='old')
    READ(iounit,NML=Bpoloidal_par,iostat=ierr)
    error_message='poloidal magnetic field'
    close(iounit)
    if(ierr/=0) call print_at_screen_nml_error
end subroutine read_Bpoloidal_nml



!--- --- --- --- --- --- ---!
!--- Ionisation
!--- --- --- --- --- --- ---!
subroutine preset_ionisation
  ionisation%L_ionisation=.false.
  ionisation%particle_per_cell=4
  ionisation%atomic_number=18.
  ionisation%mass_number=40.
end subroutine preset_ionisation

subroutine read_ionisation_nml
  NAMELIST / ionisation_parameters / ionisation
  open(iounit,file='architect.nml',status='old')
  READ(iounit,NML=ionisation_parameters,iostat=ierr)
  error_message='ionisation'
  close(iounit)
  if(ierr/=0) call print_at_screen_nml_error
end subroutine read_ionisation_nml
!--- ---!




!OSystem
subroutine preset_OS_nml
    OSys%macwin=0 !0=mac
    OSys%output_dir='.'
end subroutine preset_OS_nml


subroutine read_OS_nml
    NAMELIST / OS / OSys
    open(iounit,file='architect.nml',status='old')
    READ(iounit,NML=OS,iostat=ierr)
    error_message='OS'
    close(iounit)
    if(ierr/=0) call print_at_screen_nml_error
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
    bunch_initialization%shape(:)= 1 !1: gaussian, 2: cylinder
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
    !--- for triangular shape ---!
    bunch_initialization%Charge_right(:)=0.d0
    bunch_initialization%Charge_left(:)=1.d0
end subroutine preset_bunch_initialization


subroutine read_nml_bunch_initialization
    NAMELIST / BUNCHINIT / bunch_initialization
    open(iounit,file='architect.nml',status='old')
    READ(iounit,NML=BUNCHINIT,iostat=ierr)
    error_message='bunch parameters'
    close(iounit)
    if(ierr/=0) call print_at_screen_nml_error
end subroutine read_nml_bunch_initialization



subroutine print_at_screen_nml_error
  character(100) :: line
  backspace(iounit)
  read(iounit,fmt='(A)') line
  write(*,'(A)')    '*** *** *** *** *** *** *** *** *** *** *** ***'
  write(*,'(A)')    'Error in namelist:      '//trim(error_message)
  write(*,'(A)')    'Invalid namelist entry: '//trim(line)
  write(*,'(A,I5)') 'iostat type of error:   ',ierr
  write(*,'(A)')    '*** *** *** *** *** *** *** *** *** *** *** ***'
  stop
end subroutine print_at_screen_nml_error




 !--- --- ---!
 end module read_input_module
 !--- --- ---!
