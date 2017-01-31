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

MODULE my_types

IMPLICIT NONE



   TYPE :: simul_param	 ! sim_parameters

    REAL(8) :: sim_time,dt
    INTEGER :: Np,iter,NelectronsB(6),Np_cut,Npc,Npcz,Nbunches
    INTEGER :: jfirst_cut,j_cut,diff_order,deposition,Ndr
    REAL(8) :: Charge,CFL,scalesx,scalesz
    REAL(8) :: r0,lbunch(6),zg,zg_old,z0_first_driver,output_distance,ramp_length
    REAL(8) :: xb0(6),yb0(6),zb0(6),rB0(6)
    REAL(8) :: sigma_cut
    INTEGER :: FDTD_version
    INTEGER :: Fluid_Scheme
    INTEGER :: Output_format,reduced_PS
    INTEGER :: window_shifted_cells,I_distance_capillary
	INTEGER :: window_mode, diagnostics_with_dcut
	REAL(8) :: moving_window_speed
	real(8) :: velocity_background

    REAL(8)    :: start_n_peak_in_window,start_ramp_length_in_window,start_ramp_peak_position
	  INTEGER :: I_start_ramp_length_out_window,I_start_ramp_length_in_window,I_start_ramp_length,I_start_ramp_peak_position

	  INTEGER :: I_parabola_length_in_window,I_parabola_start,I_parabola_end,I_origin_z_axis
	  REAL(8)    :: I_parabola_normalization_factor

	  CHARACTER :: path_PS*255,path_grid*255,path_dumprestart*255,path_integrated_diagnostics*255
    CHARACTER :: out_dir*100,path_1D*100,out_root*100

	  INTEGER :: init_width_r,init_width_z
	  INTEGER :: dim_Laplacian
	  INTEGER :: jump_grid=1,jump_PS=1,sys=1

	  INTEGER :: I_exit_ramp_length,I_start_exit_ramp,I_end_exit_ramp
	  REAL(8) :: start_exit_ramp,end_exit_ramp,L_exit_ramp

	  INTEGER :: output_Integrated_params_nstep, output_grid_nstep, output_PS_nstep
	  REAL(8) :: output_grid_dist=99999999.,gridDeltaOutput=0.0,gridLastOutput=0.0
	  REAL(8) :: output_PS_dist=99999999.,PSDeltaOutput=0.0,PSLastOutput=0.0
	  REAL(8) :: output_Integrated_params_dist=99999999.,IntDeltaOutput=0.0,IntLastOutput=0.0

    LOGICAL :: L_SW_test,L_lineout
    INTEGER :: SW_sample_dimension,SW_sub_iter

	  REAL(8) :: lastInWindowCheck=0.,InWindoWCheckDelta=25.0

    LOGICAL :: L_BunchREinit=.false.,L_plasma_evolution=.true.,L_Bunch_evolve
    REAL(8) :: lastBunchREinit=0.,bunch_reinit_distance_um=1000000.
   END TYPE

   TYPE :: mesh_param   ! mesh_par
      INTEGER :: NRmax_plasma
      INTEGER :: Nxm,Nzm,Nsample_r,Nsample_z,Nxm_old,Nzm_old
      REAL(8) :: Left_Domain_boundary,Right_Domain_boundary
      REAL(8) :: ScaleX,ScaleZ,ScaleR,L_mesh,R_mesh,R_mesh_plasma,dzm,dxm,dzm_old,dxm_old
      REAL(8) :: z_min,z_max,z_min_um,z_max_um
      REAL(8) :: z_min_moving,z_max_moving,z_min_moving_um,z_max_moving_um
      REAL(8) :: Left_mesh,Right_mesh !input-parameters
      REAL(8) :: Rmax,Rmax_plasma
   END TYPE

   !--- define the background plasma profile ---!
   TYPE :: background_plasma_profile   ! bck_plasma
      INTEGER :: order_logitudinal(8),order_radial(8)
      REAL(8) :: radius_um(8),radius_internal_um(8),perturbation_amplitude(8)
      REAL(8) :: z_coordinate_um(8),z_coordinate(8),n_over_n0(8)
   END TYPE

   TYPE :: mesh_phys   ! mesh
      REAL(8) :: n_plasma_e,n_plasma_i,Zstar
      REAL(8) :: rho,n0,Jz,Jr,ux,uz,Ez,Ex,Bphi,Bphi_old,Jpe_r,Jpe_z
      REAL(8) :: Ez_bunch,Ex_bunch,Bphi_bunch,Bphi_old_bunch
      REAL(8) :: B_ex_poloidal
   END TYPE

    TYPE :: mesh_utility
      real(8), DIMENSION(:), ALLOCATABLE :: Bphi_BC_Left
    END TYPE


   TYPE :: plasma_param	 ! plasma
      REAL(8) :: n0
      REAL(8) :: l_acc
	  !~REAL :: lambda,omega_0,k_0,beta_g1,Ldeph1
	  REAL(8) :: lambda_p,omega_p,k_p,epsilon,beta_g,gamma_g,omega
   END TYPE

  !ionisation
  TYPE :: ionisation_parameters
    LOGICAL :: L_ionisation
    INTEGER :: particle_per_cell,np_per_r,np_per_z
    INTEGER :: tot_ions
    REAL(8) :: atomic_number,mass_number
  END TYPE

	!input twiss parameters
	TYPE :: twiss_param
		LOGICAL :: L_TWISS(7)
		REAL(8) :: alpha_new_factor(7),beta_new_factor(7)
	END TYPE

  !input Bpoloidal parameters
  !capillary_radius    is the DIMENSIONLESS capillary radius
  !capillary_radius_um is the DIMENSIONAL   capillary radius in microns
	TYPE :: Bpoloidal_param
		LOGICAL :: L_Bpoloidal,L_BfieldfromV
    INTEGER :: Bprofile(5)
		REAL(8) :: capillary_radius,capillary_radius_um,background_current_A(5)
    REAL(8) :: a_shape(5),z_coordinate_um(6)
	END TYPE

	!Operative System Switch
	TYPE :: OSys_param
		INTEGER :: macwin
		CHARACTER*50 :: output_dir
	END TYPE

	!Bunch-initialisation-parameters
	TYPE :: bunch_inside_initialization
		LOGICAL :: l_bunch_internal_init
		INTEGER :: n_total_bunches=1,n_particles(7)=0,self_consistent_field_bunch
		REAL(8)    :: bunch_s_x(7)=0.,bunch_s_y(7)=0.,bunch_s_z(7)=0.
		REAL(8)    :: bunch_gamma_m(7)=0.,bunch_eps_x(7)=0.,bunch_eps_y(7)=0.,bunch_dgamma(7)=0.
    real(8)    :: Charge_right(7),Charge_left(7)
		REAL(8)    :: maxnorm,wsor
		INTEGER    :: init_width_r, init_width_z, iter_max, shape(7)
		REAL(8)    :: ChargeB(7)=0., db(7)=0.
    CHARACTER  :: inbunch(6)*30
	END TYPE

  !--- dump and restart variables ---!
	TYPE :: dump_and_restart
		LOGICAL :: L_onoff,restart
    INTEGER :: nstep=9999
    REAL(8) :: distance_um=99999999.,LastOutput_um
	END TYPE


	TYPE(mesh_phys), DIMENSION(:,:), ALLOCATABLE :: mesh
	REAL(8) Zmin,Xmin,Zmin_shifted,Xmin_shifted
	INTEGER :: jcount,jfc,jcount1
	INTEGER, DIMENSION(:), ALLOCATABLE :: cut_num
	REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: savedata_dr
	REAL(8), DIMENSION(:), ALLOCATABLE :: x_mesh,z_mesh,x_mesh_shifted,z_mesh_shifted,inv_R,inv_R_shifted
	REAL(8), DIMENSION(:), ALLOCATABLE :: r_mesh,Sr,Sz,Vol
	INTEGER :: Nz,Nr,Node_min_z,Node_max_z,Node_min_r,Node_max_r,Node_end_z,Node_end_r
	INTEGER :: Node_min_lo_z,Node_max_lo_z,Node_min_lo_r,Node_max_lo_r,Node_end_lo_z,Node_end_lo_r
	INTEGER :: Node_min_ho_z,Node_max_ho_z,Node_min_ho_r,Node_max_ho_r,Node_end_ho_z,Node_end_ho_r
	REAL(8) :: DeltaR,DeltaZ,Dt,one_over_dx,one_over_dz



	REAL(8), DIMENSION(5,5) :: Db1,Db2
	DOUBLE PRECISION c,alpha,pi,c_SI
  DOUBLE PRECISION electron_mass,electron_charge
  DOUBLE PRECISION mu0
	PARAMETER (c=0.299792458D0,alpha=0.00729927d0,pi=3.141592653589793d0,c_SI=299792458D0)
  PARAMETER (electron_mass=9.10938291e-31, electron_charge=1.602176565e-19)
  PARAMETER (mu0=1.25663706D-6)
	CHARACTER :: path*90,command*200

END MODULE my_types
