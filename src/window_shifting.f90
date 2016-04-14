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

MODULE Move_Window_FDTD


USE my_types
USE use_my_types
USE pstruct_data
USE architect_class_structure
USE moments

IMPLICIT NONE

CONTAINS


!Move window FDTD version

	SUBROUTINE init_window
	   ! Initialize Moving window position - center in first driver center of mass
	   sim_parameters%zg  						= 0.
	   sim_parameters%zg_old 					= sim_parameters%zg
	   sim_parameters%window_shifted_cells      = 0

   END SUBROUTINE init_window


   SUBROUTINE control_window
		real :: dz_eff,zg_delta

   		dz_eff=mesh_par%dzm/plasma%k_p
		if (sim_parameters%window_mode.eq.0) then ! window moves with first bunch center
			sim_parameters%zg = calculate_nth_moment_bunch(1,1,3)
		else if (sim_parameters%window_mode.eq.1) then ! window moves with constant speed
			sim_parameters%zg = - sim_parameters%moving_window_speed*c*sim_parameters%sim_time
		endif

		zg_delta              = abs(sim_parameters%zg-sim_parameters%zg_old)

		if ( zg_delta>dz_eff .and. zg_delta<2.*dz_eff ) then
			sim_parameters%zg_old = sim_parameters%zg
			sim_parameters%window_shifted_cells = sim_parameters%window_shifted_cells + 1
			call Move_Window_FDTD_version_COMB

		elseif( zg_delta>2.*dz_eff .and. zg_delta<3.*dz_eff ) then
			sim_parameters%zg_old = sim_parameters%zg
			sim_parameters%window_shifted_cells = sim_parameters%window_shifted_cells + 2
			call Move_Window_FDTD_version_COMB
			write(*,*) 'Error: Moving window moved of 2 cells'
		elseif( zg_delta>=3.*dz_eff ) then
			write(*,*) 'Error: Moving window'
				stop
		endif

   END SUBROUTINE control_window



   SUBROUTINE Move_Window_FDTD_version_COMB
   IMPLICIT NONE
   INTEGER cells_advanced,i,im	,iter,nn,j
   REAL, DIMENSION(mesh_par%Nzm,mesh_par%Nxm) :: Ez,Er,Bphi,Ez_new,Er_new,Bphi_new,Bphi_old,Bphi_old_new
   REAL, DIMENSION(mesh_par%Nzm,mesh_par%Nxm) :: ux,uz,gam,ux_new,uz_new,gam_new,ne,ne_new
	 REAL, DIMENSION(mesh_par%Nzm,mesh_par%Nxm) :: Ez_bunch,Er_bunch,Bphi_bunch,Bphi_old_bunch
	 REAL, DIMENSION(mesh_par%Nzm,mesh_par%Nxm) :: Ez_bunch_new,Er_bunch_new,Bphi_bunch_new,Bphi_old_bunch_new
   INTEGER :: Nz,Nr,Node_min_z,Node_max_z,Node_min_r,Node_max_r,Node_end_z,Node_end_r
   REAL :: DeltaR,DeltaZ, test_ne

	Er         	= 0.
	Ez         	= 0.
	Bphi 	   	  = 0.
	Bphi_old   	= 0

	!--->Er_bunch   	= 0.
	!--->Ez_bunch   	= 0.
	!--->Bphi_bunch  = 0.
	!--->Bphi_old_bunch 	= 0

	Er_new     	= 0.
	Ez_new     	= 0.
	Bphi_new   	= 0.
	Bphi_old_new   	= 0.

	!--->Er_bunch_new     	= 0.
	!--->Ez_bunch_new     	= 0.
	!--->Bphi_bunch_new   	= 0.
	!--->	Bphi_old_bunch_new   	= 0.

	ux  	   	= 0.
	uz  	   	= 0.
	ne  		= 0.

	ux_new 	   	= 0.
	uz_new 	   	= 0.
	ne_new		= 0.

	Nz         	= mesh_par%Nzm
	Nr         	= mesh_par%Nxm

	Node_min_z 	= 2
	Node_max_z 	= Nz-1

	Node_min_r 	= 2
	Node_max_r 	= Nr-1

	Node_end_z 	= Nz
	Node_end_r 	= Nr

	DeltaR     	= mesh_par%dxm
	DeltaZ     	= mesh_par%dzm

	!------------------------------------------------------------------------!
	!                           Z axis
	!
	!    ghost cell        physical domain         ghost cell
	!  |______________|__________________________|______________|
	!  |              |                          |              |
	!  1         Node_min_z                 Node_max_z       Node_end_z
	!
	!------------------------------------------------------------------------!

	!------------------------------------------------------------------------!
	!                           R axis
	!
	!      Axis    physical domain                   ghost cell
	!  |________|________________________________|______________|
	!  |        |                                |              |
	!  1     Node_min_r=2                    Node_max_r       Node_end_r
	!
	!  Bphi at j=1 is physical (it is on the r axis), Ez at j=1 is at -DeltaR/2
	!  (ghost space); Ez at j=Node_end_r is physical, Ez at j=j=Node_end_r is
	!  in ghost space
	!
	!------------------------------------------------------------------------!


	!----------------------------------------------!
	!  Initial conditions from old moving window   !
	!----------------------------------------------!

	Bphi    (1:Node_end_z,1:Node_end_r)  =	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%Bphi
	Ez      (1:Node_end_z,1:Node_end_r)  =	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%Ez
	Er      (1:Node_end_z,1:Node_end_r)  =	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%Ex
	Bphi_old(1:Node_end_z,1:Node_end_r)  =	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%Bphi_old

	!--->Bphi_bunch    (1:Node_end_z,1:Node_end_r)  =	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%Bphi_bunch
	!--->Ez_bunch      (1:Node_end_z,1:Node_end_r)  =	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%Ez_bunch
	!--->Er_bunch      (1:Node_end_z,1:Node_end_r)  =	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%Ex_bunch
	!--->Bphi_old_bunch(1:Node_end_z,1:Node_end_r)  =	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%Bphi_old_bunch

	uz  (1:Node_end_z,1:Node_end_r)	     =	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%uz
	ux  (1:Node_end_z,1:Node_end_r)      =	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%ux
	ne  (1:Node_end_z,1:Node_end_r)	     =	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%n_plasma_e

	!----------------------------------------------!
	!       Longitudinal shift by one cell         !
	!----------------------------------------------!

	do i=2,Node_max_z
		Er_new      (i+1,:) = Er      (i,:)
		Ez_new      (i+1,:) = Ez      (i,:)
		Bphi_new    (i+1,:) = Bphi    (i,:)
		Bphi_old_new(i+1,:) = Bphi_old(i,:)

		!--->Er_bunch_new(i+1,:) 			= Er_bunch(i,:)
		!--->Ez_bunch_new(i+1,:) 			= Ez_bunch(i,:)
		!--->Bphi_bunch_new(i+1,:) 		= Bphi_bunch(i,:)
		!--->Bphi_old_bunch_new(i+1,:) = Bphi_old_bunch(i,:)

		ux_new      (i+1,:) = ux      (i,:)
		uz_new      (i+1,:) = uz      (i,:)
		ne_new	    (i+1,:) = ne      (i,:)
		ne_new      (i+1,mesh_par%NRmax_plasma:Node_end_r) = 0.D0
	enddo

	Er_new      (1:2,:) = 0.
	Ez_new      (1:2,:) = 0.
	Bphi_new    (1:2,:) = 0.
  Bphi_old_new(1:2,:) = 0.

	!--->Er_bunch_new	(1:2,:) = 0.
	!--->Ez_bunch_new	(1:2,:) = 0.
	!--->Bphi_bunch_new(1:2,:) = 0.
	!--->Bphi_old_bunch_new(1:2,:) = 0.


	ux_new      (1:2,:) = 0.
	uz_new      (1:2,:) = 0.


	if (sim_parameters%I_start_ramp_peak_position.gt.2) then ! whole ramp initially inside window

		if (abs(sim_parameters%zg).ge.sim_parameters%start_exit_ramp) then ! window has reached the beginning of the exit ramp

			if (sim_parameters%ramps_order.eq.1) then ! linear ramp

				if (sim_parameters%end_ramp_length.gt.0.) then ! there is a finite length exit ramp
					test_ne = 1.- 1.*(abs(sim_parameters%zg)-sim_parameters%start_exit_ramp )/(sim_parameters%end_ramp_length*plasma%lambda_p)
				else ! there is a zero length exit ramp
					test_ne = 0.
				endif

				if (test_ne .le. 0.) then
					test_ne = 0.
				endif

			else if (sim_parameters%ramps_order.eq.2) then ! parabolic ramp

			endif

			ne_new      (1:2,:) = test_ne

		else

			if (sim_parameters%order_capillary_density_z.eq.0) then ! uniform longitudinal density profile
					ne_new      (1:2,:) = 1.
			else if (sim_parameters%order_capillary_density_z.eq.2) then ! parabolic longitudinal density profile

			endif

		endif


	else ! ramp initially in part inside window, in part outside window


		if (sim_parameters%window_shifted_cells.le.(sim_parameters%I_start_ramp_length_out_window+2)) then
			! window boundary still in ramp
			ne_new      (1:2,:) = sim_parameters%start_n_peak_in_window + &
									(1.-sim_parameters%start_n_peak_in_window) &
									*(1.*(sim_parameters%window_shifted_cells-2)/(1.*sim_parameters%I_start_ramp_length_out_window) )
		else

		if (abs(sim_parameters%zg).ge.sim_parameters%start_exit_ramp) then  ! window has reached the beginning of the exit ramp

			if (sim_parameters%ramps_order.eq.1) then ! linear ramp


				if (sim_parameters%end_ramp_length.gt.0.) then
					test_ne = 1.- 1.*(abs(sim_parameters%zg)-sim_parameters%start_exit_ramp )/(sim_parameters%end_ramp_length*plasma%lambda_p)
				else
					test_ne = 0.
				endif

				if (test_ne .le. 0.) then
					test_ne = 0.
				endif

				else if (sim_parameters%ramps_order.eq.2) then ! parabolic ramp

				endif

				ne_new      (1:2,:) = test_ne
			else

				if (sim_parameters%order_capillary_density_z.eq.0) then ! uniform longitudinal density profile
					ne_new      (1:2,:) = 1.
				else if (sim_parameters%order_capillary_density_z.eq.2) then ! parabolic longitudinal density profile

				endif
			endif
		endif
	endif


	! multiplication factor for radial profile factor
	if (sim_parameters%order_capillary_density_r.eq.2) then ! parabolic radial density profile
		do j= 2,Node_max_r
			ne_new(1:2,j) 	   = ne_new(1:2,j)*radial_factor(j)
		enddo

		! BC
		ne_new(1:2,Node_end_r) = ne_new(1:2,Node_max_r) ! upper boundary
		ne_new(1:2,1         ) = ne_new(1:2,2         ) ! lower boundary
	endif

	!----------------------------!
	! Substitution of new fields !
	!----------------------------!
	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%Bphi_old   =      Bphi_old_new(1:Node_end_z,1:Node_end_r)
	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%Bphi       =      Bphi_new    (1:Node_end_z,1:Node_end_r)
	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%Ez         =      Ez_new      (1:Node_end_z,1:Node_end_r)
	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%Ex         =      Er_new      (1:Node_end_z,1:Node_end_r)

	!--->mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%Ex_bunch				= Er_bunch_new			(1:Node_end_z,1:Node_end_r)
	!--->mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%Ez_bunch				= Ez_bunch_new			(1:Node_end_z,1:Node_end_r)
	!--->mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%Bphi_bunch			= Bphi_bunch_new		(1:Node_end_z,1:Node_end_r)
	!--->mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%Bphi_old_bunch	= Bphi_old_bunch_new(1:Node_end_z,1:Node_end_r)

	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%uz         =      uz_new      (1:Node_end_z,1:Node_end_r)
	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%ux         =      ux_new      (1:Node_end_z,1:Node_end_r)
	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%n_plasma_e =      ne_new      (1:Node_end_z,1:Node_end_r)

   return

   END SUBROUTINE


END MODULE




! DO NOT DELETE UNTIL PARABOLIC Z PROFILE IS IMPLEMENTED!!!
!~if (sim_parameters%longitudinal_density_profile.eq.0) then ! linear ramp at capillary entrance
!~		if (sim_parameters%I_ramp_peak_position.gt.2) then
!~			!! whole ramp initially inside window
!~			!ne_new      (1:2,:) = 1.
!~
!~			if (abs(sim_parameters%zg).ge.sim_parameters%start_exit_ramp) then
!~
!~				! window has reached the beginning of the exit ramp
!~				if (sim_parameters%L_exit_ramp.gt.0.) then
!~					test_ne = 1.- 1.*(abs(sim_parameters%zg)-sim_parameters%start_exit_ramp )/(sim_parameters%L_exit_ramp*plasma%lambda_p)
!~				else
!~					test_ne = 0.
!~				endif
!~
!~				if (test_ne .le. 0.) then
!~					test_ne = 0.
!~				endif
!~				ne_new      (1:2,:) = test_ne
!~
!~			else
!~				ne_new      (1:2,:) = 1.
!~			endif
!~
!~
!~		else ! ramp initially in part inside window, in part outside window
!~
!~			!~if (sim_parameters%window_shifted_cells.le.sim_parameters%I_ramp_length_out_window) then
!~			!~
!~			!~	! window boundary still in ramp
!~			!~	ne_new      (1:2,:) = sim_parameters%n_peak_in_window + &
!~			!~							(1.-sim_parameters%n_peak_in_window) &
!~			!~		   *(1.*sim_parameters%window_shifted_cells/(1.*sim_parameters%I_ramp_length_out_window) )
!~			!~
!~			!~else
!~			!~	! window boundary in the uniform density zone
!~			!~	ne_new      (1:2,:) = 1.
!~			!~endif
!~
!~
!~			if (sim_parameters%window_shifted_cells.le.sim_parameters%I_ramp_length_out_window) then
!~				! window boundary still in ramp
!~				ne_new      (1:2,:) = sim_parameters%n_peak_in_window + &
!~										(1.-sim_parameters%n_peak_in_window) &
!~										*(1.*sim_parameters%window_shifted_cells/(1.*sim_parameters%I_ramp_length_out_window) )
!~
!~			else
!~
!~				if (abs(sim_parameters%zg).ge.sim_parameters%start_exit_ramp) then
!~					! window has reached the beginning of the exit ramp
!~
!~					if (sim_parameters%L_exit_ramp.gt.0.) then
!~						test_ne = 1.- 1.*(abs(sim_parameters%zg)-sim_parameters%start_exit_ramp )/(sim_parameters%L_exit_ramp*plasma%lambda_p)
!~					else
!~						test_ne = 0.
!~					endif
!~
!~					if (test_ne .le. 0.) then
!~						test_ne = 0.
!~					endif
!~					ne_new      (1:2,:) = test_ne
!~				else
!~					ne_new      (1:2,:) = 1.
!~				endif
!~
!~			endif
!~
!~
!~		endif
!~
!~	else if (sim_parameters%longitudinal_density_profile.eq.1) then ! longitudinal parabola
!~
!~	    test_ne = - sim_parameters%I_parabola_normalization_factor* &
!~				 (sim_parameters%window_shifted_cells  + sim_parameters%I_origin_z_axis - 1 - sim_parameters%I_parabola_start) &
!~				*(sim_parameters%window_shifted_cells  + sim_parameters%I_origin_z_axis - 1 - sim_parameters%I_parabola_end)
!~
!~		if (test_ne .le. 0.) then
!~			test_ne = 0.
!~		endif
!~
!~		ne_new (1:2,1:Node_end_r)  = test_ne



	! DO NOT DELETE UNTIL PARABOLIC Z PROFILE IS IMPLEMENTED!!! RAMP AT CAPILLARY END MUST BE INSERTED

	 !~   test_ne = - sim_parameters%I_parabola_normalization_factor* &
		!~		 (sim_parameters%window_shifted_cells  + sim_parameters%I_origin_z_axis - 1 - sim_parameters%I_parabola_start) &
		!~		*(sim_parameters%window_shifted_cells  + sim_parameters%I_origin_z_axis - 1 - sim_parameters%I_parabola_end)
		!~
		!~if (test_ne .le. 0.) then
		!~	test_ne = 0.
		!~endif
		!~
		!~ne_new (1:2,1:Node_end_r)  = test_ne
