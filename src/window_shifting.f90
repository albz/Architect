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
USE Compute_plasma_current
USE ion_background

IMPLICIT NONE

CONTAINS


!Move window FDTD version

	SUBROUTINE init_window
	   ! Initialize Moving window position - center in first driver center of mass
	   sim_parameters%zg  						= 0.D0
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

		zg_delta = abs(sim_parameters%zg-sim_parameters%zg_old)

		!--- begin to shift---!
		if ( zg_delta>dz_eff .and. zg_delta<2.*dz_eff ) then
			sim_parameters%zg_old = sim_parameters%zg
			sim_parameters%window_shifted_cells = sim_parameters%window_shifted_cells + 1

			mesh_par%z_min_moving_um = sim_parameters%zg + mesh_par%z_min/plasma%k_p
			mesh_par%z_max_moving_um = sim_parameters%zg + mesh_par%z_max/plasma%k_p
			mesh_par%z_min_moving = mesh_par%z_min_moving_um*plasma%k_p
			mesh_par%z_max_moving = mesh_par%z_max_moving_um*plasma%k_p

			call Move_Window_FDTD_version_COMB
			if(Bpoloidal%L_Bpoloidal)   call Move_Window_Bexternal_field
			if(ionisation%L_ionisation) call Move_Windows_ion_background

		elseif( zg_delta>2.*dz_eff .and. zg_delta<3.*dz_eff ) then
			write(*,'(A)') 'two steps moving window - stop'
			stop
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
   REAL, DIMENSION(mesh_par%Nzm,mesh_par%Nxm) :: ux,uz,gam,ux_new,uz_new,gam_new,ne,ne_new,ni,ni_new
	 REAL, DIMENSION(mesh_par%Nzm,mesh_par%Nxm) :: Ez_bunch,Er_bunch,Bphi_bunch,Bphi_old_bunch
	 REAL, DIMENSION(mesh_par%Nzm,mesh_par%Nxm) :: Ez_bunch_new,Er_bunch_new,Bphi_bunch_new,Bphi_old_bunch_new
   INTEGER :: Nz,Nr,Node_min_z,Node_max_z,Node_min_r,Node_max_r,Node_end_z,Node_end_r
   REAL :: DeltaR,DeltaZ, test_ne

	Er         	= 0.D0
	Ez         	= 0.D0
	Bphi 	   	  = 0.D0
	Bphi_old   	= 0.D0

	Er_bunch   	= 0.
	Ez_bunch   	= 0.
	Bphi_bunch  = 0.
	Bphi_old_bunch 	= 0

	Er_new     	= 0.D0
	Ez_new     	= 0.D0
	Bphi_new   	= 0.D0
	Bphi_old_new   	= 0.D0

	Er_bunch_new     	= 0.
	Ez_bunch_new     	= 0.
	Bphi_bunch_new   	= 0.
	Bphi_old_bunch_new   	= 0.

	ux  	   	= 0.D0
	uz  	   	= 0.D0
	ne  		  = 0.D0
	ni  		  = 0.D0

	ux_new 	   	= 0.D0
	uz_new 	   	= 0.D0
	ne_new		  = 0.D0
	ni_new		  = 0.D0

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

	if(sim_parameters%L_Bunch_evolve) then
		Bphi_bunch    (1:Node_end_z,1:Node_end_r)  =	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%Bphi_bunch
		Ez_bunch      (1:Node_end_z,1:Node_end_r)  =	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%Ez_bunch
		Er_bunch      (1:Node_end_z,1:Node_end_r)  =	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%Ex_bunch
		Bphi_old_bunch(1:Node_end_z,1:Node_end_r)  =	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%Bphi_old_bunch
endif

	uz  (1:Node_end_z,1:Node_end_r)	     =	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%uz
	ux  (1:Node_end_z,1:Node_end_r)      =	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%ux
	ne  (1:Node_end_z,1:Node_end_r)	     =	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%n_plasma_e
	ni  (1:Node_end_z,1:Node_end_r)	     =	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%n_plasma_i

	!----------------------------------------------!
	!       Longitudinal shift by one cell         !
	!----------------------------------------------!

	do i=2,Node_max_z
		Er_new      (i+1,:) = Er      (i,:)
		Ez_new      (i+1,:) = Ez      (i,:)
		Bphi_new    (i+1,:) = Bphi    (i,:)
		Bphi_old_new(i+1,:) = Bphi_old(i,:)

		if(sim_parameters%L_Bunch_evolve) then
			Er_bunch_new(i+1,:) 			= Er_bunch(i,:)
			Ez_bunch_new(i+1,:) 			= Ez_bunch(i,:)
			Bphi_bunch_new(i+1,:) 		= Bphi_bunch(i,:)
			Bphi_old_bunch_new(i+1,:) = Bphi_old_bunch(i,:)
		endif

		ux_new      (i+1,:) = ux      (i,:)
		uz_new      (i+1,:) = uz      (i,:)
		ne_new	    (i+1,:) = ne      (i,:)
		ni_new	    (i+1,:) = ni      (i,:)
		!---ne_new      (i+1,mesh_par%NRmax_plasma:Node_end_r) = 0.D0
	enddo

	Er_new      (1:2,:) = 0.D0
	Ez_new      (1:2,:) = 0.D0
	Bphi_new    (1:2,:) = 0.D0
	Bphi_old_new(1:2,:) = 0.D0
	ux_new      (1:2,:) = 0.D0
	uz_new      (1:2,:) = 0.D0

	if(sim_parameters%L_Bunch_evolve) then
		Er_bunch_new	(1:2,:) = 0.
		Ez_bunch_new	(1:2,:) = 0.
		Bphi_bunch_new(1:2,:) = 0.
		Bphi_old_bunch_new(1:2,:) = 0.
	endif


	!---background density---!
 	do j= 2,Node_max_r
 	 		ne_new(2,j)  = background_density_value(2,j)
			if(.not.ionisation%L_ionisation) ni_new(2,j)=ne_new(2,j)
 	enddo
 	!BC
 	ne_new(2,Node_end_r) = ne_new(2,Node_max_r) ! upper boundary
 	ne_new(2,1         ) = ne_new(2,2         ) ! lower boundary
 	ne_new(1,:         ) = ne_new(2,:         ) ! left  boundary
	if(.not.ionisation%L_ionisation) then
	 	ni_new(2,1         ) = ni_new(2,2         ) ! lower boundary
	 	ni_new(1,:         ) = ni_new(2,:         ) ! left  boundary
	endif
 	!--- ---!


	!--- Moving Background ---!
	if(Bpoloidal%L_BfieldfromV) then
		uz_new (1:2,1:mesh_par%NRmax_plasma-1) = sim_parameters%velocity_background
		Bphi_new(1,:) = mesh_util%Bphi_BC_Left(:)
		Bphi_old_new(1,:) = mesh_util%Bphi_BC_Left(:)
	end if


	!----------------------------!
	! Substitution of new fields !
	!----------------------------!
	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%Bphi_old   =      Bphi_old_new(1:Node_end_z,1:Node_end_r)
	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%Bphi       =      Bphi_new    (1:Node_end_z,1:Node_end_r)
	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%Ez         =      Ez_new      (1:Node_end_z,1:Node_end_r)
	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%Ex         =      Er_new      (1:Node_end_z,1:Node_end_r)

	if(sim_parameters%L_Bunch_evolve) then
		mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%Ex_bunch				    = Er_bunch_new			(1:Node_end_z,1:Node_end_r)
		mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%Ez_bunch				    = Ez_bunch_new			(1:Node_end_z,1:Node_end_r)
		mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%Bphi_bunch			   = Bphi_bunch_new		(1:Node_end_z,1:Node_end_r)
		mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%Bphi_old_bunch	= Bphi_old_bunch_new(1:Node_end_z,1:Node_end_r)
	endif

	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%uz         =      uz_new      (1:Node_end_z,1:Node_end_r)
	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%ux         =      ux_new      (1:Node_end_z,1:Node_end_r)
	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%n_plasma_e =      ne_new      (1:Node_end_z,1:Node_end_r)
	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%n_plasma_i =      ni_new      (1:Node_end_z,1:Node_end_r)

   return
   END SUBROUTINE


	SUBROUTINE Move_Window_Bexternal_field
		REAL(8) :: Radius,a,Bfield
		INTEGER :: i,j,Node_max_z

		!--- Rigid shifting---!
		Node_max_z 	= mesh_par%Nzm-1
		do i=Node_max_z,2,-1
			mesh(i+1,:)%B_ex_poloidal = mesh(i,:)%B_ex_poloidal
		enddo

		!--- Add Bfield on the LEFT-Boundary ---!
		Do i=1,5
			if(mesh_par%z_min_moving_um <= Bpoloidal%z_coordinate_um(i) .and. &
			   mesh_par%z_min_moving_um >  Bpoloidal%z_coordinate_um(i+1) ) then
						!--- Bfield ---!
						Bfield   = mu0*Bpoloidal%background_current_A(i)/(2.D0*pi*Bpoloidal%capillary_radius_um*1e-6) !Poloidal field from current
						Bfield   = Bfield / (96.*sqrt(plasma%n0)/3e8) !from Dimensional to DimensionLESS
						do j=2,mesh_par%Nxm

							Select Case (Bpoloidal%Bprofile(i))
								case (1) !---Linear+Cubic
									Radius=r_mesh(j)/Bpoloidal%capillary_radius
				 					a=Bpoloidal%a_shape(i)
				 					mesh(2,j)%B_ex_poloidal = -Bfield*((1.D0-a)*Radius+a*Radius**3)

								case (2) !---Linear+FlatSaturation
									a=Bpoloidal%a_shape(i)*plasma%k_p
									mesh(2,j)%B_ex_poloidal = -Bfield/a*r_mesh(j)
									if(r_mesh(j)>a) then
										mesh(2,j)%B_ex_poloidal = -Bfield
									endif

								case (3) !---Exponential Profile
									a=Bpoloidal%a_shape(i)*plasma%k_p
									mesh(i,j)%B_ex_poloidal = -Bfield*( 1.D0-exp(-r_mesh(j)/a) )/( 1.D0-exp(-r_mesh(mesh_par%Nxm)/a) )

								case(4) !--- x^a
									a=Bpoloidal%a_shape(i)
									Radius=r_mesh(j)/Bpoloidal%capillary_radius
									mesh(i,j)%B_ex_poloidal = -Bfield*( Radius )**a

								case(5) !--- a^-1 + (1-a^-1) * 1/x^a
									a=Bpoloidal%a_shape(i)
									Radius=r_mesh(j)/Bpoloidal%capillary_radius
									mesh(i,j)%B_ex_poloidal = -Bfield* (1.D0/a + (1.D0-1.D0/a)*Radius**a)

								end select

						enddo
			endif
		enddo

		!--- Boundary Conditions ---!
		mesh(1,:)%B_ex_poloidal = mesh(2,:)%B_ex_poloidal

	end subroutine Move_Window_Bexternal_field

	!--------------------------------------!
	!--- move background ions ---!
	!--------------------------------------!
	subroutine Move_Windows_ion_background
		call remove_outofboundaries_ions
		call inject_ions
	end subroutine Move_Windows_ion_background

END MODULE
