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
			bunch(1)%part(:)%cmp(14)=1.
			sim_parameters%zg = calculate_nth_moment(1,1,3,'nocentral')
		else if (sim_parameters%window_mode.eq.1) then ! window moves with constant speed
			sim_parameters%zg = - sim_parameters%moving_window_speed*c*(sim_parameters%iter+1.)*sim_parameters%dt
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
   INTEGER cells_advanced,i,j
   REAL(8) :: DeltaR,DeltaZ

	!----------------------------------------------!
	!       Longitudinal shift by one cell         !
	!----------------------------------------------!

	do i=mesh_par%Nzm,2,-1
		mesh(i,:)%Ex=mesh(i-1,:)%Ex
		mesh(i,:)%Ez=mesh(i-1,:)%Ez
		mesh(i,:)%Bphi=mesh(i-1,:)%Bphi
		mesh(i,:)%Bphi_old=mesh(i-1,:)%Bphi_old

		mesh(i,:)%uz         = mesh(i-1,:)%uz
		mesh(i,:)%ux         = mesh(i-1,:)%ux
		mesh(i,:)%n_plasma_e = mesh(i-1,:)%n_plasma_e
		mesh(i,:)%n_plasma_i = mesh(i-1,:)%n_plasma_i

		if(sim_parameters%L_Bunch_evolve) then
			mesh(i,:)%Ex_bunch			  = mesh(i-1,:)%Ex_bunch
			mesh(i,:)%Ez_bunch			  = mesh(i-1,:)%Ez_bunch
			mesh(i,:)%Bphi_bunch		  = mesh(i-1,:)%Bphi_bunch
			mesh(i,:)%Bphi_old_bunch	= mesh(i-1,:)%Bphi_old_bunch
		endif
	enddo

	mesh(1:2,:)%Ex = 0.D0
	mesh(1:2,:)%Ez = 0.D0
	mesh(1:2,:)%Bphi = 0.D0
  mesh(1:2,:)%Bphi_old= 0.D0
	mesh(1:2,:)%uz= 0.D0
	mesh(1:2,:)%ux= 0.D0

	if(sim_parameters%L_Bunch_evolve) then
		mesh(1:2,:)%Ex_bunch			  = 0.D0
		mesh(1:2,:)%Ez_bunch			  = 0.D0
		mesh(1:2,:)%Bphi_bunch		  = 0.D0
		mesh(1:2,:)%Bphi_old_bunch	= 0.D0
	endif


	!---background density---!
 	do j= 2,mesh_par%Nxm-1
 	 		mesh(2,j)%n_plasma_e = background_density_value(2,j)
			mesh(1,j)%n_plasma_e = background_density_value(1,j)
			if(.not.ionisation%L_ionisation) mesh(2,j)%n_plasma_i=mesh(2,j)%n_plasma_e
 	enddo
 	!BC
 	mesh(2,mesh_par%Nxm-1)%n_plasma_e = mesh(2,mesh_par%Nxm-1)%n_plasma_e ! upper boundary
 	mesh(2,1             )%n_plasma_e = mesh(2,2             )%n_plasma_e ! lower boundary
 	mesh(1,:             )%n_plasma_e = mesh(2,:             )%n_plasma_e ! left  boundary
	if(.not.ionisation%L_ionisation) then
	 	mesh(2,1)%n_plasma_i = mesh(2,2)%n_plasma_i ! lower boundary
	 	mesh(1,:)%n_plasma_i = mesh(2,:)%n_plasma_i ! left  boundary
	endif
 	!--- ---!

	!--- Moving Background ---!
	if(Bpoloidal%L_BfieldfromV) then
		mesh(1,:)%uz = sim_parameters%velocity_background
		mesh(1,:)%Bphi     = mesh_util%Bphi_BC_Left(:)
		mesh(1,:)%Bphi_old = mesh_util%Bphi_BC_Left(:)
	end if
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
