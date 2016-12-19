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


!--- --- --- --- --- --- --- --- --- ---!
    subroutine set_initial_background_condition
      call set_initial_plasma_density
      call set_initial_velocity
    end subroutine set_initial_background_condition
!--- --- --- --- --- --- --- --- --- ---!


  SUBROUTINE set_initial_plasma_density
    ! define capillary density at first iteration
    REAL(8) :: test_ne,Radius,Zposition,bck_density_value
		INTEGER j,i,k

    !--- *** ---!
    bck_plasma%z_coordinate=bck_plasma%z_coordinate_um*plasma%k_p
    mesh(:,:)%n_plasma_e 	= 0.d0
    !--- *** ---!

    do i=2,mesh_par%Nzm
          do j=2,mesh_par%Nxm
            mesh(i,j)%n_plasma_e=background_density_value(i,j)
      enddo
    enddo
	END SUBROUTINE set_initial_plasma_density




  SUBROUTINE set_initial_velocity
  ! define capillary velocity at first iteration
  REAL(8) :: test_ne
  REAL(8) :: ne_m3,R_current_m
  INTEGER j,i

  mesh(:,:)%ux 			    = 0.
  mesh(:,:)%uz 			    = 0.


  !---***---!
  if(Bpoloidal%L_BfieldfromV) then
    allocate(mesh_util%Bphi_BC_Left(Node_max_r))

    !--- *** ---! calculate background velocity
    ne_m3=plasma%n0*1e6
    R_current_m=Bpoloidal%capillary_radius_um*1d-6
    sim_parameters%velocity_background = Bpoloidal%background_current_A(1) / (R_current_m**2*pi*ne_m3*electron_charge*c_SI)
    sim_parameters%velocity_background = -1.D0*sim_parameters%velocity_background
    write(*,*) ''
    write(*,'(A)') 'Case :: with Moving Background >'
    write(*,'(A,1e11.3,A)') 'The background velocity is',sim_parameters%velocity_background,' * c'
    write(*,*)


    !--- *** ---! assign velocity
    do i=1,Node_max_z
      do j=1,Node_max_r
        if(mesh(i,j)%n_plasma_e>0.) mesh(i,j)%uz=sim_parameters%velocity_background
      enddo
    enddo

    !--- *** ---! caluclate Bphi
    i=2
    do j=1,Node_max_r
      mesh(i,j)%Bphi=(mesh(i,j)%uz*mesh(i,j)%n_plasma_e)*j*mesh_par%dxm/2.
      if(j>mesh_par%NRmax_plasma) then
        mesh(i,j)%Bphi=(mesh(i,mesh_par%NRmax_plasma-10)%uz*mesh(i,mesh_par%NRmax_plasma-10)%n_plasma_e) &
        *(mesh_par%NRmax_plasma*mesh_par%dxm)**2/2./j/mesh_par%dxm
      endif
      !---
      mesh_util%Bphi_BC_Left(j)=mesh(i,j)%Bphi
      !---
    enddo
    !--- now forcing the copy of i=1 to all-i
    do i=2,Node_max_z
      mesh(i,:)%Bphi=mesh_util%Bphi_BC_Left(:)
    end do
  endif
END SUBROUTINE set_initial_velocity

  FUNCTION background_density_value(i,j)
    real(8) :: background_density_value,Zposition,slope,den
    real(8) :: i_eff,j_eff, radius
    real(8) :: weightR, weightZ
    integer :: i,j,k

    i_eff=real(i-2)+.5d0
    j_eff=real(j-2)+.5d0

    radius = j_eff*mesh_par%dxm/plasma%k_p

    background_density_value=0.d0
    Zposition=mesh_par%z_min_moving_um+i_eff*mesh_par%dzm/plasma%k_p

    !--vacuum layer---!
    !if(j>=mesh_par%NRmax_plasma) return
    !----------------------!

    weightR=1.d0
    weightZ=1.d0

    do k=1,8
      if(Zposition<=bck_plasma%z_coordinate_um(k) .and. Zposition>bck_plasma%z_coordinate_um(k+1)) then

        !---*** Longitudinal profile ***---!
        if(bck_plasma%order_logitudinal(k)==0) weightZ=bck_plasma%n_over_n0(k)
        if(bck_plasma%order_logitudinal(k)==1) then
          slope=(bck_plasma%n_over_n0(k-1)-bck_plasma%n_over_n0(k))/(bck_plasma%z_coordinate_um(k)-bck_plasma%z_coordinate_um(k+1))
          weightZ=bck_plasma%n_over_n0(k-1)+slope*(Zposition-bck_plasma%z_coordinate_um(k))
          if(k==1) then
            slope=(0.d0-bck_plasma%n_over_n0(k+1))/(bck_plasma%z_coordinate_um(k)-bck_plasma%z_coordinate_um(k+1))
            weightZ=0.d0+slope*(Zposition-bck_plasma%z_coordinate_um(k))
          endif
        endif

        !---*** TRansveres-Radial profile ***---!
        if(bck_plasma%order_radial(k)==0) weightR=bck_plasma%n_over_n0(k)
        if(bck_plasma%order_radial(k)==3) weightR = cos(radius/bck_plasma%radius_um(k)*pi/2.0)**2 !COS^2 profile
        if(bck_plasma%order_radial(k)==4) weightR=1.d0+bck_plasma%perturbation_amplitude(k)*(1.d0-2.d0*cos(radius/bck_plasma%radius_um(k)*pi/2.0)**2) !1 + A (1 - 2 Cos[r/R \Pi/2]^2)
        if(radius >   bck_plasma%radius_um(k)) weightR=0.0

        background_density_value=weightZ*weightR

      endif !Zposition
    enddo

  END FUNCTION background_density_value

END MODULE
