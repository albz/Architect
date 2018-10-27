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
USE external_background_density

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

    if ( mesh(i,j)%ne_bck.le.threshold_factor) then
      Jpe_r(i,j) = 0.
      Jpe_z(i,j) = 0.
    else
      beta_x(i,j) = mesh(i,j)%ux/sqrt( 1. + mesh(i,j)%ux**2+ mesh(i,j)%uz**2 + 1e-10 )
      beta_z(i,j) = mesh(i,j)%uz/sqrt( 1. + mesh(i,j)%ux**2+ mesh(i,j)%uz**2 + 1e-10 )

      Jpe_r(i,j) = mesh(i,j)%ne_bck*beta_x(i,j) !face centered
      Jpe_z(i,j) = mesh(i,j)%ne_bck*beta_z(i,j) !face centered

      !--- with smoothing ---!
      ! if(i.ge.2 .and. j.ge.3 .and. i.le.Node_end_z-1 .and. j.le.Node_end_r-1) then
      ! Jpe_r(i,j) = &
      !              ( mesh(i  ,j  )%ne_bck*beta_x(i  ,j  )*((j-1)**2-(j-2)**2) &
      !              ! + mesh(i-1,j  )%ne_bck*beta_x(i-1,j  )*((j-1)**2-(j-2)**2) &
      !              ! + mesh(i+1,j  )%ne_bck*beta_x(i+1,j  )*((j-1)**2-(j-2)**2) &
      !              + mesh(i  ,j+1)%ne_bck*beta_x(i  ,j+1)*((j+0)**2-(j-1)**2) &
      !              ! + mesh(i-1,j+1)%ne_bck*beta_x(i-1,j+1)*((j+0)**2-(j-1)**2) &
      !              ! + mesh(i+1,j+1)%ne_bck*beta_x(i+1,j+1)*((j+0)**2-(j-1)**2) &
      !              + mesh(i  ,j-1)%ne_bck*beta_x(i  ,j-1)*((j-2)**2-(j-3)**2) &
      !              ! + mesh(i+1,j-1)%ne_bck*beta_x(i+1,j-1)*((j-2)**2-(j-3)**2) &
      !              ! + mesh(i-1,j-1)%ne_bck*beta_x(i-1,j-1)*((j-2)**2-(j-3)**2) ) &
      !              ) / (((j+0)**2-(j-3)**2)*1.)
      ! endif
      ! if(j.eq.2 .and. i.ge.2 .and. i.le.Node_end_z-1) then
      !   Jpe_r(i,j) = &
      !             ( mesh(i  ,j  )%ne_bck*beta_x(i  ,j  )*((j-1)**2-(j-2)**2) &
      !             ! + mesh(i-1,j  )%ne_bck*beta_x(i-1,j  )*((j-1)**2-(j-2)**2) &
      !             ! + mesh(i+1,j  )%ne_bck*beta_x(i+1,j  )*((j-1)**2-(j-2)**2) &
      !             + mesh(i  ,j+1)%ne_bck*beta_x(i  ,j+1)*((j+0)**2-(j-1)**2) &
      !             ! + mesh(i-1,j+1)%ne_bck*beta_x(i-1,j+1)*((j+0)**2-(j-1)**2) &
      !             ! + mesh(i+1,j+1)%ne_bck*beta_x(i+1,j+1)*((j+0)**2-(j-1)**2) ) &
      !             ) / (((j+0)**2-(j-2)**2)*1.) ! (((2+0)**2-(0)**2)*3.)
      ! endif
      !--- end smoothing ---!
    endif

  enddo
  enddo

  !--------------------------------------------------------!
	!  Centers the currents as needed by FDTD field solver   !
  !--------------------------------------------------------!

	! Jr, Jz for FDTD are centered as in the Yee cell
    do i= 2,Node_max_z
    do j= 2,Node_max_r
        !--- Jpe_r centering with 1/4 of cell weighting ---!
        !Jpe_r_for_FDTD(i,j) = &
        !         ( Jpe_r(i-1,j)*( (j-1)**2-(j-1.5)**2 )+Jpe_r(i-1,j+1)*( (j-0.5)**2-(j-1)**2 ) ) / (4.*(j-1)) &
        !  +      ( Jpe_r(i+0,j)*( (j-1)**2-(j-1.5)**2 )+Jpe_r(i+0,j+1)*( (j-0.5)**2-(j-1)**2 ) ) / (4.*(j-1))
        !--- Jpe_r centering with full cells weighting ---!
      Jpe_r_for_FDTD(i,j) = &
                  ( Jpe_r(i-1,j)*( (j-1)**2-(j-2)**2 )+Jpe_r(i-1,j+1)*( (j+0)**2-(j-1)**2 ) ) / (8.*(j-1)) &
             +    ( Jpe_r(i+0,j)*( (j-1)**2-(j-2)**2 )+Jpe_r(i+0,j+1)*( (j+0)**2-(j-1)**2 ) ) / (8.*(j-1))

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
    !--- ---!
    ! Jpe_r_for_FDTD(i,2)  = 0.5*(Jpe_r_for_FDTD(i,1)+Jpe_r_for_FDTD(i,3))
    !--- ---!
		mesh(i,1)%ne_bck = mesh(i,Node_min_r)%ne_bck
  enddo
  ! upper boundary
  Jpe_r_for_FDTD(:,Node_end_r)  = Jpe_r_for_FDTD(:,Node_max_r)
  Jpe_z_for_FDTD(:,Node_end_r)  = Jpe_z_for_FDTD(:,Node_max_r)
  mesh(:,Node_end_r)%ne_bck = mesh(:,Node_max_r)%ne_bck

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
    mesh(:,:)%ne_bck 	= 0.d0
    !--- *** ---!

    do i=2,mesh_par%Nzm
          do j=2,mesh_par%Nxm
            mesh(i,j)%ne_bck=background_density_value(i,j)
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
        if(mesh(i,j)%ne_bck>0.) mesh(i,j)%uz=sim_parameters%velocity_background
      enddo
    enddo

    !--- *** ---! caluclate Bphi
    i=2
    do j=1,Node_max_r
      mesh(i,j)%Bphi=(mesh(i,j)%uz*mesh(i,j)%ne_bck)*j*mesh_par%dr/2.
      if(j>mesh_par%Nr_plasma) then
        mesh(i,j)%Bphi=(mesh(i,mesh_par%Nr_plasma-10)%uz*mesh(i,mesh_par%Nr_plasma-10)%ne_bck) &
        *(mesh_par%Nr_plasma*mesh_par%dr)**2/2./j/mesh_par%dr
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


!--- master function that select the background generation strategy ---!
real(8) FUNCTION background_density_value(i,j)
  integer, intent(in) :: i,j
  if(bck_plasma%external_density) then
    background_density_value = background_density_value_external(i,j)
  else
    background_density_value = background_density_value_analytical(i,j)
  endif
END FUNCTION background_density_value
!--- *** ---!

  ! subroutine background_density_value_analytical(i,j,background_density_value_analyticalAAA)
  real(8) FUNCTION background_density_value_analytical(i,j)
    integer, intent(IN) :: i,j
    real(8) :: Zposition,slope,den
    real(8) :: i_eff,j_eff, radius
    real(8) :: weightR, weightZ
    real(8) :: delta_alpha,kappa_z,A,B,C
    integer :: k

    i_eff=real(i-2)+.5d0
    j_eff=real(j-2)+.5d0

    radius = j_eff*mesh_par%dr/plasma%k_p

    background_density_value_analytical=0.d0
    Zposition=mesh_par%z_min_moving_um+i_eff*mesh_par%dz/plasma%k_p

    !--vacuum layer---!
    !if(j>=mesh_par%Nr_plasma) return
    !----------------------!

    weightR=1.d0
    weightZ=1.d0

    do k=1,8
      if(Zposition<=bck_plasma%z_coordinate_um(k) .and. Zposition>bck_plasma%z_coordinate_um(k+1)) then

        !---*** Longitudinal profile ***---!
        if(  bck_plasma%order_logitudinal(k)==0 .or. bck_plasma%order_logitudinal(k)==1 ) then !linear or flat-top
            slope=(bck_plasma%n_over_n0(k)-bck_plasma%n_over_n0(k+1))/(bck_plasma%z_coordinate_um(k)-bck_plasma%z_coordinate_um(k+1))
            weightZ=bck_plasma%n_over_n0(k)+slope*(Zposition-bck_plasma%z_coordinate_um(k))
            !---------------------------------------------------------------------------------------------------------!

        else if( bck_plasma%order_logitudinal(k)==2 ) then !parabolic 'Concave' profile
          if( bck_plasma%n_over_n0(k+1)>bck_plasma%n_over_n0(k) ) then
            A=(bck_plasma%n_over_n0(k)-bck_plasma%n_over_n0(k+1))/(bck_plasma%z_coordinate_um(k+1)-bck_plasma%z_coordinate_um(k))**2
            B=2.0*(bck_plasma%n_over_n0(k+1)-bck_plasma%n_over_n0(k))*bck_plasma%z_coordinate_um(k+1)
            B=B/(bck_plasma%z_coordinate_um(k+1)-bck_plasma%z_coordinate_um(k))**2
            C=bck_plasma%n_over_n0(k+1)*bck_plasma%z_coordinate_um(k)**2
            C=C-2.0*bck_plasma%n_over_n0(k+1)*bck_plasma%z_coordinate_um(k)*bck_plasma%z_coordinate_um(k+1)
            C=C+bck_plasma%n_over_n0(k)*bck_plasma%z_coordinate_um(k+1)**2
            C=C/(bck_plasma%z_coordinate_um(k+1)-bck_plasma%z_coordinate_um(k))**2
            weightZ=A*Zposition**2+B*Zposition+C
          else
            A=(bck_plasma%n_over_n0(k+1)-bck_plasma%n_over_n0(k))/(bck_plasma%z_coordinate_um(k)-bck_plasma%z_coordinate_um(k+1))**2
            B=2.0*(bck_plasma%n_over_n0(k)-bck_plasma%n_over_n0(k+1))*bck_plasma%z_coordinate_um(k)
            B=B/(bck_plasma%z_coordinate_um(k)-bck_plasma%z_coordinate_um(k+1))**2
            C=bck_plasma%n_over_n0(k+1)*bck_plasma%z_coordinate_um(k)**2
            C=C-2.0*bck_plasma%n_over_n0(k)*bck_plasma%z_coordinate_um(k)*bck_plasma%z_coordinate_um(k+1)
            C=C+bck_plasma%n_over_n0(k)*bck_plasma%z_coordinate_um(k+1)**2
            C=C/(bck_plasma%z_coordinate_um(k)-bck_plasma%z_coordinate_um(k+1))**2
            weightZ=A*Zposition**2+B*Zposition+C
          endif
          !---------------------------------------------------------------------------------------------------------!

        else if( bck_plasma%order_logitudinal(k)==-2 ) then !parabolic 'Convex' profile
          if( bck_plasma%n_over_n0(k+1)>bck_plasma%n_over_n0(k) ) then
            A=(bck_plasma%n_over_n0(k+1)-bck_plasma%n_over_n0(k))/(bck_plasma%z_coordinate_um(k)-bck_plasma%z_coordinate_um(k+1))**2
            B=2.0*(bck_plasma%n_over_n0(k)-bck_plasma%n_over_n0(k+1))*bck_plasma%z_coordinate_um(k)
            B=B/(bck_plasma%z_coordinate_um(k)-bck_plasma%z_coordinate_um(k+1))**2
            C=bck_plasma%n_over_n0(k+1)*bck_plasma%z_coordinate_um(k)**2
            C=C-2.0*bck_plasma%n_over_n0(k)*bck_plasma%z_coordinate_um(k)*bck_plasma%z_coordinate_um(k+1)
            C=C+bck_plasma%n_over_n0(k)*bck_plasma%z_coordinate_um(k+1)**2
            C=C/(bck_plasma%z_coordinate_um(k)-bck_plasma%z_coordinate_um(k+1))**2
            weightZ=A*Zposition**2+B*Zposition+C
          else
            A=(bck_plasma%n_over_n0(k)-bck_plasma%n_over_n0(k+1))/(bck_plasma%z_coordinate_um(k+1)-bck_plasma%z_coordinate_um(k))**2
            B=2.0*(bck_plasma%n_over_n0(k+1)-bck_plasma%n_over_n0(k))*bck_plasma%z_coordinate_um(k+1)
            B=B/(bck_plasma%z_coordinate_um(k+1)-bck_plasma%z_coordinate_um(k))**2
            C=bck_plasma%n_over_n0(k+1)*bck_plasma%z_coordinate_um(k)**2
            C=C-2.0*bck_plasma%n_over_n0(k+1)*bck_plasma%z_coordinate_um(k)*bck_plasma%z_coordinate_um(k+1)
            C=C+bck_plasma%n_over_n0(k)*bck_plasma%z_coordinate_um(k+1)**2
            C=C/(bck_plasma%z_coordinate_um(k+1)-bck_plasma%z_coordinate_um(k))**2
            weightZ=A*Zposition**2+B*Zposition+C
          endif


          !---------------------------------------------------------------------------------------------------------!
        else if( bck_plasma%order_logitudinal(k)==3 ) then !COS^2 profile
          if( bck_plasma%n_over_n0(k+1)>bck_plasma%n_over_n0(k) ) then
            delta_alpha=bck_plasma%n_over_n0(k+1)-bck_plasma%n_over_n0(k)
            kappa_z= ( Zposition-bck_plasma%z_coordinate_um(k) ) / ( bck_plasma%z_coordinate_um(k)-bck_plasma%z_coordinate_um(k+1) )
            weightZ=delta_alpha*(1.d0-cos(kappa_z*pi/2.0)**2)+bck_plasma%n_over_n0(k)
          else
            delta_alpha=bck_plasma%n_over_n0(k)-bck_plasma%n_over_n0(k+1)
            kappa_z= ( Zposition-bck_plasma%z_coordinate_um(k) ) / ( bck_plasma%z_coordinate_um(k)-bck_plasma%z_coordinate_um(k+1) )
            weightZ=delta_alpha*(cos(kappa_z*pi/2.0)**2)+bck_plasma%n_over_n0(k+1)
          endif
          !---------------------------------------------------------------------------------------------------------!
        endif !endif logitudinal shaping

        !---*** TRansveres-Radial profile ***---!
        if(bck_plasma%order_radial(k)==0) weightR=1.0!bck_plasma%n_over_n0(k)
        if(bck_plasma%order_radial(k)==3) weightR = cos(radius/bck_plasma%radius_um(k)*pi/2.0)**2 !COS^2 profile
        if(bck_plasma%order_radial(k)==4) weightR=1.d0+bck_plasma%perturbation_amplitude(k)*(1.d0-2.d0*cos(radius/bck_plasma%radius_um(k)*pi/2.0)**2) !1 + A (1 - 2 Cos[r/R \Pi/2]^2)
        if(bck_plasma%order_radial(k)==5) then !Linear up to radius_internal_um, decay in COS^2
                                              if( radius < bck_plasma%radius_internal_um(k) ) then
                                                                weightR=1.0
                                              else
                                                                weightR = cos(pi/2.0*(radius-bck_plasma%radius_internal_um(k))/(bck_plasma%radius_um(k)-bck_plasma%radius_internal_um(k)))**2
                                              endif
        endif
        if(radius >   bck_plasma%radius_um(k)) weightR=0.0

        background_density_value_analytical=weightZ*weightR

      endif !Zposition
    enddo

  END FUNCTION background_density_value_analytical

END MODULE
