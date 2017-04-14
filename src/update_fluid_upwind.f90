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

MODULE Advance_Fluid_Compact_Upwind_New


USE my_types
USE use_my_types

IMPLICIT NONE

CONTAINS

! Upwind, Sprang operator splitting (half advection, electromagnetic advance, half advection)


!---! !---! !---!
   SUBROUTINE AdvanceFluid_CompactUpwind_New


   USE my_types

   IMPLICIT NONE

   INTEGER :: i,j,iter,q
   REAL(8), DIMENSION(mesh_par%Nzm,mesh_par%Nxm) :: ux,uz,ne,ux_halfDt,uz_halfDt,ux_new,uz_new,ne_new
   REAL(8), DIMENSION(mesh_par%Nzm,mesh_par%Nxm) :: Ex_f,Ez_f,Bphi_f
   REAL(8), DIMENSION(mesh_par%Nzm,mesh_par%Nxm) :: ne_ux_halfDt,ne_uz_halfDt,ne_ux_new,ne_uz_new

   REAL(8) threshold_factor

   REAL(8) beta_x_halfDt,beta_z_halfDt
   REAL(8) beta_z_left,beta_z_right,beta_x_down,beta_x_up
   REAL(8) ux_temp,uz_temp
   REAL(8) q_min,q_max
   REAL(8) A_c_left,A_c_right,A_c_up,A_c_down

   REAL(8), DIMENSION(mesh_par%Nzm,mesh_par%Nxm) :: beta_z_i_minus_halfDz,beta_x_j_minus_halfDr

   REAL(8), DIMENSION(mesh_par%Nzm,mesh_par%Nxm) :: quantity_td,quantity
   REAL(8), DIMENSION(mesh_par%Nzm,mesh_par%Nxm) :: flux_lo_z,flux_lo_right,flux_lo_up,flux_lo_r
   REAL(8), DIMENSION(mesh_par%Nzm,mesh_par%Nxm) :: flux_ho_z,flux_ho_right,flux_ho_up,flux_ho_r
   REAL(8), DIMENSION(mesh_par%Nzm,mesh_par%Nxm) :: A_left,A_down
   REAL(8), DIMENSION(mesh_par%Nzm,mesh_par%Nxm) :: C_left,C_right,C_up,C_down
   REAL(8), DIMENSION(mesh_par%Nzm,mesh_par%Nxm) :: q__plus,q__minus
   REAL(8), DIMENSION(mesh_par%Nzm,mesh_par%Nxm) :: Q_plus,Q_minus,P_plus,P_minus
   REAL(8), DIMENSION(mesh_par%Nzm,mesh_par%Nxm) :: R_plus,R_minus
   !REAL, DIMENSION(mesh_par%Nzm,mesh_par%Nxm) :: beta_x,beta_z

   !--------------!
   !     Mesh     !
   !--------------!

   threshold_factor = 1d-3

   !---------------------------------------------------------------------------------!
   !                                     Z axis                                      !
   !         									     !
   !     ghost cell                  physical domain                ghost cell       !
   !  |______________|___________________________________________|____________|      !
   !  |              |            |                |             |            |      !
   !       1           Node_min_z                     Node_max_z   Node_end_z        !
   !---------------------------------------------------------------------------------!
   !                                    Fluid Mesh				     !
   !---------------------------------------------------------------------------------!
   !                                      R axis                                     !
   !										     !
   !     ghost cell                  physical domain                ghost cell       !
   !  |______________|___________________________________________|____________|      !
   !  |              |            |                |             |            |      !
   !       1           Node_min_z                     Node_max_z   Node_end_z        !
   !                   Axis cell                                                     !
   !---------------------------------------------------------------------------------!

   !-----------------------!
   !   Initial Conditions  !
   !-----------------------!
   ux  =	mesh(:,:)%ux
   uz  =	mesh(:,:)%uz
   ne  =	mesh(:,:)%n_plasma_e

   ! Initialize auxiliary vectors and reals
   ne_ux_halfDt     = 0.D0
   ne_uz_halfDt     = 0.D0

   ux_new           = 0.D0
   uz_new           = 0.D0

   ne_new	    	= 0.D0

   beta_x_halfDt    = 0.D0
   beta_z_halfDt    = 0.D0

   Ez_f		    	= 0.D0
   Ex_f	        = 0.D0
   Bphi_f	    	= 0.D0

	!------------------------------!
	! Advection begins here 	 !
	!------------------------------!

	!--> beta-face-centered
	do i= Node_min_lo_z,Node_end_lo_z
	do j= Node_min_lo_r,Node_end_lo_r
		beta_z_left   = uz(i-1,j)/sqrt( 1.D0 + ux(i-1,j)**2 + uz(i-1,j)**2 + threshold_factor)
		beta_z_right  = uz(i  ,j)/sqrt( 1.D0 + ux(i  ,j)**2 + uz(i  ,j)**2 + threshold_factor)
		beta_x_down   = ux(i,j-1)/sqrt( 1.D0 + ux(i,j-1)**2 + uz(i,j-1)**2 + threshold_factor)
		beta_x_up     = ux(i,j  )/sqrt( 1.D0 + ux(i,j  )**2 + uz(i,j  )**2 + threshold_factor)

		beta_z_i_minus_halfDz(i,j) = 0.5D0*( beta_z_left + beta_z_right )
		beta_x_j_minus_halfDr(i,j) = 0.5D0*( beta_x_down + beta_x_up    )
	end do
	end do
	!--->Remove<---!beta_z_i_minus_halfDz(2,:)=beta_z_i_minus_halfDz(3,:)

!to be structured to be converted into a function
do q=1,3
	quantity    = 0.D0
	quantity_td = 0.D0
	if (q.eq.1) quantity = ne
	if (q.eq.2) quantity = ne*ux
	if (q.eq.3) quantity = ne*uz

	flux_lo_z=0.D0
	flux_lo_r=0.D0
	!--- donor-cell
	do i= Node_min_lo_z,Node_max_lo_z
	do j= Node_min_lo_r,Node_max_lo_r
    flux_lo_z (i,j) = quantity(i-1,j  ) * max(0.D0,beta_z_i_minus_halfDz(i,j)) + quantity(i  ,j  ) * min(0.D0,beta_z_i_minus_halfDz(i,j))
    flux_lo_r (i,j) = quantity(i  ,j-1) * max(0.D0,beta_x_j_minus_halfDr(i,j)) + quantity(i  ,j  ) * min(0.D0,beta_x_j_minus_halfDr(i,j))
	end do
	end do
	!--- left border ---!
	flux_lo_z (1,Node_min_lo_r:Node_max_lo_r) = flux_lo_z (Node_min_lo_z,Node_min_lo_r:Node_max_lo_r)
	!--- upper border ---!
	flux_lo_r (Node_min_lo_z:Node_max_lo_z,Node_max_lo_r+1) = flux_lo_r (Node_min_lo_z:Node_max_lo_z,Node_max_lo_r)


	!--- compute the low order solution
    do i= Node_min_lo_z,Node_max_lo_z
	do j= Node_min_lo_r,Node_max_lo_r
		quantity_td(i,j)= quantity(i,j) - Dt/Vol(j) * (Sr(j+1)*flux_lo_r(i,j+1)-Sr(j)*flux_lo_r(i,j))
	end do
	end do

	do i= Node_min_lo_z,Node_max_lo_z
	do j= Node_min_lo_r,Node_max_lo_r
		quantity_td(i,j)= quantity_td(i,j) - Dt*Sz(j)/Vol(j) * (flux_lo_z(i+1,j)-flux_lo_z(i,j))
	end do
	end do


	!---------------------------------------------!
	!  Boundary conditions for low order solution !
	!---------------------------------------------!

	if (q.eq.1) then
		!--- upper boundary ---!
		quantity_td(Node_min_lo_z:Node_max_lo_z,Node_end_lo_r) = quantity_td(Node_min_lo_z:Node_max_lo_z,Node_max_lo_r)
 		!--- lower boundary ---!
		quantity_td(Node_min_lo_z:Node_max_lo_z,1) = quantity_td(Node_min_lo_z:Node_max_lo_z,Node_min_lo_r)
 		!--- left boundary ---!
		!quantity_td(1,Node_min_lo_r:Node_max_lo_r) = quantity_td(Node_min_lo_z,Node_min_lo_r:Node_max_lo_r)
 		!--- right boundary ---!
        quantity_td (Node_end_lo_z,Node_min_lo_r:Node_max_lo_r) = quantity_td(Node_max_lo_z,Node_min_lo_r:Node_max_lo_r)

	else if (q.eq.2) then
		!--- upper boundary ---!
  		quantity_td(Node_min_lo_z:Node_max_lo_z,Node_end_lo_r)  =  quantity_td(Node_min_lo_z:Node_max_lo_z,Node_max_lo_r)
   		!--- lower boundary ---!
  		quantity_td(Node_min_lo_z:Node_max_lo_z,1) = -quantity_td(Node_min_lo_z:Node_max_lo_z,Node_min_lo_r)
   		!--- left boundary ---!
  		!quantity_td(1,Node_min_lo_r:Node_max_lo_r) =  0.D0
   		!--- right boundary ---!
  		quantity_td(Node_end_lo_z,Node_min_lo_r:Node_max_lo_r) = 0.D0

	else if (q.eq.3) then
		!--- upper boundary ---!
  		quantity_td(Node_min_lo_z:Node_max_lo_z,Node_end_lo_r) = quantity_td(Node_min_lo_z:Node_max_lo_z,Node_max_lo_r)
   		!--- lower boundary ---!
  		quantity_td(Node_min_lo_z:Node_max_lo_z,1) = quantity_td(Node_min_lo_z:Node_max_lo_z,Node_min_lo_r)
   		!--- left boundary ---!
  		!quantity_td(1,Node_min_lo_r:Node_max_lo_r) =  0.D0
   		!--- right boundary ---!
  		quantity_td(Node_end_lo_z,Node_min_lo_r:Node_max_lo_r)  = 0.D0

  end if



 if (quantity(i,j).ne.quantity(i,j)) write(*,*) 'Error, Nan in q,',q,'--- i,j',i,j


  do i= Node_min_lo_z,Node_max_lo_z
	do j= Node_min_lo_r,Node_max_lo_r
		quantity(i,j)= quantity_td(i,j)
	end do
	end do

	! low order scheme on axis
	do i= Node_min_ho_z,Node_max_ho_z-1
		quantity(i,2) = quantity_td(i,2)
		quantity(i,3) = quantity_td(i,3)
	end do
	! low order scheme on the left (not ghost cell)
	do j = Node_min_lo_r,Node_max_lo_r
        quantity (Node_min_ho_z,j)     = quantity_td(Node_min_ho_z,j)
		quantity (Node_min_lo_z,j)     = quantity_td(Node_min_lo_z,j)
 	enddo

    !   ------- backward substitution
	if (q.eq.1) ne_new       = quantity
	if (q.eq.2) ne_ux_halfDt = quantity
	if (q.eq.3) ne_uz_halfDt = quantity
enddo



!-------------------------- Electromagnetic Advance ------------------------------
        do i= Node_min_lo_z,Node_max_lo_z
        do j= Node_min_lo_r,Node_max_lo_r
      Ez_f   (i,j) =        mesh(i,j  )%Ez + mesh(i,j  )%Ez_bunch           ! properly centered in space
      Ex_f   (i,j) = 0.25D0*( mesh(i,j  )%Ex       + mesh(i+1,j  )%Ex   &
                          + mesh(i,j-1)%Ex       + mesh(i+1,j-1)%Ex   ) &   ! properly centered in space
                    +0.25D0*( mesh(i,j  )%Ex_bunch + mesh(i+1,j  )%Ex_bunch &
                          + mesh(i,j-1)%Ex_bunch + mesh(i+1,j-1)%Ex_bunch ) ! properly centered in space
      Bphi_f (i,j) = 0.25D0*( mesh(i,j  )%Bphi     + mesh(i  ,j-1)%Bphi &
                          + mesh(i,j  )%Bphi_old + mesh(i  ,j-1)%Bphi_old ) & ! properly centered in time
                    +0.25D0*( mesh(i,j  )%Bphi_bunch     + mesh(i  ,j-1)%Bphi_bunch &
                          + mesh(i,j  )%Bphi_old_bunch + mesh(i  ,j-1)%Bphi_old_bunch ) ! properly centered in time
        enddo
        enddo



		do i= Node_min_lo_z,Node_max_lo_z
		do j= Node_min_lo_r,Node_max_lo_r

			if (ne_new(i,j).le.threshold_factor) then
					beta_x_halfDt = 0.D0
					beta_z_halfDt = 0.D0
					!ne_new(i,j)   = threshold_factor
			else


	ux_temp = ne_ux_halfDt(i,j)/max(ne_new(i,j),threshold_factor)
	uz_temp = ne_uz_halfDt(i,j)/max(ne_new(i,j),threshold_factor)
	beta_x_halfDt = ux_temp/ sqrt(1.D0+ux_temp**2+uz_temp**2)
	beta_z_halfDt = uz_temp/ sqrt(1.D0+ux_temp**2+uz_temp**2)


			endif

  ux_new(i,j) = ux_temp + Dt*( Ex_f(i,j) - beta_z_halfDt*Bphi_f(i,j) )
	uz_new(i,j) = uz_temp + Dt*( Ez_f(i,j) + beta_x_halfDt*Bphi_f(i,j) )

			if ((ux_new(i,j).ne.ux_new(i,j)).or.(ux_new(i,j).ne.ux_new(i,j))) then


				write(*,*) 'Error in iteration ',sim_parameters%iter,'after fluid EM advance'
				write(*,*) 'localized on the mesh in row, col ',i,j
				write(*,*) 'threshold factor',threshold_factor
				write(*,*) 'ne*ux',ne_ux_halfDt(i,j)
				write(*,*) 'ne*uz',ne_uz_halfDt(i,j)
				write(*,*) 'n_new',ne_new(i,j)
				write(*,*) 'ux',ux(i,j)
				write(*,*) 'uz',uz(i,j)
				write(*,*) 'ux_new (after EM advance) ',ux_new(i,j)
				write(*,*) 'uz_new (after EM advance) ',ux_new(i,j)
				write(*,*) 'ux_temp (before EM advance)',ux_temp
				write(*,*) 'uz_temp (before EM advance) ',ux_temp

				write(*,*) 'beta_x_halfDt (for Lorentz force)',beta_x_halfDt
				write(*,*) 'beta_z_halfDt (for Lorentz force)',beta_z_halfDt
				write(*,*) 'Ex',Ex_f(i,j)
				write(*,*) 'Ez',Ez_f(i,j)
				write(*,*) 'B',Bphi_f(i,j)

				stop
			endif
		enddo
enddo



   !-----------------------!
   !  Boundary Conditions  !
   !-----------------------!
   !------------- on low order domain -------------------!
   ! lower boundary - already set in high order domain
   do i = Node_min_lo_z,Node_max_lo_z
		uz_new(i,1)             =  uz_new(i,Node_min_lo_r)
    ux_new(i,1)             = -ux_new(i,Node_min_lo_r)
   enddo

   ! left boundary
!	do j = Node_min_lo_r,Node_max_lo_r
!		ux_new(1,j)   = 0.D0
!		uz_new(1,j)   = 0.D0
!	enddo

   ! right boundary
   do j = Node_min_lo_r,Node_max_lo_r
        ux_new(Node_end_lo_z,j) = 0.D0
        uz_new(Node_end_lo_z,j) = 0.D0
   enddo

    ! upper boundary
   do i = Node_min_lo_z,Node_max_lo_z
        ux_new(i,Node_end_lo_r) = ux_new(i,Node_end_lo_r)
        uz_new(i,Node_end_lo_r) = uz_new(i,Node_end_lo_r)
   enddo


   !-----------------------------------------------------!
   ! Backward Substitution of the new fields in the mesh !
   !-----------------------------------------------------!
   mesh(:,:)%n_plasma_e = ne_new(:,:)
   mesh(:,:)%ux         = ux_new(:,:)
   mesh(:,:)%uz         = uz_new(:,:)

   END SUBROUTINE


	!--- *** ----!
	Subroutine fluid_UpWind
		call Upwind(Dt)
		call EB_forces(Dt)
		!call Upwind(Dt/2.)
	end Subroutine fluid_UpWind

	!--- *** ----!
    SUBROUTINE Upwind(DeltaT)
      INTEGER :: i,j,q
      real(8), intent(in) :: DeltaT
      REAL(8) threshold_factor, flux_down,flux_up,beta_mean
      REAL(8), DIMENSION(mesh_par%Nzm,mesh_par%Nxm) :: beta_r_upwind,beta_z_upwind,quantity,quantity_1,quantity_2
      threshold_factor=1e-10

      do i= 1,Node_end_lo_z
      do j= 1,Node_end_lo_r
        beta_r_upwind(i,j)  = mesh(i,j)%uz/sqrt( 1. + mesh(i,j)%ux**2 + mesh(i,j)%uz**2 + threshold_factor)
        beta_z_upwind(i,j)  = mesh(i,j)%uz/sqrt( 1. + mesh(i,j)%ux**2 + mesh(i,j)%uz**2 + threshold_factor)
      end do
      end do

      do q=1,3
      	quantity    = 0.
      	if (q.eq.1) quantity = mesh%n_plasma_e
      	if (q.eq.2) quantity = mesh%n_plasma_e*mesh%ux
      	if (q.eq.3) quantity = mesh%n_plasma_e*mesh%uz

        !---> along r
      	do i= Node_min_lo_z,Node_max_lo_z
      	do j= Node_min_lo_r,Node_max_lo_r
          beta_mean = 0.5*(beta_r_upwind(i,j-1)+beta_r_upwind(i,j))
          flux_down = quantity(i,j-1) * max(0.,beta_mean) + quantity(i,j) * min(0.,beta_mean)

          beta_mean = 0.5*(beta_r_upwind(i,j)+beta_r_upwind(i,j+1))
          flux_up   = quantity(i,j) * max(0.,beta_mean) + quantity(i,j+1) * min(0.,beta_mean)

          quantity_1(i,j)= quantity(i,j) - DeltaT/Vol(j) * (Sr(j+1)*flux_up-Sr(j)*flux_down)
        enddo
        enddo

        !---> along z
      	do i= Node_min_lo_z,Node_max_lo_z
      	do j= Node_min_lo_r,Node_max_lo_r
          beta_mean = 0.5*(beta_z_upwind(i-1,j)+beta_z_upwind(i,j))
          flux_down = quantity_1(i-1,j) * max(0.,beta_mean) + quantity_1(i,j) * min(0.,beta_mean)

          beta_mean = 0.5*(beta_z_upwind(i,j)+beta_z_upwind(i+1,j))
          flux_up   = quantity_1(i,j) * max(0.,beta_mean) + quantity_1(i+1,j) * min(0.,beta_mean)

          quantity_2(i,j)= quantity_1(i,j) - DeltaT*Sz(j)/Vol(j) * (flux_up-flux_down)
        enddo
        enddo
        quantity=quantity_2
      	!---------------------------------------------!
      	!  Boundary conditions for low order solution !
      	!---------------------------------------------!

      	if (q.eq.1) then
      		! upper boundary
       		do i = Node_min_lo_z,Node_max_lo_z
        		quantity(i,Node_end_lo_r) = quantity(i,Node_max_lo_r)
       		enddo
       		! lower boundary
       		do i = Node_min_lo_z,Node_max_lo_z
        		quantity(i,1)          = quantity(i,Node_min_lo_r)
       		enddo
       		! left boundary
       		do j = Node_min_lo_r,Node_max_lo_r
      			quantity(1,j)              = quantity(Node_min_lo_z,j)
       		enddo
       		! right boundary
      		do j = Node_min_lo_r,Node_max_lo_r
           	quantity (Node_end_lo_z,j)     = quantity(Node_max_lo_z,j)
      		enddo

      	else if (q.eq.2) then
      		! upper boundary
       		do i = Node_min_lo_z,Node_max_lo_z
        		quantity(i,Node_end_lo_r)  =  quantity(i,Node_max_lo_r)
       		enddo
       		! lower boundary
       		do i = Node_min_lo_z,Node_max_lo_z
        		quantity(i,1         )     = -quantity(i,Node_min_lo_r)
       		enddo
     		! left boundary
     		do j = Node_min_lo_r,Node_max_lo_r
        		quantity(1,j)              =  quantity(Node_min_lo_z,j)
     		enddo
     		! right boundary
    		do j = Node_min_lo_r,Node_max_lo_r
       		quantity(Node_end_lo_z,j)  =  0.
    		enddo

    	else if (q.eq.3) then
    		! upper boundary
    		do i = Node_min_lo_z,Node_max_lo_z
        		quantity(i,Node_end_lo_r)  =  quantity(i,Node_max_lo_r)
     		enddo
   		! lower boundary
     		do i = Node_min_lo_z,Node_max_lo_z
        		quantity(i,1         )  =  quantity(i,Node_min_lo_r)
     		enddo
   		! left boundary
     		do j = Node_min_lo_r,Node_max_lo_r
        		quantity(1,j)           =  quantity(Node_min_lo_z,j)
     		enddo
   		! right boundary
    		do j = Node_min_lo_r,Node_max_lo_r
       		quantity(Node_end_lo_z,j)  = 0.
    		enddo
  	endif

    !   ------- backward substitution
  	if (q.eq.1) mesh%n_plasma_e  = quantity
  	if (q.eq.2) then
      do i= 1,Node_end_lo_z
      do j= 1,Node_end_lo_r
        mesh(i,j)%ux = quantity(i,j)/max(mesh(i,j)%n_plasma_e,threshold_factor)
      end do
      end do
    endif
    if (q.eq.3) then
      do i= 1,Node_end_lo_z
      do j= 1,Node_end_lo_r
        mesh(i,j)%uz = quantity(i,j)/max(mesh(i,j)%n_plasma_e,threshold_factor)
      end do
      end do
    endif

  enddo

end subroutine


!-------------------------- Electromagnetic Advance ------------------------------
SUBROUTINE EB_forces(DeltaT)

  INTEGER :: i,j
  real(8), intent(in) :: DeltaT
  REAL(8) :: Ez_fc,Er_fc,Bphi_fc,ur_local,uz_local,beta_r,beta_z,threshold_factor
  threshold_factor=1e-10

    do i= 1,Node_end_lo_z
    do j= 1,Node_end_lo_r
      Ez_fc = mesh(i,j)%Ez + mesh(i,j)%Ez_bunch
      Er_fc = 0.25*( mesh(i,j)%Ex + mesh(i+1,j  )%Ex   &
                   + mesh(i,j-1)%Ex       + mesh(i+1,j-1)%Ex   ) &   ! properly centered in space
             +0.25*( mesh(i,j  )%Ex_bunch + mesh(i+1,j  )%Ex_bunch &
                          + mesh(i,j-1)%Ex_bunch + mesh(i+1,j-1)%Ex_bunch ) ! properly centered in space
      Bphi_fc = 0.25*( mesh(i,j  )%Bphi     + mesh(i  ,j-1)%Bphi &
                          + mesh(i,j  )%Bphi_old + mesh(i  ,j-1)%Bphi_old ) & ! properly centered in time
               +0.25*( mesh(i,j  )%Bphi_bunch     + mesh(i  ,j-1)%Bphi_bunch &
                          + mesh(i,j  )%Bphi_old_bunch + mesh(i  ,j-1)%Bphi_old_bunch ) ! properly centered in time

      beta_r=mesh(i,j)%uz/sqrt( 1. + mesh(i,j)%ux**2 + mesh(i,j)%uz**2 + threshold_factor)
      beta_z=mesh(i,j)%uz/sqrt( 1. + mesh(i,j)%ux**2 + mesh(i,j)%uz**2 + threshold_factor)

      ! ur_local = mesh(i,j)%ux*mesh(i,j)%n_plasma_e + DeltaT*( Er_fc - beta_z*Bphi_fc )
      ! uz_local = mesh(i,j)%uz*mesh(i,j)%n_plasma_e + DeltaT*( Ez_fc + beta_r*Bphi_fc )
      !
      ! mesh(i,j)%ux = ur_local/max( mesh(i,j)%n_plasma_e, threshold_factor )
      ! mesh(i,j)%uz = uz_local/max( mesh(i,j)%n_plasma_e, threshold_factor )

      mesh(i,j)%ux = mesh(i,j)%ux + DeltaT*( Er_fc - beta_z*Bphi_fc )
      mesh(i,j)%uz = mesh(i,j)%uz + DeltaT*( Ez_fc + beta_r*Bphi_fc )


    enddo
    enddo

   !-----------------------!
   !  Boundary Conditions  !
   !-----------------------!
   ! lower boundary - already set in high order domain
   do i = Node_min_lo_z,Node_max_lo_z
      mesh(i,1)%uz = mesh(i,Node_min_lo_r)%uz
      mesh(i,1)%ux = -mesh(i,Node_min_lo_r)%ux
   enddo
  ! left boundary
	!not necessary
   !do j = Node_min_lo_r,Node_max_lo_r
   !    mesh(1,j)%ux = mesh(1,j)%ux
   !    mesh(1,j)%uz = mesh(1,j)%uz
   !enddo

   ! right boundary
   do j = Node_min_lo_r,Node_max_lo_r
        mesh(Node_end_lo_z,j)%ux = 0.
        mesh(Node_end_lo_z,j)%uz = 0.
   enddo
  ! upper boundary
   do i = Node_min_lo_z,Node_max_lo_z
        mesh(i,Node_end_lo_r)%ux = mesh(i,Node_end_lo_r)%ux
        mesh(i,Node_end_lo_r)%uz = mesh(i,Node_end_lo_r)%uz
   enddo

  end subroutine


END MODULE
