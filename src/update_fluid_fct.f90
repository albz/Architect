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

MODULE Advance_Fluid_FCT


USE my_types
USE use_my_types

IMPLICIT NONE

CONTAINS

!---! !---! !---!

   SUBROUTINE AdvanceFluid_FCT


   IMPLICIT NONE



   INTEGER :: i,j,iter,q
   !~INTEGER :: Nz,Nr,Node_min_lo_z,Node_max_lo_z,Node_min_lo_r,Node_max_lo_r,Node_end_lo_z,Node_end_lo_r
   !~INTEGER ::       Node_min_ho_z,Node_max_ho_z,Node_min_ho_r,Node_max_ho_r,Node_end_ho_z,Node_end_ho_r
   REAL(8), DIMENSION(mesh_par%Nzm,mesh_par%Nxm) :: ux,uz,ne,ux_halfDt,uz_halfDt,ux_new,uz_new,ne_new
   REAL(8), DIMENSION(mesh_par%Nzm,mesh_par%Nxm) :: Ex_f,Ez_f,Bphi_f
   REAL(8), DIMENSION(mesh_par%Nzm,mesh_par%Nxm) :: ne_ux_halfDt,ne_uz_halfDt,ne_ux_new,ne_uz_new

   !~REAL DeltaR,DeltaZ,Dt,threshold_factor
   REAL(8) threshold_factor
   REAL(8) one_over_three,one_over_six,five_over_six

   REAL(8) beta_x_halfDt,beta_z_halfDt
   REAL(8) beta_z_left,beta_z_right,beta_x_down,beta_x_up
   REAL(8) ux_temp,uz_temp
   REAL(8) q_min,q_max
   REAL(8) A_c_left,A_c_right,A_c_up,A_c_down
   !~REAL, DIMENSION(mesh_par%Nxm) :: r_mesh,Sr,Sz,Vol

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

   !~DeltaR  = mesh_par%dxm
   !~DeltaZ  = mesh_par%dzm
   !~Dt      = sim_parameters%dt*plasma%omega_p

   threshold_factor = 1e-12
   one_over_three   = 1.D0/3.D0
   one_over_six     = 1.D0/6.D0
   five_over_six    = 5.D0/6.D0

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

!~   Nz             = mesh_par%Nzm
!~   Nr             = mesh_par%Nxm
!~
!~   !--- low order boundaries ---!
!~   Node_min_lo_z    = 2
!~   Node_max_lo_z    = Nz-1
!~
!~   Node_min_lo_r    = 2
!~   Node_max_lo_r    = Nr-1
!~
!~   Node_end_lo_z    = Nz
!~   Node_end_lo_r    = Nr
!~
!~   !--- high order boundaries ---!
!~   Node_min_ho_z    = 4
!~   Node_max_ho_z    = (Nz-1)-1
!~
!~   Node_min_ho_r    = 3
!~   Node_max_ho_r    = (Nr-1)-1
!~
!~   Node_end_ho_z    = Nz
!~   Node_end_ho_r    = Nr
!~
!~
!~   !--- --- ---!
!~   r_mesh(1)     = -DeltaR/2.
!~   r_mesh(2)     =  DeltaR/2.
!~   Sr    (1)     = 0. !2.*pi * DeltaZ * DeltaR
!~   Sr    (2)     = 0. !2.*pi * DeltaZ * (DeltaR/2.)
!~   Sz    (1)     = pi             * DeltaR**2
!~   Sz    (2)     = pi             * DeltaR**2
!~   Vol   (1)     = pi    * DeltaZ * DeltaR**2
!~   Vol   (2)     = pi    * DeltaZ * DeltaR**2
!~
!~   do j=(Node_min_lo_r+1),Node_end_lo_r
!~            r_mesh(j)   = DeltaR*(j-2)+DeltaR/2.
!~            Sr    (j)   = 2.*pi * DeltaZ * (r_mesh(j)-DeltaR/2.) !half cell down
!~            Sz    (j)   = 2.*pi *          DeltaR * r_mesh(j)
!~            Vol   (j)   = 2.*pi * DeltaZ * DeltaR * r_mesh(j)
!~   enddo

   !-----------------------!
   !   Initial Conditions  !
   !-----------------------!
   ux               =	mesh(:,:)%ux
   uz               =	mesh(:,:)%uz
   ne               =	mesh(:,:)%n_plasma_e

   ! Initialize auxiliary vectors and reals


   ne_ux_halfDt     = 0.
   ne_uz_halfDt     = 0.

   ux_new           = 0.
   uz_new           = 0.

   ne_new	    	    = 0.

   beta_x_halfDt    = 0.
   beta_z_halfDt    = 0.

   Ez_f		    	    = 0.
   Ex_f	            = 0.
   Bphi_f	    	    = 0.

!------------------------ Advection FCT -----------------------------!


	!-------------------------------------------!
	! Advection begins here 			        !
	!-------------------------------------------!

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


  !to be structured to be converted into a function
do q=1,3
	quantity    = 0.
	quantity_td = 0.
	if (q.eq.1) quantity = ne
	if (q.eq.2) quantity = ne*ux
	if (q.eq.3) quantity = ne*uz

	!--- FCT-fluxes
	!--- low-order
  flux_lo_z=0.D0
	flux_lo_r=0.D0
	!--- donor-cell
	do i= Node_min_lo_z,Node_max_lo_z
	do j= Node_min_lo_r,Node_max_lo_r
    flux_lo_z (i,j) = quantity(i-1,j  ) * max(0.D0,beta_z_i_minus_halfDz(i,j)) + quantity(i  ,j  ) * min(0.D0,beta_z_i_minus_halfDz(i,j))
    flux_lo_r (i,j) = quantity(i  ,j-1) * max(0.D0,beta_x_j_minus_halfDr(i,j)) + quantity(i  ,j  ) * min(0.D0,beta_x_j_minus_halfDr(i,j))
	end do
	end do
	!--- right border ---!
	flux_lo_z (Node_max_lo_z+1,Node_min_lo_r:Node_max_lo_r) = flux_lo_z (Node_max_lo_z,Node_min_lo_r:Node_max_lo_r)
	!--- upper border ---!
	flux_lo_r (Node_min_lo_z:Node_max_lo_z,Node_max_lo_r+1) = flux_lo_r (Node_min_lo_z:Node_max_lo_z,Node_max_lo_r)

	! right border
	do j= Node_min_lo_r,Node_max_lo_r
		flux_lo_z (Node_max_lo_z+1,j) = flux_lo_z (Node_max_lo_z,j)
	end do

	! upper border
	do i= Node_min_lo_z,Node_max_lo_z
		flux_lo_r (i,Node_max_lo_r+1) = flux_lo_r (i,Node_max_lo_r)
	end do

  flux_ho_z=0.D0
	flux_ho_r=0.D0
	!--- high-order
	do i= Node_min_ho_z,Node_max_ho_z
	do j= Node_min_ho_r,Node_max_ho_r
    !TUCKMANTEL STENCIL
    !flux_ho_z (i,j)  = ( -one_over_six *quantity(i-2,j  ) + five_over_six*quantity(i-1,j  ) + one_over_three*quantity(i  ,j  )) &
    !                   * max(0.D0,beta_z_i_minus_halfDz(i,j)) + &   ! if flux_z is positive 
    !                   ( one_over_three*quantity(i-1,j  ) + five_over_six*quantity(i  ,j  ) - one_over_six  *quantity(i+1,j  )) & ! 
    !                   * min(0.D0,beta_z_i_minus_halfDz(i,j)) ! if flux_z is negative
    !flux_ho_r (i,j)  = ( -one_over_six *quantity(i  ,j-2) + five_over_six*quantity(i  ,j-1) + one_over_three*quantity(i  ,j  )) & ! if flux_r is positive 
    !                   * max(0.D0,beta_x_j_minus_halfDr(i,j)) + &
    !                   ( one_over_three*quantity(i  ,j-1) + five_over_six*quantity(i  ,j  ) - one_over_six  *quantity(i  ,j+1)) & ! if flux_r is negative
    !                   * min(0.D0,beta_x_j_minus_halfDr(i,j))
    !ZALESAK STENCIL
    flux_ho_z (i,j)  = ( -1.D0/12.D0*(quantity(i-2,j  )+quantity(i+1,j  )) + 7.D0/12.D0*(quantity(i-1,j  ) + quantity(i  ,j  ))) &
                       * max(0.D0,beta_z_i_minus_halfDz(i,j)) + &   ! if flux_z is positive 
                       ( -1.D0/12.D0*(quantity(i-1,j  )+quantity(i+2,j  )) + 7.D0/12.D0*(quantity(i  ,j  ) + quantity(i+1,j  ))) & ! 
                       * min(0.D0,beta_z_i_minus_halfDz(i,j)) ! if flux_z is negative
    flux_ho_r (i,j)  = ( -1.D0/12.D0*(quantity(i  ,j-2)+quantity(i  ,j+1)) + 7.D0/12.D0*(quantity(i  ,j-1) + quantity(i  ,j  ))) & ! if flux_r is positive 
                       * max(0.D0,beta_x_j_minus_halfDr(i,j)) + &
                       ( -1.D0/12.D0*(quantity(i  ,j-1)+quantity(i  ,j+2)) + 7.D0/12.D0*(quantity(i  ,j  ) + quantity(i  ,j+1))) & ! ! if flux_r is negative
                       * min(0.D0,beta_x_j_minus_halfDr(i,j))


		!flux_ho_z (i,j) =  0.58333333*(quantity(i  ,j  )+quantity(i-1,j  ))- 0.08333333*(quantity(i+1,j  )+quantity(i-2,j  ))
		!flux_ho_r (i,j) =  0.58333333*(quantity(i  ,j  )+quantity(i  ,j-1))- 0.08333333*(quantity(i  ,j+1)+quantity(i  ,j-2))

		!--- correction ---!
		! antidiffusive fluxes
		A_left (i,j) = flux_ho_z (i,j) - flux_lo_z (i,j)
		A_down (i,j) = flux_ho_r (i,j) - flux_lo_r (i,j)
	end do
	end do

  !--- compute the low order solution
  do i= Node_min_lo_z,Node_max_lo_z
	do j= Node_min_lo_r,Node_max_lo_r
		quantity_td(i,j)= quantity(i,j)    - Dt/Vol(j) * (Sr(j+1)*flux_lo_r(i,j+1)-Sr(j)*flux_lo_r(i,j))
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

		! upper boundary
   		do i = Node_min_lo_z,Node_max_lo_z
        		quantity_td(i,Node_end_lo_r) = quantity_td(i,Node_max_lo_r)
   		enddo

   		! lower boundary
   		do i = Node_min_lo_z,Node_max_lo_z
        		quantity_td(i,1)             = quantity_td(i,Node_min_lo_r)
   		enddo
   		! left boundary
   		do j = Node_min_lo_r,Node_max_lo_r
			quantity_td(1,j)                   = quantity_td(Node_min_lo_z,j)
   		enddo
   		! right boundary
!   		do j = Node_min_lo_r,Node_max_lo_r
!        	quantity_td (Node_end_lo_z,j)     = quantity_td(Node_max_lo_z,j)
!   		enddo
	else if (q.eq.2) then

		! upper boundary
		! High order domain
   		do i = Node_min_lo_z,Node_max_lo_z
        		quantity_td(i,Node_end_lo_r)  =  quantity_td(i,Node_max_lo_r)
   		enddo

   		! lower boundary
   		do i = Node_min_lo_z,Node_max_lo_z
        		quantity_td(i,1         )     = -quantity_td(i,Node_min_lo_r)
   		enddo
   		! left boundary
   		do j = Node_min_lo_r,Node_max_lo_r
        		quantity_td(1,j)              =  0.D0
   		enddo
   		! right boundary
!   		do j = Node_min_lo_r,Node_max_lo_r
!        		quantity_td(Node_end_lo_z,j)  =  quantity_td(Node_max_lo_z,j)
!   		enddo
	else if (q.eq.3) then

		! upper boundary
		do i = Node_min_lo_z,Node_max_lo_z
        		quantity_td(i,Node_end_lo_r)  =  quantity_td(i,Node_max_lo_r)
   		enddo

   		! lower boundary
   		do i = Node_min_lo_z,Node_max_lo_z
        		quantity_td(i,1         )     =  quantity_td(i,Node_min_lo_r)
   		enddo
   		! left boundary
   		do j = Node_min_lo_r,Node_max_lo_r
        		quantity_td(1,j)              =  0.D0
   		enddo
   		! right boundary
!   		do j = Node_min_lo_r,Node_max_lo_r
!        		quantity_td(Node_end_lo_z,j)  = quantity_td(Node_max_lo_z,j)
!   		enddo
	endif



	!---------------------------------------------!
	!            Correction by FCT scheme         !
	!---------------------------------------------!

	A_left   = 0.D0
	A_down   = 0.D0
  C_left   = 0.D0
	C_right  = 0.D0
	C_up     = 0.D0
	C_down   = 0.D0
  q__plus  = 0.D0
	q__minus = 0.D0
  Q_plus   = 0.D0
	Q_minus  = 0.D0
	P_plus   = 0.D0
	P_minus  = 0.D0
	R_plus   = 0.D0
	R_minus  = 0.D0


	! practical correction
	! errors might be nested in here, pay attention
  do i= Node_min_ho_z,Node_max_ho_z
	do j= Node_min_ho_r,Node_max_ho_r
!	    if( A_left(i+1,j) * ( quantity_td(i+1,j)-quantity_td(i,j) ) < 0. .and. &
!	        A_left(i+1,j) * ( quantity_td(i+2,j)-quantity_td(i+1,j) ) < 0. .or. &
!	        A_left(i+1,j) * ( quantity_td(i,j)-quantity_td(i-1,j) ) < 0. ) A_left(i+1,j) = 0.
!
!	    if( A_down(i,j+1) * ( quantity_td(i,j+1)-quantity_td(i,j) ) < 0. .and. &
!	        A_down(i,j+1) * ( quantity_td(i,j+2)-quantity_td(i,j+1) ) < 0. .or. &
!	        A_down(i,j+1) * ( quantity_td(i,j)-quantity_td(i,j-1) ) < 0. ) A_down(i,j+1) = 0.
	    if( A_left(i,j) * ( quantity_td(i,j)-quantity_td(i-1,j  ) ) .le. 0.D0 ) A_left(i,j) = 0.D0
	    if( A_down(i,j) * ( quantity_td(i,j)-quantity_td(i  ,j-1) ) .le. 0.D0 ) A_down(i,j) = 0.D0
	enddo
	enddo

    do i= Node_min_ho_z,Node_max_ho_z
	do j= Node_min_ho_r,Node_max_ho_r
		P_plus (i,j) = max(A_left (i  ,j),0.D0)-min(A_left(i+1,j),0.D0)+max(A_down(i,j  ),0.D0)-min(A_down(i,j+1),0.D0)
		P_minus(i,j) = max(A_left (i+1,j),0.D0)-min(A_left(i  ,j),0.D0)+max(A_down(i,j+1),0.D0)-min(A_down(i,j  ),0.D0)

		q__plus (i,j) = max(quantity(i,j),quantity_td(i,j))
		q__minus(i,j) = min(quantity(i,j),quantity_td(i,j))
	end do
	end do

	do i= Node_min_ho_z,Node_max_ho_z
	do j= Node_min_ho_r,Node_max_ho_r

		!q_max = max(q__plus (i-1,j),q__plus (i,j),q__plus (i+1,j  ))
		!q_min = min(q__minus(i-1,j),q__minus(i,j),q__minus(i+1,j  ))
		q_max = max(q__plus (i-1,j),q__plus (i,j),q__plus (i+1,j  ),q__plus (i,j-1),q__plus (i,j+1))
		q_min = min(q__minus(i-1,j),q__minus(i,j),q__minus(i+1,j  ),q__minus(i,j-1),q__minus(i,j+1))

		Q_plus  (i,j) = (q_max           -quantity_td(i,j) )*Vol(j)
		Q_minus (i,j) = (quantity_td(i,j)-q_min            )*Vol(j)

		R_plus (i,j) = 0.D0
		R_minus(i,j) = 0.D0
		if (P_plus (i,j) >0.D0 ) R_plus (i,j) = min( 1.D0,(Q_plus (i,j)/P_plus (i,j)) )
		if (P_minus(i,j) >0.D0 ) R_minus(i,j) = min( 1.D0,(Q_minus(i,j)/P_minus(i,j)) )
	end do
	end do

	do i= Node_min_ho_z,Node_max_ho_z
	do j= Node_min_ho_r,Node_max_ho_r
		if(A_left(i,j).le.0.D0) then
			C_left(i,j) = min(R_plus(i-1,j  ),R_minus(i  ,j))
		else
			C_left(i,j) = min(R_plus(i  ,j  ),R_minus(i-1,j))
		endif

		if(A_down(i,j).le.0.) then
			C_down(i,j) = min(R_plus(i  ,j-1),R_minus(i,j  ))
		else
			C_down(i,j) = min(R_plus(i  ,j  ),R_minus(i,j-1))
		endif
	end do
	end do

	do i= Node_min_ho_z,Node_max_ho_z-1
	do j= Node_min_ho_r,Node_max_ho_r-1
		A_c_left       = C_left (i  ,j  )*A_left (i  ,j  ) * Sz(j  )
		A_c_right      = C_left (i+1,j  )*A_left (i+1,j  ) * Sz(j  )
		A_c_down       = C_down (i  ,j  )*A_down (i  ,j  ) * Sr(j  )
		A_c_up         = C_down (i  ,j+1)*A_down (i  ,j+1) * Sr(j+1)

		quantity(i,j) = quantity_td(i,j) - Dt/Vol(j) * &
		               (A_c_right-A_c_left+A_c_up-A_c_down)
	end do
	end do


	! low order scheme on axis
	do i= Node_min_ho_z,Node_max_ho_z-1
		quantity(i,2) = quantity_td(i,2)
		quantity(i,3) = quantity_td(i,3)
	end do
	! low order scheme on the right (not ghost cell)
	do j = Node_min_lo_r,Node_max_lo_r
    quantity (Node_max_ho_z,j)     = quantity_td(Node_max_ho_z,j)
		quantity (Node_max_lo_z,j)     = quantity_td(Node_max_lo_z,j)
 	enddo






	!---------------------------------------------!
	!  Boundary conditions for FCT                !
	!---------------------------------------------!
	if (q.eq.1) then

		! upper boundary
		! High order domain
		do i = Node_min_ho_z,Node_max_ho_z
        		quantity(i,Node_max_ho_r)   = quantity(i,Node_max_ho_r-1)
   		enddo

   		do i = Node_min_ho_z,Node_max_ho_z
        		quantity(i,Node_max_ho_r+1) = quantity(i,Node_max_ho_r)
   		enddo
		! Low order domain
   		do i = Node_min_lo_z,Node_max_lo_z
        		quantity(i,Node_end_lo_r)   = quantity(i,Node_max_lo_r)
   		enddo

   		! lower boundary
   		do i = Node_min_lo_z,Node_max_lo_z
        		quantity(i,1)               = quantity(i,Node_min_lo_r)
   		enddo
   		! left boundary
   		do j = Node_min_lo_r,Node_max_lo_r
			quantity(1,j)                     = quantity(Node_min_lo_z,j)
   		enddo
   		! right boundary
!   		do j = Node_min_lo_r,Node_max_lo_r
!        	quantity (Node_end_lo_z,j)     = quantity(Node_max_lo_z,j)
!   		enddo
	else if (q.eq.2) then

		! upper boundary
		! High order domain
		do i = Node_min_ho_z,Node_max_ho_z
        		quantity(i,Node_max_ho_r)    = quantity(i,Node_max_ho_r-1)
   		enddo

   		do i = Node_min_ho_z,Node_max_ho_z
        		quantity(i,Node_max_ho_r+1)  =  quantity(i,Node_max_ho_r)
   		enddo
		! Low order domain
   		do i = Node_min_lo_z,Node_max_lo_z
        		quantity(i,Node_end_lo_r)    =  quantity(i,Node_max_lo_r)
   		enddo


   		! lower boundary
   		do i = Node_min_lo_z,Node_max_lo_z
        		quantity(i,1         )       = -quantity(i,Node_min_lo_r)
   		enddo
   		! left boundary
   		do j = Node_min_lo_r,Node_max_lo_r
        		quantity(1,j)                =  0.D0
   		enddo
   		! right boundary
!   		do j = Node_min_lo_r,Node_max_lo_r
!        		quantity(Node_end_lo_z,j)  =  quantity(Node_max_lo_z,j)
!   		enddo
	else if (q.eq.3) then

		! upper boundary
   		! High order domain
		do i = Node_min_ho_z,Node_max_ho_z
        		quantity(i,Node_max_ho_r)    = quantity(i,Node_max_ho_r-1)
   		enddo

		do i = Node_min_ho_z,Node_max_ho_z
        		quantity(i,Node_max_ho_r+1)  =  quantity(i,Node_max_ho_r)
   		enddo
		! Low order domain
		do i = Node_min_lo_z,Node_max_lo_z
        		quantity(i,Node_end_lo_r)    =  quantity(i,Node_max_lo_r)
   		enddo


   		! lower boundary
   		do i = Node_min_lo_z,Node_max_lo_z
        		quantity(i,1         )       =  quantity(i,Node_min_lo_r)
   		enddo
   		! left boundary
   		do j = Node_min_lo_r,Node_max_lo_r
        		quantity(1,j)                =  0.D0
   		enddo
   		! right boundary
!   		do j = Node_min_lo_r,Node_max_lo_r
!        		quantity(Node_end_lo_z,j)  = quantity(Node_max_lo_z,j)
!   		enddo
	endif

    !   ------- backward substitution
	if (q.eq.1) ne_new       = quantity
	if (q.eq.2) ne_ux_halfDt = quantity
	if (q.eq.3) ne_uz_halfDt = quantity

enddo


!-------------------------- Electromagnetic Advance ------------------------------

        do i= Node_min_lo_z,Node_max_lo_z
        do j= Node_min_lo_r,Node_max_lo_r
      Ez_f   (i,j) =        mesh(i,j  )%Ez               + mesh(i  ,j  )%Ez_bunch         ! properly centered in space
      Ex_f   (i,j) = 0.25D0*( mesh(i,j  )%Ex             + mesh(i+1,j  )%Ex        &
                          +   mesh(i,j-1)%Ex             + mesh(i+1,j-1)%Ex      ) &      ! properly centered in space
                    +0.25D0*( mesh(i,j  )%Ex_bunch       + mesh(i+1,j  )%Ex_bunch  &
                            + mesh(i,j-1)%Ex_bunch       + mesh(i+1,j-1)%Ex_bunch)        ! properly centered in space
      Bphi_f (i,j) = 0.25D0*( mesh(i,j  )%Bphi           + mesh(i  ,j-1)%Bphi       &
                            + mesh(i,j  )%Bphi_old       + mesh(i  ,j-1)%Bphi_old ) &     ! properly centered in time
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
	beta_x_halfDt = ux_temp/ sqrt(1.+ux_temp**2+uz_temp**2)
	beta_z_halfDt = uz_temp/ sqrt(1.+ux_temp**2+uz_temp**2)

			endif

    ux_new       (i,j) = ux_temp + Dt*( Ex_f(i,j) - beta_z_halfDt*Bphi_f(i,j) )
  	uz_new       (i,j) = uz_temp + Dt*( Ez_f(i,j) + beta_x_halfDt*Bphi_f(i,j) )

			if ((ux_new(i,j).ne.ux_new(i,j)).or.(ux_new(i,j).ne.ux_new(i,j))) then
				write(*,*) 'bla'
				write(*,*) i,j
				write(*,*) 'n_new',ne_new(i,j)
				write(*,*) 'ux',ux(i,j)
				write(*,*) 'uz',uz(i,j)
				write(*,*) 'n_ux_new',ux_new(i,j)
				write(*,*) 'n_uz_new',ux_new(i,j)
				write(*,*) 'ux_temp',ux_temp
				write(*,*) 'uz_temp',ux_temp

				write(*,*) 'beta_x_halfDt',beta_x_halfDt
				write(*,*) 'beta_z_halfDt',beta_z_halfDt
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
   do j = Node_min_lo_r,Node_max_lo_r
		ux_new(1,j)             = 0.D0
    uz_new(1,j)             = 0.D0
   enddo

   ! right boundary
!   do j = Node_min_lo_r,Node_max_lo_r
!        ux_new(Node_end_lo_z,j) = ux_new(Node_end_ho_z,j)
!        uz_new(Node_end_lo_z,j) = uz_new(Node_end_ho_z,j)
!   enddo



    ! upper boundary
   do i = Node_min_lo_z,Node_max_lo_z
        ux_new(i,Node_end_lo_r) = ux_new(i,Node_end_lo_r)
        uz_new(i,Node_end_lo_r) = uz_new(i,Node_end_lo_r)
   enddo


   !--------------------------------------------!
   ! Substitution of the new fields in the mesh !
   !--------------------------------------------!

   mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%n_plasma_e =      ne_new     (1:Node_end_lo_z,1:Node_end_lo_r)
   mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%ux         =      ux_new     (1:Node_end_lo_z,1:Node_end_lo_r)
   mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%uz         =      uz_new     (1:Node_end_lo_z,1:Node_end_lo_r)



   END SUBROUTINE

END MODULE
