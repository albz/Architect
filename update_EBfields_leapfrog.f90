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


MODULE Fields_FDTD


USE my_types
USE use_my_types

IMPLICIT NONE

CONTAINS


   SUBROUTINE Kernel_Fields_FDTD_COMB ! The new subroutine

   IMPLICIT NONE

   INTEGER i,iter,j



   REAL(8), DIMENSION(mesh_par%Nzm,mesh_par%Nxm) :: Ez,Er,Bphi,Ez_new,Er_new,Bphi_new,Jr,Jz,Jpe_r,Jpe_z
   !~INTEGER Nz,Nr,Node_min_z,Node_max_z,Node_min_r,Node_max_r,Node_end_z,Node_end_r
   REAL(8), DIMENSION(mesh_par%Nxm) :: r_mesh_EM
   REAL(8) :: r_factor1,r_factor2
   REAL(8) :: threshold_factor=1e-3
   !~REAL :: DeltaR,DeltaZ,r_factor1,r_factor2,Dt

	!~Nz         = mesh_par%Nzm
	!~Nr         = mesh_par%Nxm

!~	Node_min_z = 2
!~	Node_max_z = Nz-1
!~
!~	Node_min_r = 2
!~	Node_max_r = Nr-1
!~
!~	Node_end_z = Nz
!~	Node_end_r = Nr

	!~! Adimensional integration step and cell sizes
	!~DeltaR     = mesh_par%dxm
	!~DeltaZ     = mesh_par%dzm
	!~Dt         = sim_parameters%dt*plasma%omega_p

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




	!----------------------!
	!    Mesh creation
	!----------------------!


	r_factor1     		= 0.
	r_factor2     		= 0.
	r_mesh_EM        	= 0.
	do j=1,Node_max_r
            r_mesh_EM(j) =  DeltaR*(j-1)
	enddo

	!----------------------!
	!  Initial conditions
	!----------------------!

	Er   = 0.
	Ez   = 0.
	Bphi = 0.

	Er_new   = 0.
	Ez_new   = 0.
	Bphi_new = 0.

	Bphi(:,:)  =	mesh(:,:)%Bphi
	Ez  (:,:)  =	mesh(:,:)%Ez
	Er  (:,:)  =	mesh(:,:)%Ex


	!-----------------------!
	!  Boundary Conditions  !
	!-----------------------!

	! upper boundary
	do i = Node_min_z,Node_max_z
		Bphi(i         ,Node_end_r) = Bphi(i         ,Node_max_r)
!     	Bphi_bunch(i   ,Node_end_r) = Bphi_bunch(i   ,Node_max_r) HERE HERE HERE
	enddo

	! lower boundary
	do i = Node_min_z,Node_max_z
		Bphi(i         ,1         ) = 0.
! 		Bphi_bunch(i   ,1         ) = 0.
	enddo

	! left boundary
	do j = Node_min_r,Node_max_r
		Bphi(1         ,j         ) = Bphi(Node_min_z,j         )
! 		Bphi_bunch(1         ,j         ) = Bphi_bunch(Node_min_z,j         ) HERE HERE HERE
	enddo

	! right boundary
	do j = Node_min_r,Node_max_r
		Bphi(Node_end_z,j         ) = Bphi(Node_max_z,j         )
! 		Bphi_bunch(Node_end_z,j   ) = Bphi_bunch(Node_max_z,j   ) HERE HERE HERE
	enddo



	!------------------------!
	!        Sources          !
	!------------------------!

	Jz       =   0.
	Jr       =   0.

	Jpe_z	 =   0.
	Jpe_r    =   0.

	! Reads the currents in physical space




!~    ! beam current
!~	Jz   (Node_min_z:Node_max_z,Node_min_r:Node_max_r)  =   mesh(2:(mesh_par%Nzm-1),2:(mesh_par%Nxm-1))%Jz
!~	Jr   (Node_min_z:Node_max_z,Node_min_r:Node_max_r)  =   mesh(2:(mesh_par%Nzm-1),2:(mesh_par%Nxm-1))%Jr
!~
!~	! plasma electron current
!~    Jpe_z(Node_min_z:Node_max_z,Node_min_r:Node_max_r)  =   mesh(2:(mesh_par%Nzm-1),2:(mesh_par%Nxm-1))%Jpe_z
!~	Jpe_r(Node_min_z:Node_max_z,Node_min_r:Node_max_r)  =   mesh(2:(mesh_par%Nzm-1),2:(mesh_par%Nxm-1))%Jpe_r


	! beam current
	Jz   (1:Node_end_z,1:Node_end_r)  =   mesh(1:Node_end_z,1:Node_end_r)%Jz
	Jr   (1:Node_end_z,1:Node_end_r)  =   mesh(1:Node_end_z,1:Node_end_r)%Jr

	! plasma electron current
    Jpe_z(1:Node_end_z,1:Node_end_r)  =   mesh(1:Node_end_z,1:Node_end_r)%Jpe_z
	Jpe_r(1:Node_end_z,1:Node_end_r)  =   mesh(1:Node_end_z,1:Node_end_r)%Jpe_r



	!~if (sim_parameters%iter.ne.1) then

	!------------------------!
	!     E Fields Advance
	!------------------------!

	do i= Node_min_z,Node_max_z
    do j= Node_min_r,Node_max_r

      Ez_new (i,j) = Ez (i,j) + Dt/DeltaR*( Bphi(i,j)*r_mesh_EM(j) - Bphi(i,j-1)*r_mesh_EM(j-1) )*2./(r_mesh_EM(j)+r_mesh_EM(j-1)) - Dt * Jpe_z(i,j)
      Er_new (i,j) = Er (i,j) - Dt/DeltaZ*( Bphi  (i,j)            - Bphi(i-1,j  )              )                                  - Dt * Jpe_r(i,j)

    enddo
	enddo

	!-----------------------!
	!  Boundary conditions  !
	!         for E         !
	!-----------------------!

	! upper boundary
    do i=Node_min_z,Node_max_z
		Er_new(i         ,Node_end_r) = Er_new(i,Node_max_r)
		Ez_new(i         ,Node_end_r) = Ez_new(i,Node_max_r)
    enddo

	! lower boundary
    do i=Node_min_z,Node_max_z
		Er_new(i         ,1         ) = 0.
		Ez_new(i         ,1         ) = Ez_new(i,Node_min_r)
!         j=2
!         Ez_new(i         ,j         ) = Ez (i,j) + Dt/DeltaR*( Bphi(i,j) )*8. - Dt * Jpe_z(i,j)
    enddo

	! left boundary
    do j=Node_min_r,Node_max_r
		Er_new(1         ,j         ) = Er_new(Node_min_z,j)
		Ez_new(1         ,j         ) = Ez_new(Node_min_z,j)
    enddo

	! right boundary
    do j=Node_min_r,Node_max_r
  	Er_new(Node_end_z,j         ) = Er_new(Node_max_z,j)
  	Ez_new(Node_end_z,j         ) = Ez_new(Node_max_z,j)
    enddo


	!------------------------!
	!     B Field Advance
	!------------------------!

	do i= Node_min_z,Node_max_z
    do j= Node_min_r,Node_max_r

      Bphi_new(i,j) = Bphi(i,j) + Dt/DeltaR*( Ez_new(i,j+1) - Ez_new(i,j)  ) - Dt/DeltaZ*( Er_new(i+1,j) - Er_new(i,j) )

    enddo
	enddo

	!-----------------------!
	!  Boundary conditions  !
	!         for B         !
	!-----------------------!

	! upper boundary
    do i=Node_min_z,Node_max_z
      Bphi_new(i         ,Node_end_r) = Bphi_new(i,Node_max_r)
    enddo

	! lower boundary
    do i=Node_min_z,Node_max_z
      Bphi_new(i         ,1         ) = 0.
    enddo

	! left boundary
    do j=Node_min_r,Node_max_r
      Bphi_new(1         ,j         ) = Bphi_new(Node_min_z,j)
    enddo


	! right boundary
    do j=Node_min_r,Node_max_r
      Bphi_new(Node_end_z,j         ) = Bphi_new(Node_max_z,j)
    enddo

	!~endif






	!--------------------------------!
	! Substitution of the new fields !
	!--------------------------------!

	!~if (sim_parameters%iter.eq.1) then
	!~	Bphi_new(1:Node_end_z,1:Node_end_r)		= Bphi(1:Node_end_z,1:Node_end_r)
	!~	Ez_new  (1:Node_end_z,1:Node_end_r)		= Ez  (1:Node_end_z,1:Node_end_r)
	!~	Er_new  (1:Node_end_z,1:Node_end_r)		= Er  (1:Node_end_z,1:Node_end_r)
	!~endif



	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%Bphi_old =   Bphi(1:Node_end_z,1:Node_end_r) ! Needed to center Bphi in time in Boris Pusher

	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%Bphi     =   Bphi_new(1:Node_end_z,1:Node_end_r)
	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%Ez       =   Ez_new  (1:Node_end_z,1:Node_end_r)
	mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%Ex       =   Er_new  (1:Node_end_z,1:Node_end_r)



  !--->
!   mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%Bphi_old = 0.0
!   mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%Bphi     = 0.0
!   mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%Ez       = 0.0
!   mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)%Ex       = 0.0






! Filtering fields:
!   if(filter) then

!      mesh%Ez=matmul(matmul(Fl,mesh%Ez),Fr)

!      mesh%Ex=matmul(matmul(Fl,mesh%Ex),Fr)

!      mesh%Bphi=matmul(matmul(Fl,mesh%Bphi),Fr)

!   endif

   return

   END SUBROUTINE


END MODULE
