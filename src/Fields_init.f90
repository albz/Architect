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


MODULE init_fields
! initializes electromagnetic fields in cold approximation (all the bunches must have the same average initial gammma)
USE my_types
USE use_my_types
USE linear_algebra
USE Compute_beam_current_FDTD
USE moments

INTEGER :: dim_Lapl_r, dim_Lapl_half_z, dim_Lapl
real(8), dimension(:,:), allocatable :: a0,a1,a2,a3,a4,rhoR
!implicit none

contains

	subroutine init_EM_fields
	! all the bunches must have the same average initial gammma

	REAL(8), DIMENSION(:,:), ALLOCATABLE:: Phi_matrix,Laplacian_matrix,Phi_whole_domain,Laplacian_matrix_sparse
	REAL(8), DIMENSION(:), ALLOCATABLE::	Laplacian_matrix_sparse_vector
	INTEGER, DIMENSION(:), ALLOCATABLE::Laplacian_matrix_sparse_vector_row,Laplacian_matrix_sparse_vector_column
	REAL(8), DIMENSION(:), ALLOCATABLE  :: rho_vector,Phi_vector
	INTEGER, DIMENSION(7) :: bunches_position
	INTEGER:: b,dim,i,j,ierr
	INTEGER :: left,right,bottom,up
	INTEGER :: dim_Laplacian, count_non_null_elements
	REAL(8), DIMENSION(3,3) :: A
	REAL(8), DIMENSION(3) :: c
	REAL(8) :: distance


		dim_Lapl_r			=	bunch_initialization%init_width_r*mesh_par%Nsample_r/2
		dim_Lapl_r			= min(dim_Lapl_r,Node_max_r) !force maximum matrix dimension
		dim_Lapl_half_z	=	bunch_initialization%init_width_z*mesh_par%Nsample_z
		dim_Lapl			  = 2*dim_Lapl_half_z*dim_Lapl_r
		sim_parameters%dim_Laplacian = dim_Lapl

		!--- allocate ---!
		allocate(a0(mesh_par%Nzm,mesh_par%Nxm))
		allocate(a1(mesh_par%Nzm,mesh_par%Nxm))
		allocate(a2(mesh_par%Nzm,mesh_par%Nxm))
		allocate(a3(mesh_par%Nzm,mesh_par%Nxm))
		allocate(a4(mesh_par%Nzm,mesh_par%Nxm))
		allocate(rhoR(mesh_par%Nzm,mesh_par%Nxm))
		!--- *** ---!


		if(bunch_initialization%self_consistent_field_bunch==2) then
			if(.not.allocated(Laplacian_matrix)) ALLOCATE(Laplacian_matrix( dim_Lapl, dim_Lapl 	) )
		else if(bunch_initialization%self_consistent_field_bunch==3) then
			if(.not.allocated(Laplacian_matrix_sparse_vector)) ALLOCATE(Laplacian_matrix_sparse_vector			( 5*dim_Lapl ) )
			if(.not.allocated(Laplacian_matrix_sparse_vector_row)) ALLOCATE(Laplacian_matrix_sparse_vector_row		( 5*dim_Lapl ) )
			if(.not.allocated(Laplacian_matrix_sparse_vector_column)) ALLOCATE(Laplacian_matrix_sparse_vector_column	( 5*dim_Lapl ) )
		else if(bunch_initialization%self_consistent_field_bunch==4) then
			dim_Lapl_r			=	bunch_initialization%init_width_r*mesh_par%Nsample_r/2
			dim_Lapl_half_z	=	bunch_initialization%init_width_z*mesh_par%Nsample_z
			dim_Lapl			  = 2*dim_Lapl_half_z*dim_Lapl_r
			sim_parameters%dim_Laplacian = dim_Lapl
			if(.not.allocated(Laplacian_matrix_sparse_vector)) ALLOCATE(Laplacian_matrix_sparse_vector			( 5*dim_Lapl ) )
			if(.not.allocated(Laplacian_matrix_sparse_vector_row)) ALLOCATE(Laplacian_matrix_sparse_vector_row		( 5*dim_Lapl ) )
			if(.not.allocated(Laplacian_matrix_sparse_vector_column)) ALLOCATE(Laplacian_matrix_sparse_vector_column	( 5*dim_Lapl ) )
		endif

		ALLOCATE(Phi_whole_domain(mesh_par%Nzm,mesh_par%Nxm))
		Phi_whole_domain=0.0

		! -- finds bunches position in order to know where to cut the partial domain
		do b=1,bunch_initialization%n_total_bunches
			bunches_position(b) = 1 + abs(mesh_par%z_min)/mesh_par%dzm
			! distance = calculate_nth_moment_bunch(b,1,3)-calculate_nth_moment_bunch(1,1,3)
			bunch(b)%part(:)%cmp(14)=1.
			distance = calculate_nth_moment(b,1,3,'nocentral') &
			         - calculate_nth_moment(1,1,3,'nocentral')
			distance = distance/mesh_par%dzm*plasma%k_p
			bunches_position(b) = bunches_position(b) + int(distance)
		enddo



		write(*,'(A)') 'Initializing EM fields:'

		do b=1,bunch_initialization%n_total_bunches
			write(*,'(A,I2)') 'Computing Potential for bunch:',b


			!---project charge on the grid: bunch 'b'
			call Compute_beam_3current_bunch_FDTD(b)
			!call analytic_rho_only_first_bunch

			if(bunch_initialization%self_consistent_field_bunch==2) then
				! ---- Laplacian matrix for the Poisson problem
				Laplacian_matrix(:,:)=0.
				write(*,*) 'generating Laplacian matrix'
				call Laplacian(Laplacian_matrix) !domain values
				write(*,*) 'Laplacian matrix generated'
			else if(bunch_initialization%self_consistent_field_bunch==3) then
				! ---- Laplacian matrix for the Poisson problem
				write(*,'(A)') 'generating sparse Laplacian matrix'
				call Laplacian_sparse(Laplacian_matrix_sparse_vector,Laplacian_matrix_sparse_vector_row, &
				Laplacian_matrix_sparse_vector_column,count_non_null_elements)
				write(*,'(A)') 'Sparse Laplacian matrix generated'
			else if(bunch_initialization%self_consistent_field_bunch==4) then
				! ---- Laplacian matrix for the Poisson problem
				write(*,'(A)') 'generating sparse Laplacian matrix'
				call Laplacian_sparse(Laplacian_matrix_sparse_vector,Laplacian_matrix_sparse_vector_row, &
				Laplacian_matrix_sparse_vector_column,count_non_null_elements)
				write(*,'(A)') 'Sparse Laplacian matrix generated'
			else if(bunch_initialization%self_consistent_field_bunch==5) then
				! ---- Laplacian matrix for the Poisson problem
				write(*,'(A)') 'generating sparse Laplacian matrix for the SOR cleaner version'
				call Laplacian_sparse_SOR_new(ierr)
				write(*,'(A)') 'Sparse Laplacian matrix generated'
			endif

			! ---- The part of the domain involved is around the bunch
			left  =  max(2,(bunches_position(b)-dim_Lapl_half_z))
			right =  min(mesh_par%Nzm-1, (bunches_position(b)+dim_Lapl_half_z-1))
			bottom = 2
			if(bunch_initialization%self_consistent_field_bunch==5) bottom=1
			up     = (dim_Lapl_r+1)

			if(.not.allocated(Phi_matrix)) ALLOCATE(Phi_matrix((right-left+1),(up-bottom+1)) )
			if(.not.allocated(rho_vector)) ALLOCATE(rho_vector((right-left+1)*(up-bottom+1)) )
			if(.not.allocated(Phi_vector)) ALLOCATE(Phi_vector((right-left+1)*(up-bottom+1)) )
			Phi_matrix(:,:)	= 0.
			Phi_vector(:)	= 0.
			rho_vector(:) 	= 0.



			if(bunch_initialization%self_consistent_field_bunch==3) then
					call from_M_to_v(mesh(left:right,bottom:up)%rho,rho_vector)
			else if(bunch_initialization%self_consistent_field_bunch==4) then
					call from_M_to_v2(mesh(left:right,bottom:Node_max_r)%rho,rho_vector)
			endif

			! ---- Source for Poisson Problem
			Phi_vector = -1.*rho_vector

			! ---- Solution of the Poisson problem
			if(bunch_initialization%self_consistent_field_bunch==2) then
				write(*,*) 'solving Laplacian matrix with LU technique'
				call GaussEliminationLU(Laplacian_matrix,Phi_vector)
				write(*,*) 'Laplacian matrix solved with LU'
			else if(bunch_initialization%self_consistent_field_bunch==3) then
				write(*,*) 'solving Laplacian matrix with SOR technique'
				call SOR_sparse_Matrix_sparse(Laplacian_matrix_sparse_vector,Phi_vector, &
					     Laplacian_matrix_sparse_vector_row,Laplacian_matrix_sparse_vector_column, &
					     bunch_initialization%maxnorm,bunch_initialization%iter_max, &
					     bunch_initialization%wsor,count_non_null_elements,(right-left+1)*(up-bottom+1))
				write(*,'(A)') 'Laplacian matrix solved with SOR'
			else if(bunch_initialization%self_consistent_field_bunch==4) then
				write(*,*) 'solving Laplacian matrix with SOR technique'
				call SOR_sparse_Matrix_sparse(Laplacian_matrix_sparse_vector,Phi_vector, &
					     Laplacian_matrix_sparse_vector_row,Laplacian_matrix_sparse_vector_column, &
					     bunch_initialization%maxnorm,bunch_initialization%iter_max, &
					     bunch_initialization%wsor,count_non_null_elements,(right-left+1)*(up-bottom+1))
				write(*,'(A)') 'Laplacian matrix solved with SOR'
			! else if(bunch_initialization%self_consistent_field_bunch==3) then
			! 	write(*,*) 'solving Laplacian with Conjugate-Gradient technique'
			! 	call CG_sparse(Laplacian_matrix_sparse_vector,Phi_vector, &
			! 		     Laplacian_matrix_sparse_vector_row,Laplacian_matrix_sparse_vector_column, &
			! 		     bunch_initialization%maxnorm,bunch_initialization%iter_max, &
			! 		     bunch_initialization%wsor,count_non_null_elements,(right-left+1)*(up-bottom+1))
			! 	write(*,'(A)') 'Laplacian matrix solved with CG'
		else if(bunch_initialization%self_consistent_field_bunch==5) then
			write(*,*) 'solving Laplacian matrix with the SOR technique (clean code version)'
			write(*,*) shape(a1(left:right,bottom:up))
			call sor_homo( &
			a1(left:right,bottom:up), &
			a3(left:right,bottom:up), &
			a2(left:right,bottom:up), &
			a4(left:right,bottom:up), &
			a0(left:right,bottom:up), &
			rhoR(left:right,bottom:up), & !right-hand-side
			Phi_matrix, &
			2, bunch_initialization%wsor, 0.8d0, bunch_initialization%maxnorm, bunch_initialization%iter_max, &
			100, .TRUE., .FALSE.)
			write(*,'(A)') 'Laplacian matrix solved with SOR'
			endif

			if(bunch_initialization%self_consistent_field_bunch<5) call from_v_to_M(Phi_vector,Phi_matrix)
			! ----- Put partial domain potential in its place in the whole domain, summing to the computed potentials

			Phi_whole_domain(left:right,bottom:up) = Phi_whole_domain(left:right,bottom:up) + Phi_matrix(:,:)

			write(*,'(A,I2)') 'Computed Potential for bunch:',b
		enddo


		! -- Left boundary
		Phi_whole_domain(1,:) 			= Phi_whole_domain(2,:)
		! -- Right boundary
		Phi_whole_domain(Node_end_z,:) 	= Phi_whole_domain(Node_max_z,:)
		! -- Lower boundary
		Phi_whole_domain(:,1) 			= Phi_whole_domain(:,2)
		! -- Upper boundary
		Phi_whole_domain(:,Node_end_r) 	= Phi_whole_domain(:,Node_max_r)



		call init_E_field(Phi_whole_domain)
		write(*,'(A)') 'E field initialized'

		call init_B_field
		write(*,'(A)') 'B field initialized'


		if(allocated(Laplacian_matrix))                    DEALLOCATE(Laplacian_matrix)
		if(allocated(Laplacian_matrix_sparse))             DEALLOCATE(Laplacian_matrix_sparse)
		if(allocated(Laplacian_matrix_sparse_vector))      DEALLOCATE(Laplacian_matrix_sparse_vector		)
		if(allocated(Laplacian_matrix_sparse_vector_row))  DEALLOCATE(Laplacian_matrix_sparse_vector_row	)
		if(allocated(Laplacian_matrix_sparse_vector_column))  DEALLOCATE(Laplacian_matrix_sparse_vector_column)
		if(allocated(rho_vector))                          DEALLOCATE(rho_vector)
		if(allocated(Phi_vector))                          DEALLOCATE(Phi_vector)
		if(allocated(Phi_whole_domain))                    DEALLOCATE(Phi_whole_domain)
		if(allocated(Phi_matrix))                          DEALLOCATE(Phi_matrix)

		write(*,'(A)') '--- EM fields initialized ---'
	end subroutine init_EM_fields


	subroutine init_E_field(Phi)
		! computes E field from electrostatic potential, all in laboratory frame
		! Ez = -1/gamma_0^2 d(Phi)/dz
		! Er = - d(Phi)/dr
		! all the bunches must have the same average initial gammma

		INTEGER i,j
		REAL(8) gamma_0
		REAL(8), DIMENSION(mesh_par%Nzm,mesh_par%Nxm), intent(in) :: Phi

		!--- old version 2 points stencil ---!
		!gamma_0 = bunch_initialization%bunch_gamma_m(1)
		!do i=2,(mesh_par%Nzm-2)
		!	do j=2,(mesh_par%Nxm-1)
		!		mesh(i,j)%Ex_bunch = -1.*(Phi(i,j)-Phi(i,j-1))/mesh_par%dxm
		!		mesh(i,j)%Ez_bunch = -1./gamma_0**2*(Phi(i,j)-Phi(i-1,j))/mesh_par%dzm
		!	enddo
		!enddo

		!---> using 3points stencil
		gamma_0 = bunch_initialization%bunch_gamma_m(1)
		do i=2,(mesh_par%Nzm-2)
			do j=2,(mesh_par%Nxm-1)
				mesh(i,j)%Ex_bunch = -1.*(Phi(i,j+1)-Phi(i,j-1))/2./mesh_par%dxm
				mesh(i,j)%Ez_bunch = -1./gamma_0**2*(Phi(i+1,j)-Phi(i-1,j))/2./mesh_par%dzm
			enddo
		enddo

		!--> Extra BC on axis for Ez on lower boundary (trust me it's necessary)
		!--->mesh(:,1)%Ez = mesh(:,2)%Ez
		mesh(:,1)%Ez_bunch = mesh(:,2)%Ez_bunch
		mesh(1:2,:)%Ex_bunch = 0.D0
		mesh(1:2,:)%Ez_bunch = 0.D0
		mesh(mesh_par%Nzm-2:mesh_par%Nzm,:)%Ex_bunch = 0.D0
		mesh(mesh_par%Nzm-2:mesh_par%Nzm,:)%Ez_bunch = 0.D0
		mesh(:,mesh_par%Nxm-2:mesh_par%Nxm)%Ex_bunch = 0.D0
		mesh(:,mesh_par%Nxm-2:mesh_par%Nxm)%Ez_bunch = 0.D0
		mesh(:,:)%Bphi = 0.D0
		mesh(:,:)%Bphi_bunch = 0.D0

		! ----- center E fields for FDTD
		call center_E_fields
		! ----- Boundary conditions
		call BC_Efield
		call BC_Bfield

	end subroutine init_E_field


	subroutine Laplacian(Lapl_matrix)
		REAL(8), DIMENSION(dim_Lapl, dim_Lapl ), intent(inout):: Lapl_matrix
		REAL(8), DIMENSION(2*dim_Lapl_half_z*(dim_Lapl_r+1)):: radiusj
		REAL(8) gamma_0
		INTEGER i,j,k,nn,b
		INTEGER bla

			nn=2*dim_Lapl_half_z
			gamma_0 = bunch_initialization%bunch_gamma_m(1)

			k=1
			do j=2,(dim_Lapl_r+1)
			do i=1,(2*dim_Lapl_half_z)
				radiusj(k)=x_mesh_shifted(j)
				k=k+1
			enddo
			enddo

			do k=1,dim_Lapl

				Lapl_matrix(k,k) 		= 					- 2./mesh_par%dzm**2/gamma_0**2
				Lapl_matrix(k,k) 		= Lapl_matrix(k,k) 	- 2./mesh_par%dxm**2
				Lapl_matrix(k,k-1)  	= 1./mesh_par%dzm**2/gamma_0**2

				if(k.lt.(dim_Lapl-1)) then
					Lapl_matrix(k,k+1)  = 1./mesh_par%dzm**2/gamma_0**2
				endif

				if(radiusj(k).ne.0.) then

					if (k.ge.(nn+1)) then
						Lapl_matrix(k,k-nn) 	= 1./mesh_par%dxm**2 * (1.-mesh_par%dxm/2./radiusj(k))
					endif

					if ((k+nn).le.dim_Lapl) then
						Lapl_matrix(k,k+nn) = 1./mesh_par%dxm**2 * (1.+mesh_par%dxm/2./radiusj(k))
					endif

				endif

			enddo

			! ------ Implicit boundary conditions in Laplacian

			! -- Left boundary
			k=1
			do j=2,(dim_Lapl_r+1)
			do i=1,(2*dim_Lapl_half_z)
				if(i.eq.1) then
					!------------->>> Lapl_matrix(k,k-1) = 0.
				endif
				k=k+1
			enddo
			enddo

			! -- Right boundary
			k=1
			do j=2,(dim_Lapl_r+1)
			do i=1,(2*dim_Lapl_half_z)
				if(   (i.eq.(2*dim_Lapl_half_z)) .and. &
				(k.lt.(2*dim_Lapl_half_z*(dim_Lapl_r)))    ) then
					!------------->>> Lapl_matrix(k,k+1) = 0.
				endif
				k=k+1
			enddo
			enddo

			! -- Upper boundary
			k=1
			do j=2,(dim_Lapl_r+1)
			do i=1,(2*dim_Lapl_half_z)
				if( j.eq.(dim_Lapl_r) ) then
					if ( k.lt.(dim_Lapl-nn) ) then
					!if ( k.lt.(dim_Lapl_r+1) ) then
						!------------->>> Lapl_matrix(k,k+nn) = 0.
					endif
				endif
				k=k+1
			enddo
			enddo

			! -- Lower Boundary
			k=1
			do j=2,(dim_Lapl_r+1)
			do i=1,(2*dim_Lapl_half_z)
					if (j.eq.2) then
						Lapl_matrix(k,k) 	= -(1. + mesh_par%dxm/2./radiusj(k))/mesh_par%dxm**2 - 2./mesh_par%dzm**2/gamma_0**2
						if (k.ge.(nn + 1)) then
							Lapl_matrix(k,k-nn) = 0.
						endif
						if ( k.lt. (dim_Lapl-nn   )   ) then
							Lapl_matrix(k,k+nn) =  (1. + mesh_par%dxm/2./radiusj(k))/mesh_par%dxm**2
						endif
					endif
					k=k+1
			enddo
			enddo

	end subroutine Laplacian

	subroutine Laplacian_sparse(Lapl_matrix_sparse_vector,Lapl_matrix_sparse_vector_row,Lapl_matrix_sparse_vector_column,count_non_null_elements)
		REAL(8), DIMENSION(2*dim_Lapl_half_z*(dim_Lapl_r+1)):: radiusj
		REAL(8) gamma_0
		INTEGER i,j,k,nn,b,count
		REAL(8), DIMENSION(5*dim_Lapl), intent(inout):: Lapl_matrix_sparse_vector
		INTEGER, DIMENSION(5*dim_Lapl), intent(inout):: Lapl_matrix_sparse_vector_row
		INTEGER, DIMENSION(5*dim_Lapl), intent(inout):: Lapl_matrix_sparse_vector_column
		INTEGER, intent(inout):: count_non_null_elements
		INTEGER bla

			Lapl_matrix_sparse_vector(:) 		= 0.
			Lapl_matrix_sparse_vector_row(:) 	= 0
			Lapl_matrix_sparse_vector_column(:) = 0
			count_non_null_elements = 0

			nn=2*dim_Lapl_half_z
			gamma_0 = bunch_initialization%bunch_gamma_m(1)
			write(*,'(A,I10)') 'Laplacian matrix size: n x n, with n:',dim_Lapl
			k=1
			count=1
			do j=2,(dim_Lapl_r+1)
			do i=1,(2*dim_Lapl_half_z)


				radiusj(k)=x_mesh_shifted(j)

				! ------------ !
				! Index k,k-nn !
				! ------------ !
				if(radiusj(k).ne.0.) then
					if (k.ge.(nn+1)) then
						Lapl_matrix_sparse_vector		(count)	= 1./mesh_par%dxm**2 * (1.-mesh_par%dxm/2./radiusj(k))
						Lapl_matrix_sparse_vector_row	(count)	= k
						Lapl_matrix_sparse_vector_column(count)	= k-nn

						if (j.eq.2) then ! -- Lower Boundary BC
							Lapl_matrix_sparse_vector		(count)	= 0.
							Lapl_matrix_sparse_vector_row	(count)	= k
							Lapl_matrix_sparse_vector_column(count)	= k-nn
						endif

						count = count+1
					endif
				endif

				! ------------ !
				! Index k,k-1  !
				! ------------ !
				if (k.gt.1) then
					Lapl_matrix_sparse_vector		(count) = 1./mesh_par%dzm**2/gamma_0**2
					Lapl_matrix_sparse_vector_row	(count) = k
					Lapl_matrix_sparse_vector_column(count) = k-1


					if(i.eq.1) then ! Left boundary BC
						Lapl_matrix_sparse_vector		(count) = 0.
						Lapl_matrix_sparse_vector_row	(count) = k
						Lapl_matrix_sparse_vector_column(count) = k-1
					endif

					count = count+1
				endif


				! ------------ !
				! Index k,k    !
				! ------------ !
				Lapl_matrix_sparse_vector		(count) = - 2./mesh_par%dzm**2/gamma_0**2
				Lapl_matrix_sparse_vector		(count) = Lapl_matrix_sparse_vector(count) - 2./mesh_par%dxm**2
				Lapl_matrix_sparse_vector_row	(count) = k
				Lapl_matrix_sparse_vector_column(count) = k


				if (j.eq.2) then ! -- Lower Boundary BC
					Lapl_matrix_sparse_vector		(count) = -(1. + mesh_par%dxm/2./radiusj(k))/mesh_par%dxm**2 - 2./mesh_par%dzm**2/gamma_0**2
					Lapl_matrix_sparse_vector_row	(count) = k
					Lapl_matrix_sparse_vector_column(count) = k
				endif

				count = count+1

				! ------------ !
				! Index k,k+1  !
				! ------------ !
				if(k.lt.(dim_Lapl)) then
					Lapl_matrix_sparse_vector		(count) = 1./mesh_par%dzm**2/gamma_0**2
					Lapl_matrix_sparse_vector_row	(count) = k
					Lapl_matrix_sparse_vector_column(count) = k+1

					if(   i.eq.(2*dim_Lapl_half_z)    ) then ! -- Right boundary BC
						Lapl_matrix_sparse_vector		(count) = 0.
						Lapl_matrix_sparse_vector_row	(count) = k
						Lapl_matrix_sparse_vector_column(count) = k+1
					endif

					count = count+1
				endif

				! ------------ !
				! Index k,k+nn !
				! ------------ !
				if(radiusj(k).ne.0.) then
					if (k.lt.(dim_Lapl-nn)) then
						Lapl_matrix_sparse_vector		(count) = 1./mesh_par%dxm**2 * (1.+mesh_par%dxm/2./radiusj(k))
						Lapl_matrix_sparse_vector_row	(count) = k
						Lapl_matrix_sparse_vector_column(count) = k+nn

						if( j.eq.(dim_Lapl_r) ) then ! -- Upper boundary BC
						!if ( k.lt.(dim_Lapl_r+1) ) then
							Lapl_matrix_sparse_vector		(count) = 0. !***ALBZ***
							Lapl_matrix_sparse_vector_row	(count) = k
							Lapl_matrix_sparse_vector_column(count) = k+nn
						endif

						if (j.eq.2) then ! -- Lower Boundary BC
							Lapl_matrix_sparse_vector		(count) = (1. + mesh_par%dxm/2./radiusj(k))/mesh_par%dxm**2
							Lapl_matrix_sparse_vector_row	(count) = k
							Lapl_matrix_sparse_vector_column(count) = k+nn
						endif

						count = count+1
					endif


				endif

				k=k+1


			enddo
			enddo

			count_non_null_elements=count-1

	end subroutine Laplacian_sparse


	subroutine Laplacian_sparse_SOR_new(ierr)
		integer, intent(out) :: ierr
		integer :: i,j
		real :: radius1,radius2
	! 	real(8), DIMENSION(mesh_par%Nzm,mesh_par%Nxm), intent(inout) :: a0,a1,a2,a3,a4,rhoR
  !
		!--- allocate ---!
		a0=0.d0
		a1=0.d0
		a2=0.d0
		a3=0.d0
		a4=0.d0
		rhoR=0.d0

		!--- matrix filling  ---!
		do i=2,mesh_par%Nzm
			do j=2,mesh_par%Nxm
				a1(i,j)=mesh_par%dxm/mesh_par%dzm * (j-.5)*mesh_par%dxm
				a2(i,j)=mesh_par%dzm/mesh_par%dxm * (j-1)*mesh_par%dxm
				a3(i,j)=mesh_par%dxm/mesh_par%dzm * (j-.5)*mesh_par%dxm
				a4(i,j)=mesh_par%dzm/mesh_par%dxm * (j-2)*mesh_par%dxm
				a0(i,j)=-(a1(i,j)+a2(i,j)+a3(i,j)+a4(i,j))
				radius1 = (j-1.5) * mesh_par%dxm
				rhoR(i,j)=-mesh(i,j)%rho * radius1 * mesh_par%dzm*mesh_par%dxm
			enddo
			!--- adding a ghost cell at the bottom layer as a boundary condition ---!
			j=1
			a1(i,j)=a1(i,j+1)
			a2(i,j)=a2(i,j+1)
			a3(i,j)=a3(i,j+1)
			a4(i,j)=a4(i,j+1)
			a0(i,j)=a0(i,j+1)
			rhoR(i,j)=rhoR(i,j+1)
		enddo

		!--- boundary conditions ---!
		!set directly inside the SOR method
		ierr=0
	end subroutine Laplacian_sparse_SOR_new


	subroutine from_M_to_v(rho_matrix,rho_vector)
		REAL(8), DIMENSION(2*dim_Lapl_half_z,dim_Lapl_r):: rho_matrix
		REAL(8), DIMENSION(dim_Lapl), intent(out) :: rho_vector
		INTEGER k,i,j

		k=1
		do j=1,dim_Lapl_r
		do i=1,(2*dim_Lapl_half_z)
				rho_vector(k)= rho_matrix(i,j)
				k=k+1
		enddo
		enddo

	end subroutine from_M_to_v

	subroutine from_M_to_v2(rho_matrix,rho_vector)
		REAL(8), DIMENSION(2*dim_Lapl_half_z,Node_max_r-1):: rho_matrix
		REAL(8), DIMENSION(dim_Lapl), intent(out) :: rho_vector
		INTEGER k,i,j
		k=1
		do j=1,dim_Lapl_r
		do i=1,(2*dim_Lapl_half_z)
				if(j<Node_max_r) then
					rho_vector(k)= rho_matrix(i,j)
				else
					rho_vector(k)=0.0
				endif
				k=k+1
		enddo
		enddo
	end subroutine from_M_to_v2



	subroutine from_v_to_M(Phi_vector,Phi_matrix)
		REAL(8), DIMENSION(2*dim_Lapl_half_z,dim_Lapl_r), intent(out) :: Phi_matrix
		REAL(8), DIMENSION(dim_Lapl), intent(in) 	:: Phi_vector
		INTEGER k,i,j

		k=1
		do j=1,dim_Lapl_r
		do i=1,(2*dim_Lapl_half_z)
				Phi_matrix(i,j)=Phi_vector(k)
				k=k+1
		enddo
		enddo

	end subroutine from_v_to_M

	subroutine init_B_field
		REAL(8) beta_0
		INTEGER i,j

		! computes B field from E field, all in laboratory frame
		! B = || beta_0 || unit vector_z x Efield
		! all the bunches must have the same average initial gammma

		beta_0 = sqrt( 1.-1./(bunch_initialization%bunch_gamma_m(1))**2 ) ! all the bunches must have the same average initial gammma

		do i=2,(mesh_par%Nzm-2)
			do j=2,(mesh_par%Nxm-1)
				mesh(i,j)%Bphi_bunch = beta_0*mesh(i,j)%Ex_bunch
			enddo
		enddo

		! ----- Center properly B field for FDTD
		call center_B_fields

		! ----- Boundary conditions for FDTD
		call BC_Bfield

	end subroutine init_B_field

	subroutine center_E_fields
	INTEGER i,j
	REAL(8), DIMENSION(mesh_par%Nzm,mesh_par%Nxm) :: Ex_not_centered,Ez_not_centered

	Ex_not_centered(:,:)=mesh(:,:)%Ex_bunch
	Ez_not_centered(:,:)=mesh(:,:)%Ez_bunch

		!--- old 2 points stencil ---!
		!do i=2,(mesh_par%Nzm-2)
		!	do j=2,(mesh_par%Nxm-2)
		!		mesh(i,j)%Ex_bunch = 0.25*(Ex_not_centered(i,j) + Ex_not_centered(i,j+1) + Ex_not_centered(i-1,j) + Ex_not_centered(i-1,j+1))
		!		mesh(i,j)%Ez_bunch = 0.25*(Ez_not_centered(i,j) + Ez_not_centered(i,j-1) + Ez_not_centered(i+1,j) + Ez_not_centered(i+1,j-1))
		!	enddo
		!enddo

		!---> using 3points stencil
		do i=1,(mesh_par%Nzm-1)
			do j=1,(mesh_par%Nxm-1)
				mesh(i,j)%Ex_bunch = 0.5*(Ex_not_centered(i,j) + Ex_not_centered(i,j+1) )
				mesh(i,j)%Ez_bunch = 0.5*(Ez_not_centered(i,j) + Ez_not_centered(i,j+1) )
			enddo
		enddo
	end subroutine center_E_fields



	subroutine center_B_fields
		INTEGER i
		REAL(8), DIMENSION(mesh_par%Nzm,mesh_par%Nxm) :: B_not_centered

		B_not_centered(:,:) = mesh(:,:)%Bphi_bunch

		do i=2,(mesh_par%Nzm-2)
			mesh(i,:)%Bphi_bunch = 0.5*(B_not_centered(i,:) + B_not_centered(i+1,:))
		enddo
	end subroutine center_B_fields


	subroutine BC_Efield
		! -- Ez
		! Left Boundary
		mesh(1         ,:         )%Ez = mesh(Node_min_z,:         )%Ez
		mesh(1         ,:         )%Ez_bunch = mesh(Node_min_z,:   )%Ez_bunch
		! Right Boundary
		mesh(Node_end_z,:         )%Ez = mesh(Node_max_z,:         )%Ez
		mesh(Node_end_z,:         )%Ez_bunch = mesh(Node_max_z,:   )%Ez_bunch
		! Lower Boundary
		mesh(:         ,1         )%Ez = mesh(:         ,2         )%Ez
		mesh(:         ,1         )%Ez_bunch = mesh(:   ,2         )%Ez_bunch
		! Upper Boundary
		mesh(:         ,Node_end_r)%Ez = mesh(:         ,Node_max_r)%Ez
		mesh(:         ,Node_end_r)%Ez_bunch = mesh(:   ,Node_max_r)%Ez_bunch

		! -- Er
		! Left Boundary
		mesh(1         ,:         )%Ex = mesh(Node_min_z,:         )%Ex
		mesh(1         ,:         )%Ex_bunch = mesh(Node_min_z,:   )%Ex_bunch
		! Right Boundary
		mesh(Node_end_z,:         )%Ex = mesh(Node_max_z,:         )%Ex
		mesh(Node_end_z,:         )%Ex_bunch = mesh(Node_max_z,:   )%Ex_bunch
		! Lower Boundary
		mesh(:         ,1         )%Ex = 0.
		mesh(:         ,1         )%Ex_bunch = 0.
		! Upper Boundary
		mesh(:         ,Node_end_r)%Ex = mesh(:         ,Node_max_r)%Ex
		mesh(:         ,Node_end_r)%Ex_bunch    = mesh(:   ,Node_end_r-3)%Ex_bunch !***ALBZ***
		mesh(:         ,Node_end_r-1)%Ex_bunch = mesh(:   ,Node_end_r-3)%Ex_bunch
		mesh(:         ,Node_end_r-2)%Ex_bunch = mesh(:   ,Node_end_r-3)%Ex_bunch
	end subroutine BC_Efield



	subroutine BC_Bfield
		! ----- Boundary conditions
		! Left Boundary
		mesh(1         ,:         )%Bphi = mesh(Node_min_z,:         )%Bphi
		mesh(1         ,:         )%Bphi_bunch = mesh(Node_min_z,:         )%Bphi_bunch
		! Right Boundary
		mesh(Node_end_z,:         )%Bphi = mesh(Node_max_z,:         )%Bphi
		mesh(Node_end_z,:         )%Bphi_bunch = mesh(Node_max_z,:         )%Bphi_bunch
		! Lower Boundary
		mesh(:         ,1         )%Bphi = 0.
		mesh(:         ,1         )%Bphi_bunch = 0.
		! Upper Boundary
		mesh(:         ,Node_end_r)%Bphi = mesh(:         ,Node_max_r)%Bphi
		mesh(:         ,Node_end_r)%Bphi_bunch = mesh(:         ,Node_max_r)%Bphi_bunch
	end subroutine BC_Bfield



	subroutine init_EM_fields_coax_shells
	! all the bunches must have the same average initial gammma
		REAL(8) gamma_0
		INTEGER j, j_prime
		write(*,'(A)') 'Initialization bunch self-consistent field :: Initializing EM fields, coaxial shells method '

		!---project charge on the grid

		mesh(:,:)%Ez_bunch=0. ! with this initialization Ez is completely null
		mesh(:,:)%Ex_bunch=0.
		gamma_0 = bunch_initialization%bunch_gamma_m(1)

		do j=2,(mesh_par%Nxm-2)
			do j_prime=1,j

				!Ex on axis is zero, thus cycle starts from j=2

				if (j_prime.eq.1) then 		!contribution by axis cell

					mesh(2:(mesh_par%Nzm-2),j)%Ex_bunch = mesh(2:(mesh_par%Nzm-2),j)%Ex_bunch + &
					(0.25*mesh_par%dxm**2)*0.5*mesh(2:(mesh_par%Nzm-2),j_prime)%rho/gamma_0/(x_mesh_shifted(j))

				else !contribution by cells not on axis

					if(j_prime.eq.j) then 	!contribution by the cell itself
						mesh(2:(mesh_par%Nzm-2),j)%Ex_bunch = mesh(2:(mesh_par%Nzm-2),j)%Ex_bunch + &
						0.5*mesh(2:(mesh_par%Nzm-2),j_prime)%rho/gamma_0*(mesh_par%dxm-1./x_mesh_shifted(j_prime)*0.25*mesh_par%dxm**2)
					else 					! contribution by cells with smaller radius
						mesh(2:(mesh_par%Nzm-2),j)%Ex_bunch = mesh(2:(mesh_par%Nzm-2),j)%Ex_bunch + &
						2.*(x_mesh_shifted(j_prime)*mesh_par%dxm)*0.5*mesh(2:(mesh_par%Nzm-2),j_prime)%rho/gamma_0/(x_mesh_shifted(j)-x_mesh_shifted(j_prime))
					endif

				endif

			enddo
		enddo

		call center_E_field_coax_shells
		mesh(:,:)%Ex_bunch = mesh(:,:)%Ex_bunch*gamma_0*0.5 ! gamma factor needed to go back to the lab frame
		write(*,'(A)') 'E field initialized'

		call init_B_field
		write(*,'(A)') 'B field initialized'
	end subroutine init_EM_fields_coax_shells


	subroutine center_E_field_coax_shells
		! centers Ex field for FDTD from the one found by coaxial shells method
		! all the bunches must have the same average initial gammma
		INTEGER i,j
		REAL(8), DIMENSION(mesh_par%Nzm,mesh_par%Nxm) :: Ex_not_centered

		Ex_not_centered(:,:)=mesh(:,:)%Ex_bunch
		mesh(:,:)%Ex_bunch = 0.

		do i=2,(mesh_par%Nzm-2)
			do j=1,(mesh_par%Nxm-2)
				mesh(i,j)%Ex_bunch = 0.5*(Ex_not_centered(i,j) + Ex_not_centered(i+1,j))
			enddo
		enddo

		! ----- Boundary conditions
		call BC_Efield

	end subroutine center_E_field_coax_shells

!-----------------------------------------------!
SUBROUTINE ReInitialise_EB_Bunch_fields
	if ( sim_parameters%L_BunchREinit ) then
		bunch(1)%part(:)%cmp(14)=1.
		if ( abs(calculate_nth_moment(1,1,3,'nocentral') &
		   -sim_parameters%lastBunchREinit) > sim_parameters%bunch_reinit_distance_um ) then
		! if ( abs(calculate_nth_moment_bunch(1,1,3)-sim_parameters%lastBunchREinit) > sim_parameters%bunch_reinit_distance_um ) then
			sim_parameters%lastBunchREinit = calculate_nth_moment(1,1,3,'nocentral')
			call Compute_bunch_density
			call init_EM_fields
		end if
	end if
END SUBROUTINE ReInitialise_EB_Bunch_fields


subroutine init_null_EM_fields
	! Nonphysical option: null initial EM fields, only for comparisons with initialization
	write(*,*) 'Initializing EM fields to zero'
	mesh(:,:)%Ex = 0.D0
	mesh(:,:)%Ez = 0.D0
	mesh(:,:)%Bphi = 0.D0
	mesh(:,:)%Bphi_old 	= 0.D0

	mesh(:,:)%Ex_bunch 	= 0.D0
	mesh(:,:)%Ez_bunch 	= 0.D0
	mesh(:,:)%Bphi_bunch = 0.D0
	mesh(:,:)%Bphi_old_bunch = 0.D0
	write(*,*) 'EM fields initialized to zero'
end subroutine init_null_EM_fields


!--- *** ---!
subroutine init_external_Bfields
	if(Bpoloidal%L_Bpoloidal) then
		write(*,'(A)') 'B-field :: Externally imposed'
		call set_external_Bfield
	endif
end subroutine init_external_Bfields
!--- *** ---!

subroutine set_external_Bfield
	!*** Poloidal External Field ***!
	!--- Controlled by input file:
	!--- Activate with: Bpoloidal%L_Bpoloidal=.TRUE.
	real(8) :: Bfield,Radius,a,Zposition
	Bfield   = mu0*Bpoloidal%background_current_A(1)/(2.D0*pi*Bpoloidal%capillary_radius_um*1e-6) !Poloidal field from current
	write(*,'(A,1p1e14.5,A)') 'B_external activated :: Peak B-field > ',Bfield,' (T)'
	Bfield   = Bfield / (96.*sqrt(plasma%n0)/3e8) !from Dimensional to DimensionLESS
	Bpoloidal%capillary_radius = Bpoloidal%capillary_radius_um * plasma%k_p !from Dimensional to DimensionLESS
	!Bpoloidal%B_ex_poloidal   = Bpoloidal%B_ex_poloidal / (electron_mass*(plasma%omega_p*1e15)/electron_charge)

	do i=2,(mesh_par%Nzm) !---*** ***---!
		Zposition=(mesh_par%z_min+i*mesh_par%dzm)/plasma%k_p
		if(Zposition < Bpoloidal%z_coordinate_um(1)) then
			do j=2,(mesh_par%Nxm)

				  Select Case (Bpoloidal%Bprofile(1))
						case (1) !---Linear+Cubic
							Radius=r_mesh(j)/Bpoloidal%capillary_radius
							a=Bpoloidal%a_shape(1)
							mesh(i,j)%B_ex_poloidal = -Bfield*((1.D0-a)*Radius+a*Radius**3)

						case (2) !---Linear+FlatSaturation
							a=Bpoloidal%a_shape(1)*plasma%k_p
							mesh(i,j)%B_ex_poloidal = -Bfield/a*r_mesh(j)
							if(r_mesh(j)>a) then
								mesh(i,j)%B_ex_poloidal = -Bfield
							endif

						case (3) !---Exponential Profile
							a=Bpoloidal%a_shape(1)*plasma%k_p
							mesh(i,j)%B_ex_poloidal = -Bfield*( 1.D0-exp(-r_mesh(j)/a) )/( 1.D0-exp(-r_mesh(mesh_par%Nxm)/a) )

						case(4) !--- x^a
							a=Bpoloidal%a_shape(1)
							Radius=r_mesh(j)/Bpoloidal%capillary_radius
							mesh(i,j)%B_ex_poloidal = -Bfield*( Radius )**a

						case(5) !--- a^-1 + (1-a^-1) * 1/x^a
							a=Bpoloidal%a_shape(1)
							Radius=r_mesh(j)/Bpoloidal%capillary_radius
							mesh(i,j)%B_ex_poloidal = -Bfield* (1.D0/a + (1.D0-1.D0/a)*Radius**a)

					end select !---and Shape Select

			enddo
		endif
	enddo
	!--->BC
  mesh(1,:)%B_ex_poloidal						=-mesh(2,:)%B_ex_poloidal
	mesh(mesh_par%Nzm,:)%B_ex_poloidal= mesh(mesh_par%Nzm-1,:)%B_ex_poloidal
	mesh(:,1)%B_ex_poloidal						= mesh(:,2)%B_ex_poloidal
	mesh(:,mesh_par%Nxm)%B_ex_poloidal= mesh(:,mesh_par%Nxm-1)%B_ex_poloidal
end subroutine set_external_Bfield



END MODULE init_fields
