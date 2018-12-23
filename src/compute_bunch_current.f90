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

MODULE Compute_beam_current_FDTD

USE my_types
USE use_my_types
USE class_species
USE class_particle
USE utilities


IMPLICIT NONE

CONTAINS

! Compute_beam_3current_FDTD
! Computes the density and current distribution of a bunch in the
! longitudinal (bunchip%Z) and transverse (bunchip%X) direction.
	SUBROUTINE Compute_beam_3current_FDTD

   	IMPLICIT NONE

		LOGICAL :: L_dump_particle_ongrid,L_CheckParticle
		INTEGER ip,i,j
		REAL(8) :: nel
		REAL(8) :: beta_x,beta_y,beta_z,beta_r
		INTEGER:: indx,indx_s,indz,indz_s
		REAL(8) :: Wr,Wz,Wr_s,Wz_s,fraz,fraz_s,pos_r,pos_z
		REAL(8) :: w00,w10,w01,w11



  ! Please remember: time and space centering of rho and J are different!

  ! -------------------------------------------------------------------------!
  !                Beam charge and current deposition              !
  ! -------------------------------------------------------------------------!



   	mesh%Jz          = 0.
   	mesh%Jr          = 0.

	L_dump_particle_ongrid=.false.
	sim_parameters%gridDeltaOutput=abs(sim_parameters%sim_time*c-sim_parameters%gridLastOutput)
    if ((mod(sim_parameters%iter,sim_parameters%output_grid_nstep).eq.0 ).or.(sim_parameters%iter.eq.1) &
        .or. sim_parameters%gridDeltaOutput>sim_parameters%output_grid_dist) then
		mesh%ne_b         = 0.
		L_dump_particle_ongrid=.true.
	endif

	L_CheckParticle=.false.
	if(abs(sim_parameters%zg-sim_parameters%lastInWindowCheck)>sim_parameters%InWindoWCheckDelta) then
		L_CheckParticle=.true.
		sim_parameters%lastInWindowCheck=sim_parameters%zg
	endif

	if(L_CheckParticle) then
		do j  = 1,sim_parameters%Nbunches
	    do ip = 1,bunchip%n_particles(j)
				call inwindow(j,ip)
			enddo
		enddo
	endif


  if(L_dump_particle_ongrid .or. L_dump_particle_ongrid .or. sim_parameters%L_Bunch_evolve) then
	do j  = 1,sim_parameters%Nbunches
    do ip = 1,bunchip%n_particles(j)

		nel	  = bunch(j)%part(ip)%cmp(13)


		if(L_CheckParticle) call inwindow(j,ip) ! checks if particle is inside moving window


		if (bunch(j)%part(ip)%cmp(7).eq. 1.) then !particle inside window, deposition of its current and charge


			call weights_particle_ongrid(0,bunch(j)%part(ip)%cmp(1),bunch(j)%part(ip)%cmp(2) ,bunch(j)%part(ip)%cmp(3) , &
										 bunch(j)%part(ip)%cmp(9),bunch(j)%part(ip)%cmp(10),bunch(j)%part(ip)%cmp(11), &
										 Wr,Wz,Wr_s,Wz_s,fraz,fraz_s,pos_r,pos_z,indx,indx_s,indz,indz_s)
			! ----- Jr ----
			call particle_weights(0,Wr,Wz,Wr_s,Wz_s,fraz,fraz_s,indx,indx_s,pos_r,w00,w10,w01,w11)

			beta_x = bunch(j)%part(ip)%cmp(4) &
			/sqrt(1. + (bunch(j)%part(ip)%cmp(4))**2 + (bunch(j)%part(ip)%cmp(5))**2 + (bunch(j)%part(ip)%cmp(6))**2 )
			beta_y = bunch(j)%part(ip)%cmp(5) &
			/sqrt(1. + (bunch(j)%part(ip)%cmp(4))**2 + (bunch(j)%part(ip)%cmp(5))**2 + (bunch(j)%part(ip)%cmp(6))**2 )

			beta_r = 0.5*( beta_x*(bunch(j)%part(ip)%cmp(1)+bunch(j)%part(ip)%cmp(9 )) &
			              +beta_y*(bunch(j)%part(ip)%cmp(2)+bunch(j)%part(ip)%cmp(10))   )
			beta_r = beta_r*plasma%k_p/pos_r

			mesh(indz  ,indx_s  )%Jr     = mesh(indz  ,indx_s  )%Jr+nel*w00*beta_r
      mesh(indz+1,indx_s+1)%Jr     = mesh(indz+1,indx_s+1)%Jr+nel*w11*beta_r
      mesh(indz+1,indx_s  )%Jr     = mesh(indz+1,indx_s  )%Jr+nel*w10*beta_r
      mesh(indz  ,indx_s+1)%Jr     = mesh(indz  ,indx_s+1)%Jr+nel*w01*beta_r



			! ----- Jz ----
			call particle_weights(1,Wr,Wz,Wr_s,Wz_s,fraz,fraz_s,indx,indx_s,pos_r,w00,w10,w01,w11)

			beta_z = bunch(j)%part(ip)%cmp(6) &
			/sqrt(1. + (bunch(j)%part(ip)%cmp(4))**2 +(bunch(j)%part(ip)%cmp(5))**2 + (bunch(j)%part(ip)%cmp(6))**2 )

	mesh(indz_s  ,indx  )%Jz     = mesh(indz_s  ,indx  )%Jz+nel*w00*beta_z
	mesh(indz_s+1,indx+1)%Jz     = mesh(indz_s+1,indx+1)%Jz+nel*w11*beta_z
	mesh(indz_s+1,indx  )%Jz     = mesh(indz_s+1,indx  )%Jz+nel*w10*beta_z
	mesh(indz_s  ,indx+1)%Jz     = mesh(indz_s  ,indx+1)%Jz+nel*w01*beta_z


			! ---- rho ---
			! being unnecessary to compute the fields, it is computed only in the iterations in which it is written on file
			if(L_dump_particle_ongrid) then
 				call particle_weights(2,Wr,Wz,Wr_s,Wz_s,fraz,fraz_s,indx,indx_s,pos_r,w00,w10,w01,w11)
				mesh(indz_s  ,indx_s  )%ne_b    = mesh(indz_s  ,indx_s  )%ne_b+nel*w00
				mesh(indz_s+1,indx_s+1)%ne_b    = mesh(indz_s+1,indx_s+1)%ne_b+nel*w11
				mesh(indz_s+1,indx_s  )%ne_b    = mesh(indz_s+1,indx_s  )%ne_b+nel*w10
				mesh(indz_s  ,indx_s+1)%ne_b    = mesh(indz_s  ,indx_s+1)%ne_b+nel*w01
			endif

		endif

	enddo
	enddo
endif

   return

   END SUBROUTINE













! CURRENT ON THE GRID BUNCH by BUNCH
! we use it, for now, to initialize the bunch
	SUBROUTINE Compute_beam_3current_bunch_FDTD(nbunch)

   	IMPLICIT NONE

	INTEGER, INTENT(IN) :: nbunch
   	INTEGER ip,i,j
	REAL(8) nel
  REAL(8) :: beta_x,beta_y,beta_z,beta_r
	INTEGER:: indx,indx_s,indz,indz_s
	REAL(8) :: Wr,Wz,Wr_s,Wz_s,fraz,fraz_s,pos_r,pos_z
	REAL(8) :: w00,w10,w01,w11
	REAL(8) :: conv


  ! Please remember: time and space centering of rho and J are different!

  ! -------------------------------------------------------------------------!
  !                   Beam charge and current deposition                     !
  ! -------------------------------------------------------------------------!



   	mesh%Jz          = 0.
   	mesh%Jr          = 0.
	mesh%ne_b         = 0.

	j  = nbunch !nth bunch
    do ip = 1,bunchip%n_particles(j)


		nel	  = bunch(j)%part(ip)%cmp(13)

		call inwindow(j,ip) ! checks if particle is inside moving window

		if (bunch(j)%part(ip)%cmp(7).eq. 1.) then !particle inside window, deposition of its current and charge


			call weights_particle_ongrid(0,bunch(j)%part(ip)%cmp(1),bunch(j)%part(ip)%cmp(2) ,bunch(j)%part(ip)%cmp(3) , &
										bunch(j)%part(ip)%cmp(9),bunch(j)%part(ip)%cmp(10),bunch(j)%part(ip)%cmp(11), &
										Wr,Wz,Wr_s,Wz_s,fraz,fraz_s,pos_r,pos_z,indx,indx_s,indz,indz_s)

			! ---- rho ---
			! being unnecessary to compute the fields, it is computed only in the iterations in which it is written on file
			call particle_weights(2,Wr,Wz,Wr_s,Wz_s,fraz,fraz_s,indx,indx_s,pos_r,w00,w10,w01,w11)

			mesh(indz_s  ,indx_s  )%ne_b    = mesh(indz_s  ,indx_s  )%ne_b+nel*w00
			mesh(indz_s+1,indx_s+1)%ne_b    = mesh(indz_s+1,indx_s+1)%ne_b+nel*w11
			mesh(indz_s+1,indx_s  )%ne_b    = mesh(indz_s+1,indx_s  )%ne_b+nel*w10
			mesh(indz_s  ,indx_s+1)%ne_b    = mesh(indz_s  ,indx_s+1)%ne_b+nel*w01

		endif !inside window

	enddo

   conv     = (1./(2.*pi*mesh_par%dz*mesh_par%dr))*(plasma%k_p*1.e4)**3   	!Divide by cell volume (1/r factor included in subroutine Compute_beam_3current_FDTD)
   conv     = conv/plasma%n0  						!dimentionless
   mesh(:,:)%ne_b = conv*mesh(:,:)%ne_b


   return

   END SUBROUTINE




	SUBROUTINE Compute_bunch_density

   	IMPLICIT NONE

   	INTEGER ip,i,j
	REAL(8) nel
  REAL(8) :: beta_x,beta_y,beta_z,beta_r
	INTEGER:: indx,indx_s,indz,indz_s
	REAL(8) :: Wr,Wz,Wr_s,Wz_s,fraz,fraz_s,pos_r,pos_z
	REAL(8) :: w00,w10,w01,w11
	REAL(8) :: conv



  ! --------------------------------------------------------------!
  !                   Beam charge on grid                         !
  ! --------------------------------------------------------------!

	mesh%ne_b         = 0.

	do j  = 1,sim_parameters%Nbunches
    do ip = 1,bunchip%n_particles(j)


		nel	  = bunch(j)%part(ip)%cmp(13)

		call inwindow(j,ip) ! checks if particle is inside moving window

		if (bunch(j)%part(ip)%cmp(7).eq. 1.) then !particle inside window, deposition of its current and charge


			call weights_particle_ongrid(0,bunch(j)%part(ip)%cmp(1),bunch(j)%part(ip)%cmp(2) ,bunch(j)%part(ip)%cmp(3) , &
										 bunch(j)%part(ip)%cmp(9),bunch(j)%part(ip)%cmp(10),bunch(j)%part(ip)%cmp(11), &
										 Wr,Wz,Wr_s,Wz_s,fraz,fraz_s,pos_r,pos_z,indx,indx_s,indz,indz_s)

			! ---- rho ---
			call particle_weights(2,Wr,Wz,Wr_s,Wz_s,fraz,fraz_s,indx,indx_s,pos_r,w00,w10,w01,w11)
			mesh(indz_s  ,indx_s  )%ne_b    = mesh(indz_s  ,indx_s  )%ne_b+nel*w00
			mesh(indz_s+1,indx_s+1)%ne_b    = mesh(indz_s+1,indx_s+1)%ne_b+nel*w11
			mesh(indz_s+1,indx_s  )%ne_b    = mesh(indz_s+1,indx_s  )%ne_b+nel*w10
			mesh(indz_s  ,indx_s+1)%ne_b    = mesh(indz_s  ,indx_s+1)%ne_b+nel*w01

		endif

	enddo
	enddo

   conv     = (1./(2.*pi*mesh_par%dz*mesh_par%dr))*(plasma%k_p*1.e4)**3   	!Divide by cell volume (1/r factor included in subroutine Compute_beam_3current_FDTD)
   conv     = conv/plasma%n0  						!dimensionless
   mesh(:,:)%ne_b = conv*mesh(:,:)%ne_b



   return
   END SUBROUTINE







	SUBROUTINE analytic_rho_only_first_bunch

	IMPLICIT NONE

	REAL(8)    :: gamma_0,beta_0,alpha,sigma_r,sigma_z
	INTEGER :: i,j

	gamma_0 = bunchip%gamma(1)
	beta_0  = sqrt(1-1/gamma_0**2.)

	! bunch normalized peak density
	alpha   = bunchip%ChargeB(1)*1e-9/(plasma%n0*1e6)/1.602e-19
	alpha   = alpha/(2.*pi)**1.5/bunchip%sx_um(1)/bunchip%sy_um(1)/bunchip%sz_um(1)/1e-18

	! adimensional bunch dimensions
	sigma_r = bunchip%sx_um(1)*plasma%k_p
	sigma_z = bunchip%sz_um(1)*plasma%k_p

	! force rho to analytic value of multivariate gaussian
	mesh(:,:)%ne_b = 0.
	mesh(:,:)%Jz  = 0.
	mesh(:,:)%Jr  = 0.

	do i= Node_min_z,Node_max_z
    	do j= Node_min_r,Node_max_r
      	mesh(i,j)%ne_b    = alpha*exp(-x_mesh(j)**2./2./sigma_r**2)*exp(-z_mesh(i)**2./2./sigma_z**2)
      	mesh(i,j)%Jz     = mesh(i,j)%ne_b * beta_0
    	enddo
	enddo

	! BC on lower boundary
	do i=Node_min_z,Node_max_z
		mesh(i,1)%ne_b = mesh(i,2)%ne_b
		mesh(i,1)%Jz  = mesh(i,2)%Jz
	enddo

	return
	END SUBROUTINE analytic_rho_only_first_bunch





END MODULE
