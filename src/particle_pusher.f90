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

MODULE MoveParticle_FDTD

 USE my_types
 USE use_my_types
 USE pstruct_data
 USE architect_class_structure
 USE utilities


IMPLICIT NONE



CONTAINS

!Kernel_MoveParticle_FDTD
!Evoluzione del beam in coordinate cilindriche


  SUBROUTINE Kernel_MoveParticle_FDTD



   IMPLICIT NONE


   INTEGER ip,i,j
   REAL(8) :: ai,Exi,Ezi,Bphii,gammaend
   REAL(8) :: xacc,yacc,zacc,racc
   REAL(8) :: Beta_beam_x,Beta_beam_y,Beta_beam_z
   INTEGER :: indz,indx,indz_s,indx_s
   REAL(8) :: Wr,Wz,Wr_s,Wz_s,fraz,fraz_s,pos_r,pos_z
   REAL(8) :: w00,w10,w01,w11
   REAL(8) :: Ex,Ey,Ez,Er,Bx,By,sin_theta,cos_theta,ux,uy,uz,ux_minus,uy_minus,uz_minus,ux_plus,uy_plus,uz_plus
   REAL(8) :: gamma_half_timestep,ux_prime,uy_prime,uz_prime,ux_new,uy_new,uz_new
   REAL(8) :: tx,ty,tz,sx,sy,sz


   !--------------------------------------------------------------------------------!
   !                             Cycle on every particle                            !
   !--------------------------------------------------------------------------------!
   do j =1,sim_parameters%Nbunches
   do ip=1,bunch_initialization%n_particles(j)

		! ----- Field weighting on grid

		if (bunch(j)%part(ip)%cmp(7).eq. 1.) then ! particle in moving window --> weight force from the grid
			call weights_particle_ongrid(1,bunch(j)%part(ip)%cmp(1),bunch(j)%part(ip)%cmp(2) ,bunch(j)%part(ip)%cmp(3) , &
											 bunch(j)%part(ip)%cmp(9),bunch(j)%part(ip)%cmp(10),bunch(j)%part(ip)%cmp(11), &
											 Wr,Wz,Wr_s,Wz_s,fraz,fraz_s,pos_r,pos_z,indx,indx_s,indz,indz_s)
		! ----- Exi -------

			call particle_weights(3,Wr,Wz,Wr_s,Wz_s,fraz,fraz_s,indx,indx_s,pos_r,w00,w10,w01,w11)

			Exi   =mesh(indz    ,indx_s  )%Ex  		* w00 &
						+mesh(indz    ,indx_s+1)%Ex  		* w01 &
						+mesh(indz+1  ,indx_s  )%Ex  		* w10 &
						+mesh(indz+1  ,indx_s+1)%Ex  		* w11
            !<--> +mesh(indz    ,indx_s  )%Ex_bunch	* w00 &
      			!<--> +mesh(indz    ,indx_s+1)%Ex_bunch	* w01 &
      			!<--> +mesh(indz+1  ,indx_s  )%Ex_bunch	* w10 &
      			!<--> +mesh(indz+1  ,indx_s+1)%Ex_bunch	* w11


		! ----- Ezi -------
			call particle_weights(4,Wr,Wz,Wr_s,Wz_s,fraz,fraz_s,indx,indx_s,pos_r,w00,w10,w01,w11)

			Ezi   =   mesh(indz_s    ,indx  )%Ez  		* w00 &
						   +mesh(indz_s    ,indx+1)%Ez  		* w01 &
						   +mesh(indz_s+1  ,indx  )%Ez   		* w10 &
						   +mesh(indz_s+1  ,indx+1)%Ez      * w11
               !<--> +mesh(indz_s    ,indx  )%Ez_bunch * w00 &
         			 !<--> +mesh(indz_s    ,indx+1)%Ez_bunch * w01 &
         			 !<--> +mesh(indz_s+1  ,indx  )%Ez_bunch * w10 &
         			 !<--> +mesh(indz_s+1  ,indx+1)%Ez_bunch * w11

		! ----- Bphii ----- (please remember, B must be interpolated between new and old value )
			call particle_weights(5,Wr,Wz,Wr_s,Wz_s,fraz,fraz_s,indx,indx_s,pos_r,w00,w10,w01,w11)

			Bphii =	0.5*(mesh(indz_s  ,indx_s  )%Bphi     	* w00 &
						      +mesh(indz_s  ,indx_s+1)%Bphi     	* w01 &
					      	+mesh(indz_s+1,indx_s  )%Bphi     	* w10 &
				      		+mesh(indz_s+1,indx_s+1)%Bphi     	* w11	)
             !<--> +0.5*(mesh(indz_s  ,indx_s  )%Bphi_bunch  * w00 &
      	      		!<--> +mesh(indz_s  ,indx_s+1)%Bphi_bunch  * w01 &
      		      	!<--> +mesh(indz_s+1,indx_s  )%Bphi_bunch  * w10 &
      	      		!<--> +mesh(indz_s+1,indx_s+1)%Bphi_bunch  * w11	)

			Bphii = Bphii + &
					0.5*(mesh(indz_s  ,indx_s  )%Bphi_old 	* w00 &
						+mesh(indz_s  ,indx_s+1)%Bphi_old 	* w01 &
						+mesh(indz_s+1,indx_s  )%Bphi_old 	* w10 &
						+mesh(indz_s+1,indx_s+1)%Bphi_old 	* w11	)
            !<--> +0.5*(mesh(indz_s  ,indx_s  )%Bphi_old_bunch 	* w00 &
  					!<-->+mesh(indz_s  ,indx_s+1)%Bphi_old_bunch	* w01 &
  					!<-->+mesh(indz_s+1,indx_s  )%Bphi_old_bunch	* w10 &
  					!<-->+mesh(indz_s+1,indx_s+1)%Bphi_old_bunch	* w11	)

      !External poloidal effect
      if(Bpoloidal%L_Bpoloidal) then
        Bphii = Bphii + &
  					  (mesh(indz_s  ,indx_s  )%B_ex_poloidal 	* w00 &
  						+mesh(indz_s  ,indx_s+1)%B_ex_poloidal 	* w01 &
  						+mesh(indz_s+1,indx_s  )%B_ex_poloidal 	* w10 &
  						+mesh(indz_s+1,indx_s+1)%B_ex_poloidal 	* w11	)
      end if


		else ! particle out of window, no dynamical evolution
			Ezi   = 0.
			Exi   = 0.
			Bphii = 0.
		endif

		! ----- Leap-Frog with Boris rotation
		cos_theta = plasma%k_p* bunch(j)%part(ip)%cmp(1)/pos_r  ! x/r
		sin_theta = plasma%k_p* bunch(j)%part(ip)%cmp(2)/pos_r	! y/r

		Ex        = Exi * cos_theta
		Ey        = Exi * sin_theta

		Bx        = -Bphii*sin_theta
		By        =  Bphii*cos_theta
		!Bz=0

		ux        = bunch(j)%part(ip)%cmp(4)
		uy        = bunch(j)%part(ip)%cmp(5)
		uz        = bunch(j)%part(ip)%cmp(6)

		ux_minus  = ux + 0.5*Ex * plasma%omega_p * sim_parameters%dt
		uy_minus  = uy + 0.5*Ey * plasma%omega_p * sim_parameters%dt
		uz_minus  = uz + 0.5*Ezi* plasma%omega_p * sim_parameters%dt

		gamma_half_timestep = sqrt(1. + ux_minus**2 + uy_minus**2 + uz_minus**2 )

		tx        = 0.5*Bx/gamma_half_timestep * plasma%omega_p * sim_parameters%dt
		ty        = 0.5*By/gamma_half_timestep * plasma%omega_p * sim_parameters%dt
		!tz=0

		ux_prime  = ux_minus - uz_minus*ty !+uy_minus*tz
		uy_prime  = uy_minus + uz_minus*tx !-ux_minus*tz
		uz_prime  = uz_minus + ux_minus*ty - uy_minus*tx

		sx        = 2.* tx / (1. + tx**2 + ty**2 )
		sy        = 2.* ty / (1. + tx**2 + ty**2 )
		!sz=0

		ux_plus   = ux_minus - uz_prime*sy !+uy_prime*sz
		uy_plus   = uy_minus + uz_prime*sx !-ux_prime*sz
		uz_plus   = uz_minus + ux_prime*sy - uy_prime*sx

		ux_new    = ux_plus + 0.5*Ex * plasma%omega_p * sim_parameters%dt
		uy_new    = uy_plus + 0.5*Ey * plasma%omega_p * sim_parameters%dt
		uz_new    = uz_plus + 0.5*Ezi* plasma%omega_p * sim_parameters%dt

		gammaend  = sqrt(1. + ux_new**2 + uy_new**2 + uz_new**2)

		Beta_beam_x = ux_new/gammaend !updated velocity
		Beta_beam_y = uy_new/gammaend !updated velocity
		Beta_beam_z = uz_new/gammaend !updated velocity

		xacc  = bunch(j)%part(ip)%cmp(1)+c*Beta_beam_x*sim_parameters%dt
		yacc  = bunch(j)%part(ip)%cmp(2)+c*Beta_beam_y*sim_parameters%dt
		zacc  = bunch(j)%part(ip)%cmp(3)+c*Beta_beam_z*sim_parameters%dt


		! ------------------ Diagnostics for Nan and Infinity after advancement -------- !
		if ((xacc.ne.xacc).or.(yacc.ne.yacc).or.(zacc.ne.zacc)) then
			write(*,*) '-------------------------'
			write(*,*) 'Nan in accelerated position'
			write(*,*) 'bunch number    = ',j
			write(*,*) 'particle number = ',ip
			write(*,*) 'Iteration = ',sim_parameters%iter
			write(*,*) 'x old =', bunch(j)%part(ip)%cmp(1)
			write(*,*) 'y old =', bunch(j)%part(ip)%cmp(2)
			write(*,*) 'z old =', bunch(j)%part(ip)%cmp(3)
			write(*,*) 'x acc =', xacc
			write(*,*) 'y acc =', yacc
			write(*,*) 'z acc =', zacc
			write(*,*) ' '
			write(*,*) 'P_x start = ',bunch(j)%part(ip)%cmp(4)
			write(*,*) 'P_y start = ',bunch(j)%part(ip)%cmp(5)
			write(*,*) 'P_z start = ',bunch(j)%part(ip)%cmp(6)
			write(*,*) '   '
			write(*,*) 'rlocal = ',pos_r
			write(*,*) 'Ezi = ',Ezi
			write(*,*) 'fraz = ',fraz
			write(*,*) 'fraz_s = ',fraz_s
			write(*,*) 'indx = ',indx
			write(*,*) 'indx_s = ',indx_s
			write(*,*) 'Exi = ',Exi
			write(*,*) 'Bphii = ',Bphii
			write(*,*) '   '
			write(*,*) 'Beta_x = ',Beta_beam_x
			write(*,*) 'Beta_y = ',Beta_beam_y
			write(*,*) 'Beta_z = ',Beta_beam_x
			write(*,*) 'gamma end = ',gammaend
			write(*,*) '--------------------------'
			stop
		endif

		if ((xacc.eq.(xacc+ 1.0)).or.(yacc.eq.(yacc+ 1.0)).or.(zacc.eq.(zacc+ 1.0)) ) then
			write(*,*) '--------------------------------------'
			write(*,*) 'bunch number    = ',j
			write(*,*) 'particle number = ',ip
			write(*,*) 'x old =', bunch(j)%part(ip)%cmp(1)
			write(*,*) 'y old =', bunch(j)%part(ip)%cmp(2)
			write(*,*) 'z old =', bunch(j)%part(ip)%cmp(3)
			write(*,*) 'x acc =', xacc
			write(*,*) 'y acc =', yacc
			write(*,*) 'z acc =', zacc
			write(*,*) 'Iteration = ',sim_parameters%iter
			write(*,*) 'position after pusher is Infinity'
			write(*,*) '--------------------------------------'
			stop
		endif

 ! ---------------------    Updates beam momenta    ----------------------------- !

		bunch(j)%part(ip)%cmp(4)     = ux_new
		bunch(j)%part(ip)%cmp(5)     = uy_new
		bunch(j)%part(ip)%cmp(6)     = uz_new

 ! ---------  Stores the old positions, needed for correct J computation -------- !

		bunch(j)%part(ip)%cmp(9 )=bunch(j)%part(ip)%cmp(1)
		bunch(j)%part(ip)%cmp(10)=bunch(j)%part(ip)%cmp(2)
		bunch(j)%part(ip)%cmp(11)=bunch(j)%part(ip)%cmp(3)

 ! ----------------    Updates beam positions and velocities   ------------------ !
		bunch(j)%part(ip)%cmp(1) = xacc
		bunch(j)%part(ip)%cmp(2) = yacc
		bunch(j)%part(ip)%cmp(3) = zacc

	enddo
	enddo

   return

   END SUBROUTINE

END MODULE
