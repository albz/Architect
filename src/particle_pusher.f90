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

use digit_precision
use my_types
use use_my_types
use pstruct_data
use architect_class_structure
use utilities

IMPLICIT NONE



CONTAINS

!Kernel_MoveParticle_FDTD
!Evoluzione del beam in coordinate cilindriche


  SUBROUTINE Kernel_MoveParticle_FDTD



   IMPLICIT NONE


   INTEGER ip,i,j
   REAL(8) :: gamma,gammaend
   REAL(8) :: Er_onparticle,Ez_onparticle,Bphi_onparticle 				!original fields interpolated on particle position
   REAL(8) :: Ex_onparticle,Ey_onparticle,Bx_onparticle,By_onparticle	!Er-Bphi on particle position projected on {x,y} components
   REAL(8) :: q,m
   REAL(8) :: xacc,yacc,zacc,racc
   REAL(8) :: Beta_beam_x,Beta_beam_y,beta_beam_r,Beta_beam_z
   INTEGER :: indz,indx,indz_s,indx_s
   REAL(8) :: Wr,Wz,Wr_s,Wz_s,fraz,fraz_s,pos_r,pos_z
   REAL(8) :: w00,w10,w01,w11
   REAL(8) :: Ex,Ey,Ez,Er,Bx,By,sin_theta,cos_theta,ux,uy,uz,ux_minus,uy_minus,uz_minus,ux_plus,uy_plus,uz_plus
   REAL(8) :: gamma_half_timestep,ux_prime,uy_prime,uz_prime,ux_new,uy_new,uz_new
   REAL(8) :: tx,ty,tz,sx,sy,sz,alpha,t2
   REAL(8) :: inv_gamma,pxsm,pysm,pzsm,us2,s,upx,upy,upz


   !--------------------------------------------------------------------------------!
   !                             Cycle on every particle                            !
   !--------------------------------------------------------------------------------!
   do j =1,sim_parameters%Nbunches
   do ip=1,bunchip%n_particles(j)

		! ----- Field weighting on grid

		if (bunch(j)%part(ip)%cmp(7).eq. 1.) then ! particle inside the moving window --> weight force from the grid
			call weights_particle_ongrid(1,bunch(j)%part(ip)%cmp(1),bunch(j)%part(ip)%cmp(2) ,bunch(j)%part(ip)%cmp(3) , &
											 bunch(j)%part(ip)%cmp(9),bunch(j)%part(ip)%cmp(10),bunch(j)%part(ip)%cmp(11), &
											 Wr,Wz,Wr_s,Wz_s,fraz,fraz_s,pos_r,pos_z,indx,indx_s,indz,indz_s)
		

			!--- Er on particle ---!
			call particle_weights(3,Wr,Wz,Wr_s,Wz_s,fraz,fraz_s,indx,indx_s,pos_r,w00,w10,w01,w11)

			Er_onparticle = mesh(indz    ,indx_s  )%Ex  		* w00 &
						  + mesh(indz    ,indx_s+1)%Ex  		* w01 &
						  + mesh(indz+1  ,indx_s  )%Ex  		* w10 &
						  + mesh(indz+1  ,indx_s+1)%Ex  		* w11
			if(sim_parameters%L_selffield) Er_onparticle = Er_onparticle &
                          + mesh(indz    ,indx_s  )%Ex_bunch	* w00 &
      			          + mesh(indz    ,indx_s+1)%Ex_bunch	* w01 &
      			          + mesh(indz+1  ,indx_s  )%Ex_bunch	* w10 &
      			          + mesh(indz+1  ,indx_s+1)%Ex_bunch	* w11

			!--- Ez on particle ---!
			call particle_weights(4,Wr,Wz,Wr_s,Wz_s,fraz,fraz_s,indx,indx_s,pos_r,w00,w10,w01,w11)

			Ez_onparticle  = mesh(indz_s    ,indx  )%Ez  		* w00 &
						   + mesh(indz_s    ,indx+1)%Ez  		* w01 &
						   + mesh(indz_s+1  ,indx  )%Ez   		* w10 &
						   + mesh(indz_s+1  ,indx+1)%Ez         * w11
			if(sim_parameters%L_selffield)  Ez_onparticle = Ez_onparticle &
			               + mesh(indz_s    ,indx  )%Ez_bunch   * w00 &
         			       + mesh(indz_s    ,indx+1)%Ez_bunch   * w01 &
         			       + mesh(indz_s+1  ,indx  )%Ez_bunch   * w10 &
         			       + mesh(indz_s+1  ,indx+1)%Ez_bunch   * w11

			!--- Bphi on particle ---!
			!--- Bphi must be interpolated between new and old value ---!
			call particle_weights(5,Wr,Wz,Wr_s,Wz_s,fraz,fraz_s,indx,indx_s,pos_r,w00,w10,w01,w11)

			Bphi_onparticle = 0.5*(   mesh(indz_s  ,indx_s  )%Bphi     	* w00 &
						            + mesh(indz_s  ,indx_s+1)%Bphi     	* w01 &
					      	        + mesh(indz_s+1,indx_s  )%Bphi     	* w10 &
							        + mesh(indz_s+1,indx_s+1)%Bphi     	* w11 )
			if(sim_parameters%L_selffield) Bphi_onparticle = Bphi_onparticle &
                            + 0.5*(   mesh(indz_s  ,indx_s  )%Bphi_bunch  * w00 &
      	                            + mesh(indz_s  ,indx_s+1)%Bphi_bunch  * w01 &
      		                        + mesh(indz_s+1,indx_s  )%Bphi_bunch  * w10 &
      	                            + mesh(indz_s+1,indx_s+1)%Bphi_bunch  * w11	)
			
			Bphi_onparticle = Bphi_onparticle + &
					          0.5*(   mesh(indz_s  ,indx_s  )%Bphi_old 	* w00 &
						    +         mesh(indz_s  ,indx_s+1)%Bphi_old 	* w01 &
						    +         mesh(indz_s+1,indx_s  )%Bphi_old 	* w10 &
						    +         mesh(indz_s+1,indx_s+1)%Bphi_old 	* w11 )
			if(sim_parameters%L_selffield) Bphi_onparticle = Bphi_onparticle &
            		        + 0.5*(   mesh(indz_s  ,indx_s  )%Bphi_old_bunch 	* w00 &
  									+ mesh(indz_s  ,indx_s+1)%Bphi_old_bunch	* w01 &
  									+ mesh(indz_s+1,indx_s  )%Bphi_old_bunch	* w10 &
  									+ mesh(indz_s+1,indx_s+1)%Bphi_old_bunch	* w11 )

      !--- add External poloidal effect ---!
      		if(Bpoloidal%L_Bpoloidal) then
        		Bphi_onparticle = Bphi_onparticle + &
  					    (mesh(indz_s  ,indx_s  )%B_ex_poloidal 	* w00 &
  						+mesh(indz_s  ,indx_s+1)%B_ex_poloidal 	* w01 &
  						+mesh(indz_s+1,indx_s  )%B_ex_poloidal 	* w10 &
  						+mesh(indz_s+1,indx_s+1)%B_ex_poloidal 	* w11	)
      		end if


		else ! particle out of window, no dynamical evolution
			Ez_onparticle   = 0.
			Er_onparticle   = 0.
			Bphi_onparticle = 0.
      		pos_z = plasma%k_p*     ( bunch(j)%part(ip)%cmp(3)-sim_parameters%zg )
      		pos_r = plasma%k_p* sqrt( bunch(j)%part(ip)%cmp(1)**2+bunch(j)%part(ip)%cmp(2)**2 )
		endif

	! ----- Leap-Frog with Boris rotation for 'equal' particle ----- !
	if( trim(bunchip%PWeights(j))=='equal' .and. .not. sim_parameters%L_vay) then

		cos_theta      = plasma%k_p* bunch(j)%part(ip)%cmp(1)/pos_r  ! x/r
		sin_theta      = plasma%k_p* bunch(j)%part(ip)%cmp(2)/pos_r	! y/r

		Ex_onparticle  = Er_onparticle * cos_theta
		Ey_onparticle  = Er_onparticle * sin_theta

		Bx_onparticle  = -Bphi_onparticle * sin_theta
		By_onparticle  =  Bphi_onparticle * cos_theta
		!Bz=0

		ux             = bunch(j)%part(ip)%cmp(4) !beta-x * gamma
		uy             = bunch(j)%part(ip)%cmp(5) !beta-y * gamma
		uz             = bunch(j)%part(ip)%cmp(6) !beta-z * gamma
		q              = bunch(j)%part(ip)%cmp(15)
		m              = bunch(j)%part(ip)%cmp(16)

		ux_minus  = ux + q/m * 0.5*Ex_onparticle * sim_parameters%dt
		uy_minus  = uy + q/m * 0.5*Ey_onparticle * sim_parameters%dt
		uz_minus  = uz + q/m * 0.5*Ez_onparticle * sim_parameters%dt

		gamma_half_timestep = sqrt(1. + ux_minus**2 + uy_minus**2 + uz_minus**2 )

		tx        = q/m * 0.5*Bx_onparticle/gamma_half_timestep * sim_parameters%dt
		ty        = q/m * 0.5*By_onparticle/gamma_half_timestep * sim_parameters%dt
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

		ux_new    = ux_plus + q/m * 0.5*Ex_onparticle * sim_parameters%dt
		uy_new    = uy_plus + q/m * 0.5*Ey_onparticle * sim_parameters%dt
		uz_new    = uz_plus + q/m * 0.5*Ez_onparticle * sim_parameters%dt

		gammaend  = sqrt(1. + ux_new**2 + uy_new**2 + uz_new**2)

		Beta_beam_x = ux_new/gammaend !updated velocity
		Beta_beam_y = uy_new/gammaend !updated velocity
		Beta_beam_z = uz_new/gammaend !updated velocity

		xacc  = bunch(j)%part(ip)%cmp(1)+c*Beta_beam_x*sim_parameters%dt_fs
		yacc  = bunch(j)%part(ip)%cmp(2)+c*Beta_beam_y*sim_parameters%dt_fs
		zacc  = bunch(j)%part(ip)%cmp(3)+c*Beta_beam_z*sim_parameters%dt_fs


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
			write(*,*) 'Ez_onparticle = ',Ez_onparticle
			write(*,*) 'fraz = ',fraz
			write(*,*) 'fraz_s = ',fraz_s
			write(*,*) 'indx = ',indx
			write(*,*) 'indx_s = ',indx_s
			write(*,*) 'Er_onparticle = ',Er_onparticle
			write(*,*) 'Bphi_onparticle = ',Bphi_onparticle
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

	endif !casse for 'equal particles' and with no vay pusher

	!--- Leap-Frog with VAY pusher ---!
	if( trim(bunchip%PWeights(j))=='equal' .and. sim_parameters%L_vay) then
		cos_theta = plasma%k_p* bunch(j)%part(ip)%cmp(1)/pos_r  ! x/r
		sin_theta = plasma%k_p* bunch(j)%part(ip)%cmp(2)/pos_r	! y/r

		Ex_onparticle = Er_onparticle * cos_theta
		Ey_onparticle = Er_onparticle * sin_theta

		Bx_onparticle = -Bphi_onparticle * sin_theta
		By_onparticle =  Bphi_onparticle * cos_theta
		!Bz=0

		ux             = bunch(j)%part(ip)%cmp(4) !beta-x * gamma
		uy             = bunch(j)%part(ip)%cmp(5) !beta-y * gamma
		uz             = bunch(j)%part(ip)%cmp(6) !beta-z * gamma
		q              = bunch(j)%part(ip)%cmp(15)
		m              = bunch(j)%part(ip)%cmp(16)
		
		inv_gamma = one_dp/sqrt(one_dp + ux**2 + uy**2 + uz**2 )

		upx  = ux + q/m * sim_parameters%dt * Ex_onparticle
		upy  = uy + q/m * sim_parameters%dt * Ey_onparticle
		upz  = uz + q/m * sim_parameters%dt * Ez_onparticle

        Tx   = q/m * sim_parameters%dt/2.* Bx_onparticle
        Ty   = q/m * sim_parameters%dt/2.* By_onparticle
        Tz   = q/m * sim_parameters%dt/2.* zero_dp

        upx = upx + inv_gamma*(uy*Tz - uz*Ty) 
        upy = upy + inv_gamma*(uz*Tx - ux*Tz)
        upz = upz + inv_gamma*(ux*Ty - uy*Tx)

        alpha = one_dp + upx*upx + upy*upy + upz*upz
        T2    = Tx**2 + Ty**2 + Tz**2

        s     = alpha - T2
        us2   = (upx*Tx + upy*Ty + upz*Tz)**2

		alpha = one_dp/sqrt(0.5*(s + sqrt(s**2 + 4.0*( T2 + us2 ))))

        Tx = alpha *Tx
        Ty = alpha *Ty
        Tz = alpha *Tz

        s = one_dp/(one_dp+Tx**2+Ty**2+Tz**2)
        alpha   = upx*Tx + upy*Ty + upz*Tz
        
        pxsm = s*(upx + alpha*Tx + Tz*upy - Ty*upz)
        pysm = s*(upy + alpha*Ty + Tx*upz - Tz*upx)
		pzsm = s*(upz + alpha*Tz + Ty*upx - Tx*upy)

 		!--- Updates beam momenta ---!
		bunch(j)%part(ip)%cmp(4)  = pxsm
		bunch(j)%part(ip)%cmp(5)  = pysm
		bunch(j)%part(ip)%cmp(6)  = pzsm
		
		!--- Stores the old positions, needed for correct J computation ---!
		bunch(j)%part(ip)%cmp(9 )=bunch(j)%part(ip)%cmp(1)
		bunch(j)%part(ip)%cmp(10)=bunch(j)%part(ip)%cmp(2)
		bunch(j)%part(ip)%cmp(11)=bunch(j)%part(ip)%cmp(3)

		!--- Updates beam positions and velocities ---!
		gammaend  = sqrt(1. + pxsm**2 + pysm**2 + pzsm**2)

		Beta_beam_x = pxsm/gammaend !updated velocity
		Beta_beam_y = pysm/gammaend !updated velocity
		Beta_beam_z = pzsm/gammaend !updated velocity

		bunch(j)%part(ip)%cmp(1)  = bunch(j)%part(ip)%cmp(1) + c*Beta_beam_x*sim_parameters%dt_fs
		bunch(j)%part(ip)%cmp(2)  = bunch(j)%part(ip)%cmp(2) + c*Beta_beam_y*sim_parameters%dt_fs
		bunch(j)%part(ip)%cmp(3)  = bunch(j)%part(ip)%cmp(3) + c*Beta_beam_z*sim_parameters%dt_fs

	endif !casse for 'equal particles' with VAY pusher


	! ----- Leap-Frog with Boris rotation for 'weighted' particle ----- !
	if(trim(bunchip%PWeights(j))=='weighted') then
		gamma       = sqrt(1.D0 + bunch(j)%part(ip)%cmp(4)**2 + bunch(j)%part(ip)%cmp(6)**2 )
		beta_beam_r = bunch(j)%part(ip)%cmp(4)/gamma ! r-component
		beta_beam_z = bunch(j)%part(ip)%cmp(6)/gamma ! z-component

		bunch(j)%part(ip)%cmp(4) = bunch(j)%part(ip)%cmp(4) + (Er_onparticle-beta_beam_z*Bphi_onparticle) * sim_parameters%dt
		bunch(j)%part(ip)%cmp(6) = bunch(j)%part(ip)%cmp(6) + (Ez_onparticle+beta_beam_r*Bphi_onparticle) * sim_parameters%dt

		!--- storing old position ---!
		bunch(j)%part(ip)%cmp(9 )=bunch(j)%part(ip)%cmp(1)
		bunch(j)%part(ip)%cmp(10)=bunch(j)%part(ip)%cmp(2)
		bunch(j)%part(ip)%cmp(11)=bunch(j)%part(ip)%cmp(3)

		gamma       = sqrt(1.D0 + bunch(j)%part(ip)%cmp(4)**2 + bunch(j)%part(ip)%cmp(6)**2 )
		beta_beam_r = bunch(j)%part(ip)%cmp(4)/gamma ! r-component
		beta_beam_z = bunch(j)%part(ip)%cmp(6)/gamma ! z-component

		bunch(j)%part(ip)%cmp(1) = bunch(j)%part(ip)%cmp(1) + c*beta_beam_r*sim_parameters%dt_fs
		bunch(j)%part(ip)%cmp(3) = bunch(j)%part(ip)%cmp(3) + c*beta_beam_z*sim_parameters%dt_fs
		if(bunch(j)%part(ip)%cmp(1) < 0.) then
			bunch(j)%part(ip)%cmp(1) = - bunch(j)%part(ip)%cmp(1)
			bunch(j)%part(ip)%cmp(4) = - bunch(j)%part(ip)%cmp(4)
		endif
	endif !weighted case

	enddo
	enddo

   return

   END SUBROUTINE

END MODULE
