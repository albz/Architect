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

MODULE ionisation_module

 USE my_types
 USE use_my_types
 USE pstruct_data
 USE architect_class_structure
 USE utilities
 USE ion_background

IMPLICIT NONE

CONTAINS


!--- kernel ---!
subroutine ionise
	integer, allocatable :: coeff_lm(:)
	real, allocatable :: V(:),coeff_nlstar(:)
	integer :: N,i,j,l,m, k, ion_id
	integer :: material, count
	real(8) :: E_max,E_field_max ,Ea = 5.14d0, VH = 13.5984, omega_a = 41.3, dt_fs = 1.
	real(8) :: Vn, nstar, w_adk,cstar,W_one_lev
	real(8) :: xnew,ynew,znew
	real(8) :: xold,yold,zold
	real(8) :: Wr,Wz,Wr_s,Wz_s,fraz,fraz_s,pos_z,pos_r
	integer :: indx,indx_s,indz,indz_s
	real(8) :: f_pos_r,f_pos_rs,f_pos_z,f_pos_zs
	real(8) :: w00,w10,w01,w11,Er_on_p,Ez_on_p
	real(8) :: rnd

	material = ionisation%atomic_number

	allocate(coeff_lm(material))
	allocate(coeff_nlstar(material))
	allocate(V(material))

	call fillvector_V(material,V)
	call fillvector_coeff_lm(material,coeff_lm)
	call fillvector_coeff_nlstar(material,coeff_nlstar,V)


	!---calculate E-field---!
	do ion_id=1,ionisation%tot_ions
		if(static_ion(ion_id)%cmp(7)>0.d0 .and. static_ion(ion_id)%cmp(5)<static_ion(ion_id)%cmp(3) ) then

			xnew=static_ion(ion_id)%cmp(2)/sqrt(2.)/plasma%k_p
			ynew=static_ion(ion_id)%cmp(2)/sqrt(2.)/plasma%k_p
			znew=static_ion(ion_id)%cmp(1)/plasma%k_p
			xold=xnew
			yold=ynew
			zold=znew
			call weights_particle_ongrid(1,xnew,ynew,znew,xold,yold,zold,Wr,Wz,Wr_s,Wz_s,fraz,fraz_s,pos_r,pos_z,indx,indx_s,indz,indz_s)

		  ! ----- Er -------
			call particle_weights(3,Wr,Wz,Wr_s,Wz_s,fraz,fraz_s,indx,indx_s,pos_r,w00,w10,w01,w11)
			Er_on_p  = mesh(indz    ,indx_s  )%Ex  		* w00 &
								+mesh(indz    ,indx_s+1)%Ex  		* w01 &
								+mesh(indz+1  ,indx_s  )%Ex  		* w10 &
								+mesh(indz+1  ,indx_s+1)%Ex  		* w11 &
		            +mesh(indz    ,indx_s  )%Ex_bunch	* w00 &
		      			+mesh(indz    ,indx_s+1)%Ex_bunch	* w01 &
		      			+mesh(indz+1  ,indx_s  )%Ex_bunch	* w10 &
		      			+mesh(indz+1  ,indx_s+1)%Ex_bunch	* w11

			! ----- Ez -------
			call particle_weights(4,Wr,Wz,Wr_s,Wz_s,fraz,fraz_s,indx,indx_s,pos_r,w00,w10,w01,w11)
			Ez_on_p   =   mesh(indz_s    ,indx  )%Ez  		* w00 &
								   +mesh(indz_s    ,indx+1)%Ez  		* w01 &
								   +mesh(indz_s+1  ,indx  )%Ez   		* w10 &
								   +mesh(indz_s+1  ,indx+1)%Ez      * w11 &
		               +mesh(indz_s    ,indx  )%Ez_bunch * w00 &
		         			 +mesh(indz_s    ,indx+1)%Ez_bunch * w01 &
		         			 +mesh(indz_s+1  ,indx  )%Ez_bunch * w10 &
		         			 +mesh(indz_s+1  ,indx+1)%Ez_bunch * w11


			E_max = (from_dimlessE_to_dimE(Er_on_p)+from_dimlessE_to_dimE(Ez_on_p))/1d9/Ea
			E_max = max(E_max,1d-10)

			!E_field_max = 5.3 !maximum field in eV
			k=static_ion(ion_id)%cmp(5)+1 !ionising shell

			Vn = sqrt(V(k)/VH)
			nstar=float(k)/Vn
			cstar=coeff_nlstar(k)*coeff_lm(k)*0.5*Vn**2
			w_adk=omega_a*coeff_nlstar(k)*coeff_lm(k)*0.5*(Vn)**2
			w_adk=w_adk*((2./E_max)*Vn**3)**(2.*nstar-1.)
			w_adk=max(w_adk,10.0)
			w_adk=w_adk*exp((-2./3./E_max)*Vn**3)
			W_one_lev=1.-exp(-sim_parameters%dt*w_adk)

			call random_number(rnd)
			if(rnd <= W_one_lev) then
				!print *,'yes'
				static_ion(ion_id)%cmp(5) = static_ion(ion_id)%cmp(5)+1

				i=int((static_ion(ion_id)%cmp(1)-mesh_par%z_min_moving)/mesh_par%dzm)+2
		    j=int(static_ion(ion_id)%cmp(2)/mesh_par%dxm)+2
				mesh(i,j)%n_plasma_e = mesh(i,j)%n_plasma_e + 1d0/real(ionisation%particle_per_cell) * static_ion(ion_id)%cmp(6)
				mesh(i,j)%n_plasma_i = mesh(i,j)%n_plasma_i + 1d0/real(ionisation%particle_per_cell) * static_ion(ion_id)%cmp(6)
			endif


		endif
	enddo

	deallocate(coeff_lm,V,coeff_nlstar)

end subroutine ionise


!--- *** ---!
	FUNCTION C_nlstar2(nstar,lstar)
		real, intent(in) :: lstar,nstar
		real :: C_nlstar2
		C_nlstar2=2.**(2.*nstar) / (nstar*gamma(nstar+lstar+1.)*gamma(nstar-lstar))
	END FUNCTION C_nlstar2
!---***

	SUBROUTINE fillvector_coeff_nlstar(material,coeff_nlstar,V)
	integer, intent(in) :: material
	real, intent(out) :: coeff_nlstar(material)
	real, intent(in) :: V(material)
	integer :: i
	real :: lstar,nstar,VH=13.5984

	lstar=sqrt(VH/V(1))-1.
	DO i=material,1,-1
		nstar=real(i)*sqrt(VH/V(i))
		coeff_nlstar(material-i+1)=C_nlstar2(nstar,lstar)
	ENDDO
	END SUBROUTINE fillvector_coeff_nlstar
!---***

	FUNCTION calculate_coeff_lm(l,m)
		integer, intent(in) :: l,m
		integer :: calculate_coeff_lm
		calculate_coeff_lm= ((2*l+1)*Factorial(l+abs(m)))/(2**abs(m)*Factorial(abs(m))*Factorial(l-abs(m)))

	END FUNCTION calculate_coeff_lm
!---***

	SUBROUTINE fillvector_coeff_lm(material,coeff_lm)
	integer, intent(in) :: material
	integer, intent(out) :: coeff_lm(material)
	integer, allocatable :: coeff_lm_reversed(:)
	integer :: count,n,l,m,spin

	allocate(coeff_lm_reversed(material))

	count=0
	DO n=1,6
	DO l=0,n-1
	DO m=-l,l
	DO spin=0,1
		count=count+1
		coeff_lm_reversed(count)=calculate_coeff_lm(l,m)
		if(count==material) goto 100
	ENDDO
	ENDDO
	ENDDO
	ENDDO
100 continue
	DO n=1,material
		coeff_lm(material+1-n)=coeff_lm_reversed(n)
	ENDDO

	deallocate(coeff_lm_reversed)
	END SUBROUTINE fillvector_coeff_lm
!---***
	SUBROUTINE fillvector_V(material,V)
	integer, intent(in) :: material
	real, intent(out) :: V(material)

	select case(material)
		case(1)
			V(1) = 13.5984

		case(2) !     [He]2s
		    V(1)=24.5874
   			V(2)=54.41776


		case(3) !Li : [He]2s
			V(1)= 5.3917
			V(2)= 75.6400
			V(3)= 122.4543

		case(6) !C : [He] 2s^2 2p^2
			V(1)= 11.2603
			V(2)= 24.3833
			V(3)= 47.8878
			V(4)= 64.4939
			V(5)= 392.087
			V(6)= 489.993

		case(7) !N : [He] 2s^2 2p^3
			V(1)= 14.53
			V(2)= 29.601
			V(3)= 47.449
			V(4)= 77.473
			V(5)= 97.89
			V(6)= 552.07
			V(7)= 667.046

		case(10) !Ne : [He]2s^2 2p^6
			V(1)= 21.5646
			V(2)= 40.9633
			V(3)= 63.45
			V(4)= 97.12
			V(5)= 126.21
			V(6)= 157.93
			V(7)= 207.276
			V(8)= 239.099
			V(9)= 1195.8286
			V(10)= 1362.199

		case(18) !Ar : [Ne] 3s^2 3p^6
			V(1) = 15.75962
			V(2) = 27.62967
			V(3) = 40.74
			V(4) = 59.81
			V(5) = 75.02
			V(6) = 91.009
			V(7) = 124.323
			V(8) = 143.460
			V(9) = 422.45
			V(10)= 478.69
			V(11)= 538.96
			V(12)= 618.26
			V(13)= 686.10
			V(14)= 755.74
			V(15)= 854.77
			V(16)= 918.03
			V(17)= 4120.8857
			V(18)= 4426.2296
		case(29) !Cu =[Ar] 3d^10 4s^1
			V(1)  = 7.726
			V(2)  = 20.29
  			V(3)  = 36.84
  			V(4)  = 57.38
			V(5)  = 79.8
			V(6)  = 103.0
  			V(7)  = 139.0
  			V(8)  = 166.0
  			V(9)  = 199.0
  			V(10) = 232.0
  			V(11) = 265.0
  			V(12) = 369.0
  			V(13) = 384.0
  			V(14) = 401.0
  			V(15) = 435.0
  			V(16) = 484.0
  			V(17) = 520.0
  			V(18) = 557.0
  			V(19) = 670.6
  			V(20) = 1697.0
  			V(21) = 1804.0
  			V(22) = 1916.0
  			V(23) = 2060.0
  			V(24) = 2182.0
  			V(25) = 2308.0
  			V(26) = 2478.0
  			V(27) = 2587.8
  			V(28) = 11062.0
  			V(29) = 11567.0

		case default
			write(*,*) 'No material:',material
			stop
	end select
	END SUBROUTINE fillvector_V


END MODULE
