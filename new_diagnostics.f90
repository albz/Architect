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

MODULE Diagnostics_on_Bunches


IMPLICIT NONE

CONTAINS

!--- --- --- --- --- --- ---!
SUBROUTINE apply_Sigma_cut(sim_par,plasma_k_0)

USE my_types
USE pstruct_data
USE architect_class_structure

   TYPE(simul_param), INTENT(IN) :: sim_par
   INTEGER :: ss, npc, ni, nf, ii
   REAL(8) avgz,plasma_k_0,radius,drelZ,max_ellipsoidal



	! --- --- --- !
	! bunch(es) diagnostic
	! --- --- --- !
	do ii=1,sim_par%Nbunches
		npc=sum(bunch(ii)%part(:)%cmp(8))
		avgz=0.
		avgz=sum(bunch(ii)%part(:)%cmp(3)*bunch(ii)%part(:)%cmp(8))
		avgz=avgz/npc


		max_ellipsoidal = max( sim_par%rB0(ii) , sim_par%lbunch(ii) )
		max_ellipsoidal = max_ellipsoidal * sim_par%sigma_cut

		do ss=1,size(bunch(ii)%part(:))
			if( bunch(ii)%part(ss)%cmp(8) .eq. 1.) then
				radius = sqrt(bunch(ii)%part(ss)%cmp(1)**2+bunch(ii)%part(ss)%cmp(2)**2+(bunch(ii)%part(ss)%cmp(3)-avgz)**2)
				if( radius > max_ellipsoidal ) bunch(ii)%part(ss)%cmp(8)=0.
			endif
		enddo
	enddo



END SUBROUTINE




!--- --- --- --- --- --- ---!
SUBROUTINE write_integrated_diagnostics(sim_par,plasma_k_0)


USE my_types
USE pstruct_data
USE architect_class_structure


   TYPE(simul_param), INTENT(IN) :: sim_par
   INTEGER :: ss, npc, npl, nplost, ni, nf, ii
   REAL(8) :: avgz,plasma_k_0,radius,drelZ
   REAL(8) :: mu_x,mu_y,mu_z,mu_px,mu_py,mu_pz,sigma_x,sigma_y,sigma_z
   REAL(8) :: sigma_px,sigma_py,sigma_pz,Corr_x_px,Corr_y_py
   REAL(8) :: eps_x,eps_y,gamma_mean,sigma_g,dgsug
   REAL(8), DIMENSION(:), ALLOCATABLE :: sigma__x,sigma__y,sigma__z,sigma__pz,sigma__py,sigma__px,gamma_sp
   CHARACTER(1) :: num2str

	! --- --- --- !
	! bunch(es) diagnostic
	! --- --- --- !


	do ii=1,sim_par%Nbunches

!---
		npc=sum(bunch(ii)%part(:)%cmp(8))
		npl=size(bunch(ii)%part(:))
		nplost=npl-npc
		allocate( sigma__x(1:npl),sigma__y(1:npl),sigma__z(1:npl),sigma__pz(1:npl),sigma__py(1:npl),sigma__px(1:npl),gamma_sp(1:npl))
!---
		mu_x  = sum(bunch(ii)%part(:)%cmp(1)*bunch(ii)%part(:)%cmp(8))/npc
		mu_y  = sum(bunch(ii)%part(:)%cmp(2)*bunch(ii)%part(:)%cmp(8))/npc
		mu_z  = sum(bunch(ii)%part(:)%cmp(3)*bunch(ii)%part(:)%cmp(8))/npc
		mu_px = sum(bunch(ii)%part(:)%cmp(4)*bunch(ii)%part(:)%cmp(8))/npc
		mu_py = sum(bunch(ii)%part(:)%cmp(5)*bunch(ii)%part(:)%cmp(8))/npc
		mu_pz = sum(bunch(ii)%part(:)%cmp(6)*bunch(ii)%part(:)%cmp(8))/npc

		sigma__x  = (bunch(ii)%part(:)%cmp(1)-mu_x)*bunch(ii)%part(:)%cmp(8)
		sigma_x   = sqrt( dot_product(dble(sigma__x),dble(sigma__x))/dble(npc) )
		sigma__y  = (bunch(ii)%part(:)%cmp(2)-mu_y)*bunch(ii)%part(:)%cmp(8)
		sigma_y   = sqrt( dot_product(dble(sigma__y),dble(sigma__y))/dble(npc) )
		sigma__z  = (bunch(ii)%part(:)%cmp(3)-mu_z)*bunch(ii)%part(:)%cmp(8)
		sigma_z   = sqrt( dot_product(dble(sigma__z),dble(sigma__z))/dble(npc) )
		sigma__px = (bunch(ii)%part(:)%cmp(4)-mu_px)*bunch(ii)%part(:)%cmp(8)
		sigma_px  = sqrt( dot_product(dble(sigma__px),dble(sigma__px))/dble(npc) )
		sigma__py = (bunch(ii)%part(:)%cmp(5)-mu_py)*bunch(ii)%part(:)%cmp(8)
		sigma_py  = sqrt( dot_product(dble(sigma__py),dble(sigma__py))/dble(npc) )
		sigma__pz = (bunch(ii)%part(:)%cmp(6)-mu_pz)*bunch(ii)%part(:)%cmp(8)
		sigma_pz  = sqrt( dot_product(dble(sigma__pz),dble(sigma__pz))/dble(npc) )

		Corr_x_px = dot_product(dble(sigma__x),dble(sigma__px))/dble(npc)
		Corr_y_py = dot_product(dble(sigma__y),dble(sigma__py))/dble(npc)

		eps_x = sqrt(sigma_x**2 * sigma_px**2 -Corr_x_px**2)
		eps_y = sqrt(sigma_y**2 * sigma_py**2 -Corr_y_py**2)

		gamma_sp=real(bunch(ii)%part(:)%cmp(8))* &
		sqrt(1.+(bunch(ii)%part(:)%cmp(4))**2+(bunch(ii)%part(:)%cmp(5))**2+(bunch(ii)%part(:)%cmp(6))**2)

		gamma_mean=sum(gamma_sp)/real(npc)
		gamma_sp=(gamma_sp-gamma_mean)*real(bunch(ii)%part(:)%cmp(8))
		sigma_g=sqrt(dot_product(gamma_sp,gamma_sp)/real(npc))
		dgsug=sigma_g/gamma_mean

		write(num2str,'(I1)') ii
		open(15,file='bunch_integrated_diagnostic_'//num2str//'.dat',status='unknown',position='append')
			write(15,'(1p100e14.5)') sim_par%sim_time,-mu_z,mu_x,mu_y,mu_z,mu_px,mu_py,mu_pz,sigma_x,sigma_y,sigma_z,Corr_x_px,Corr_y_py,eps_x,eps_y,sigma_g,gamma_mean,dgsug,real(npc),real(nplost)
		close(15)

		deallocate( sigma__x,sigma__y,sigma__z,sigma__pz,sigma__py,sigma__px,gamma_sp)

	ENDDO

END SUBROUTINE


































END MODULE
