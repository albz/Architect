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

MODULE ComputeCurrentFDTD

USE my_types
USE use_my_types
USE Compute_beam_current_FDTD
USE Compute_plasma_current
USE pstruct_data
USE architect_class_structure
USE utilities



IMPLICIT NONE

CONTAINS




   SUBROUTINE Kernel_ComputeCurrent_FDTD
   IMPLICIT NONE
   INTEGER inc,tt
   REAL(8) :: k_0,conv

   k_0 = plasma%k_p

   call Compute_beam_3current_FDTD

   ! Conversion factor
   conv     = (1./(2.*pi*mesh_par%dzm*mesh_par%dxm))*(k_0*1.e4)**3   	!Divide by cell volume (1/r factor included in subroutine Compute_beam_3current_FDTD)
   conv     = conv/plasma%n0           						!dimensionless

   mesh%rho =   mesh%rho*conv
   mesh%Jz  =   mesh%Jz*conv
   mesh%Jr  =   mesh%Jr*conv

   if ( sim_parameters%Fluid_Scheme==0 ) then
	     write(*,*) 'error: old upwind centering is not supported in this version'
	     stop
   !call Compute_plasma_electron_current_FDTD
   else if (sim_parameters%Fluid_Scheme==2   .or.   sim_parameters%Fluid_Scheme==1) then
       call Compute_plasma_electron_current
   endif

   return
   END SUBROUTINE


END MODULE
