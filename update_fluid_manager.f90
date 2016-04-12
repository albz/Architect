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

MODULE FluidAdvance_FDTD

USE my_types
USE use_my_types


IMPLICIT NONE



CONTAINS

    ! Advances the fluid quantities with cold fluid plasma model.
    SUBROUTINE Kernel_FluidAdvance_FDTD

      USE Advance_Fluid_FCT
      !USE Advance_Fluid_Compact_Upwind
      USE Advance_Fluid_Compact_Upwind_New
      !USE Advance_Fluid_FCT_Alternative

      IMPLICIT NONE

      if  (sim_parameters%Fluid_Scheme==0) then
        write(*,*) 'fluid routine not valid - old upwind is not used in this version'
        stop
        !call AdvanceFluid_CompactUpwind
      else if (sim_parameters%Fluid_Scheme==1)  then
        !call AdvanceFluid_FCT_Alternative(mesh_par,sim_parameters,plasma,Dt_fs)
        call AdvanceFluid_FCT
      else if (sim_parameters%Fluid_Scheme==2)  then
        call AdvanceFluid_CompactUpwind_New
        !call fluid_UpWind
      endif

    return
    END SUBROUTINE
END MODULE
