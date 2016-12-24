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

module use_my_types
USE my_types

IMPLICIT NONE

TYPE(simul_param) :: sim_parameters
TYPE(mesh_param) :: mesh_par
TYPE(plasma_param) :: plasma
TYPE(bunch_inside_initialization) :: bunch_initialization
TYPE(OSys_param) :: OSys
TYPE(twiss_param) :: twiss
TYPE(Bpoloidal_param) :: Bpoloidal
  TYPE(mesh_utility) :: mesh_util
TYPE(background_plasma_profile) :: bck_plasma
TYPE(ionisation_parameters) :: ionisation
TYPE(dump_and_restart) :: dump_restart

end module use_my_types
