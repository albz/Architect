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


module architect_class_structure
 implicit none

 integer,parameter :: P_ncmp=16

 type particle
  real(8) :: cmp(P_ncmp)

	! ---- cmp - particle components ---------- !

	! 1  : X
	! 2  : Y
	! 3  : Z
	! 4  : Px
	! 5  : Py
	! 6  : Pz

	! 7  : cut
	! 8  : dcut, cut from integrated parameters diagnostics

	! 9  : X_old (in previous iteration, necessary for correct time centering in charge deposition)
	! 10 : Y_old (in previous iteration, necessary for correct time centering in charge deposition)
	! 11 : Z_old (in previous iteration, necessary for correct time centering in charge deposition)

	! 12 : part_charge (macroparticle charge)
	! 13 : elperpart (electrons per macroparticle)

	! 14 : select particle for diagnostics
	!      originally this part was called maskbunch a logical variable
	!      used to select the particle, it comes to be more conveninet
	!      to create a specific particle-component and then to convert
	!      it to a logical variable
	!      1 :: use it for diagnostics --- 0 :: exclude it

	! 15 : q, particle charge in unit of elementary charge (eg electron=-1.)
	! 16 : m, particle mass in unit of electron mass (eg electron=1.)

 end type particle

 type static_background_ion
  real(8) :: cmp(7)
	! ---- cmp - background ion components ---------- !
	! 1  : Z_coordinate
	! 2  : R_coordinate
	! 3  : A - Mass number
	! 4  : Z - Atomic number
	! 5  : Zstar - ionisation value
  ! 6  : n_plasma : original ion density
  ! 7  : inout - check particle Active(1) or inactive(0)
 end type static_background_ion


 type species
  type(particle),allocatable :: part(:)
 end type species

end module architect_class_structure


  !--- *** ---!
module pstruct_data
	use architect_class_structure
	implicit none

	type(species) :: bunch(6)

	type(static_background_ion), allocatable :: static_ion(:)
end module pstruct_data
