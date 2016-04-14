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

 integer,parameter :: P_ncmp=13

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

	!~! 14 : indx    --> to be converted in int
	!~! 15 : indx_s  --> to be converted in int
	!~! 16 : indz    --> to be converted in int
	!~! 17 : indz_s  --> to be converted in int
	!~! 18 : Wr
	!~! 19 : Wz
	!~! 20 : Wr_s
	!~! 21 : Wz_s
	!~! 22 : fraz
	!~! 23 : fraz_s

	! ---------------------------------------- !

 end type particle

 type species
  type(particle),allocatable :: part(:)
 end type species


end module architect_class_structure


 !--------------------------


module pstruct_data
 use architect_class_structure
 implicit none

 type(species) :: bunch(6)

end module pstruct_data
