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

MODULE dump_status

USE my_types
USE use_my_types
USE pstruct_data
USE architect_class_structure

IMPLICIT NONE

CONTAINS


!--- *** ---!
! SUBROUTINE dump_status
! 	call dump_status_particles
! END SUBROUTINE dump_status
!--- *** ---!




!---------------------------------!
!---------------------------------!
SUBROUTINE dump_status_particles
	integer :: i,j,k
	CHARACTER(255) :: filename

	filename=TRIM(sim_parameters%path_dumprestart)//'particle_dump.arch'
	open(15,file=filename,status='unknown',access='stream')

	write(15) 1 !output-identifier-format
	write(15) bunch_initialization%n_total_bunches
	write(15) (bunch_initialization%n_particles(I),I=1,bunch_initialization%n_total_bunches)
	write(15) (((bunch(i)%part(j)%cmp(k),k=1,13),J=1,size(bunch(i)%part(:))),I=1,sim_parameters%Nbunches)
	close(15)

END SUBROUTINE dump_status_particles




END MODULE
