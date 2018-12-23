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
	USE class_species
	USE class_particle

IMPLICIT NONE

CONTAINS

!--- *** ---!
SUBROUTINE print_dump_and_restart_files
			if(dump_restart%L_onoff) then
						if(mod(sim_parameters%iter,dump_restart%nstep).eq.0 &
								 .or. dump_restart%distance_um <	abs(sim_parameters%sim_time*c-dump_restart%LastOutput_um)   ) then
											 write(*,'(A,I10,A,f12.3)') 'Printing files for Dump and Restart :: at Iteration =',sim_parameters%iter,'  -  at run distance =',sim_parameters%zg
											 call dump_whole_status
											 dump_restart%LastOutput_um=sim_parameters%sim_time*c
						endif
			endif
END SUBROUTINE print_dump_and_restart_files
!--- *** ---!

!---------------------------------!
!---------------------------------!
SUBROUTINE dump_whole_status
			CHARACTER(255) :: filename
			INTEGER :: i,j,k
						filename=TRIM(sim_parameters%path_dumprestart)//'sim_parameters.arch'
						open(15,file=filename,status='replace',access='stream')
						write(15) sim_parameters
						close(15)

						filename=TRIM(sim_parameters%path_dumprestart)//'mesh_par.arch'
						open(15,file=filename,status='replace',access='stream')
						write(15) mesh_par
						close(15)

						filename=TRIM(sim_parameters%path_dumprestart)//'mesh.arch'
						open(15,file=filename,status='replace',access='stream')
						write(15) mesh
						close(15)

						filename=TRIM(sim_parameters%path_dumprestart)//'plasma.arch'
						open(15,file=filename,status='replace',access='stream')
						write(15) plasma
						close(15)

						filename=TRIM(sim_parameters%path_dumprestart)//'bunch.arch'
						open(15,file=filename,status='replace',access='stream')
						write(15) bunchip
						close(15)

						filename=TRIM(sim_parameters%path_dumprestart)//'OSys.arch'
						open(15,file=filename,status='replace',access='stream')
						write(15) OSys
						close(15)

						filename=TRIM(sim_parameters%path_dumprestart)//'twiss.arch'
						open(15,file=filename,status='replace',access='stream')
						write(15) twiss
						close(15)

						filename=TRIM(sim_parameters%path_dumprestart)//'Bpoloidal.arch'
						open(15,file=filename,status='replace',access='stream')
						write(15) Bpoloidal
						close(15)

						! filename=TRIM(sim_parameters%path_dumprestart)//'mesh_util.arch'
						! open(15,file=filename,status='replace',access='stream')
						! write(15) mesh_util
						! close(15)

						filename=TRIM(sim_parameters%path_dumprestart)//'bck_plasma.arch'
						open(15,file=filename,status='replace',access='stream')
						write(15) bck_plasma
						close(15)

						filename=TRIM(sim_parameters%path_dumprestart)//'ionisation.arch'
						open(15,file=filename,status='replace',access='stream')
						write(15) ionisation
						close(15)

						filename=TRIM(sim_parameters%path_dumprestart)//'particles.arch'
						open(15,file=filename,status='replace',access='stream')
						write(15) (((bunch(i)%part(j)%cmp(k),k=1,13),J=1,size(bunch(i)%part(:))),I=1,sim_parameters%Nbunches)
						close(15)

						! filename=TRIM(sim_parameters%path_dumprestart)//'bunch.arch'
						! open(15,file=filename,status='replace',access='stream')
						! write(15) bunch
						! close(15)
						!
						! filename=TRIM(sim_parameters%path_dumprestart)//'static_ion.arch'
						! open(15,file=filename,status='replace',access='stream')
						! write(15) static_ion
						! close(15)
END SUBROUTINE dump_whole_status

SUBROUTINE read_whole_dumped_status
			call read_dumped_sim_parameters
			call read_dumped_mesh_par
			call read_dumped_mesh
			call read_dumped_plasma
			call read_dumped_bunch
			call read_dumped_OSys
			call read_dumped_twiss
			call read_dumped_Bpoloidal
			call read_dumped_bck_plasma
			call read_dumped_ionisation
			call read_dumped_particles
END SUBROUTINE read_whole_dumped_status

SUBROUTINE read_dumped_sim_parameters
			CHARACTER(255) :: filename
			filename='dumprestart/sim_parameters.arch'
			open(15,file=filename,status='old',access='stream')
			read(15) sim_parameters
			close(15)
END SUBROUTINE read_dumped_sim_parameters

SUBROUTINE read_dumped_mesh_par
			CHARACTER(255) :: filename
			filename='dumprestart/mesh_par.arch'
			open(15,file=filename,status='old',access='stream')
			read(15) mesh_par
			close(15)
END SUBROUTINE read_dumped_mesh_par

SUBROUTINE read_dumped_mesh
			CHARACTER(255) :: filename
			filename='dumprestart/mesh.arch'
			open(15,file=filename,status='old',access='stream')
			read(15) mesh
			close(15)
END SUBROUTINE read_dumped_mesh

SUBROUTINE read_dumped_plasma
			CHARACTER(255) :: filename
			filename='dumprestart/plasma.arch'
			open(15,file=filename,status='old',access='stream')
			read(15) plasma
			close(15)
END SUBROUTINE read_dumped_plasma

SUBROUTINE read_dumped_bunch
			CHARACTER(255) :: filename
			filename='dumprestart/bunch.arch'
			open(15,file=filename,status='old',access='stream')
			read(15) bunchip
			close(15)
END SUBROUTINE read_dumped_bunch

SUBROUTINE read_dumped_OSys
			CHARACTER(255) :: filename
			filename='dumprestart/OSys.arch'
			open(15,file=filename,status='old',access='stream')
			read(15) OSys
			close(15)
END SUBROUTINE read_dumped_OSys

SUBROUTINE read_dumped_twiss
			CHARACTER(255) :: filename
			filename='dumprestart/twiss.arch'
			open(15,file=filename,status='old',access='stream')
			read(15) twiss
			close(15)
END SUBROUTINE read_dumped_twiss

SUBROUTINE read_dumped_Bpoloidal
			CHARACTER(255) :: filename
			filename='dumprestart/Bpoloidal.arch'
			open(15,file=filename,status='old',access='stream')
			read(15) Bpoloidal
			close(15)
END SUBROUTINE read_dumped_Bpoloidal

SUBROUTINE read_dumped_bck_plasma
			CHARACTER(255) :: filename
			filename='dumprestart/bck_plasma.arch'
			open(15,file=filename,status='old',access='stream')
			read(15) bck_plasma
			close(15)
END SUBROUTINE read_dumped_bck_plasma

SUBROUTINE read_dumped_ionisation
			CHARACTER(255) :: filename
			filename='dumprestart/ionisation.arch'
			open(15,file=filename,status='old',access='stream')
			read(15) ionisation
			close(15)
END SUBROUTINE read_dumped_ionisation

SUBROUTINE read_dumped_particles
			CHARACTER(255) :: filename
			INTEGER :: i,j,k
			filename='dumprestart/particles.arch'
			open(15,file=filename,status='old',access='stream')
			read(15) (((bunch(i)%part(j)%cmp(k),k=1,13),J=1,size(bunch(i)%part(:))),I=1,sim_parameters%Nbunches)
			close(15)
END SUBROUTINE read_dumped_particles

END MODULE
