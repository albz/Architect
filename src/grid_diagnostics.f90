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

MODULE grid_diagnostics

  ! USE pstruct_data
  USE my_types
  USE use_my_types
  USE moments

  IMPLICIT NONE

CONTAINS

!---controller---!
  subroutine lineout

    call density_lineout_2bunch_cm
    call Er_lineout_2bunch_cm
    call Ez_onaxis

  end subroutine lineout

!------------------------------ *** -----------------------------------!
!------------------------------ *** -----------------------------------!
!------------------------------ *** -----------------------------------!
!------------------------------ *** -----------------------------------!
!------------------------------ *** -----------------------------------!
!------------------------------ *** -----------------------------------!
!------------------------------ *** -----------------------------------!
!------------------------------ *** -----------------------------------!



subroutine density_lineout_2bunch_cm
  real(8) :: mu_z
  integer :: idx_mu_z,bunch_number,ss
  character*90 :: filename

  bunch_number=2
  mu_z = calculate_nth_moment_bunch(bunch_number,1,3)
  idx_mu_z=1 + int( (mu_z*plasma%k_p-mesh_par%z_min_moving) * one_over_dz )

  filename=TRIM(sim_parameters%path_integrated_diagnostics)//'density_lineout_2bunch_cm.dat'
  call open_file(OSys%macwin,filename)
    write(11,'(1p1000e14.5)') ( mesh(idx_mu_z,ss)%n_plasma_e , ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid  )
  close(11)
end subroutine density_lineout_2bunch_cm




subroutine Er_lineout_2bunch_cm
  real(8) :: mu_z
  integer :: idx_mu_z,bunch_number,ss
  character*90 :: filename

  bunch_number=2
  mu_z = calculate_nth_moment_bunch(bunch_number,1,3)
  idx_mu_z=1 + int( (mu_z*plasma%k_p-mesh_par%z_min_moving) * one_over_dz )

  filename=TRIM(sim_parameters%path_integrated_diagnostics)//'Er_lineout_2bunch_cm.dat'
  call open_file(OSys%macwin,filename)
    write(11,'(1p1000e14.5)') ( mesh(idx_mu_z,ss)%Ex , ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid  )
  close(11)
end subroutine Er_lineout_2bunch_cm

subroutine Ez_onaxis
  integer :: ss
  character*90 :: filename

  filename=TRIM(sim_parameters%path_integrated_diagnostics)//'Ez_tot_onaxis_lineout.dat'
  call open_file(OSys%macwin,filename)
    !write(11,'(1p1000e14.5)') ( mesh(2,ss)%Ez+mesh(2,ss)%Ez_bunch, ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid  )
    write(11,'(1p1000e14.5)') ( mesh(ss,2)%Ez+mesh(ss,2)%Ez_bunch, ss=2,(mesh_par%Nzm-1),sim_parameters%jump_grid  )
  close(11)

  filename=TRIM(sim_parameters%path_integrated_diagnostics)//'Ez_bck_onaxis_lineout.dat'
  call open_file(OSys%macwin,filename)
    write(11,'(1p1000e14.5)') ( mesh(ss,2)%Ez, ss=2,(mesh_par%Nzm-1),sim_parameters%jump_grid  )
  close(11)
end subroutine Ez_onaxis


END MODULE grid_diagnostics
