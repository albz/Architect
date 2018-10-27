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


MODULE external_background_density

USE my_types
USE use_my_types

IMPLICIT NONE

CONTAINS

  subroutine load_background_external_density_profile()
    integer :: Nz,Nr,i,j
    real(8) :: Dz,Dr,resolution
    real(8), dimension(:,:), allocatable :: flipped
    character(199) :: string

    write(*,'(A,A)') 'loading background density from external file: ',TRIM(bck_plasma%filename)
    open(9,file=TRIM(bck_plasma%filename),status='old')

   !--- read: header ---!
   read(9,'(A)') string !first header comment
   read(9,'(A)') string !Nz-header
   read(9,*) Nz
   read(9,'(A)') string !Nr-header
   read(9,*) Nr
   read(9,'(A)') string !Dz-header in um
   read(9,*) Dz
   read(9,'(A)') string !Dr-header un um
   read(9,*) Dr
   read(9,'(A)') string !ne-header un cm^-3


    allocate(mesh_util%ne_ext(Nz,Nr))
    Do j=1,Nr
      read(9,*) mesh_util%ne_ext(:,j)
    ENDDO
    close(9)

    !--- from dimensional to dimensionless        ---!
    !--- External file provided in units of cm^-3 ---!
    Do i=1,Nz
      Do j=1,Nr
        mesh_util%ne_ext(i,j)=mesh_util%ne_ext(i,j)/plasma%n0
        if(mesh_util%ne_ext(i,j) .le. bck_plasma%threshold_suppression) mesh_util%ne_ext(i,j)=0.0d0
      Enddo
    Enddo

    !--- FLIP-Left-Right        ---!
    !--- The matrix extends on the right, should be on the left, I am flipping it ---!
    allocate(flipped(Nz,Nr))
    Do i=1,Nz
        flipped(i,:)=mesh_util%ne_ext(Nz+1-i,:)
    Enddo
    mesh_util%ne_ext=flipped !copy back
    deallocate(flipped)


    !--- cross check mesh resolution ---!
    resolution=1000.d0
    if(Nint(Dz*resolution)==Nint(mesh_par%dz/plasma%k_p*resolution)) then
      write(*,*) 'the same resolution has been used for Dz mesh importing'
    else
      write(*,*) 'ERROR :: NOT the same resolution has been used for Dz mesh importing'
      stop
    endif
    if(Nint(Dr*resolution)==Nint(mesh_par%dr/plasma%k_p*resolution)) then
      write(*,*) 'the same resolution has been used for Dr mesh importing'
    else
      write(*,*) 'ERROR :: NOT the same resolution has been used for Dr mesh importing'
      stop
    endif

  end subroutine load_background_external_density_profile


  real(8) FUNCTION background_density_value_external(i,j)
    integer, intent(in) :: i,j
    integer :: i_eff,j_eff,dim1,dim2
    real(8) :: dzm_um

    !--- *** ---!
    dzm_um=mesh_par%dz/plasma%k_p
    background_density_value_external=0.d0
    !--- *** ---!
    dim1=size(mesh_util%ne_ext,1)
    dim2=size(mesh_util%ne_ext,2)

    i_eff=Node_end_z-i+1
    ! whatch out :: you need index swappint
    ! in Architect i=1 is the first cell on the Left (z_max_moving_um)
    ! in the matrix: i=1 is the Right index, so i=1 has to become Nzm
    j_eff=j

    !--- adding moving window idx-shifting ---!
    i_eff = i_eff - Nint(mesh_par%z_max_moving_um/dzm_um) + 1 + Nint(bck_plasma%shift_um/dzm_um)

    !--- check matrix dimensions and return values ---!
    if(i_eff>dim1) return
    if(j_eff>dim2) return
    if(i_eff>0) background_density_value_external=mesh_util%ne_ext(i_eff,j_eff)
  END FUNCTION background_density_value_external





END MODULE
