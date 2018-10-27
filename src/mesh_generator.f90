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


MODULE Make_a_Mesh

use digit_precision
use my_types
use use_my_types
use external_background_density

IMPLICIT NONE

CONTAINS

    SUBROUTINE Kernel_Make_a_mesh
        integer :: i,j,idx
        real :: dzm_min ,dummy1(8,8),dummy2(8,8)

        write(*,'(A)') ' --- Mesh Generation --- '
        sim_parameters%r0     = sim_parameters%rB0(1) ! Transverse size of mesh is based on transverse size of first bunch

        if(bck_plasma%external_density) call load_background_external_density_profile()

        !--- mesh cells > dimensionless
        mesh_par%dz               = mesh_par%dz_um * plasma%k_p
        mesh_par%dr               = mesh_par%dr_um * plasma%k_p

        !--- left box
        mesh_par%z_min_um          = - mesh_par%Nzm * mesh_par%dz_um
        mesh_par%z_min             = - mesh_par%Nzm * mesh_par%dz
        mesh_par%z_min_moving_um   =   mesh_par%z_min_um
        mesh_par%z_min_moving      =   mesh_par%z_min

        !--- right box
        mesh_par%z_max_um           =  zero_dp
        mesh_par%z_max_moving_um    =  zero_dp
        mesh_par%z_max              =  zero_dp
        mesh_par%z_max_moving       =  zero_dp

        !--- upper side *R*
        mesh_par%R_mesh_um          = mesh_par%Nxm * mesh_par%dr_um
        mesh_par%R_mesh             = mesh_par%Nxm * mesh_par%dr

        !--- plasma channel extension ---!
        mesh_par%R_plasma_um       = maxval(bck_plasma%radius_um(:)) !Plasma extension
        mesh_par%R_plasma          = mesh_par%R_plasma_um * plasma%k_p
        mesh_par%Nr_plasma         = int( mesh_par%R_plasma_um /mesh_par%dr_um)+1 !Plasma extension (number of cells)


        !--- nsample_r
        mesh_par%Nsample_r = sim_parameters%rB0(1)/mesh_par%dr


        allocate(   z_mesh(1:mesh_par%Nzm)                  )
        allocate(   x_mesh(1:mesh_par%Nxm)                  )
        allocate(   z_mesh_shifted(1:mesh_par%Nzm)          )
        allocate(   x_mesh_shifted(1:mesh_par%Nxm)          )
        allocate(   inv_R(1:mesh_par%Nxm)                   )
        allocate(   inv_R_shifted(1:mesh_par%Nxm)           )
        allocate(   mesh(1:mesh_par%Nzm,1:mesh_par%Nxm)     )
        allocate(   r_mesh(1:mesh_par%Nxm)                  )
        allocate(   Sr(1:mesh_par%Nxm)                      )
        allocate(   Sz(1:mesh_par%Nxm)                      )
        allocate(   Vol(1:mesh_par%Nxm)                     )

        mesh%ne_b	= zero_dp
        mesh%Jz		= zero_dp
        mesh%Jr		= zero_dp
        mesh%ux		= zero_dp
        mesh%uz		= zero_dp
        mesh%Ez		= zero_dp
        mesh%Ex		= zero_dp
        mesh%Bphi	= zero_dp


        !--- Longitudinal mesh ---!
        do idx=1,mesh_par%Nzm
            z_mesh(idx) = mesh_par%dz*real(idx-1) + mesh_par%z_min
        enddo

        !--- Radial mesh ---!
        do idx=1,mesh_par%Nxm
            x_mesh(idx) = mesh_par%dr*(real(idx)-1.5)	                    ! mesh transverse lattice (positive semiaxis)
        enddo

        ! Shifted mesh for weighting purposes
        z_mesh_shifted   = z_mesh+0.5*mesh_par%dz
        x_mesh_shifted   = x_mesh+0.5*mesh_par%dr

        ! Minimum values of mesh points
        Zmin             = minval(z_mesh)		  		! Moving window coordinates
        Xmin             = minval(x_mesh)
        Zmin_shifted     = minval(z_mesh_shifted)
        Xmin_shifted     = minval(x_mesh_shifted)

        ! Inverse of distance from axis, for current deposition
        inv_R            = 1./abs(x_mesh)  				! first element defined at Dr/2 from axis
        inv_R_shifted(1) = 2./mesh_par%dr
        do i=2,mesh_par%Nxm
            inv_R_shifted(i) = 1./abs(x_mesh_shifted(i))
        enddo

        ! --- auxiliary variables
        ! DeltaR = mesh_par%dr
        ! DeltaZ = mesh_par%dz
        !Dt     = sim_parameters%dt_fs*plasma%omega_p

        one_over_dx = 1./mesh_par%dr
        one_over_dz = 1./mesh_par%dz

        ! ---------------------------------
        ! Boundaries
        Nz             = mesh_par%Nzm
        Nr             = mesh_par%Nxm

        !--- low order for upwind ---!
        Node_min_z    = 2
        Node_max_z    = Nz-1

        Node_min_r    = 2
        Node_max_r    = Nr-1

        Node_end_z    = Nz
        Node_end_r    = Nr

        !--- low order boundaries for fct
        Node_min_lo_z    = 2
        Node_max_lo_z    = Nz-1

        Node_min_lo_r    = 2
        Node_max_lo_r    = Nr-1

        Node_end_lo_z    = Nz
        Node_end_lo_r    = Nr

        !--- high order boundaries for fct ---!
        Node_min_ho_z    = 4
        Node_max_ho_z    = (Nz-1)-1

        Node_min_ho_r    = 3
        Node_max_ho_r    = (Nr-1)-1

        Node_end_ho_z    = Nz
        Node_end_ho_r    = Nr

        ! --------------------------------

        !--- --- ---!
        r_mesh(1)     = -mesh_par%dr/2.
        r_mesh(2)     =  mesh_par%dr/2.
        Sr    (1)     = 0.
        Sr    (2)     = 0.
        Sz    (1)     = pi                  * mesh_par%dr**2
        Sz    (2)     = pi                  * mesh_par%dr**2
        Vol   (1)     = pi    * mesh_par%dz * mesh_par%dr**2
        Vol   (2)     = pi    * mesh_par%dz * mesh_par%dr**2

        do j=(Node_min_lo_r+1),Node_end_lo_r
            r_mesh(j)   = mesh_par%dr*(j-2)+mesh_par%dr/2.
            Sr    (j)   = 2.*pi * mesh_par%dz * mesh_par%dr*(j-2) !half cell down
            Sz    (j)   = 2.*pi *               mesh_par%dr * r_mesh(j)
            Vol   (j)   = 2.*pi * mesh_par%dz * mesh_par%dr * r_mesh(j)
        enddo


        mesh_par%Nxm_old=mesh_par%Nxm
        mesh_par%Nzm_old=mesh_par%Nzm
        mesh_par%dxm_old=mesh_par%dr
        mesh_par%dzm_old=mesh_par%dz

        write(*,'(A)') 'mesh generated'
        write(*,*)
    END SUBROUTINE Kernel_Make_a_mesh


    SUBROUTINE set_dt
        !--- from CFL to Dt in (omega_p) and in (fs) ---!
        sim_parameters%dt_fs = 0.5D0 * sim_parameters%CFL * 1.D0/sqrt(mesh_par%dr**(-2)+mesh_par%dz**(-2))/plasma%k_p *1e-6 / 3e8 / 1e-15
        sim_parameters%dt    = sim_parameters%dt_fs*plasma%omega_p

        write(*,*)
        write(*,'(A)') ' --- CFL condition'
        write(*,'(A,f6.3)') 'CFL          >',sim_parameters%CFL
        write(*,'(A,f6.3)') 'CFL (fs)     >',sim_parameters%dt_fs
        write(*,'(A,f6.3)') 'CFL (um)     >',sim_parameters%dt_fs*c
        write(*,'(A,f6.3)') 'CFL (omega_p)>',sim_parameters%dt
        write(*,*)
    END SUBROUTINE set_dt

END MODULE
