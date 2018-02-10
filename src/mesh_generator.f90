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

USE my_types
USE use_my_types
USE external_background_density

IMPLICIT NONE

CONTAINS
   SUBROUTINE Kernel_Make_a_mesh

   !USE my_types

   IMPLICIT NONE

   INTEGER :: ss,i,j
   REAL :: dzm_min ,dummy1(8,8),dummy2(8,8)

   write(*,'(A)') ' --- Mesh Generation --- '
   sim_parameters%r0     = sim_parameters%rB0(1) ! Transverse size of mesh is based on transverse size of first bunch

   if(.not.allocated(z_mesh).and..not.allocated(x_mesh)) then
         	call define_mesh_par

	     	ALLOCATE(	z_mesh(1:mesh_par%Nzm),x_mesh(1:mesh_par%Nxm), &
						z_mesh_shifted(1:mesh_par%Nzm),x_mesh_shifted(1:mesh_par%Nxm), &
						inv_R(1:mesh_par%Nxm),inv_R_shifted(1:mesh_par%Nxm), &
						mesh(1:mesh_par%Nzm,1:mesh_par%Nxm),&
						r_mesh(1:mesh_par%Nxm),Sr(1:mesh_par%Nxm),Sz(1:mesh_par%Nxm),Vol(1:mesh_par%Nxm) )



			call init_mesh(mesh_par%Nzm,mesh_par%Nxm)

			! Longitudinal mesh of moving window
         	do ss=1,mesh_par%Nzm
            		z_mesh(ss) =(mesh_par%z_min-mesh_par%dzm)+mesh_par%dzm*real(ss-1)   ! mesh longitudinal lattice
         	enddo

			! Transverse mesh of moving window (first mesh point at -Dx/2, second at Dx/2 and so on)
         	do ss=1,mesh_par%Nxm
            		x_mesh(ss) = mesh_par%dxm*(real(ss)-1.5)	                    ! mesh transverse lattice (positive semiaxis)
         	enddo

			! Shifted mesh for weighting purposes
			z_mesh_shifted   = z_mesh+0.5*mesh_par%dzm
			x_mesh_shifted   = x_mesh+0.5*mesh_par%dxm

			! Minimum values of mesh points
			Zmin             = minval(z_mesh)		  		! Moving window coordinates
			Xmin             = minval(x_mesh)
			Zmin_shifted     = minval(z_mesh_shifted)
			Xmin_shifted     = minval(x_mesh_shifted)

			! Inverse of distance from axis, for current deposition


			inv_R            = 1./abs(x_mesh)  				! first element defined at Dr/2 from axis
			inv_R_shifted(1) = 2./mesh_par%dxm
			do i=2,mesh_par%Nxm
				inv_R_shifted(i) = 1./abs(x_mesh_shifted(i))
			enddo


			! --- auxiliary variables
			DeltaR = mesh_par%dxm
			DeltaZ = mesh_par%dzm
			!Dt     = sim_parameters%dt*plasma%omega_p

			one_over_dx = 1./mesh_par%dxm
			one_over_dz = 1./mesh_par%dzm


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
		   r_mesh(1)     = -DeltaR/2.
		   r_mesh(2)     =  DeltaR/2.
		   Sr    (1)     = 0.
		   Sr    (2)     = 0.
		   Sz    (1)     = pi             * DeltaR**2
		   Sz    (2)     = pi             * DeltaR**2
		   Vol   (1)     = pi    * DeltaZ * DeltaR**2
		   Vol   (2)     = pi    * DeltaZ * DeltaR**2

		   do j=(Node_min_lo_r+1),Node_end_lo_r
					r_mesh(j)   = DeltaR*(j-2)+DeltaR/2.
					Sr    (j)   = 2.*pi * DeltaZ * DeltaR*(j-2) !half cell down
					Sz    (j)   = 2.*pi *          DeltaR * r_mesh(j)
					Vol   (j)   = 2.*pi * DeltaZ * DeltaR * r_mesh(j)
		   enddo


  endif

          !-------------------------------------------------------------------------------

	  mesh_par%Nxm_old=mesh_par%Nxm
	  mesh_par%Nzm_old=mesh_par%Nzm
	  mesh_par%dxm_old=mesh_par%dxm
	  mesh_par%dzm_old=mesh_par%dzm

	  write(*,'(A)') 'mesh generated'
	  write(*,*)
  return
  END SUBROUTINE Kernel_Make_a_mesh


  SUBROUTINE define_mesh_par
  USE my_types
  IMPLICIT NONE
  REAL :: dzm_min
  INTEGER :: k

  !--- Set Mesh parameters --- !
         sim_parameters%z0_first_driver = 0.
         !--- *** Transverse *** ---!
         if(mesh_par%R_mesh<0.D0) then
             mesh_par%R_mesh=mesh_par%Rmax*plasma%k_p*sim_parameters%r0
             write(*,'(A)')
             write(*,'(A)') 'From INPUT :: you have selected Rmax :: maximum extension in unit of sigma_r'
             write(*,'(A,f8.3,A,f8.3,A,f8.3,A)') 'Rmax:        ',mesh_par%Rmax,' x sigma_r (first bunch) --- R_mesh (dimensionless):        ',mesh_par%R_mesh,' x Kp --- R_mesh:       ',mesh_par%R_mesh/plasma%k_p,' (um)'
         else
             mesh_par%R_mesh=mesh_par%R_mesh*plasma%k_p
             mesh_par%Rmax  =mesh_par%R_mesh / (plasma%k_p*sim_parameters%r0)
             write(*,'(A)')
             write(*,'(A)') 'From INPUT :: you have selected R_mesh :: maximum extension in microns'
             write(*,'(A,f8.3,A,f8.3,A,f8.3,A)') 'R_mesh:        ',mesh_par%R_mesh/plasma%k_p,' (um) --- R_mesh (dimensionless):        ',mesh_par%R_mesh, 'x Kp --- Rmax (derived): ', mesh_par%Rmax,' x sigma_r (first bunch)'
         end if
         !--- *** plasma radial extension *** ---!
         if(mesh_par%R_mesh_plasma<0.D0) then
             mesh_par%R_mesh_plasma=mesh_par%Rmax_plasma*plasma%k_p*sim_parameters%r0
             write(*,'(A)')
             write(*,'(A)') 'From INPUT :: you have selected R_mesh_plasma :: maximum plasma extension in unit of sigma_r'
             write(*,'(A,f8.3,A,f8.3,A,f8.3,A)') 'Rmax_plasma: ',mesh_par%Rmax_plasma,' x sigma_r (first bunch) --- R_mesh_plasma (dimensionless): ',mesh_par%R_mesh_plasma,' x Kp --- R_mesh_plasma:',mesh_par%R_mesh_plasma/plasma%k_p,' (um)'
         else
             mesh_par%R_mesh_plasma=mesh_par%R_mesh_plasma*plasma%k_p
             mesh_par%Rmax_plasma  =mesh_par%R_mesh_plasma / (plasma%k_p*sim_parameters%r0)
             write(*,'(A)')
             write(*,'(A)') 'From INPUT :: you have selected R_mesh_plasma :: maximum plasma extension in microns'
             write(*,'(A,f8.3,A,f8.3,A,f8.3,A)') 'R_mesh_plasma: ',mesh_par%R_mesh_plasma/plasma%k_p,' (um) --- R_mesh_plasma (dimensionless): ',mesh_par%R_mesh_plasma, 'x Kp --- Rmax (derived): ', mesh_par%Rmax_plasma,' x sigma_r (first bunch)'
         end if



         mesh_par%ScaleX=2.*sim_parameters%r0*plasma%k_p 	! Transverse size of mesh is based on transverse size of first bunch
         !--- *** DELTA_X *** ---! Mesh cell size (transverse) in adimensional units
         if(mesh_par%dxm<0.D0) then      !---> from number of points
             mesh_par%dxm =mesh_par%ScaleX/real(mesh_par%Nsample_r)
             write(*,'(A)')
             write(*,'(A)') 'From INPUT :: transverse direction r : resolution in number of points per sigma_r'
             !write(*,'(A,f8.3,A,I3,A)') 'Delta_R : ',mesh_par%dxm/plasma%k_p,' [um]  ::  sigma_r/Delta_R :',mesh_par%Nsample_r/2,' points (first bunch)'
             write(*,'(A,f8.3,A)') 'Delta_R (um)             > ',mesh_par%dxm/plasma%k_p,'   --- derived ---'
             write(*,'(A,f8.3,A)') 'Delta_R (dimensionless)  > ',mesh_par%dxm           ,'   --- derived ---'
             write(*,'(A,I4,A)')   'Delta_R (points)         > ',mesh_par%Nsample_r/2   ,'       sigma_r/Delta_R first bunch  --- from nml ---'
         else                            !---> from resolution
             mesh_par%dxm = mesh_par%dxm * plasma%k_p
             mesh_par%Nsample_r = 2.*sim_parameters%rB0(1)/mesh_par%dxm*plasma%k_p
             write(*,'(A)')
             write(*,'(A)') 'From INPUT :: transverse direction r : resolution in microns'
             !write(*,'(A,f8.3,A,I3,A)') 'Delta_R : ',mesh_par%dxm/plasma%k_p,' [um]  ::  sigma_r/Delta_R :',mesh_par%Nsample_r/2,' points (first bunch)'
             write(*,'(A,f8.3,A)') 'Delta_R (um)             > ',mesh_par%dxm/plasma%k_p,'   --- from nml ---'
             write(*,'(A,f8.3,A)') 'Delta_R (dimensionless)  > ',mesh_par%dxm           ,'   --- from nml ---'
             write(*,'(A,I4,A)')   'Delta_R (points)         > ',mesh_par%Nsample_r/2   ,'       sigma_r/Delta_R first bunch  --- derived ---'
         endif
         ! Domain radial boundary
         mesh_par%R_mesh=mesh_par%R_mesh 			! Mesh transv. size
         ! mesh nodes (Transverse)
         mesh_par%Nxm         = int(mesh_par%R_mesh/mesh_par%dxm)+1 	! Mesh transv. dimension
         mesh_par%NRmax_plasma= int(mesh_par%R_mesh_plasma/mesh_par%dxm)+1 !Plasma extension (number of cells)
         if(mod(mesh_par%Nxm,2).eq.1) mesh_par%Nxm=mesh_par%Nxm+1


         !--- *** LONGITUDINAL *** ---!
         mesh_par%ScaleZ=sim_parameters%lbunch(1)*plasma%k_p
         !---DELTA_Z---! Mesh cell size (longitudinal) in adimensional units
         if(mesh_par%dzm<0.D0) then      !---> from number of points
             mesh_par%dzm =mesh_par%ScaleZ/real(mesh_par%Nsample_z)	! Mesh longitudinal step
             write(*,'(A)')
             write(*,'(A)') 'From INPUT :: longitudinal direction z : resolution in number of points per sigma_z'
             write(*,'(A,f8.3,A)') 'Delta_Z (um)             > ',mesh_par%dzm/plasma%k_p,'   --- derived ---'
             write(*,'(A,f8.3,A)') 'Delta_Z (dimensionless)  > ',mesh_par%dzm           ,'   --- derived ---'
             write(*,'(A,I4,A)')   'Delta_Z (points)         > ',mesh_par%Nsample_z     ,'       sigma_z/Delta_Z first bunch  --- from nml ---'
         else                            !---> from resolution
             mesh_par%dzm = mesh_par%dzm * plasma%k_p
             mesh_par%Nsample_z = sim_parameters%lbunch(1)/mesh_par%dzm*plasma%k_p
             write(*,'(A)')
             write(*,'(A)') 'From INPUT :: longitudinal direction z : resolution in microns'
             write(*,'(A,f8.3,A)') 'Delta_Z (um)             > ',mesh_par%dzm/plasma%k_p,'   --- from nml ---'
             write(*,'(A,f8.3,A)') 'Delta_Z (dimensionless)  > ',mesh_par%dzm           ,'   --- from nml ---'
             write(*,'(A,I4,A)')   'Delta_Z (points)         > ',mesh_par%Nsample_z     ,'       sigma_z/Delta_Z first bunch  --- derived ---'
         endif


         !---*** Domain boundaries ***---!
        if(mesh_par%Left_mesh<0.D0 .or. mesh_par%Right_mesh<0.D0) then
          mesh_par%z_min=plasma%k_p*(sim_parameters%z0_first_driver-mesh_par%Left_Domain_boundary*plasma%lambda_p)
          mesh_par%z_max=plasma%k_p*(sim_parameters%z0_first_driver+mesh_par%Right_Domain_boundary*plasma%lambda_p)
          mesh_par%z_min_um=mesh_par%z_min/plasma%k_p
          mesh_par%z_max_um=mesh_par%z_max/plasma%k_p

          mesh_par%z_min_moving=mesh_par%z_min
          mesh_par%z_max_moving=mesh_par%z_max
          mesh_par%z_min_moving_um=mesh_par%z_min_um
          mesh_par%z_max_moving_um=mesh_par%z_max_um

          write(*,'(A)')
          write(*,'(A)') 'From INPUT :: you have selected Longitudinal Mesh Boundaries :: in unit of lambda_p'
          write(*,'(A,f8.3,A,f8.3,A)') 'z_min: ',mesh_par%Left_Domain_boundary,' x lambda_p --- z_min (derived): ',mesh_par%z_min_um,' (um)'
          write(*,'(A,f8.3,A,f8.3,A)') 'z_max: ',mesh_par%Right_Domain_boundary,' x lambda_p --- z_max (derived): ',mesh_par%z_max_um,' (um)'
        else
          mesh_par%z_min=plasma%k_p*(sim_parameters%z0_first_driver-mesh_par%Left_mesh)
          mesh_par%z_max=plasma%k_p*(sim_parameters%z0_first_driver+mesh_par%Right_mesh)
          mesh_par%z_min_um=mesh_par%z_min/plasma%k_p
          mesh_par%z_max_um=mesh_par%z_max/plasma%k_p

          mesh_par%z_min_moving=mesh_par%z_min
          mesh_par%z_max_moving=mesh_par%z_max
          mesh_par%z_min_moving_um=mesh_par%z_min_um
          mesh_par%z_max_moving_um=mesh_par%z_max_um

          write(*,'(A)')
          write(*,'(A)') 'From INPUT :: you have selected Longitudinal Mesh Boundaries :: in microns'
          write(*,'(A,f8.3,A,f8.3,A)') 'z_min: ',mesh_par%z_min_um,' (um) --- z_min (derived): ',mesh_par%Left_mesh/plasma%lambda_p,' x lambda_p'
          write(*,'(A,f8.3,A,f8.3,A)') 'z_max: ',mesh_par%z_max_um,' (um) --- z_max (derived): ',mesh_par%Right_mesh/plasma%lambda_p,' x lambda_p'
        endif
         mesh_par%L_mesh = mesh_par%z_max-mesh_par%z_min 	! Mesh long. size
         ! mesh nodes (longitudinal)
         mesh_par%Nzm=int(mesh_par%L_mesh/mesh_par%dzm)+1   	! Mesh long. dimension
	 !! In FDTD version, uses only x>0 axis, add ghost cells to mesh
         if (sim_parameters%FDTD_version.eq.1) then
            		 mesh_par%Nxm=mesh_par%Nxm/2
                 mesh_par%Nzm=mesh_par%Nzm+2
            		 mesh_par%Nzm=mesh_par%Nzm+2
         endif


        !---*** number of point per bunch dimensions ***---!
        write(*,'(A)')
        write(*,'(A)') 'number of discretising points for each bunch:'
        do k=1,sim_parameters%Nbunches
          write(*,'(A,I1,A,I5)') 'Bunch(',k,') :: sigma_r/Delta_R (points) =',int(bunch_initialization%bunch_s_x(k)/mesh_par%dxm*plasma%k_p)
          write(*,'(A,I1,A,I5)') 'Bunch(',k,') :: sigma_z/Delta_Z (points) =',int(bunch_initialization%bunch_s_z(k)/mesh_par%dzm*plasma%k_p)
          write(*,'(A)')
         enddo


         !--- if an external density file is given and required the code goes through this piece of code ---!
         if(bck_plasma%external_density) call load_background_external_density_profile()

      return
   END SUBROUTINE define_mesh_par


   SUBROUTINE init_mesh(n,m)
   USE my_types
   IMPLICIT NONE
   INTEGER n,m

      mesh%rho		= 0.
      mesh%Jz		= 0.
      mesh%Jr		= 0.
      mesh%ux		= 0.
      mesh%uz		= 0.
      mesh%Ez		= 0.
      mesh%Ex		= 0.
      mesh%Bphi		= 0.

	  n=n
	  m=m

      return

   END SUBROUTINE init_mesh


   SUBROUTINE dt_calculation
     !-from CFL to Dt in (fs)-!
     !sim_parameters%dt = sim_parameters%CFL * min(mesh_par%dxm,mesh_par%dzm)/plasma%k_p *1e-6 / 3e8 / 1e-15
     sim_parameters%dt = 0.5D0 * sim_parameters%CFL * 1.D0/sqrt(mesh_par%dxm**(-2)+mesh_par%dzm**(-2))/plasma%k_p *1e-6 / 3e8 / 1e-15
     Dt = sim_parameters%dt*plasma%omega_p

     write(*,*)
     write(*,'(A)') ' --- CFL condition'
     write(*,'(A,f6.3)') 'CFL     >',sim_parameters%CFL
     write(*,'(A,f6.3)') 'CFL (fs)>',sim_parameters%dt
     write(*,'(A,f6.3)') 'CFL (um)>',sim_parameters%dt*c
     write(*,*)
   END SUBROUTINE dt_calculation

END MODULE
