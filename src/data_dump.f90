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

MODULE Data_dumping

USE my_types
USE use_my_types
USE pstruct_data
USE architect_class_structure
USE ion_background
USE Diagnostics_on_Bunches
USE grid_diagnostics
USE moments
USE ISO_C_BINDING


IMPLICIT NONE

CONTAINS

	SUBROUTINE first_print_at_screen
		write(*,'(A,f15.8)') 'Time(fs)    = ',sim_parameters%sim_time
		write(*,'(A,f15.8)') 'dt(fs)      = ',sim_parameters%dt_fs
		write(*,'(A,f8.3)')  'l_p (um)    = ',plasma%lambda_p
		write(*,*)
		write(*,'(A)') ' --- Mesh number of Points:'
		write(*,'(A,I5)') 'Nz_mesh > ',(mesh_par%Nzm-2)
		write(*,'(A,I5)') 'Nx_mesh > ',(mesh_par%Nxm-2)
		write(*,*)
		write(*,'(A)') ' --- Mesh Physical Size:'
		write(*,'(A,f8.3)') 'Delta_Z(um) > ',mesh_par%dz/plasma%k_p
		write(*,'(A,f8.3)') 'Delta_R(um) > ',mesh_par%dr/plasma%k_p
		write(*,*)
		write(*,'(A)') ' --- Box (Domain) physical extension:'
		write(*,'(A,f11.3,A)')       'Z_min_mesh(um) =',maxval(z_mesh)/plasma%k_p, ' ___ R_min_mesh(um) = 0.000'
		write(*,'(A,f11.3,A,f11.3)') 'Z_max_mesh(um) =',minval(z_mesh)/plasma%k_p, ' ___ R_max_mesh(um) = ',maxval(x_mesh)/plasma%k_p
		write(*,*)
		write(*,*)
	END SUBROUTINE first_print_at_screen


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                         !
!                      Diagnostics Subroutines - Binary Format                            !
!                                                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE data_dump
	sim_parameters%gridDeltaOutput=abs(sim_parameters%sim_time*c-sim_parameters%gridLastOutput)
	sim_parameters%PSDeltaOutput=abs(sim_parameters%sim_time*c-sim_parameters%PSLastOutput)

	!--- dump of grid quantitites ---!
	if(mod(sim_parameters%iter,sim_parameters%output_grid_nstep).eq.0 &
		.or. sim_parameters%gridDeltaOutput>sim_parameters%output_grid_dist)then
			if(sim_parameters%iter>30) then
				write(*,'(A,I10,A,f12.3)') 'Printing Grid quantities :: at Iteration =',sim_parameters%iter,'  -  after run distance =',sim_parameters%zg
				call from_particleZstart_to_meshZstar
				if (sim_parameters%Output_format.eq.0) call savedata
				if (sim_parameters%Output_format.eq.1) call savedata_bin
				if (sim_parameters%Output_format.eq.2) then
					call savedata_bin
					call savedata_vtk
				endif
				sim_parameters%gridLastOutput=sim_parameters%sim_time*c
			endif
	endif

	!--- dump of Phase Space ---!
		if(mod(sim_parameters%iter,sim_parameters%output_PS_nstep).eq.0 &
				 .or. sim_parameters%PSDeltaOutput>sim_parameters%output_PS_dist ) then
				 if(sim_parameters%iter>30) then
						write(*,'(A,I10,A,f12.3)') 'Printing PS   quantities :: at Iteration =',sim_parameters%iter,'  -  at run distance =',sim_parameters%zg
						if (sim_parameters%Output_format.eq.0) call save_beam
						if (sim_parameters%Output_format.eq.1) call save_beam_bin
						if (sim_parameters%Output_format.eq.2) then
																									 call save_beam_bin
																									! call save_beam_vtk
					  endif
						sim_parameters%PSLastOutput=sim_parameters%sim_time*c
		endif
	endif
END SUBROUTINE data_dump


	SUBROUTINE final_data_dump

		call from_particleZstart_to_meshZstar

		if (sim_parameters%Output_format.eq.0) then
			call savedata
			call save_beam
		else if (sim_parameters%Output_format.eq.1) then
			call savedata_bin
			call save_beam_bin
		endif

	END SUBROUTINE final_data_dump

	SUBROUTINE write_read_nml
		NAMELIST / write_nml / plasma, sim_parameters, bck_plasma, bunchip, &
													 ionisation, Bpoloidal, twiss, mesh_par, OSys

		open(9,file='datain.arch')!,status='unknown')
			write(9,NML=write_nml,err=30)
30  continue
		close(9)

	END SUBROUTINE write_read_nml




   ! Beam quantities
   SUBROUTINE save_beam_bin

   IMPLICIT NONE

   INTEGER ss,delta_z,n_act,n_trasv,n_long,Output_version,i,n
   CHARACTER :: variable*150,filename*255,positions*6,frmt*4,pas*1,position*10
   REAL(8) avgz,k0plasma,sigma_x,sigma_y,sigma_w,mu_z_witness,zmed,np
   !SAVE n0

!--------------------------- Define positions string for filenames ------------------------------------------------------
!    avgz = sum( bunch(1)%part(:)%cmp(3) * bunch(1)%part(:)%cmp(8) ) / sum(bunch(1)%part(:)%cmp(8))
   avgz = calculate_nth_moment(1,1,3,'nocentral') + mesh_par%z_max_moving_um
   write(position,'(I7.7)') int( -avgz )


!--------------------------- Saving Phase Space ------------------------------------------------------------------------
	filename=TRIM(sim_parameters%path_PS)//TRIM(ADJUSTL(position))//'.arch'
	Output_version = 3

	open(15,file=filename,status='unknown',access='stream')

	!--------------------------- Writes Output version -----------------------------------------------------------------
    	write(15) Output_version

	!--------------------------- Writes number of Bunches and Bunch Charges ----------------------------------------------------------------
			write(15) bunchip%n_total_bunches
			write(15) bunchip%ChargeB(1:bunchip%n_total_bunches)

	!--------------------------- Writes Z position -----------------------------------------------------------------
    	write(15) int(-avgz)

	!--------------------------- Writes Phase Space -----------------------------------------------------------------
	if (sim_parameters%reduced_PS.eq.0) then
		! writes whole phase space
		n=0
		do i=1,sim_parameters%Nbunches
			n=n+size(bunch(i)%part(:))
		enddo
		write(15) n !size(bunchip%X)

		do i =1,sim_parameters%Nbunches
		do ss=1,size(bunch(i)%part(:)),sim_parameters%jump_PS
			write(15) bunch(i)%part(ss)%cmp(1), bunch(i)%part(ss)%cmp(2),bunch(i)%part(ss)%cmp(3), &
			          bunch(i)%part(ss)%cmp(4), bunch(i)%part(ss)%cmp(5),bunch(i)%part(ss)%cmp(6), &
			          1.D0*i, bunch(i)%part(ss)%cmp(7), bunch(i)%part(ss)%cmp(8), &
								bunch(i)%part(ss)%cmp(12), bunch(i)%part(ss)%cmp(13)
		enddo
		enddo
	else if (sim_parameters%reduced_PS.eq.1) then
		! writes the phase space of the particles not cut by the diagnostics

		n=0
		do i=1,sim_parameters%Nbunches
			n=n+sum(bunch(i)%part(:)%cmp(8))
		enddo
		write(15) int(n)

		do i =1,sim_parameters%Nbunches
		do ss=1,size(bunch(i)%part(:)),sim_parameters%jump_PS
	 		if( bunch(i)%part(ss)%cmp(8) .eq. 1. ) then
				write(15) bunch(i)%part(ss)%cmp(1), bunch(i)%part(ss)%cmp(2),bunch(i)%part(ss)%cmp(3), &
						  bunch(i)%part(ss)%cmp(4), bunch(i)%part(ss)%cmp(5),bunch(i)%part(ss)%cmp(6), &
						  1.D0*i, bunch(i)%part(ss)%cmp(7), bunch(i)%part(ss)%cmp(8), &
							bunch(i)%part(ss)%cmp(12), bunch(i)%part(ss)%cmp(13)
			endif
		enddo
		enddo

	endif

	close(15)

   return

   END SUBROUTINE




   ! Grid quantities
	SUBROUTINE savedata_bin
	IMPLICIT NONE
	INTEGER   :: idx,delta_z,n0,n0_w,nf_w,n_w,n_act,n_trasv,n_long,Output_version,i
	CHARACTER :: variable*150,filename*255,positions*6,frmt*4,pas*1,position*10
	REAL(dp) avgz,k0plasma,sigma_x,sigma_y,sigma_w,mu_z_witness,zmed,np
	SAVE n0

	!--- output version ---!
	!Output_version = 3 :prior to ionisatio
	Output_version = 5

	!--- Define positions string for filenames ---!
	avgz = calculate_nth_moment(1,1,3,'nocentral') + mesh_par%z_max_moving_um
	write(position, '(I7.7)' ) int( -avgz )

    filename=TRIM(sim_parameters%path_grid)//TRIM(ADJUSTL(position))//'.arch'
    open(15,file=filename,status='unknown',access='stream')

	!--- Writes Output version ---!
    write(15) Output_version

	!--- Writes header info ---!
	write(15) plasma%k_p
	write(15) plasma%omega_p
	write(15) -avgz*plasma%k_p
	write(15) plasma%n0
	write(15) -avgz

	!--- Writes matrix dimensions (Half plane, no ghost cells) ---!
    write(15) 2*(mesh_par%Nxm-2)/sim_parameters%jump_grid
    write(15) (mesh_par%Nzm-2)/sim_parameters%jump_grid

	!--- Writes grid coordinates ---!
    write(15) (-x_mesh(idx), idx=mesh_par%Nxm-1,2,-sim_parameters%jump_grid)
    write(15) ( x_mesh(idx), idx=2,mesh_par%Nxm-1, sim_parameters%jump_grid)
	write(15) ( z_mesh(idx), idx=2,(mesh_par%Nzm-1),sim_parameters%jump_grid)

	!--- rho_bunch ---!
	write(15) ( mesh(2:(mesh_par%Nzm-1),idx)%ne_b , idx=mesh_par%Nxm-1,2,  -sim_parameters%jump_grid )
	write(15) ( mesh(2:(mesh_par%Nzm-1),idx)%ne_b , idx=2,(mesh_par%Nxm-1),sim_parameters%jump_grid  )

	!--- rho_background ---!
	write(15) ( mesh(2:(mesh_par%Nzm-1),idx)%ne_bck , idx=mesh_par%Nxm-1,2,  -sim_parameters%jump_grid )
	write(15) ( mesh(2:(mesh_par%Nzm-1),idx)%ne_bck , idx=2,(mesh_par%Nxm-1),sim_parameters%jump_grid  )

	!--- Er ---!
	write(15) (-mesh(2:(mesh_par%Nzm-1),idx)%Ex-mesh(2:(mesh_par%Nzm-1),idx)%Ex_bunch , idx=mesh_par%Nxm-2,1,  -sim_parameters%jump_grid )
	write(15) ( mesh(2:(mesh_par%Nzm-1),idx)%Ex+mesh(2:(mesh_par%Nzm-1),idx)%Ex_bunch , idx=2,(mesh_par%Nxm-1),sim_parameters%jump_grid  )
	write(15) (-mesh(2:(mesh_par%Nzm-1),idx)%Ex                                      , idx=mesh_par%Nxm-2,1,  -sim_parameters%jump_grid )
	write(15) ( mesh(2:(mesh_par%Nzm-1),idx)%Ex                                      , idx=2,(mesh_par%Nxm-1),sim_parameters%jump_grid  )
	write(15) (                               -mesh(2:(mesh_par%Nzm-1),idx)%Ex_bunch , idx=mesh_par%Nxm-2,1,  -sim_parameters%jump_grid )
	write(15) (                               +mesh(2:(mesh_par%Nzm-1),idx)%Ex_bunch , idx=2,(mesh_par%Nxm-1),sim_parameters%jump_grid  )

	!--- Ez ---!
	write(15) ( mesh(2:(mesh_par%Nzm-1),idx)%Ez+mesh(2:(mesh_par%Nzm-1),idx)%Ez_bunch , idx=mesh_par%Nxm-1,2,  -sim_parameters%jump_grid )
	write(15) ( mesh(2:(mesh_par%Nzm-1),idx)%Ez+mesh(2:(mesh_par%Nzm-1),idx)%Ez_bunch , idx=2,(mesh_par%Nxm-1),sim_parameters%jump_grid  )
	write(15) ( mesh(2:(mesh_par%Nzm-1),idx)%Ez                                      , idx=mesh_par%Nxm-1,2,  -sim_parameters%jump_grid )
	write(15) ( mesh(2:(mesh_par%Nzm-1),idx)%Ez                                      , idx=2,(mesh_par%Nxm-1),sim_parameters%jump_grid  )
	write(15) (                               +mesh(2:(mesh_par%Nzm-1),idx)%Ez_bunch , idx=mesh_par%Nxm-1,2,  -sim_parameters%jump_grid )
	write(15) (                               +mesh(2:(mesh_par%Nzm-1),idx)%Ez_bunch , idx=2,(mesh_par%Nxm-1),sim_parameters%jump_grid  )

	!--- Bphi ---!
	write(15) (-mesh(2:(mesh_par%Nzm-1),idx)%Bphi-mesh(2:(mesh_par%Nzm-1),idx)%Bphi_bunch-mesh(2:(mesh_par%Nzm-1),idx)%B_ex_poloidal , idx=mesh_par%Nxm-1,2,  -sim_parameters%jump_grid )
	write(15) ( mesh(2:(mesh_par%Nzm-1),idx)%Bphi+mesh(2:(mesh_par%Nzm-1),idx)%Bphi_bunch+mesh(2:(mesh_par%Nzm-1),idx)%B_ex_poloidal , idx=2,(mesh_par%Nxm-1),sim_parameters%jump_grid  )
	write(15) (-mesh(2:(mesh_par%Nzm-1),idx)%Bphi                                                                                  , idx=mesh_par%Nxm-1,2,  -sim_parameters%jump_grid )
	write(15) ( mesh(2:(mesh_par%Nzm-1),idx)%Bphi                                        																					, idx=2,(mesh_par%Nxm-1),sim_parameters%jump_grid  )
	write(15) (                                 -mesh(2:(mesh_par%Nzm-1),idx)%Bphi_bunch 																					, idx=mesh_par%Nxm-1,2,  -sim_parameters%jump_grid )
	write(15) (                                 +mesh(2:(mesh_par%Nzm-1),idx)%Bphi_bunch 																					, idx=2,(mesh_par%Nxm-1),sim_parameters%jump_grid  )
	write(15) (																																				 -mesh(2:(mesh_par%Nzm-1),idx)%B_ex_poloidal , idx=mesh_par%Nxm-1,2,  -sim_parameters%jump_grid )
	write(15) (																																				 +mesh(2:(mesh_par%Nzm-1),idx)%B_ex_poloidal , idx=2,(mesh_par%Nxm-1),sim_parameters%jump_grid  )

	!--- Jr_bunch ---!
	write(15) ( mesh(1:(mesh_par%Nzm-2),idx)%Jr , idx=mesh_par%Nxm-2,1,  -sim_parameters%jump_grid )
	write(15) ( mesh(1:(mesh_par%Nzm-2),idx)%Jr , idx=1,(mesh_par%Nxm-2),sim_parameters%jump_grid  )

	!--- Jr_background ---!
	write(15) ( mesh(2:(mesh_par%Nzm-1),idx)%Jpe_r , idx=mesh_par%Nxm-2,1,  -sim_parameters%jump_grid )
	write(15) ( mesh(2:(mesh_par%Nzm-1),idx)%Jpe_r , idx=2,(mesh_par%Nxm-1),sim_parameters%jump_grid  )

	!--- Jz_bunch ---!
	write(15) ( mesh(2:(mesh_par%Nzm-1),idx)%Jz , idx=mesh_par%Nxm-1,2,  -sim_parameters%jump_grid )
	write(15) ( mesh(2:(mesh_par%Nzm-1),idx)%Jz , idx=2,(mesh_par%Nxm-1),sim_parameters%jump_grid  )

	!--- Jz_background ---!
	write(15) ( mesh(2:(mesh_par%Nzm-1),idx)%Jpe_z , idx=mesh_par%Nxm-1,2,  -sim_parameters%jump_grid )
	write(15) ( mesh(2:(mesh_par%Nzm-1),idx)%Jpe_z , idx=2,(mesh_par%Nxm-1),sim_parameters%jump_grid  )

	!--- Zstar ---!
	write(15) ( mesh(2:(mesh_par%Nzm-1),idx)%Zstar , idx=mesh_par%Nxm-1,2,  -sim_parameters%jump_grid )
	write(15) ( mesh(2:(mesh_par%Nzm-1),idx)%Zstar , idx=2,(mesh_par%Nxm-1),sim_parameters%jump_grid  )

	!--- Writes ni_bck ---!
	write(15) ( mesh(2:(mesh_par%Nzm-1),idx)%ni_bck , idx=mesh_par%Nxm-1,2,  -sim_parameters%jump_grid )
	write(15) ( mesh(2:(mesh_par%Nzm-1),idx)%ni_bck , idx=2,(mesh_par%Nxm-1),sim_parameters%jump_grid  )

    close(15)
	END SUBROUTINE


!--- Savedata_bin subroutine with cpp call to write VTK output ---!
	 SUBROUTINE savedata_vtk
   IMPLICIT NONE
   CHARACTER :: filename*255,position*7

!--- Define utils for vtk writing ---!
   write(position,'(I7.7)') NINT(-sum(bunch(1)%part(:)%cmp(3)*bunch(1)%part(:)%cmp(8))/sum(bunch(1)%part(:)%cmp(8)))

 !--- Writes bunch rho ---!
	filename=TRIM(sim_parameters%path_vtk)//'rho_bunch_um_'//position//'.vtr'//C_NULL_CHAR
	call print_matrix(1,z_mesh(2:mesh_par%Nzm-1)/plasma%k_p,x_mesh(2:mesh_par%Nxm-1)/plasma%k_p,mesh(2:mesh_par%Nzm-2,2:mesh_par%Nxm-2)%ne_b,mesh_par%Nzm-2,mesh_par%Nxm-2, filename)

	filename=TRIM(sim_parameters%path_vtk)//'rho_bck_um_'//TRIM(ADJUSTL(position))//'.vtr'//C_NULL_CHAR
	call print_matrix(2,z_mesh(2:mesh_par%Nzm-1)/plasma%k_p,x_mesh(2:mesh_par%Nxm-1)/plasma%k_p,mesh(2:mesh_par%Nzm-2,2:mesh_par%Nxm-2)%ne_bck,mesh_par%Nzm-2,mesh_par%Nxm-2, filename)

	filename=TRIM(sim_parameters%path_vtk)//'rho_tot_um_'//TRIM(ADJUSTL(position))//'.vtr'//C_NULL_CHAR
	call print_matrix(3,z_mesh(2:mesh_par%Nzm-1)/plasma%k_p,x_mesh(2:mesh_par%Nxm-1)/plasma%k_p,mesh(2:mesh_par%Nzm-2,2:mesh_par%Nxm-2)%ne_b+mesh(2:mesh_par%Nzm-2,2:mesh_par%Nxm-2)%ne_bck,mesh_par%Nzm-2,mesh_par%Nxm-2, filename)

	filename=TRIM(sim_parameters%path_vtk)//'Er_bck_um_'//TRIM(ADJUSTL(position))//'.vtr'//C_NULL_CHAR
	call print_matrix(4,z_mesh(2:mesh_par%Nzm-1)/plasma%k_p,x_mesh(2:mesh_par%Nxm-1)/plasma%k_p,mesh(2:mesh_par%Nzm-2,2:mesh_par%Nxm-2)%Ex,mesh_par%Nzm-2,mesh_par%Nxm-2, filename)

	filename=TRIM(sim_parameters%path_vtk)//'Er_bunch_um_'//TRIM(ADJUSTL(position))//'.vtr'//C_NULL_CHAR
	call print_matrix(5,z_mesh(2:mesh_par%Nzm-1)/plasma%k_p,x_mesh(2:mesh_par%Nxm-1)/plasma%k_p,mesh(2:mesh_par%Nzm-2,2:mesh_par%Nxm-2)%Ex_bunch,mesh_par%Nzm-2,mesh_par%Nxm-2, filename)

	filename=TRIM(sim_parameters%path_vtk)//'Er_tot_um_'//TRIM(ADJUSTL(position))//'.vtr'//C_NULL_CHAR
	call print_matrix(6,z_mesh(2:mesh_par%Nzm-1)/plasma%k_p,x_mesh(2:mesh_par%Nxm-1)/plasma%k_p,mesh(2:mesh_par%Nzm-2,2:mesh_par%Nxm-2)%Ex+mesh(2:mesh_par%Nzm-2,2:mesh_par%Nxm-2)%Ex_bunch,mesh_par%Nzm-2,mesh_par%Nxm-2, filename)

	filename=TRIM(sim_parameters%path_vtk)//'Ez_bck_um_'//TRIM(ADJUSTL(position))//'.vtr'//C_NULL_CHAR
	call print_matrix(7,z_mesh(2:mesh_par%Nzm-1)/plasma%k_p,x_mesh(2:mesh_par%Nxm-1)/plasma%k_p,mesh(2:mesh_par%Nzm-2,2:mesh_par%Nxm-2)%Ez,mesh_par%Nzm-2,mesh_par%Nxm-2, filename)

	filename=TRIM(sim_parameters%path_vtk)//'Ez_bunch_um_'//TRIM(ADJUSTL(position))//'.vtr'//C_NULL_CHAR
	call print_matrix(8,z_mesh(2:mesh_par%Nzm-1)/plasma%k_p,x_mesh(2:mesh_par%Nxm-1)/plasma%k_p,mesh(2:mesh_par%Nzm-2,2:mesh_par%Nxm-2)%Ez_bunch,mesh_par%Nzm-2,mesh_par%Nxm-2, filename)

	filename=TRIM(sim_parameters%path_vtk)//'Ez_tot_um_'//TRIM(ADJUSTL(position))//'.vtr'//C_NULL_CHAR
	call print_matrix(9,z_mesh(2:mesh_par%Nzm-1)/plasma%k_p,x_mesh(2:mesh_par%Nxm-1)/plasma%k_p,mesh(2:mesh_par%Nzm-2,2:mesh_par%Nxm-2)%Ez+mesh(2:mesh_par%Nzm-2,2:mesh_par%Nxm-2)%Ez_bunch,mesh_par%Nzm-2,mesh_par%Nxm-2, filename)

	!Still missing as quantities in the vtk format
	!
	!
	! mesh(i,j)%Bphi
	! mesh(i,j)%Bphi_bunch
	! mesh(i,j)%Bphi
	! mesh(i,j)%B_ex_poloidal

! !--------------------------- Writes Bphi --------------------------------------------------------------------------------
! !--------------------------- Writes bunch Jr ----------------------------------------------------------------------------
! 	write(15) ( mesh(1:(mesh_par%Nzm-2),ss)%Jr , ss=mesh_par%Nxm-2,1,  -sim_parameters%jump_grid )
! 	write(15) ( mesh(1:(mesh_par%Nzm-2),ss)%Jr , ss=1,(mesh_par%Nxm-2),sim_parameters%jump_grid  )
!
! !--------------------------- Writes plasma Jr ---------------------------------------------------------------------------
! 	write(15) ( mesh(2:(mesh_par%Nzm-1),ss)%Jpe_r , ss=mesh_par%Nxm-2,1,  -sim_parameters%jump_grid )
! 	write(15) ( mesh(2:(mesh_par%Nzm-1),ss)%Jpe_r , ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid  )
!
! !--------------------------- Writes bunch Jz ----------------------------------------------------------------------------
! 	write(15) ( mesh(2:(mesh_par%Nzm-1),ss)%Jz , ss=mesh_par%Nxm-1,2,  -sim_parameters%jump_grid )
! 	write(15) ( mesh(2:(mesh_par%Nzm-1),ss)%Jz , ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid  )
!
! !--------------------------- Writes plasma Jz ---------------------------------------------------------------------------
! 	write(15) ( mesh(2:(mesh_par%Nzm-1),ss)%Jpe_z , ss=mesh_par%Nxm-1,2,  -sim_parameters%jump_grid )
! 	write(15) ( mesh(2:(mesh_par%Nzm-1),ss)%Jpe_z , ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid  )
!
! !--------------------------- Writes plasma Zstar ---------------------------------------------------------------------------
! 	write(15) ( mesh(2:(mesh_par%Nzm-1),ss)%Zstar , ss=mesh_par%Nxm-1,2,  -sim_parameters%jump_grid )
! 	write(15) ( mesh(2:(mesh_par%Nzm-1),ss)%Zstar , ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid  )
!
! !--------------------------- Writes ni_bck ---------------------------------------------------------------------------
! 	write(15) ( mesh(2:(mesh_par%Nzm-1),ss)%ni_bck , ss=mesh_par%Nxm-1,2,  -sim_parameters%jump_grid )
! 	write(15) ( mesh(2:(mesh_par%Nzm-1),ss)%ni_bck , ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid  )
!
!     close(15)
!
   return
   END SUBROUTINE savedata_vtk



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                         !
!                      Diagnostics Subroutines - ASCII  Format                            !
!                                                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !Beam quantities
   SUBROUTINE save_beam

   IMPLICIT NONE

   INTEGER ss,zmed,n0,i,n
   CHARACTER :: variable*150,filename*255,position*10
   REAL avgz,k0plasma,np
   SAVE n0

   avgz = calculate_nth_moment(1,1,3,'nocentral') + mesh_par%z_max_moving_um
   write(position,'(I7.7)') int( -avgz )

	filename=TRIM(sim_parameters%path_PS)//'PS_'//TRIM(ADJUSTL(position))//'_um.dat'

	open(15,file=filename,status='unknown')

	if (sim_parameters%reduced_PS.eq.0) then

	!--- whole phase space  ---!
		do i =1,sim_parameters%Nbunches
		do ss=1,size(bunch(i)%part(:)),sim_parameters%jump_PS
			write(15,'(6e14.5,3I4.4)') bunch(i)%part(ss)%cmp(1),bunch(i)%part(ss)%cmp(2),bunch(i)%part(ss)%cmp(3), &
			bunch(i)%part(ss)%cmp(4),bunch(i)%part(ss)%cmp(5),bunch(i)%part(ss)%cmp(6), &
			i,int(bunch(i)%part(ss)%cmp(7)),int(bunch(i)%part(ss)%cmp(8))
		enddo
		enddo

	else if (sim_parameters%reduced_PS.eq.1) then

	!--- reduced phase space 'dcut' apply---!
		do i =1,sim_parameters%Nbunches,sim_parameters%jump_PS
		do ss=1,size(bunch(i)%part(:))
			if( bunch(i)%part(ss)%cmp(8) .eq. 1. ) then
				write(15,'(6e14.5,3I4.4)') bunch(i)%part(ss)%cmp(1),bunch(i)%part(ss)%cmp(2),bunch(i)%part(ss)%cmp(3), &
				bunch(i)%part(ss)%cmp(4),bunch(i)%part(ss)%cmp(5),bunch(i)%part(ss)%cmp(6), &
				i,int(bunch(i)%part(ss)%cmp(7)),int(bunch(i)%part(ss)%cmp(8))
			endif
		enddo
		enddo
	endif

	close(15)

	return

   END SUBROUTINE






   !Grid quantities

   SUBROUTINE savedata

   IMPLICIT NONE

   INTEGER ss,zmed,n0,i
   CHARACTER :: variable*150,filename*255,position*10
   REAL avgz,k0plasma,np
   SAVE n0

   avgz = calculate_nth_moment(1,1,3,'nocentral') + mesh_par%z_max_moving_um
   write(position,'(I7.7)') int( -avgz )


19       format(3000(1x,e14.5))


   variable='rho_'
   filename=TRIM(sim_parameters%path_grid)//TRIM(variable)//TRIM(ADJUSTL(position))//'_um.dat'
   open(15,file=filename,status='unknown')
   do ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid
	write(15,19) mesh(2:(mesh_par%Nzm-1),(mesh_par%Nxm+1-ss))%ne_b
   enddo
   do ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid
	write(15,19) mesh(2:(mesh_par%Nzm-1),ss)%ne_b
   enddo

   close(15)

   variable='Jz_'
   filename=TRIM(sim_parameters%path_grid)//TRIM(variable)//TRIM(ADJUSTL(position))//'_um.dat'
   open(15,file=filename,status='unknown')
   do ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid
	write(15,19) mesh(2:(mesh_par%Nzm-1),(mesh_par%Nxm+1-ss))%Jz
   enddo
   do ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid
	write(15,19) mesh(2:(mesh_par%Nzm-1),ss)%Jz
   enddo
   close(15)

   variable='Jr_'
   filename=TRIM(sim_parameters%path_grid)//TRIM(variable)//TRIM(ADJUSTL(position))//'_um.dat'
   open(15,file=filename,status='unknown')
   do ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid
 	write(15,19) -mesh(2:(mesh_par%Nzm-1),(mesh_par%Nxm+1-ss))%Jr
   enddo
   do ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid
	write(15,19) mesh(2:(mesh_par%Nzm-1),ss)%Jr
   enddo
   close(15)

   variable='Ex_bunch_'
   filename=TRIM(sim_parameters%path_grid)//TRIM(variable)//TRIM(ADJUSTL(position))//'_um.dat'
   open(15,file=filename,status='unknown')
   do ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid
		write(15,19) -mesh(2:(mesh_par%Nzm-1),(mesh_par%Nxm+1-ss))%Ex_bunch
   enddo
   do ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid
	write(15,19) mesh(2:(mesh_par%Nzm-1),ss)%Ex_bunch
   enddo
   close(15)

   variable='Ex_bck_'
   filename=TRIM(sim_parameters%path_grid)//TRIM(variable)//TRIM(ADJUSTL(position))//'_um.dat'
   open(15,file=filename,status='unknown')
   do ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid
		write(15,19) -mesh(2:(mesh_par%Nzm-1),(mesh_par%Nxm+1-ss))%Ex
   enddo
   do ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid
	write(15,19) mesh(2:(mesh_par%Nzm-1),ss)%Ex
   enddo
   close(15)

   variable='Ex_'
   filename=TRIM(sim_parameters%path_grid)//TRIM(variable)//TRIM(ADJUSTL(position))//'_um.dat'
   open(15,file=filename,status='unknown')
   do ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid
		write(15,19) -mesh(2:(mesh_par%Nzm-1),(mesh_par%Nxm+1-ss))%Ex-mesh(2:(mesh_par%Nzm-1),(mesh_par%Nxm+1-ss))%Ex_bunch
   enddo
   do ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid
	write(15,19) mesh(2:(mesh_par%Nzm-1),ss)%Ex+mesh(2:(mesh_par%Nzm-1),ss)%Ex_bunch
   enddo
   close(15)

   variable='Ez_bunch_'
   filename=TRIM(sim_parameters%path_grid)//TRIM(variable)//TRIM(ADJUSTL(position))//'_um.dat'
   open(15,file=filename,status='unknown')
   do ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid
		write(15,19) mesh(2:(mesh_par%Nzm-1),(mesh_par%Nxm+1-ss))%Ez_bunch
   enddo
   do ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid
	write(15,19) mesh(2:(mesh_par%Nzm-1),ss)%Ez_bunch
   enddo
   close(15)

   variable='Ez_bck_'
   filename=TRIM(sim_parameters%path_grid)//TRIM(variable)//TRIM(ADJUSTL(position))//'_um.dat'
   open(15,file=filename,status='unknown')
   do ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid
		write(15,19) mesh(2:(mesh_par%Nzm-1),(mesh_par%Nxm+1-ss))%Ez
   enddo
   do ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid
	write(15,19) mesh(2:(mesh_par%Nzm-1),ss)%Ez
   enddo
   close(15)

   variable='Ez_'
   filename=TRIM(sim_parameters%path_grid)//TRIM(variable)//TRIM(ADJUSTL(position))//'_um.dat'
   open(15,file=filename,status='unknown')
   do ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid
		write(15,19) mesh(2:(mesh_par%Nzm-1),(mesh_par%Nxm+1-ss))%Ez_bunch+mesh(2:(mesh_par%Nzm-1),(mesh_par%Nxm+1-ss))%Ez
   enddo
   do ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid
	write(15,19) mesh(2:(mesh_par%Nzm-1),ss)%Ez_bunch+mesh(2:(mesh_par%Nzm-1),ss)%Ez
   enddo
   close(15)

   variable='Bphi_bunch_'
   filename=TRIM(sim_parameters%path_grid)//TRIM(variable)//TRIM(ADJUSTL(position))//'_um.dat'
   open(15,file=filename,status='unknown')
   do ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid
		write(15,19) -mesh(2:(mesh_par%Nzm-1),(mesh_par%Nxm+1-ss))%Bphi_bunch
   enddo
   do ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid
	write(15,19) mesh(2:(mesh_par%Nzm-1),ss)%Bphi_bunch
   enddo
   close(15)

   variable='Bphi_bck_'
   filename=TRIM(sim_parameters%path_grid)//TRIM(variable)//TRIM(ADJUSTL(position))//'_um.dat'
   open(15,file=filename,status='unknown')
   do ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid
		write(15,19) -mesh(2:(mesh_par%Nzm-1),(mesh_par%Nxm+1-ss))%Bphi
   enddo
   do ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid
	write(15,19) mesh(2:(mesh_par%Nzm-1),ss)%Bphi
   enddo
   close(15)

   variable='Bphi_ex_poloidal'
   filename=TRIM(sim_parameters%path_grid)//TRIM(variable)//TRIM(ADJUSTL(position))//'_um.dat'
   open(15,file=filename,status='unknown')
   do ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid
		write(15,19) -mesh(2:(mesh_par%Nzm-1),(mesh_par%Nxm+1-ss))%B_ex_poloidal
   enddo
   do ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid
	write(15,19) mesh(2:(mesh_par%Nzm-1),ss)%B_ex_poloidal
   enddo
   close(15)

   variable='Bphi_'
   filename=TRIM(sim_parameters%path_grid)//TRIM(variable)//TRIM(ADJUSTL(position))//'_um.dat'
   open(15,file=filename,status='unknown')
   do ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid
		write(15,19) -mesh(2:(mesh_par%Nzm-1),(mesh_par%Nxm+1-ss))%Bphi-mesh(2:(mesh_par%Nzm-1),(mesh_par%Nxm+1-ss))%Bphi_bunch-mesh(2:(mesh_par%Nzm-1),(mesh_par%Nxm+1-ss))%B_ex_poloidal
   enddo
   do ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid
	write(15,19) mesh(2:(mesh_par%Nzm-1),ss)%Bphi+mesh(2:(mesh_par%Nzm-1),ss)%Bphi_bunch+mesh(2:(mesh_par%Nzm-1),ss)%B_ex_poloidal
   enddo
   close(15)

   variable='ux_'
   filename=TRIM(sim_parameters%path_grid)//TRIM(variable)//TRIM(ADJUSTL(position))//'_um.dat'
   open(15,file=filename,status='unknown')
   do ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid
        write(15,19) -mesh(2:(mesh_par%Nzm-1),(mesh_par%Nxm+1-ss))%ux
   enddo
   do ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid
 	write(15,19) mesh(2:(mesh_par%Nzm-1),ss)%ux
   enddo
   close(15)

   variable='uz_'
   filename=TRIM(sim_parameters%path_grid)//TRIM(variable)//TRIM(ADJUSTL(position))//'_um.dat'
   open(15,file=filename,status='unknown')
   do ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid
	write(15,19) mesh(2:(mesh_par%Nzm-1),(mesh_par%Nxm+1-ss))%uz
   enddo
   do ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid
	write(15,19) mesh(2:(mesh_par%Nzm-1),ss)%uz
   enddo
   close(15)

   variable='ne_'
   filename=TRIM(sim_parameters%path_grid)//TRIM(variable)//TRIM(ADJUSTL(position))//'_um.dat'
   open(15,file=filename,status='unknown')
   do ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid
 	write(15,19) mesh(2:(mesh_par%Nzm-1),(mesh_par%Nxm+1-ss))%ne_bck
   enddo
   do ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid
	write(15,19) mesh(2:(mesh_par%Nzm-1),ss)%ne_bck
   enddo
   close(15)

   variable='Jeplasma_r_'
   filename=TRIM(sim_parameters%path_grid)//TRIM(variable)//TRIM(ADJUSTL(position))//'_um.dat'
   open(15,file=filename,status='unknown')
   do ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid
	write(15,19) -mesh(2:(mesh_par%Nzm-1),(mesh_par%Nxm+1-ss))%Jpe_r
   enddo
   do ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid
	write(15,19) mesh(2:(mesh_par%Nzm-1),ss)%Jpe_r
   enddo
   close(15)

   variable='Jeplasma_z_'
   filename=TRIM(sim_parameters%path_grid)//TRIM(variable)//TRIM(ADJUSTL(position))//'_um.dat'
   open(15,file=filename,status='unknown')
   do ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid
	write(15,19) mesh(2:(mesh_par%Nzm-1),(mesh_par%Nxm+1-ss))%Jpe_z
   enddo
   do ss=2,(mesh_par%Nxm-1),sim_parameters%jump_grid
	write(15,19) mesh(2:(mesh_par%Nzm-1),ss)%Jpe_z
   enddo
   close(15)

   return

   END SUBROUTINE


	 SUBROUTINE print_integrated_diagnostics
		   sim_parameters%IntDeltaOutput=abs(sim_parameters%sim_time*c-sim_parameters%IntLastOutput)
			 if(mod(sim_parameters%iter,sim_parameters%output_Integrated_params_nstep).eq.0 &
					  .or. sim_parameters%IntDeltaOutput>sim_parameters%output_Integrated_params_dist) then
						       call bunch_diagnostics_Integrated_AllBunches
						       if(sim_parameters%L_lineout) call lineout
						       sim_parameters%IntLastOutput=sim_parameters%sim_time*c
			endif
	 END SUBROUTINE print_integrated_diagnostics


	SUBROUTINE print_at_screen
			if(mod(sim_parameters%iter,10)==0) then
						write(*,'(A,I10,A,f12.3)') 'Iteration =',sim_parameters%iter,'    -    Run distance =',sim_parameters%zg
			endif
	END SUBROUTINE print_at_screen


	SUBROUTINE print_cpu_file
			logical :: exist
			character :: filename*99

			filename=trim(sim_parameters%out_dir)//'cpu'

			if(mod(sim_parameters%iter,200)==0) then
						inquire(file=trim(filename), exist=exist)
						if (exist) then
								open(99, file=trim(filename), status="old", position="append", action="write")
						else
						    open(99, file=trim(filename), status="new", action="write")
					  end if
				  	write(99,'(I10,f15.3)') sim_parameters%iter, sim_parameters%zg
						close(99)
			endif
	END SUBROUTINE print_cpu_file


END MODULE
