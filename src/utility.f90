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

MODULE utilities

USE my_types
USE use_my_types
USE pstruct_data
USE architect_class_structure

IMPLICIT NONE

CONTAINS

!======================================================================

subroutine open_file(operating_system,filename)
	integer, intent(in) 	:: operating_system
	character, intent(in) 	:: filename*90

	select case (operating_system)

		case (0) !mac
			open(11,file=TRIM(ADJUSTL(filename)),form='formatted', position='append')

		case (1) !win
			open(11,file=TRIM(sim_parameters%out_root)//TRIM(ADJUSTL(filename)),form='formatted', position='append')
	end select

end subroutine open_file

!---!


! -------------------------------------------------------------------------!
!                   Check particle inside moving window                    !
! -------------------------------------------------------------------------!
	subroutine inwindow(j,ip)
	integer, intent(in) :: j,ip !j is the number of the bunch, ip is the number of the particle within the j-th bunch
	real(8) pos_z,pos_r

		! in or out with respect of particle pusher
		if(bunch(j)%part(ip)%cmp(7)==1.) then
			pos_z = bunch(j)%part(ip)%cmp(3)-sim_parameters%zg
			pos_r = sqrt( (bunch(j)%part(ip)%cmp(1))**2 +(bunch(j)%part(ip)%cmp(2))**2 )
			pos_z = plasma%k_p * pos_z
			pos_r = plasma%k_p * pos_r

			if( pos_z .le. (minval(z_mesh_shifted)+20.*mesh_par%dzm)       )  bunch(j)%part(ip)%cmp(7)=0.
			if( pos_z .ge. (maxval(z_mesh_shifted)-20.*mesh_par%dzm)       )  bunch(j)%part(ip)%cmp(7)=0.
			if( pos_r .ge. maxval(bck_plasma%radius_um(:))-2.*mesh_par%dxm )  bunch(j)%part(ip)%cmp(7)=0.
	    endif

	end subroutine



! -------------------------------------------------------------------------!
!       from General Weights to Jr_w, Jz_w, rho_w                          !
! -------------------------------------------------------------------------!
	subroutine particle_weights(grid_choice,Wr,Wz,Wr_s,Wz_s,fraz,fraz_s,idr,idrs,r,w00,w10,w01,w11)
	!--> grid_choice=0 - Jr
	!--> grid_choice=1 - Jz
	!--> grid_choice=2 - rho
	!--> grid_choice=3 - Er
	!--> grid_choice=4 - Ez
	!--> grid_choice=5 - Bphi
	integer, intent(in) :: grid_choice,idr,idrs
	real(8), intent(in)  :: Wr,Wz,Wr_s,Wz_s,fraz,fraz_s,r
	real(8), intent(out) :: w00,w10,w01,w11

	w00=0.
	w01=0.
	w10=0.
	w11=0.

	select case (grid_choice)
		case (0) !--> grid_choice=0 - Jr
			if (idrs.eq.1) then !onaxis
				w00=Wz     *fraz_s*Wr_s   	/r 		!*inv_R_shifted(idrs)
				w01=Wz     *(1.-fraz_s*Wr_s)/r/2. 	!*inv_R_shifted(idrs)
				w10=(1.-Wz)*fraz_s*Wr_s		/r 		!*inv_R_shifted(idrs)
				w11=(1.-Wz)*(1.-fraz_s*Wr_s)/r/2. 	!*inv_R_shifted(idrs)
			else
				w00=Wz     *fraz_s*Wr_s   	*inv_R_shifted(idrs)
				w01=Wz     *(1.-fraz_s*Wr_s)*inv_R_shifted(idrs+1)
				w10=(1.-Wz)*fraz_s*Wr_s		*inv_R_shifted(idrs)
				w11=(1.-Wz)*(1.-fraz_s*Wr_s)*inv_R_shifted(idrs+1)
			endif

		case (1) !--> grid_choice=1 - Jz
			if ((idr.eq.2).and.(r.le.(0.5*mesh_par%dxm))) then !onaxis
				w00=Wz_s				  *inv_R(idr)
				w10=(1.-Wz_s)			  *inv_R(idr)
			else
				w00=Wz_s 	 *fraz*Wr     *inv_R(idr)
				w01=Wz_s     *(1.-fraz*Wr)*inv_R(idr+1)
				w10=(1.-Wz_s)*fraz*Wr     *inv_R(idr)
				w11=(1.-Wz_s)*(1.-fraz*Wr)*inv_R(idr+1)
			endif

		case (2) !--> grid_choice=2 - rho
			if (idrs.eq.1) then !onaxis
				w00=Wz_s	 *fraz_s*Wr_s		  *2./mesh_par%dxm*4. 	!*inv_R_shifted(idrs)
				w01=Wz_s	 *(1.-fraz_s*Wr_s)*inv_R_shifted(idrs+1)
				w10=(1.-Wz_s)*fraz_s*Wr_s	  *2./mesh_par%dxm*4.	!*inv_R_shifted(idrs)
				w11=(1.-Wz_s)*(1.-fraz_s*Wr_s)*inv_R_shifted(idrs+1)
			else
				w00=Wz_s	 *fraz_s*Wr_s		  *inv_R_shifted(idrs)
				w01=Wz_s	 *(1.-fraz_s*Wr_s)*inv_R_shifted(idrs+1)
				w10=(1.-Wz_s)*fraz_s*Wr_s	  *inv_R_shifted(idrs)
				w11=(1.-Wz_s)*(1.-fraz_s*Wr_s)*inv_R_shifted(idrs+1)
			endif

		case (3) !--> grid_choice=3 - Er
			if (idrs.eq.1) then !onaxis
				w00=Wz     *fraz_s*Wr_s
				w01=Wz     *(1.-fraz_s*Wr_s)/2.
				w10=(1.-Wz)*fraz_s*Wr_s
				w11=(1.-Wz)*(1.-fraz_s*Wr_s)/2.
			else
				w00=Wz     *fraz_s*Wr_s
				w01=Wz     *(1.-fraz_s*Wr_s)
				w10=(1.-Wz)*fraz_s*Wr_s
				w11=(1.-Wz)*(1.-fraz_s*Wr_s)
			endif

		case (4) !--> grid_choice=4 - Ez
			if ((idr.eq.2).and.(r.le.(0.5*mesh_par%dxm))) then !onaxis
				w00=Wz_s
				w10=(1.-Wz_s)
			else
				w00=Wz_s 	 *fraz*Wr
				w01=Wz_s     *(1.-fraz*Wr)
				w10=(1.-Wz_s)*fraz*Wr
				w11=(1.-Wz_s)*(1.-fraz*Wr)
			endif

		case (5) !--> grid_choice=5 - Bphi
			if (idrs.eq.1) then !onaxis
				w00=Wz_s	 *fraz*Wr_s
				w01=Wz_s	 *(1.-fraz_s*Wr_s)/2.
				w10=(1.-Wz_s)*fraz_s*Wr_s
				w11=(1.-Wz_s)*(1.-fraz_s*Wr_s)/2.
			else
				w00=Wz_s	 *fraz*Wr_s
				w01=Wz_s	 *(1.-fraz_s*Wr_s)
				w10=(1.-Wz_s)*fraz_s*Wr_s
				w11=(1.-Wz_s)*(1.-fraz_s*Wr_s)
			endif

	end select
	end subroutine

! -------------------------------------------------------------------------!
!                   Square weighting, Cylindrical mesh                     !
! -------------------------------------------------------------------------!
	subroutine weights_particle_ongrid(routine_choice,xnew,ynew,znew,xold,yold,zold,Wr,Wz,Wr_s,Wz_s,fraz,fraz_s,pos_r,pos_z,indx,indx_s,indz,indz_s)
		!check the includes
		!--> routine_choice=0 - charge and current deposition, the position is interpolated between new and old position
		!--> routine_choice=1 - pusher, force weighting from grid using present positions
		integer, intent(in) :: routine_choice
		real(8), intent(in)  :: xnew,ynew,znew,xold,yold,zold
		real(8), intent(out) :: Wr,Wz,Wr_s,Wz_s,fraz,fraz_s
		real(8), intent(inout) :: pos_r,pos_z
		integer, intent(inout) :: indx,indx_s
		integer, intent(inout) :: indz, indz_s
		real(8) ::  f_pos_r,f_pos_rs,f_pos_z,f_pos_zs

		!weights cleanup
		Wr    = 0.D0
		Wr_s  = 0.D0
		Wz    = 0.D0
		Wz_s  = 0.D0
		fraz  = 0.D0
		fraz_s= 0.D0
		indx  = 0
		indx_s= 0
		pos_r = 0.D0
		indz  = 0
		indz_s= 0

  		!RZ coordinates
		pos_z = plasma%k_p*     ( 0.5D0*(znew+zold)-sim_parameters%zg )
		pos_r = plasma%k_p* sqrt( (xnew+xold)**2+(ynew+yold)**2 )*0.5D0

		select case (routine_choice)
			case (0) !--> routine_choice=0 - charge and current deposition
				pos_z = plasma%k_p*     ( 0.5D0*(znew+zold)-sim_parameters%zg )
				pos_r = plasma%k_p* sqrt( (xnew+xold)**2+(ynew+yold)**2 )*0.5D0
			case (1) !--> routine_choice=1 - pusher
				pos_z = plasma%k_p*     (  znew-sim_parameters%zg )
				pos_r = plasma%k_p* sqrt( (xnew)**2+(ynew)**2 )
		end select

		f_pos_r  = (pos_r-Xmin)        *one_over_dx
		f_pos_rs = (pos_r-Xmin_shifted)*one_over_dx
		f_pos_z  = (pos_z-Zmin)        *one_over_dz
		f_pos_zs = (pos_z-Zmin_shifted)*one_over_dz


		! Indices in moving window
		indx	 = 1 + int(f_pos_r)
		indx_s	 = 1 + int(f_pos_rs) !shifted x-grid
		indz	 = 1 + int(f_pos_z)
		indz_s	 = 1 + int(f_pos_zs) !shifted z-grid

		!BC r=0
		if (indx.eq.1) indx=2

		Wz    = 1.D0-(f_pos_z - indz+1.D0 )
		Wz_s  = 1.D0-(f_pos_zs - indz_s+1.D0 )
		Wr_s  = 1.D0-(f_pos_rs - indx_s+1.D0 )

		fraz_s = (2.D0*x_mesh_shifted(indx_s)*one_over_dx + Wr_s) &
				/ (2.D0*x_mesh_shifted(indx_s)*one_over_dx + 1.D0  )


		if ((indx.eq.2).and.(pos_r.le.(0.5D0*mesh_par%dxm))) then
			!do nothing
		else
			Wr    = 1.D0-(f_pos_r - indx+1.D0 )
			fraz   = (2.D0*x_mesh(indx)*one_over_dx + Wr ) &
					/ (2.D0*x_mesh(indx)*one_over_dx + 1.D0 )
		endif


			if ((Wz.gt.1D0).or.(Wz_s.gt.1D0)) then
				write(*,*) 'Error in Weighting along z'
				write(*,*) 'wz',Wz
				write(*,*) 'wz_s',Wz_s
				write(*,*) '-----------------'
				write(*,*) 'pos_z',pos_z
				write(*,*) 'dz',mesh_par%dzm
				write(*,*) 'indz',indz
				write(*,*) 'indz_s',indz_s
				write(*,*) '           '
				write(*,*) 'z_mesh(indz-1)',z_mesh(indz-1)
				write(*,*) 'z_mesh(indz)',z_mesh(indz)
				write(*,*) 'z_mesh(indz+1)',z_mesh(indz+1)
				write(*,*) '           '
				write(*,*) 'z_mesh_shifted(indz_s-1)',z_mesh_shifted(indz_s-1)
				write(*,*) 'z_mesh_shifted(indz_s)',z_mesh_shifted(indz_s)
				write(*,*) 'z_mesh_shifted(indz_s+1)',z_mesh_shifted(indz_s+1)
				write(*,*) '-----------------'
				!stop
				!--- forcing ---!
				if ( Wz  .gt.1D0 ) Wz  =1.D0
				if ( Wz_s.gt.1D0 ) Wz_s=1.D0
            endif

			if (Wr_s.gt.1D0) then
				write(*,*) 'Error in Force Weighting for x_shifted'
				write(*,*) 'wr_s',Wr_s
				write(*,*) '-----------------'
				write(*,*) 'pos_r',pos_r
				write(*,*) 'dx',mesh_par%dxm
				write(*,*) 'indx_s',indx_s
				write(*,*) '           '
				write(*,*) 'x_mesh_shifted(indx_s-1)',x_mesh_shifted(indx_s-1)
				write(*,*) 'x_mesh_shifted(indx_s)',x_mesh_shifted(indx_s)
				write(*,*) 'x_mesh_shifted(indx_s+1)',x_mesh_shifted(indx_s+1)
				write(*,*) '           '
				write(*,*) '-----------------'
				!stop
				!--- forcing ---!
				if ( Wr_s.gt.1D0 ) Wr_s=1.D0
			endif

			if (Wr.gt.1D0) then
				write(*,*) 'Error in Force Weighting for x'
				write(*,*) 'wr',Wr
				write(*,*) '-----------------'
				write(*,*) 'pos_r',pos_r
				write(*,*) 'dx',mesh_par%dxm
				write(*,*) 'indx',indx
				write(*,*) '           '
				write(*,*) 'x_mesh(indx-1)',x_mesh(indx-1)
				write(*,*) 'x_mesh(indx)',x_mesh(indx)
				write(*,*) 'x_mesh(indx+1)',x_mesh(indx+1)
				write(*,*) '           '
				write(*,*) '-----------------'
				!stop
				!--- forcing ---!
				if ( Wr  .gt.1D0 ) Wr  =1.D0
			endif

			if ((indx.lt.1).or.(indx_s.lt.1).or.(indz.lt.1).or.(indz_s.lt.1)) then
				write(*,*) 'Error in index'
				write(*,*) '-----------------'
				write(*,*) 'pos_z',pos_z
				write(*,*) 'dz',mesh_par%dzm
				write(*,*) 'indz',indz
				write(*,*) 'indz_s',indz_s
				write(*,*) 'pos_r',pos_r
				write(*,*) 'indz',indx
				write(*,*) 'indz_s',indx_s
				write(*,*) '           '
				write(*,*) 'x_mesh(indx-1)',x_mesh(indx-1)
				write(*,*) 'x_mesh(indx)',x_mesh(indx)
				write(*,*) 'x_mesh(indx+1)',x_mesh(indx+1)
				write(*,*) '           '
				write(*,*) 'x_mesh_shifted(indx_s-1)',x_mesh_shifted(indx_s-1)
				write(*,*) 'x_mesh_shifted(indx_s)',x_mesh_shifted(indx_s)
				write(*,*) 'x_mesh_shifted(indx_s+1)',x_mesh_shifted(indx_s+1)
				write(*,*) '           '
				write(*,*) 'z_mesh(indz-1)',z_mesh(indz-1)
				write(*,*) 'z_mesh(indz)',z_mesh(indz)
				write(*,*) 'z_mesh(indz+1)',z_mesh(indz+1)
				write(*,*) '           '
				write(*,*) 'z_mesh_shifted(indz_s-1)',z_mesh_shifted(indz_s-1)
				write(*,*) 'z_mesh_shifted(indz_s)',z_mesh_shifted(indz_s)
				write(*,*) 'z_mesh_shifted(indz_s+1)',z_mesh_shifted(indz_s+1)
				write(*,*) '-----------------'
				stop
	  endif

	end subroutine


	!--- *** ---!
		RECURSIVE FUNCTION Factorial(n)  RESULT(Fact)
			INTEGER :: Fact
			INTEGER, INTENT(IN) :: n
			IF (n == 0) THEN
			   Fact = 1
			ELSE
			   Fact = n * Factorial(n-1)
			END IF
		END FUNCTION Factorial
	!---***

	!--- *** ---!
		real(8) FUNCTION from_dimlessE_to_dimE(E)
			real(8), INTENT(IN) :: E
			from_dimlessE_to_dimE = E*96. * sqrt( plasma%n0 )
		END FUNCTION from_dimlessE_to_dimE
	!---***


END MODULE utilities
