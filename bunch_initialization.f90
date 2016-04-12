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

 module initialize_bunch


	USE my_types
	USE use_my_types
	USE pstruct_data
	USE architect_class_structure
	USE bunch_generation
	USE moments


 implicit none

 !--- --- ---!
 contains
 !--- --- ---!

 subroutine init_bunch

	INTEGER :: nb,i,j
	REAL(8) zmmedB,Zb
	REAL(8), DIMENSION(:,:), ALLOCATABLE ::  bunch_init
	CHARACTER :: name_file*255


	write(*,*) "initializing the bunches"
	! Bunch(es) memory allocation
   DO nb = 1,bunch_initialization%n_total_bunches
	allocate(   bunch(nb)%part(bunch_initialization%n_particles(nb)) )! ,STAT=AllocStatus)
   ENDDO


    !init-bunch(es)
   if(bunch_initialization%l_bunch_internal_init) then

	do i=1,bunch_initialization%n_total_bunches

			allocate( bunch_init(6,bunch_initialization%n_particles(i)) )

			sim_parameters%rB0(i)    = bunch_initialization%bunch_s_x(i)
			sim_parameters%lbunch(i) = bunch_initialization%bunch_s_z(i)

			if(.not.twiss%L_TWISS) then
				call generate_bunch( &
				0.D0 , 0.D0, 0.D0,&
				bunch_initialization%bunch_s_x(i), &
				bunch_initialization%bunch_s_y(i), &
				bunch_initialization%bunch_s_z(i), &
				bunch_initialization%bunch_gamma_m(i),&
				bunch_initialization%bunch_eps_x(i), &
				bunch_initialization%bunch_eps_y(i),bunch_initialization%bunch_dgamma(i),&
				bunch_initialization%n_particles(i),bunch_init)
			else !bunch with twiss
				call generate_bunch_twiss( &
				0.D0 , 0.D0, 0.D0,&
				bunch_initialization%bunch_s_x(i), &
				bunch_initialization%bunch_s_y(i), &
				bunch_initialization%bunch_s_z(i),&
				bunch_initialization%bunch_gamma_m(i),&
				bunch_initialization%bunch_eps_x(i), &
				bunch_initialization%bunch_eps_y(i),bunch_initialization%bunch_dgamma(i),&
				bunch_initialization%n_particles(i),bunch_init, &
				twiss%alpha_new_factor(i),  &
				twiss%beta_new_factor(i)   )
			endif

			!il nome di questa funzione e' incomprensibile
			!call modify_bunch(bunch_init,i)

			do j=1,bunch_initialization%n_particles(i)
				bunch(i)%part(j)%cmp(1)=bunch_init(1,j)
				bunch(i)%part(j)%cmp(2)=bunch_init(2,j)
				bunch(i)%part(j)%cmp(3)=bunch_init(3,j)
				bunch(i)%part(j)%cmp(4)=bunch_init(4,j)!-bunch_init(6,j)*bunch_initialization%f0(i)*cos(atan2(bunch_init(2,j),bunch_init(1,j)))
				bunch(i)%part(j)%cmp(5)=bunch_init(5,j)!-bunch_init(6,j)*bunch_initialization%f0(i)*sin(atan2(bunch_init(2,j),bunch_init(1,j)))
				bunch(i)%part(j)%cmp(6)=-bunch_init(6,j) !backward velocity
				bunch(i)%part(j)%cmp(7)=1.
				bunch(i)%part(j)%cmp(8)=1.

				bunch(i)%part(j)%cmp(9)=bunch_init(1,j)  !Xold
				bunch(i)%part(j)%cmp(10)=bunch_init(2,j) !Yold
				bunch(i)%part(j)%cmp(11)=bunch_init(3,j) !Zold
			enddo

			if (allocated(bunch_init)) deallocate(bunch_init)
	enddo

   else !read from external file

	do i=1,bunch_initialization%n_total_bunches

		allocate( bunch_init(6,bunch_initialization%n_particles(i)) )

		name_file = sim_parameters%inbunch(i)
		call read_from_external_file(name_file, bunch_initialization%n_particles(i),bunch_init)

            if (twiss%beta_new_factor(i).le.0) then

                write(*,*) 'Non positive new beta function value.'
                write(*,*) 'Doing nothing on bunch', i

            else

                call modify_bunch(bunch_init,i)

            endif

		do j=1,bunch_initialization%n_particles(i)
			bunch(i)%part(j)%cmp(1)=bunch_init(1,j)
			bunch(i)%part(j)%cmp(2)=bunch_init(2,j)
			bunch(i)%part(j)%cmp(3)=bunch_init(3,j)
			bunch(i)%part(j)%cmp(4)=bunch_init(4,j)
			bunch(i)%part(j)%cmp(5)=bunch_init(5,j)
			bunch(i)%part(j)%cmp(6)=bunch_init(6,j)
			bunch(i)%part(j)%cmp(7)=1.
			bunch(i)%part(j)%cmp(8)=1.

			bunch(i)%part(j)%cmp(9)=bunch_init(1,j)  !Xold
			bunch(i)%part(j)%cmp(10)=bunch_init(2,j) !Yold
			bunch(i)%part(j)%cmp(11)=bunch_init(3,j) !Zold
		enddo
		if (allocated(bunch_init)) deallocate(bunch_init)

		sim_parameters%rB0(i)    = calculate_nth_central_moment_bunch(i,2,1)
		sim_parameters%lbunch(i) = calculate_nth_central_moment_bunch(i,2,3)
	enddo

   endif

!------------------------------------------------------!

   Zb=0.
   !~if  (bunch_initialization%n_total_bunches.gt.1) then
   !~	do i=1,(bunch_initialization%n_total_bunches-1),-1
   !~		bunch_initialization%db(i+1)=bunch_initialization%db(i)
   !~	enddo
	!~bunch_initialization%db(1) = 0.
   !~endif
   !~bunch_initialization%db(1) = 0.

   do i=1,bunch_initialization%n_total_bunches
      zmmedB                    		= calculate_nth_moment_bunch(i,1,3)
	  if (i.ge.1) then
		Zb                          	= Zb+bunch_initialization%db(i)*plasma%lambda_p
	  else if (i.eq.1) then
		Zb = 0.
	  endif
      bunch(i)%part(:)%cmp(3)			= bunch(i)%part(:)%cmp(3)- zmmedB + Zb 							! shift Z of bunch
      bunch(i)%part(:)%cmp(11)			= bunch(i)%part(:)%cmp(3)  										! change Z_old accordingly
	  bunch(i)%part(:)%cmp(12)			= bunch_initialization%ChargeB(i)/bunch_initialization%n_particles(i) ! macroparticle charge
	  bunch(i)%part(:)%cmp(13)			= 1.e10*bunch(i)%part(:)%cmp(12)/1.6021766 						! electrons per macroparticle
   enddo

   !---initial diagnostic
      do i=1,bunch_initialization%n_total_bunches
         write(*,*) 'Bunch(',i,')  sigma_x (um) =',calculate_nth_central_moment_bunch(i,2,1)
         write(*,*) 'Bunch(',i,')  sigma_z (um) =',calculate_nth_central_moment_bunch(i,2,3)
         write(*,*)
      enddo
      write(*,*)'Plasma wavelength (um) =',plasma%lambda_p
      write(*,*)


   ! Computes total charge and number of electrons
   sim_parameters%Charge=sum(bunch_initialization%ChargeB(1:bunch_initialization%n_total_bunches))
   sim_parameters%NelectronsB(1:bunch_initialization%n_total_bunches) = &
    1.e10*bunch_initialization%ChargeB(1:bunch_initialization%n_total_bunches)/1.6021766 !Electrons of every macroparticle
	!1.e10*bunch_initialization%ChargeB(1:bunch_initialization%n_total_bunches)/1.6 !Electrons of every macroparticle



	 write(*,*) "bunches initialized"
	 write(*,*)



 end subroutine init_bunch




 subroutine modify_bunch(in_bunch,nb)
    integer nb
    real(8) :: in_bunch(6,bunch_initialization%n_particles(nb)),stats(9)
    real(8) :: Tmatrix_X(2,2),Tmatrix_Y(2,2)

    call def_stats(stats,in_bunch,nb)

    write(*,*)
    write(*,*) 'BUNCH',nb
    write(*,*)
    write(*,*) 'Sx_old=', stats(1), 'Spx_old=', stats(2)
    write(*,*) 'Sxpx_old=', stats(3), 'enx_old=',stats(7)
    write(*,*) 'Sy_old=', stats(4), 'Spy_old=', stats(5)
    write(*,*) 'Sypy_old=', stats(6), 'enx_old=',stats(8)
    write(*,*) 'Alphax_old=', -stats(3)/stats(7), 'Betax_old=', stats(9)*stats(1)**2/stats(7)
    write(*,*) 'Alphay_old=', -stats(6)/stats(8), 'Betay_old=', stats(9)*stats(4)**2/stats(8)
    write(*,*)

    !call define_transport_matrices(Tmatrix_X,Tmatrix_Y,stats,nb)

    !call set_twiss_parameters(in_bunch,Tmatrix_X,Tmatrix_Y,nb)

    !call def_stats(stats,in_bunch,nb)

    write(*,*) 'Sx_new=', stats(1), 'Spx_new=', stats(2)
    write(*,*) 'Sxpx_new=', stats(3), 'enx_new=',stats(7)
    write(*,*) 'Sy_new=', stats(4), 'Spy_new=', stats(5)
    write(*,*) 'Sypy_new=', stats(6), 'enx_new=',stats(8)
    write(*,*) 'Alphax_new=', -stats(3)/stats(7), 'Betax_new=', stats(9)*stats(1)**2/stats(7)
    write(*,*) 'Alphay_new=', -stats(6)/stats(8), 'Betay_new=', stats(9)*stats(4)**2/stats(8)
    write(*,*)

    return

 end subroutine modify_bunch







 subroutine def_stats(stats,in_bunch,nb)
! Calculates the beam parameters sigma_x,y, sigma_px,py and sigma_xpx,ypy
! Structure of vector stats is (Sx,Spx,Sxpx,Sy,Spy,Sypy,enx,eny)
    integer nb,np
    real(8) :: in_bunch(6,bunch_initialization%n_particles(nb)),stats(9)
    real(8) :: mx=0.,mpx=0.,my=0.,mpy=0.

    stats=0.
    np=bunch_initialization%n_particles(nb)


    mx=sum(in_bunch(1,:))/real(np)
    mpx=sum(in_bunch(4,:))/real(np)
    my=sum(in_bunch(2,:))/real(np)
    mpy=sum(in_bunch(5,:))/real(np)


    stats(1)=sum((in_bunch(1,:)-mx)**2)/real(np)
    stats(1)=sqrt(stats(1))
    stats(2)=sum((in_bunch(4,:)-mpx)**2)/real(np)
    stats(2)=sqrt(stats(2))
    stats(3)=sum((in_bunch(1,:)-mx)*(in_bunch(4,:)-mpx))/real(np)
    stats(4)=sum((in_bunch(2,:)-mx)**2)/real(np)
    stats(4)=sqrt(stats(4))
    stats(5)=sum((in_bunch(5,:)-mpx)**2)/real(np)
    stats(5)=sqrt(stats(5))
    stats(6)=sum((in_bunch(2,:)-mx)*(in_bunch(5,:)-mpx))/real(np)
    stats(7)=sqrt(stats(1)**2*stats(2)**2-stats(3)**2)
    stats(8)=sqrt(stats(4)**2*stats(5)**2-stats(6)**2)
    stats(9)=sqrt(1.+(dot_product(in_bunch(4,:),in_bunch(4,:))+dot_product(in_bunch(5,:),in_bunch(5,:))+dot_product(in_bunch(6,:),in_bunch(6,:)))/real(np))

    return

 end subroutine def_stats








 subroutine define_transport_matrices(Tx,Ty,stats,nb)
! Defines the two 2x2 matrices Tx and Ty. Once applied to the (x,px) and (y,py) vectors respectively,
! they realize the transformation on the Twiss parameters alpha->t*alpha and beta->r*beta at constant emittance.
! The Ti(2,1)=0 choice is customary and arbitrary.
    integer nb
    real :: Tx(2,2),Ty(2,2),stats(9)
    real :: Kx=0,Ky=0,r,t

    Tx=0
    Ty=0
    r=twiss%beta_new_factor(nb)
    t=twiss%alpha_new_factor(nb)
    Kx=stats(1)*stats(2)*sqrt(r)/sqrt(stats(1)**2*stats(2)**2+(t**2-1)*stats(3)**2)
    Ky=stats(4)*stats(5)*sqrt(r)/sqrt(stats(4)**2*stats(5)**2+(t**2-1)*stats(6)**2)

    Tx(1,1)=Kx
    Tx(1,2)=Kx*(t-1.)*stats(3)/stats(2)**2
    Tx(2,1)=0
    Tx(2,2)=1./Kx

    Ty(1,1)=Ky
    Ty(1,2)=Ky*(t-1.)*stats(6)/stats(5)**2
    Ty(2,1)=0
    Ty(2,2)=1./Ky

    return

 end subroutine define_transport_matrices




 subroutine set_twiss_parameters(in_bunch,Tx,Ty,nb)
! Applies the Tx and Ty transport matrices
   integer i,nb
   real in_bunch(6,bunch_initialization%n_particles(nb)),out_bunch(6,bunch_initialization%n_particles(nb)),Tx(2,2),Ty(2,2)

   out_bunch=in_bunch

   do i=1,bunch_initialization%n_particles(nb)
      out_bunch(1,i)=Tx(1,1)*in_bunch(1,i)+Tx(1,2)*in_bunch(4,i)
      out_bunch(4,i)=Tx(2,1)*in_bunch(1,i)+Tx(2,2)*in_bunch(4,i)
      out_bunch(2,i)=Ty(1,1)*in_bunch(2,i)+Ty(1,2)*in_bunch(5,i)
      out_bunch(5,i)=Ty(2,1)*in_bunch(2,i)+Ty(2,2)*in_bunch(5,i)
   enddo

   in_bunch=out_bunch

   return

end subroutine set_twiss_parameters






 !--- --- ---!
 end module initialize_bunch
 !--- --- ---!
