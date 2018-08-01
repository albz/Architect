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
	REAL(8) mu_z,Zb
	REAL(8), DIMENSION(:,:), ALLOCATABLE ::  bunch_init
  REAL(8) :: alphaTwiss, betaTwiss, ax11,ay11,ax12,ay12,s_x,s_y,eps_x,eps_y
	CHARACTER :: name_file*255


  write(*,'(A)') ' --- Bunch(es) Initialisation ---'

    !init-bunch(es)   --- Internal init
    if(bunch_initialization%l_bunch_internal_init) then

      !here I need a funciton to calculate the number of particles, so I can inialize the class
      do nb = 1,bunch_initialization%n_total_bunches
        bunch_initialization%n_particles(nb) = calculate_bunch_number_of_particles( &
                                                nb,bunch_initialization%shape(nb), &
                                                bunch_initialization%PWeights(nb), &
                                                bunch_initialization%n_particles(nb))
        allocate(   bunch(nb)%part(bunch_initialization%n_particles(nb)) )
        write(*,'(A,I1,A,I8,A)') 'Bunch(',nb,') :: initialised with >>> ', bunch_initialization%n_particles(nb),' particles'
      enddo


      do i=1,bunch_initialization%n_total_bunches
        allocate( bunch_init(8,bunch_initialization%n_particles(i)) )
        sim_parameters%rB0(i)    = bunch_initialization%bunch_s_x(i)
  			sim_parameters%lbunch(i) = bunch_initialization%bunch_s_z(i)

        !--- BiGaussian + Equal particles + no optimisation
        if(   trim(bunch_initialization%shape(i))=='bigaussian'  &
        .and. trim(bunch_initialization%PWeights(i))=='equal'    &
        .and. trim(bunch_initialization%optimisation(i))=='no')  &
                call generate_bunch_bigaussian_equal(i, 0.D0 , 0.D0, 0.D0,&
              				bunch_initialization%bunch_s_x(i), &
              				bunch_initialization%bunch_s_y(i), &
              				bunch_initialization%bunch_s_z(i), &
              				bunch_initialization%bunch_gamma_m(i),&
              				bunch_initialization%bunch_eps_x(i), &
              				bunch_initialization%bunch_eps_y(i),bunch_initialization%bunch_dgamma(i),&
              				bunch_initialization%n_particles(i),bunch_initialization%ChargeB(i), &
                      bunch_initialization%sigma_cut(i))

        !--- BiGaussian + Equal particles + yes optimisation
        if(   trim(bunch_initialization%shape(i))=='bigaussian'  &
        .and. trim(bunch_initialization%PWeights(i))=='equal'    &
        .and. trim(bunch_initialization%optimisation(i))=='yes')  &
                call generate_bunch_bigaussian_equal_optimised(i, 0.D0 , 0.D0, 0.D0,&
              				bunch_initialization%bunch_s_x(i), &
              				bunch_initialization%bunch_s_y(i), &
              				bunch_initialization%bunch_s_z(i), &
              				bunch_initialization%bunch_gamma_m(i),&
              				bunch_initialization%bunch_eps_x(i), &
              				bunch_initialization%bunch_eps_y(i),bunch_initialization%bunch_dgamma(i),&
              				bunch_initialization%n_particles(i),bunch_initialization%ChargeB(i), &
                      bunch_initialization%sigma_cut(i))

        !--- BiGaussian + weighted particles + no optimisation
        if(   trim(bunch_initialization%shape(i))=='bigaussian'  &
        .and. trim(bunch_initialization%PWeights(i))=='weighted'    &
        .and. trim(bunch_initialization%optimisation(i))=='no')  &
                call generate_bunch_bigaussian_weighted(i, 0.D0 , 0.D0, 0.D0,&
                      bunch_initialization%bunch_s_x(i), &
                      bunch_initialization%bunch_s_y(i), &
                      bunch_initialization%bunch_s_z(i), &
                      bunch_initialization%bunch_gamma_m(i),&
                      bunch_initialization%bunch_eps_x(i), &
                      bunch_initialization%bunch_eps_y(i),bunch_initialization%bunch_dgamma(i),&
                      bunch_initialization%n_particles(i),bunch_initialization%ChargeB(i),&
                      bunch_initialization%npZ(i), &
                      bunch_initialization%npR(i), bunch_initialization%sigma_cut(i))

        !--- BiGaussian + weighted particles + yes optimisation
        if(   trim(bunch_initialization%shape(i))=='bigaussian'  &
        .and. trim(bunch_initialization%PWeights(i))=='weighted'    &
        .and. trim(bunch_initialization%optimisation(i))=='yes')  &
                call generate_bunch_bigaussian_weighted_optimised(i, 0.D0 , 0.D0, 0.D0,&
                      bunch_initialization%bunch_s_x(i), &
                      bunch_initialization%bunch_s_y(i), &
                      bunch_initialization%bunch_s_z(i), &
                      bunch_initialization%bunch_gamma_m(i),&
                      bunch_initialization%bunch_eps_x(i), &
                      bunch_initialization%bunch_eps_y(i),bunch_initialization%bunch_dgamma(i),&
                      bunch_initialization%n_particles(i),bunch_initialization%ChargeB(i),&
                      bunch_initialization%npZ(i), &
                      bunch_initialization%npR(i), bunch_initialization%sigma_cut(i))


        !--- cylindrical + equal particles + no optimisation
        if(   trim(bunch_initialization%shape(i))=='cylindrical'  &
        .and. trim(bunch_initialization%PWeights(i))=='equal'     &
        .and. trim(bunch_initialization%optimisation(i))=='no')   &
                      call generate_bunch_trapezoidalZ_uniformR_equal( &
                          i, 0.D0 , 0.D0, 0.D0,&
                          bunch_initialization%bunch_s_x(i), &
                          bunch_initialization%bunch_s_y(i), &
                          bunch_initialization%bunch_s_z(i),&
                          bunch_initialization%bunch_gamma_m(i),&
                          bunch_initialization%bunch_eps_x(i), &
                          bunch_initialization%bunch_eps_y(i),bunch_initialization%bunch_dgamma(i),&
                          bunch_initialization%n_particles(i),&
                          bunch_initialization%Charge_right(i),bunch_initialization%Charge_left(i), &
                          bunch_initialization%ChargeB(i))

          !--- cylindrical + equal particles + yes optimisation
          if(   trim(bunch_initialization%shape(i))=='cylindrical'  &
          .and. trim(bunch_initialization%PWeights(i))=='equal'     &
          .and. trim(bunch_initialization%optimisation(i))=='yes')   &
                        call generate_bunch_trapezoidalZ_uniformR_equal_optimised( &
                            i, 0.D0 , 0.D0, 0.D0,&
                            bunch_initialization%bunch_s_x(i), &
                            bunch_initialization%bunch_s_y(i), &
                            bunch_initialization%bunch_s_z(i),&
                            bunch_initialization%bunch_gamma_m(i),&
                            bunch_initialization%bunch_eps_x(i), &
                            bunch_initialization%bunch_eps_y(i),bunch_initialization%bunch_dgamma(i),&
                            bunch_initialization%n_particles(i),&
                            bunch_initialization%Charge_right(i),bunch_initialization%Charge_left(i), &
                            bunch_initialization%ChargeB(i))

            !--- cylindrical + weighted particles + no optimisation
            if(   trim(bunch_initialization%shape(i))=='cylindrical'  &
            .and. trim(bunch_initialization%PWeights(i))=='weighted'     &
            .and. trim(bunch_initialization%optimisation(i))=='no')   &
                          call generate_bunch_trapezoidalZ_uniformR_weighted( &
                              i, 0.D0 , 0.D0, 0.D0,&
                              bunch_initialization%bunch_s_x(i), &
                              bunch_initialization%bunch_s_y(i), &
                              bunch_initialization%bunch_s_z(i),&
                              bunch_initialization%bunch_gamma_m(i),&
                              bunch_initialization%bunch_eps_x(i), &
                              bunch_initialization%bunch_eps_y(i),bunch_initialization%bunch_dgamma(i),&
                              bunch_initialization%n_particles(i),&
                              bunch_initialization%Charge_right(i),bunch_initialization%Charge_left(i), &
                              bunch_initialization%ChargeB(i), &
                              bunch_initialization%npZ(i),bunch_initialization%npR(i))

              !--- cylindrical + weighted particles + yes optimisation
              if(   trim(bunch_initialization%shape(i))=='cylindrical'  &
              .and. trim(bunch_initialization%PWeights(i))=='weighted'     &
              .and. trim(bunch_initialization%optimisation(i))=='yes')   &
                            call generate_bunch_trapezoidalZ_uniformR_weighted_optimised( &
                                i, 0.D0 , 0.D0, 0.D0,&
                                bunch_initialization%bunch_s_x(i), &
                                bunch_initialization%bunch_s_y(i), &
                                bunch_initialization%bunch_s_z(i),&
                                bunch_initialization%bunch_gamma_m(i),&
                                bunch_initialization%bunch_eps_x(i), &
                                bunch_initialization%bunch_eps_y(i),bunch_initialization%bunch_dgamma(i),&
                                bunch_initialization%n_particles(i),&
                                bunch_initialization%Charge_right(i),bunch_initialization%Charge_left(i), &
                                bunch_initialization%ChargeB(i), &
                                bunch_initialization%npZ(i),bunch_initialization%npR(i))



                !--- trapezoidal + equal particles + no optimisation
                if(   trim(bunch_initialization%shape(i))=='trapezoidal'  &
                .and. trim(bunch_initialization%PWeights(i))=='equal'     &
                .and. trim(bunch_initialization%optimisation(i))=='no')   &
                              call generate_bunch_trapezoidalZ_gaussianR_equal( &
                                    i, 0.D0 , 0.D0, 0.D0,&
                            				bunch_initialization%bunch_s_x(i), &
                            				bunch_initialization%bunch_s_y(i), &
                            				bunch_initialization%bunch_s_z(i),&
                            				bunch_initialization%bunch_gamma_m(i),&
                            				bunch_initialization%bunch_eps_x(i), &
                            				bunch_initialization%bunch_eps_y(i),bunch_initialization%bunch_dgamma(i),&
                            				bunch_initialization%n_particles(i),&
                                    bunch_initialization%Charge_right(i),bunch_initialization%Charge_left(i), &
                                    bunch_initialization%sigma_cut(i),bunch_initialization%ChargeB(i))

                  !--- trapezoidal + equal particles + yes optimisation
                  if(   trim(bunch_initialization%shape(i))=='trapezoidal'  &
                  .and. trim(bunch_initialization%PWeights(i))=='equal'     &
                  .and. trim(bunch_initialization%optimisation(i))=='yes')   &
                                call generate_bunch_trapezoidalZ_gaussianR_equal_optimised( &
                                      i, 0.D0 , 0.D0, 0.D0,&
                              				bunch_initialization%bunch_s_x(i), &
                              				bunch_initialization%bunch_s_y(i), &
                              				bunch_initialization%bunch_s_z(i),&
                              				bunch_initialization%bunch_gamma_m(i),&
                              				bunch_initialization%bunch_eps_x(i), &
                              				bunch_initialization%bunch_eps_y(i),bunch_initialization%bunch_dgamma(i),&
                              				bunch_initialization%n_particles(i),&
                                      bunch_initialization%Charge_right(i),bunch_initialization%Charge_left(i), &
                                      bunch_initialization%sigma_cut(i),bunch_initialization%ChargeB(i))

                    !--- trapezoidal + weighted particles + no optimisation
                    if(   trim(bunch_initialization%shape(i))=='trapezoidal'  &
                    .and. trim(bunch_initialization%PWeights(i))=='weighted'  &
                    .and. trim(bunch_initialization%optimisation(i))=='no')   &
                                  call generate_bunch_trapezoidalZ_gaussianR_weighted( &
                                        i, 0.D0 , 0.D0, 0.D0,&
                                        bunch_initialization%bunch_s_x(i), &
                                        bunch_initialization%bunch_s_y(i), &
                                        bunch_initialization%bunch_s_z(i),&
                                        bunch_initialization%bunch_gamma_m(i),&
                                        bunch_initialization%bunch_eps_x(i), &
                                        bunch_initialization%bunch_eps_y(i),bunch_initialization%bunch_dgamma(i),&
                                        bunch_initialization%n_particles(i),&
                                        bunch_initialization%Charge_right(i),bunch_initialization%Charge_left(i), &
                                        bunch_initialization%ChargeB(i), &
                                        bunch_initialization%npZ(i), bunch_initialization%npR(i), &
                                        bunch_initialization%sigma_cut(i))

                      !--- trapezoidal + weighted particles + yes optimisation
                      if(   trim(bunch_initialization%shape(i))=='trapezoidal'  &
                      .and. trim(bunch_initialization%PWeights(i))=='weighted'     &
                      .and. trim(bunch_initialization%optimisation(i))=='yes')   &
                                    call generate_bunch_trapezoidalZ_gaussianR_weighted_optimised( &
                                        i, 0.D0 , 0.D0, 0.D0,&
                                        bunch_initialization%bunch_s_x(i), &
                                        bunch_initialization%bunch_s_y(i), &
                                        bunch_initialization%bunch_s_z(i),&
                                        bunch_initialization%bunch_gamma_m(i),&
                                        bunch_initialization%bunch_eps_x(i), &
                                        bunch_initialization%bunch_eps_y(i),bunch_initialization%bunch_dgamma(i),&
                                        bunch_initialization%n_particles(i),&
                                        bunch_initialization%Charge_right(i),bunch_initialization%Charge_left(i), &
                                        bunch_initialization%ChargeB(i), &
                                        bunch_initialization%npZ(i), bunch_initialization%npR(i), &
                                        bunch_initialization%sigma_cut(i))

      !---!
			if (allocated(bunch_init)) deallocate(bunch_init)
	enddo

      else !read from external file
            do i=1,bunch_initialization%n_total_bunches
                name_file = bunch_initialization%inbunch(i)
                write(*,*) name_file
                call read_header_external_file(name_file,bunch_initialization%chargeB(i),bunch_initialization%n_particles(i))
                bunch_initialization%chargeB(i)=-bunch_initialization%chargeB(i)
                allocate( bunch_init(6,bunch_initialization%n_particles(i)) )
                allocate( bunch(i)%part(bunch_initialization%n_particles(i)) )! ,STAT=AllocStatus)
                call read_from_external_file(name_file,bunch_initialization%n_particles(i),bunch_init)

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

                    bunch(i)%part(j)%cmp(12) = bunch_initialization%ChargeB(i)/bunch_initialization%n_particles(i) ! macroparticle charge
                    bunch(i)%part(j)%cmp(13) = 1.e10*bunch(i)%part(j)%cmp(12)/1.6021766 						! electrons per macroparticle

                    bunch(i)%part(j)%cmp(14)=1.0 !particle diagnostic counting
                enddo
            		if (allocated(bunch_init)) deallocate(bunch_init)
                sim_parameters%rB0(i)    = sqrt( calculate_nth_moment(i,2,1,'central') )
            		sim_parameters%lbunch(i) = sqrt( calculate_nth_moment(i,2,3,'central') )
          	enddo
   endif

!--- Apply Twiss Rotation ---!
  do i=1,bunch_initialization%n_total_bunches
   if( twiss%L_TWISS(i)  ) then

     alphaTwiss = twiss%alpha_new_factor(i)
     betaTwiss  =	twiss%beta_new_factor(i)
     eps_x=bunch_initialization%bunch_eps_x(i)
     eps_y=bunch_initialization%bunch_eps_y(i)
     s_x=bunch_initialization%bunch_s_x(i)
     s_y=bunch_initialization%bunch_s_y(i)

     ax11=sqrt( eps_x*betaTwiss/(s_x**2+s_x**2*alphaTwiss**2) )
     ay11=sqrt( eps_y*betaTwiss/(s_y**2+s_y**2*alphaTwiss**2) )
     ax12=-ax11*alphaTwiss*s_x**2/eps_x
     ay12=-ay11*alphaTwiss*s_y**2/eps_y

     do j=1,bunch_initialization%n_particles(i)
        bunch(i)%part(j)%cmp(1)=ax11*bunch(i)%part(j)%cmp(1)+ax12*bunch(i)%part(j)%cmp(4)
        bunch(i)%part(j)%cmp(2)=ay11*bunch(i)%part(j)%cmp(2)+ay12*bunch(i)%part(j)%cmp(5)
        bunch(i)%part(j)%cmp(3)=bunch(i)%part(j)%cmp(3)
        bunch(i)%part(j)%cmp(4)=bunch(i)%part(j)%cmp(4)/ax11
        bunch(i)%part(j)%cmp(5)=bunch(i)%part(j)%cmp(5)/ay11
        bunch(i)%part(j)%cmp(6)=bunch(i)%part(j)%cmp(6)

        bunch(i)%part(j)%cmp(9)=bunch(i)%part(j)%cmp(1)  !Xold
        bunch(i)%part(j)%cmp(10)=bunch(i)%part(j)%cmp(2) !Yold
        bunch(i)%part(j)%cmp(11)=bunch(i)%part(j)%cmp(3) !Zold

        bunch(i)%part(j)%cmp(14)=1.0 !particle diagnostic counting
    enddo
  endif
enddo



!------------------------------------------------------!

!--- calculate bunch charge for the non-gaussian bunches and write a diagnostic---!
  do i=1,bunch_initialization%n_total_bunches
    if(trim(bunch_initialization%shape(i))=='bigaussian') write(*,'(A,I1,A,f8.4,A)') 'Charge Bunch(',i,') :: ', bunch_initialization%ChargeB(i),'[nC] --- selected from nml file'
    if(trim(bunch_initialization%shape(i))=='cylindrical') write(*,'(A,I1,A,f8.4,A)') 'Charge Bunch(',i,') :: ', bunch_initialization%ChargeB(i),'[nC] --- computed'
    if(trim(bunch_initialization%shape(i))=='trapezoidal') write(*,'(A,I1,A,f8.4,A)') 'Charge Bunch(',i,') :: ',bunch_initialization%ChargeB(i),'[nC] --- computed'
  enddo
  write(*,'(A)')


  !--- apply shifting ---!
   Zb=0.
    do i=1,bunch_initialization%n_total_bunches
      bunch(i)%part(:)%cmp(14)=1.
      ! mu_z = abs(calculate_nth_moment(i,1,3,'nocentral'))
      if (i.ge.1) Zb = Zb+int((bunch_initialization%db(i)*plasma%lambda_p)/(mesh_par%dzm/plasma%k_p))*(mesh_par%dzm/plasma%k_p)
      write(*,'(A,f11.3)') 'distance between bunches (um)',Zb
      bunch(i)%part(:)%cmp(3)  = bunch(i)%part(:)%cmp(3) + Zb		! shift Z of bunch
      bunch(i)%part(:)%cmp(11) = bunch(i)%part(:)%cmp(3)  			! change Z_old accordingly
    enddo

   !---initial diagnostic
      do i=1,bunch_initialization%n_total_bunches
        write(*,'(A,I1,A,A)') 'Bunch(',i,') :: shape > ',trim(bunch_initialization%shape(i))
        write(*,'(A,I1,A,f11.3)') 'Bunch(',i,') :: sigma_x (um) =',sqrt( calculate_nth_moment(i,2,1,'central') )
        write(*,'(A,I1,A,f11.3)') 'Bunch(',i,') :: sigma_z (um) =',sqrt( calculate_nth_moment(i,2,3,'central') )
      enddo
      write(*,'(A,f11.3)')'Plasma wavelength   (um) =',plasma%lambda_p
      write(*,'(A,f11.3)')'Plasma wavenumber (1/um) =',plasma%k_p


   ! Computes total charge and number of electrons
   sim_parameters%Charge=sum(bunch_initialization%ChargeB(1:bunch_initialization%n_total_bunches))
   sim_parameters%NelectronsB(1:bunch_initialization%n_total_bunches) = &
    1.e10*bunch_initialization%ChargeB(1:bunch_initialization%n_total_bunches)/1.6021766 !Electrons of every macroparticle
	!1.e10*bunch_initialization%ChargeB(1:bunch_initialization%n_total_bunches)/1.6 !Electrons of every macroparticle



	 write(*,'(A)') ' bunches initialized ---'
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





subroutine dimension_first_bunch
 INTEGER :: nb,i,j
 REAL(8) zmmedB,Zb
 REAL(8), DIMENSION(:,:), ALLOCATABLE ::  bunch_init
 CHARACTER :: name_file*255

 write(*,'(A)')
 write(*,'(A)') ' >>> Initialisation :: Identify dimension for the first bunch <<<'

  if(bunch_initialization%l_bunch_internal_init) then !from input file
    sim_parameters%rB0(1)    = bunch_initialization%bunch_s_x(1)
    sim_parameters%lbunch(1) = bunch_initialization%bunch_s_z(1)
  else !read from external file
    name_file = bunch_initialization%inbunch(1)
    call read_header_external_file(name_file,bunch_initialization%chargeB(1),bunch_initialization%n_particles(1))
    bunch_initialization%chargeB(1)=-bunch_initialization%chargeB(1)
    allocate( bunch_init(6,bunch_initialization%n_particles(1)) )
    allocate( bunch(1)%part(bunch_initialization%n_particles(1)) )
    call read_from_external_file(name_file,bunch_initialization%n_particles(1),bunch_init)
    do j=1,bunch_initialization%n_particles(1)
      bunch(1)%part(j)%cmp(1)=bunch_init(1,j)
      bunch(1)%part(j)%cmp(2)=bunch_init(2,j)
      bunch(1)%part(j)%cmp(3)=bunch_init(3,j)
      bunch(1)%part(j)%cmp(4)=bunch_init(4,j)
      bunch(1)%part(j)%cmp(5)=bunch_init(5,j)
      bunch(1)%part(j)%cmp(6)=bunch_init(6,j)
      bunch(1)%part(j)%cmp(7)=1.
      bunch(1)%part(j)%cmp(8)=1.
      bunch(1)%part(j)%cmp(9)=bunch_init(1,j)  !Xold
      bunch(1)%part(j)%cmp(10)=bunch_init(2,j) !Yold
      bunch(1)%part(j)%cmp(11)=bunch_init(3,j) !Zold
      bunch(1)%part(j)%cmp(12)= bunch_initialization%ChargeB(1)/bunch_initialization%n_particles(1) ! macroparticle charge
      bunch(1)%part(j)%cmp(13)= 1.e10*bunch(1)%part(j)%cmp(12)/1.6021766 						! electrons per macroparticle
      bunch(1)%part(j)%cmp(14)=1.0 !particle diagnostic counting
    enddo
    if (allocated(bunch_init)) deallocate(bunch_init)
    ! sim_parameters%rB0(1)    = sqrt( calculate_nth_central_moment_bunch(1,2,1) )
    ! sim_parameters%lbunch(1) = sqrt( calculate_nth_central_moment_bunch(1,2,3) )
    sim_parameters%rB0(1)    = sqrt( calculate_nth_moment(1,2,1,'central') )
    sim_parameters%lbunch(1) = sqrt( calculate_nth_moment(1,2,3,'central') )
    if (allocated(bunch(1)%part)) deallocate(bunch(1)%part)
  endif

  write(*,'(A)') '>>> Initialisation :: first bunch dimension :: identified <<<'
  write(*,*)
end subroutine dimension_first_bunch



 !--- --- ---!
 end module initialize_bunch
 !--- --- ---!
