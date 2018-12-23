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

 module bunch_generation

 USE random_numbers_functions
 USE my_types
 USE use_my_types
 USE class_species
 USE class_particle
 USE moments
 implicit none

 !--- --- ---!
 contains
 !--- --- ---!



 !----------------------------------------------!
 !----------------------------------------------!
 !---            BIAGAUSSIAN                 ---!
 !----------------------------------------------!
 !----------------------------------------------!

  subroutine generate_bunch_bigaussian_equal(bunch_number)
    integer,intent(in)   :: bunch_number
    integer              :: nparticles,npZ,npR
    real(8)              :: x_cm,y_cm,z_cm,s_x,s_y,s_z,gamma_m,eps_x,eps_y,dgamma,bunch_charge,s_cut,q,m
    real(8), allocatable :: rnumber(:),rnumber3D(:,:)
    character(3)         :: optimisation
      
    !--- storying input variables into shorter-variable-name ---!
    x_cm          = 0.
    y_cm          = 0.
    z_cm          = bunchip%z_cm(bunch_number)
    s_x           = bunchip%sx_um(bunch_number)
    s_y           = bunchip%sy_um(bunch_number)
    s_z           = bunchip%sz_um(bunch_number)
    gamma_m       = bunchip%gamma(bunch_number)
    eps_x         = bunchip%epsx_um(bunch_number)
    eps_x         = bunchip%epsx_um(bunch_number)
    dgamma        = bunchip%dgamma(bunch_number)
    nparticles    = bunchip%n_particles(bunch_number)
    npZ           = bunchip%npZ(bunch_number)
    npR           = bunchip%npR(bunch_number)
    bunch_charge  = bunchip%chargeB(bunch_number)
    s_cut         = bunchip%sigma_cut(bunch_number)
    q             = bunchip%particle_charge(bunch_number)
    m             = bunchip%particle_mass(bunch_number)
    optimisation  = bunchip%optimisation(bunch_number)
    !--- --------------------------------------------------- ---!
    allocate(rnumber(nparticles),rnumber3D(3,nparticles))
      
    !--- charges and weights ---!
    !--- they need to be initialised here since they are used within the diagnostic subroutine ---!
    bunch(bunch_number)%part(:)%cmp(12)=bunch_charge/nparticles ! macroparticle charge
    bunch(bunch_number)%part(:)%cmp(13)=1.e10*bunch(bunch_number)%part(:)%cmp(12)/1.6021766
    !--- flags ---!
    bunch(bunch_number)%part(:)%cmp(7)=1.0
    bunch(bunch_number)%part(:)%cmp(8)=1.0
    bunch(bunch_number)%part(:)%cmp(14)=1.0

    call boxmuller_vector_cutDimensional(rnumber3D,nparticles,s_cut,3)
    if(       trim(optimisation) == 'no'  ) then
      bunch(bunch_number)%part(:)%cmp(1)=rnumber3D(1,:)*s_x
      bunch(bunch_number)%part(:)%cmp(2)=rnumber3D(2,:)*s_y
      bunch(bunch_number)%part(:)%cmp(3)=rnumber3D(3,:)*s_z
    else if ( trim(optimisation) == 'yes' ) then
      call initialise_rmsdimension(bunch_number,rnumber3D(1,:),nparticles,s_x,'x')
      call initialise_rmsdimension(bunch_number,rnumber3D(2,:),nparticles,s_y,'y')
      call initialise_rmsdimension(bunch_number,rnumber3D(3,:),nparticles,s_z,'z')
    endif
    !--- shifting ---!
    bunch(bunch_number)%part(:)%cmp(3)=bunch(bunch_number)%part(:)%cmp(3) + z_cm

    call boxmuller_vector(rnumber,nparticles)
    if( trim(optimisation) ==  'no' ) bunch(bunch_number)%part(:)%cmp(4)=rnumber*eps_x/s_x
    if( trim(optimisation) == 'yes' ) call initialise_emittance_x_weightedbunch(bunch_number,rnumber,nparticles,eps_x,s_x)
    call boxmuller_vector(rnumber,nparticles)
    if( trim(optimisation) ==  'no' ) bunch(bunch_number)%part(:)%cmp(5)=rnumber*eps_y/s_y
    if( trim(optimisation) == 'yes' ) call initialise_emittance_y_weightedbunch(bunch_number,rnumber,nparticles,eps_y,s_y)
    call boxmuller_vector(rnumber,nparticles)
    if( trim(optimisation) ==  'no' ) bunch(bunch_number)%part(:)%cmp(6)=-rnumber*(0.01*dgamma)*gamma_m-gamma_m
    if( trim(optimisation) == 'yes' ) call initialise_energy_spread_weightedbunch(bunch_number,rnumber,nparticles,dgamma,gamma_m)
    
    !--- X-Y-Z-old ---!
    bunch(bunch_number)%part(:)%cmp(9)=bunch(bunch_number)%part(:)%cmp(1)
    bunch(bunch_number)%part(:)%cmp(10)=bunch(bunch_number)%part(:)%cmp(2)
    bunch(bunch_number)%part(:)%cmp(11)=bunch(bunch_number)%part(:)%cmp(3)
    !--- q m ---!
    bunch(bunch_number)%part(:)%cmp(15)=q
    bunch(bunch_number)%part(:)%cmp(16)=m
    !--- ---!
    deallocate(rnumber,rnumber3D)
  end subroutine generate_bunch_bigaussian_equal




  !--- *** ---!
  subroutine generate_bunch_bigaussian_weighted(bunch_number)
    integer,intent(in)   :: bunch_number
    integer              :: nparticles,npZ,npR,iz,ir,p,ppcr,ppcz
    real(8)              :: x_cm,y_cm,z_cm,s_x,s_y,s_z,gamma_m,eps_x,eps_y,dgamma,bunch_charge,s_cut,q,m
    real(8), allocatable :: rnumber(:),rnumber3D(:,:)
    real(8)              :: dr,dz,alpha,r
    real(8)              :: r_part_dim,r_part_dimless,r_cc_dim,r_cc_dimless,theta,z_part_dim,z_cc_dim
    character(3)         :: optimisation
      
    !--- storying input variables into shorter-variable-name ---!
    x_cm          = 0.
    y_cm          = 0.
    z_cm          = bunchip%z_cm(bunch_number)
    s_x           = bunchip%sx_um(bunch_number)
    s_y           = bunchip%sy_um(bunch_number)
    s_z           = bunchip%sz_um(bunch_number)
    gamma_m       = bunchip%gamma(bunch_number)
    eps_x         = bunchip%epsx_um(bunch_number)
    eps_y         = bunchip%epsy_um(bunch_number)
    dgamma        = bunchip%dgamma(bunch_number)
    nparticles    = bunchip%n_particles(bunch_number)
    npZ           = bunchip%npZ(bunch_number)
    npR           = bunchip%npR(bunch_number)
    bunch_charge  = bunchip%chargeB(bunch_number)
    s_cut         = bunchip%sigma_cut(bunch_number)
    q             = bunchip%particle_charge(bunch_number)
    m             = bunchip%particle_mass(bunch_number)
    optimisation  = bunchip%optimisation(bunch_number)
    !--- --------------------------------------------------- ---!
    allocate(rnumber(nparticles))

    dr=mesh_par%dr/plasma%k_p
    dz=mesh_par%dz/plasma%k_p

    alpha = (bunch_charge*1e-9) /(2.d0*pi)/sqrt(2.d0*pi)/(s_x*1e-6)/(s_y*1e-6)/(s_z*1e-6)
    alpha = alpha / 1.6e-19 / (plasma%n0*1e6) ! * (plasma%k_p/1e6)**3

    p=1
    do iz=-int(s_cut*s_z/dz),int(s_cut*s_z/dz)
      do ir=0,int(s_cut*s_x/dr)
        if( ((ir+.5)*dr/s_cut/s_x)**2+((iz+.5)*dz/s_cut/s_z)**2<1.) then
          do ppcr=1,npR
            do ppcz=1,npZ
              r_part_dim          =ir*dr+dr/(npR+1.)*ppcr
              r_part_dimless      =ir*mesh_par%dr+mesh_par%dr/(npR+1.)*ppcr
              r_cc_dim            =ir*dr+dr/2.
              r_cc_dimless        =ir*mesh_par%dr+mesh_par%dr/2.
              !  theta               =2.d0*pi/(npR*npZ)*(ppcz+(ppcr-1)*npZ) !random_number_range(0.d0,2.d0*pi)
              theta               =0.D0
              z_part_dim          =iz*dz+dz/(npZ+1.)*ppcz
              z_cc_dim            =iz*dz+dz/2.
              bunch(bunch_number)%part(p)%cmp(1)=r_part_dim*cos(theta)
              bunch(bunch_number)%part(p)%cmp(2)=r_part_dim*sin(theta)
              bunch(bunch_number)%part(p)%cmp(3)=z_part_dim + z_cm
              bunch(bunch_number)%part(p)%cmp(12)=alpha *r_part_dimless*mesh_par%dr*mesh_par%dz /npZ/npR
              bunch(bunch_number)%part(p)%cmp(12)=bunch(bunch_number)%part(p)%cmp(12) * 1.5101/sqrt(plasma%n0/1e16)
              bunch(bunch_number)%part(p)%cmp(12)=bunch(bunch_number)%part(p)%cmp(12) * exp(-r_part_dim**2/2./s_x**2)
              bunch(bunch_number)%part(p)%cmp(12)=bunch(bunch_number)%part(p)%cmp(12) * exp(-z_part_dim**2/2./s_z**2)
              bunch(bunch_number)%part(p)%cmp(13)=1.e10*bunch(bunch_number)%part(p)%cmp(12)/1.6021766
              p=p+1
            enddo
          enddo
        endif
      enddo
    enddo

    call boxmuller_vector(rnumber,nparticles)
    if( trim(optimisation) ==  'no' ) bunch(bunch_number)%part(:)%cmp(4)=rnumber*eps_x/s_x
    if( trim(optimisation) == 'yes' ) call initialise_emittance_x_weightedbunch(bunch_number,rnumber,nparticles,eps_x,s_x)
    bunch(bunch_number)%part(:)%cmp(5)=0.D0
    call boxmuller_vector(rnumber,nparticles)
    if( trim(optimisation) ==  'no' ) bunch(bunch_number)%part(:)%cmp(6)=-rnumber*(0.01*dgamma)*gamma_m-gamma_m
    if( trim(optimisation) == 'yes' ) call initialise_energy_spread_weightedbunch(bunch_number,rnumber,nparticles,dgamma,gamma_m)

    !--- flags ---!
    bunch(bunch_number)%part(:)%cmp(7)=1.0
    bunch(bunch_number)%part(:)%cmp(8)=1.0
    bunch(bunch_number)%part(:)%cmp(14)=1.0
    !--- X-Y-Z-old ---!
    bunch(bunch_number)%part(:)%cmp(9)=bunch(bunch_number)%part(:)%cmp(1)
    bunch(bunch_number)%part(:)%cmp(10)=bunch(bunch_number)%part(:)%cmp(2)
    bunch(bunch_number)%part(:)%cmp(11)=bunch(bunch_number)%part(:)%cmp(3)
    !--- q m ---!
    bunch(bunch_number)%part(:)%cmp(15)=q
    bunch(bunch_number)%part(:)%cmp(16)=m
    !--- ---!
    deallocate(rnumber)
  end subroutine generate_bunch_bigaussian_weighted


!----------------------------------------------!
!----------------------------------------------!
!---    TRAPEZOIDAL::Z  + UNIFROM::R        ---!
!----------------------------------------------!
!----------------------------------------------!
  subroutine generate_bunch_trapezoidalZ_uniformR_equal(bunch_number)
    integer,intent(in)   :: bunch_number
    integer              :: nparticles,npZ,npR,i
    real(8)              :: x_cm,y_cm,z_cm,s_x,s_y,s_z,gamma_m,eps_x,eps_y,dgamma,bunch_charge,s_cut,q,m,Charge_right,Charge_left
    real(8)              :: z,y,x,a
    real(8), allocatable :: rnumber(:),particle_radius(:),particle_theta(:)
    character(3)         :: optimisation
      
    !--- storying input variables into shorter-variable-name ---!
    x_cm          = 0.
    y_cm          = 0.
    z_cm          = bunchip%z_cm(bunch_number)
    s_x           = bunchip%sx_um(bunch_number)
    s_y           = bunchip%sy_um(bunch_number)
    s_z           = bunchip%sz_um(bunch_number)
    gamma_m       = bunchip%gamma(bunch_number)
    eps_x         = bunchip%epsx_um(bunch_number)
    eps_y         = bunchip%epsy_um(bunch_number)
    dgamma        = bunchip%dgamma(bunch_number)
    nparticles    = bunchip%n_particles(bunch_number)
    npZ           = bunchip%npZ(bunch_number)
    npR           = bunchip%npR(bunch_number)
    bunch_charge  = bunchip%chargeB(bunch_number)
    Charge_left   = bunchip%Charge_left(bunch_number)
    Charge_right  = bunchip%Charge_right(bunch_number)
    s_cut         = bunchip%sigma_cut(bunch_number)
    q             = bunchip%particle_charge(bunch_number)
    m             = bunchip%particle_mass(bunch_number)
    optimisation  = bunchip%optimisation(bunch_number)
    !--- --------------------------------------------------- ---!
    allocate(rnumber(nparticles),particle_radius(nparticles),particle_theta(nparticles))

    !--- bunch_charge ---!
    bunch_charge = electron_charge*pi * (s_x*1d-6) * (s_y*1d-6) &
                                      * (Charge_left+Charge_right)*s_z*1d-6/2.d0 * (plasma%n0*1d6)
    bunch_charge = bunch_charge*1e9 !converting to [nC]
    bunchip%chargeB(bunch_number) = bunch_charge !forcing the input value

    !--- charges and weights ---!
    bunch(bunch_number)%part(:)%cmp(12)=bunch_charge/nparticles ! macroparticle charge
    bunch(bunch_number)%part(:)%cmp(13)=1.e10*bunch(bunch_number)%part(:)%cmp(12)/1.6021766

    do i=1,nparticles
        call random_number(z)
        call random_number(a)
        Do while(a>Charge_left+(Charge_right-Charge_left)*z)
            call random_number(z)
            call random_number(a)
        enddo
        x=random_number_range(-1.d0,1.d0)
        y=random_number_range(-1.d0,1.d0)
        Do while(sqrt(x**2+y**2)>1.d0)
            x=random_number_range(-1.d0,1.d0)
            y=random_number_range(-1.d0,1.d0)
        enddo
        bunch(bunch_number)%part(i)%cmp(1)=x*s_x!+x_cm
        bunch(bunch_number)%part(i)%cmp(2)=y*s_y!+y_cm
        bunch(bunch_number)%part(i)%cmp(3)=z*s_z+z_cm
    enddo

    call boxmuller_vector(rnumber,nparticles)
    if( trim(optimisation) ==  'no' ) bunch(bunch_number)%part(:)%cmp(4)=rnumber*eps_x/s_x
    if( trim(optimisation) == 'yes' ) call initialise_emittance_x_weightedbunch(bunch_number,rnumber,nparticles,eps_x,s_x)
    
    call boxmuller_vector(rnumber,nparticles)
    if( trim(optimisation) ==  'no' ) bunch(bunch_number)%part(:)%cmp(5)=rnumber*eps_y/s_y
    if( trim(optimisation) == 'yes' ) call initialise_emittance_y_weightedbunch(bunch_number,rnumber,nparticles,eps_y,s_y)

    call boxmuller_vector(rnumber,nparticles)
    if( trim(optimisation) == 'no'  ) bunch(bunch_number)%part(:)%cmp(6)=-rnumber * 0.01*dgamma*gamma_m-gamma_m
    if( trim(optimisation) == 'yes' ) call initialise_energy_spread_weightedbunch(bunch_number,rnumber,nparticles,dgamma,gamma_m)    

    !--- flags ---!
    bunch(bunch_number)%part(:)%cmp(7)=1.0
    bunch(bunch_number)%part(:)%cmp(8)=1.0
    bunch(bunch_number)%part(:)%cmp(14)=1.0
    !--- X-Y-Z-old ---!
    bunch(bunch_number)%part(:)%cmp(9)=bunch(bunch_number)%part(:)%cmp(1)
    bunch(bunch_number)%part(:)%cmp(10)=bunch(bunch_number)%part(:)%cmp(2)
    bunch(bunch_number)%part(:)%cmp(11)=bunch(bunch_number)%part(:)%cmp(3)
    !--- q m ---!
    bunch(bunch_number)%part(:)%cmp(15)=q
    bunch(bunch_number)%part(:)%cmp(16)=m
    !--- ---!
    deallocate(rnumber,particle_radius,particle_theta)
  end subroutine generate_bunch_trapezoidalZ_uniformR_equal



  subroutine generate_bunch_trapezoidalZ_uniformR_weighted(bunch_number)
    integer,intent(in)   :: bunch_number
    integer              :: nparticles,npZ,npR,i,iz,ir,p,ppcr,ppcz
    real(8)              :: x_cm,y_cm,z_cm,s_x,s_y,s_z,gamma_m,eps_x,eps_y,dgamma,bunch_charge,s_cut,q,m,Charge_right,Charge_left
    real(8)              :: z,y,x,a,dr,dz
    real(8), allocatable :: rnumber(:),particle_radius(:),particle_theta(:)
    real(8)              :: r_part_dim,r_part_dimless,r_cc_dim,r_cc_dimless,z_part_dim,z_cc_dim,theta
    character(3)         :: optimisation
      
    !--- storying input variables into shorter-variable-name ---!
    x_cm          = 0.
    y_cm          = 0.
    z_cm          = bunchip%z_cm(bunch_number)
    s_x           = bunchip%sx_um(bunch_number)
    s_y           = bunchip%sy_um(bunch_number)
    s_z           = bunchip%sz_um(bunch_number)
    gamma_m       = bunchip%gamma(bunch_number)
    eps_x         = bunchip%epsx_um(bunch_number)
    eps_y         = bunchip%epsy_um(bunch_number)
    dgamma        = bunchip%dgamma(bunch_number)
    nparticles    = bunchip%n_particles(bunch_number)
    npZ           = bunchip%npZ(bunch_number)
    npR           = bunchip%npR(bunch_number)
    bunch_charge  = bunchip%chargeB(bunch_number)
    Charge_left   = bunchip%Charge_left(bunch_number)
    Charge_right  = bunchip%Charge_right(bunch_number)
    s_cut         = bunchip%sigma_cut(bunch_number)
    q             = bunchip%particle_charge(bunch_number)
    m             = bunchip%particle_mass(bunch_number)
    optimisation  = bunchip%optimisation(bunch_number)
    !--- --------------------------------------------------- ---!
    allocate(rnumber(nparticles),particle_radius(nparticles),particle_theta(nparticles))

    dr=mesh_par%dr/plasma%k_p
    dz=mesh_par%dz/plasma%k_p

    p=1
    do iz=-int(s_z/dz),0
      do ir=0,int(s_x/dr)
        do ppcr=1,npR
          do ppcz=1,npZ
            r_part_dim          =ir*dr+dr/(npR+1.)*ppcr !particle radius
            r_part_dimless      =ir*mesh_par%dr+mesh_par%dr/(npR+1.)*ppcr !particle radius dimensionless
            r_cc_dim            =ir*dr+dr/2. !cell centre :: now not used
            r_cc_dimless        =ir*mesh_par%dr+mesh_par%dr/2. !cell centre dimensionless:: now not used
            !  theta               =2.d0*pi/(npR*npZ)*(ppcz+(ppcr-1)*npZ) !random_number_range(0.d0,2.d0*pi)
            theta               =0.D0
            z_part_dim          =iz*dz+dz/(npZ+1.)*ppcz
            z_cc_dim            =iz*dz+dz/2.
            bunch(bunch_number)%part(p)%cmp(1)=r_part_dim*cos(theta)
            bunch(bunch_number)%part(p)%cmp(2)=r_part_dim*sin(theta)
            bunch(bunch_number)%part(p)%cmp(3)=z_part_dim+s_z
            bunch(bunch_number)%part(p)%cmp(12)=Charge_right + (Charge_left-Charge_right)/s_z*abs(z_part_dim)
            bunch(bunch_number)%part(p)%cmp(12)=bunch(bunch_number)%part(p)%cmp(12) * 1.5101/sqrt(plasma%n0/1e16) * r_part_dimless *mesh_par%dr*mesh_par%dz /npZ/npR
            bunch(bunch_number)%part(p)%cmp(13)=1.e10*bunch(bunch_number)%part(p)%cmp(12)/1.6021766
            p=p+1
          enddo
        enddo
      enddo
    enddo

    call boxmuller_vector(rnumber,nparticles)
    if( trim(optimisation) == 'no'  ) bunch(bunch_number)%part(:)%cmp(4)=rnumber*eps_x/s_x
    if( trim(optimisation) == 'yes' ) call initialise_emittance_x_weightedbunch(bunch_number,rnumber,nparticles,eps_x,s_x)
  
    bunch(bunch_number)%part(:)%cmp(5)=0.D0

    call boxmuller_vector(rnumber,nparticles)
    if( trim(optimisation) == 'no'  ) bunch(bunch_number)%part(:)%cmp(6)=-rnumber * 0.01*dgamma*gamma_m-gamma_m
    if( trim(optimisation) == 'yes' ) call initialise_energy_spread_weightedbunch(bunch_number,rnumber,nparticles,dgamma,gamma_m)

    !--- flags ---!
    bunch(bunch_number)%part(:)%cmp(7)=1.0
    bunch(bunch_number)%part(:)%cmp(8)=1.0
    bunch(bunch_number)%part(:)%cmp(14)=1.0
    !--- X-Y-Z-old ---!
    bunch(bunch_number)%part(:)%cmp(9)=bunch(bunch_number)%part(:)%cmp(1)
    bunch(bunch_number)%part(:)%cmp(10)=bunch(bunch_number)%part(:)%cmp(2)
    bunch(bunch_number)%part(:)%cmp(11)=bunch(bunch_number)%part(:)%cmp(3)
    !--- q m ---!
    bunch(bunch_number)%part(:)%cmp(15)=q
    bunch(bunch_number)%part(:)%cmp(16)=m
    !--- ---!
    deallocate(rnumber,particle_radius,particle_theta)
  end subroutine generate_bunch_trapezoidalZ_uniformR_weighted


!----------------------------------------------!
!----------------------------------------------!
!---    TRAPEZOIDAL::Z  + GAUSSIAN::R       ---!
!----------------------------------------------!
!----------------------------------------------!
  subroutine generate_bunch_trapezoidalZ_gaussianR_equal(bunch_number)
    integer,intent(in)   :: bunch_number
    integer              :: nparticles,npZ,npR,i
    real(8)              :: x_cm,y_cm,z_cm,s_x,s_y,s_z,gamma_m,eps_x,eps_y,dgamma,bunch_charge,s_cut,q,m,Charge_right,Charge_left
    real(8)              :: z,y,x,a
    real(8), allocatable :: rnumber(:),rnumber3D(:,:),particle_radius(:),particle_theta(:)
    character(3)         :: optimisation
      
    !--- storying input variables into shorter-variable-name ---!
    x_cm          = 0.
    y_cm          = 0.
    z_cm          = bunchip%z_cm(bunch_number)
    s_x           = bunchip%sx_um(bunch_number)
    s_y           = bunchip%sy_um(bunch_number)
    s_z           = bunchip%sz_um(bunch_number)
    gamma_m       = bunchip%gamma(bunch_number)
    eps_x         = bunchip%epsx_um(bunch_number)
    eps_y         = bunchip%epsy_um(bunch_number)
    dgamma        = bunchip%dgamma(bunch_number)
    nparticles    = bunchip%n_particles(bunch_number)
    npZ           = bunchip%npZ(bunch_number)
    npR           = bunchip%npR(bunch_number)
    bunch_charge  = bunchip%chargeB(bunch_number)
    Charge_left   = bunchip%Charge_left(bunch_number)
    Charge_right  = bunchip%Charge_right(bunch_number)
    s_cut         = bunchip%sigma_cut(bunch_number)
    q             = bunchip%particle_charge(bunch_number)
    m             = bunchip%particle_mass(bunch_number)
    optimisation  = bunchip%optimisation(bunch_number)
    !--- --------------------------------------------------- ---!
    allocate(rnumber(nparticles),rnumber3D(3,nparticles),particle_radius(nparticles),particle_theta(nparticles))

    !--- bunch_charge ---!
    bunch_charge = electron_charge* 2.d0*pi* (s_x*1d-6) * (s_y*1d-6) &
                          * (Charge_left+Charge_right)*s_z*1d-6/2.d0 * (plasma%n0*1d6)
    bunch_charge = bunch_charge*1e9 !converting to [nC]
    bunchip%chargeB(bunch_number) = bunch_charge  !forcing the input value

    !--- charges and weights ---!
    bunch(bunch_number)%part(:)%cmp(12)=bunch_charge/nparticles ! macroparticle charge
    bunch(bunch_number)%part(:)%cmp(13)=1.e10*bunch(bunch_number)%part(:)%cmp(12)/1.6021766
    !--- flags ---!
    bunch(bunch_number)%part(:)%cmp(7)=1.0
    bunch(bunch_number)%part(:)%cmp(8)=1.0
    bunch(bunch_number)%part(:)%cmp(14)=1.0

    do i=1,nparticles
        call random_number(z)
        call random_number(a)
        if(Charge_left >= Charge_right) then !triangular shape left-right
            Do while(a*Charge_left>=Charge_left+(Charge_right-Charge_left)*z)
                call random_number(z)
                call random_number(a)
            enddo
        elseif(Charge_left < Charge_right) then !triangular shape right-left
            Do while(a*Charge_right>=Charge_left+(Charge_right-Charge_left)*z)
                call random_number(z)
                call random_number(a)
            enddo
        endif
        ! bunch(bunch_number)%part(i)%cmp(3)=z*(s_z+mesh_par%dz/plasma%k_p) + .5*mesh_par%dz/plasma%k_p !+z_cm
        bunch(bunch_number)%part(i)%cmp(3)=z*s_z + z_cm
    enddo

    call boxmuller_vector_cutDimensional(rnumber3D,nparticles,s_cut,2)
    if( trim(optimisation) ==  'no' ) bunch(bunch_number)%part(:)%cmp(1)=rnumber3D(1,:)*s_x
    if( trim(optimisation) ==  'no' ) bunch(bunch_number)%part(:)%cmp(2)=rnumber3D(2,:)*s_y
    if( trim(optimisation) ==  'yes') call initialise_rmsdimension(bunch_number,rnumber3D(1,:),nparticles,s_x,'x')
    if( trim(optimisation) ==  'yes') call initialise_rmsdimension(bunch_number,rnumber3D(2,:),nparticles,s_y,'y')

    call boxmuller_vector(rnumber,nparticles)
    if( trim(optimisation) ==  'no' ) bunch(bunch_number)%part(:)%cmp(4)=rnumber*eps_x/s_x
    if( trim(optimisation) ==  'yes') call initialise_emittance_x_weightedbunch(bunch_number,rnumber,nparticles,eps_x,s_x)

    call boxmuller_vector(rnumber,nparticles)
    if( trim(optimisation) ==  'no' ) bunch(bunch_number)%part(:)%cmp(5)=rnumber*eps_y/s_y
    if( trim(optimisation) ==  'yes') call initialise_emittance_y_weightedbunch(bunch_number,rnumber,nparticles,eps_y,s_y)

    call boxmuller_vector(rnumber,nparticles)
    if( trim(optimisation) ==  'no' ) bunch(bunch_number)%part(:)%cmp(6)=-rnumber*(0.01*dgamma)*gamma_m-gamma_m
    if( trim(optimisation) ==  'yes') call initialise_energy_spread_weightedbunch(bunch_number,rnumber,nparticles,dgamma,gamma_m)

    !--- X-Y-Z-old ---!
    bunch(bunch_number)%part(:)%cmp(9)=bunch(bunch_number)%part(:)%cmp(1)
    bunch(bunch_number)%part(:)%cmp(10)=bunch(bunch_number)%part(:)%cmp(2)
    bunch(bunch_number)%part(:)%cmp(11)=bunch(bunch_number)%part(:)%cmp(3)
    !--- q m ---!
    bunch(bunch_number)%part(:)%cmp(15)=q
    bunch(bunch_number)%part(:)%cmp(16)=m
    !--- ---!
    deallocate(rnumber,rnumber3D,particle_radius,particle_theta)
  end subroutine generate_bunch_trapezoidalZ_gaussianR_equal



!--- *** ---!
  subroutine generate_bunch_trapezoidalZ_gaussianR_weighted(bunch_number)
    integer,intent(in)   :: bunch_number
    integer              :: nparticles,npZ,npR,i,iz,ir,p,ppcr,ppcz,Ns_z
    real(8)              :: x_cm,y_cm,z_cm,s_x,s_y,s_z,gamma_m,eps_x,eps_y,dgamma,bunch_charge,s_cut,q,m,Charge_right,Charge_left
    real(8)              :: z,y,x,a,dr,dz
    real(8), allocatable :: rnumber(:),particle_radius(:),particle_theta(:)
    real(8)              :: r_part_dim,r_part_dimless,r_cc_dim,r_cc_dimless,z_part_dim,z_cc_dim,theta
    character(3)         :: optimisation
      
    !--- storying input variables into shorter-variable-name ---!
    x_cm          = 0.
    y_cm          = 0.
    z_cm          = bunchip%z_cm(bunch_number)
    s_x           = bunchip%sx_um(bunch_number)
    s_y           = bunchip%sy_um(bunch_number)
    s_z           = bunchip%sz_um(bunch_number)
    gamma_m       = bunchip%gamma(bunch_number)
    eps_x         = bunchip%epsx_um(bunch_number)
    eps_y         = bunchip%epsy_um(bunch_number)
    dgamma        = bunchip%dgamma(bunch_number)
    nparticles    = bunchip%n_particles(bunch_number)
    npZ           = bunchip%npZ(bunch_number)
    npR           = bunchip%npR(bunch_number)
    bunch_charge  = bunchip%chargeB(bunch_number)
    Charge_left   = bunchip%Charge_left(bunch_number)
    Charge_right  = bunchip%Charge_right(bunch_number)
    s_cut         = bunchip%sigma_cut(bunch_number)
    q             = bunchip%particle_charge(bunch_number)
    m             = bunchip%particle_mass(bunch_number)
    optimisation  = bunchip%optimisation(bunch_number)
    !--- --------------------------------------------------- ---!
    allocate(rnumber(nparticles),particle_radius(nparticles),particle_theta(nparticles))

    !--- bunch_charge ---!
    bunch_charge= electron_charge* 2.d0*pi* (s_x*1d-6) * (s_y*1d-6) &
                        * (Charge_left+Charge_right)*s_z*1d-6/2.d0 * (plasma%n0*1d6)
    bunch_charge = bunch_charge*1e9 !converting to [nC]
    bunchip%chargeB(bunch_number) = bunch_charge !forcing the input value

    !---
    dr=mesh_par%dr/plasma%k_p
    dz=mesh_par%dz/plasma%k_p

    Ns_z = int(s_z/dz)-1
    p=1
      do iz=-Ns_z,0
        do ir=0,int(s_cut*s_x/dr)-1
          do ppcr=1,npR
            do ppcz=1,npZ
              r_part_dim          =ir*dr+dr/(npR+1.)*ppcr
              r_part_dimless      =ir*mesh_par%dr+mesh_par%dr/(npR+1.)*ppcr
              r_cc_dim            =ir*dr+dr/2. !cc:cell-centre radius cell centre
              r_cc_dimless        =ir*mesh_par%dr+mesh_par%dr/2. !cc:cell-centre radius cell centre
              ! theta               =2.d0*pi/(npR*npZ)*(ppcz+(ppcr-1)*npZ) !random_number_range(0.d0,2.d0*pi)
              theta               =0.D0
              z_part_dim          =iz*dz+dz/(npZ+1.)*ppcz
              z_cc_dim            =(iz+.5)*dz+dz/2.
              bunch(bunch_number)%part(p)%cmp(1)=r_part_dim*cos(theta)
              bunch(bunch_number)%part(p)%cmp(2)=r_part_dim*sin(theta)
              bunch(bunch_number)%part(p)%cmp(3)=z_part_dim +Ns_z*dz
              bunch(bunch_number)%part(p)%cmp(12)=Charge_right + (Charge_left-Charge_right)/(int(s_z/dz)-.5)*abs(iz)
              bunch(bunch_number)%part(p)%cmp(12)=bunch(bunch_number)%part(p)%cmp(12) * 1.5101/sqrt(plasma%n0/1e16) * r_part_dimless *mesh_par%dr*mesh_par%dz /npZ/npR
              bunch(bunch_number)%part(p)%cmp(12)=bunch(bunch_number)%part(p)%cmp(12) * exp(-r_part_dim**2/2./s_x**2)
              bunch(bunch_number)%part(p)%cmp(13)=1.e10*bunch(bunch_number)%part(p)%cmp(12)/1.6021766
              p=p+1
            enddo
          enddo
      enddo
    enddo

    call boxmuller_vector(rnumber,nparticles)
    if( trim(optimisation) ==  'no' ) bunch(bunch_number)%part(:)%cmp(4)=rnumber*eps_x/s_x
    if( trim(optimisation) == 'yes' ) call initialise_emittance_x_weightedbunch(bunch_number,rnumber,nparticles,eps_x,s_x)

    bunch(bunch_number)%part(:)%cmp(5)=0.D0
    
    call boxmuller_vector(rnumber,nparticles)
    if( trim(optimisation) ==  'no' ) bunch(bunch_number)%part(:)%cmp(6)=-rnumber*(0.01*dgamma)*gamma_m-gamma_m
    if( trim(optimisation) == 'yes' ) call initialise_energy_spread_weightedbunch(bunch_number,rnumber,nparticles,dgamma,gamma_m)

    !--- flags ---!
    bunch(bunch_number)%part(:)%cmp(7)=1.0
    bunch(bunch_number)%part(:)%cmp(8)=1.0
    bunch(bunch_number)%part(:)%cmp(14)=1.0
    !--- X-Y-Z-old ---!
    bunch(bunch_number)%part(:)%cmp(9)=bunch(bunch_number)%part(:)%cmp(1)
    bunch(bunch_number)%part(:)%cmp(10)=bunch(bunch_number)%part(:)%cmp(2)
    bunch(bunch_number)%part(:)%cmp(11)=bunch(bunch_number)%part(:)%cmp(3)
    !--- q m ---!
    bunch(bunch_number)%part(:)%cmp(15)=q
    bunch(bunch_number)%part(:)%cmp(16)=m
    !--- ---!
    deallocate(rnumber,particle_radius,particle_theta)
  end subroutine generate_bunch_trapezoidalZ_gaussianR_weighted

!----------------------------------------------!
!----------------------------------------------!
!---    READ FROM EXTERNAL FILE             ---!
!----------------------------------------------!
!----------------------------------------------!
 subroutine read_from_external_file(name_file,number_of_particles,generated_bunch)
  character,intent(in) :: name_file*255
  integer,intent(in)   :: number_of_particles
  real(8),intent(out)   :: generated_bunch(:,:)
  real(8) :: dummy_charge
  integer :: i,dummy_nparticles
  character(199) :: string

  open(9,file=trim(name_file),status='old')

  write(*,'(A,A)') 'loading particle from external file: ',trim(name_file)

  !--- read: header ---!
  read(9,*) string
  read(9,*) dummy_charge
  read(9,*) string
  read(9,*) dummy_nparticles
  read(9,*) string

  DO i=1,number_of_particles
    read(9,*) generated_bunch(1,i), generated_bunch(2,i), generated_bunch(3,i), generated_bunch(4,i), generated_bunch(5,i), generated_bunch(6,i)
  ENDDO

  close(9)

  write(*,'(A,A)') trim(name_file),'   loaded'
  end subroutine read_from_external_file


 !--- *** ---!
  subroutine read_header_external_file(name_file,charge,number_of_particles)
    character,intent(in) :: name_file*255
    real(8),intent(out)   :: charge
    integer,intent(out)   :: number_of_particles
    character(199) :: string

    open(9,file=trim(name_file),status='unknown')
    !--- read: header ---!
    read(9,*) string
    read(9,*) charge
    read(9,*) string
    read(9,*) number_of_particles
    write(*,'(A,f6.3,A)') 'reading header :: Bunch Charge::  ',charge,' (nC)'
    write(*,'(A,1I8)') 'Reading header :: Bunch number of particles:: ',number_of_particles
    close(9)
  end subroutine read_header_external_file


  real(8) function calculate_bunch_number_of_particles(nb,shape,PWeights,nparticles)
    integer, intent(in) :: nb,nparticles
    character(len=*), intent(in) :: shape,PWeights
    integer :: iz,ir
    real(8) :: dz,dr, s_x, s_y, s_z, s_cut

    calculate_bunch_number_of_particles=0
    dr=mesh_par%dr/plasma%k_p
    dz=mesh_par%dz/plasma%k_p
    s_x=bunchip%sx_um(nb)
    s_z=bunchip%sz_um(nb)
    s_cut = bunchip%sigma_cut(nb)

    if( trim(PWeights)=='equal') calculate_bunch_number_of_particles=nparticles

    if( trim(shape)=='bigaussian' .and. trim(PWeights)=='weighted') then
      do iz=-int(s_cut*s_z/dz),int(s_cut*s_z/dz)
        do ir=0,int(s_cut*s_x/dr)+1
          if( ((ir+.5)*dr/s_cut/s_x)**2+((iz+.5)*dz/s_cut/s_z)**2<1.) then
            calculate_bunch_number_of_particles=calculate_bunch_number_of_particles+ &
                bunchip%npZ(nb)*bunchip%npR(nb)
          endif
        enddo
      enddo
    endif


    if( trim(shape)=='cylindrical' .and. trim(PWeights)=='weighted') then
      do iz=-int(s_z/dz),0
        do ir=0,int(s_x/dr)
            calculate_bunch_number_of_particles=calculate_bunch_number_of_particles+ &
                bunchip%npZ(nb)*bunchip%npR(nb)
        enddo
      enddo
    endif


    if ( trim(shape)=='trapezoidal' .and. trim(PWeights)=='weighted') then
      do iz=-int(s_z/dz)+1,0
        do ir=0,int(s_cut*s_x/dr)-1
            calculate_bunch_number_of_particles=calculate_bunch_number_of_particles+ &
                bunchip%npZ(nb)*bunchip%npR(nb)
        enddo
      enddo
    endif
  end function calculate_bunch_number_of_particles


 subroutine initialise_emittance_x_weightedbunch(bunch_number,rnumber,nparticles,eps_x,s_x)
   integer, intent(in)    :: bunch_number,nparticles
   real(8), intent(in)    :: rnumber(nparticles),s_x,eps_x
   real(8) :: err,eps_local

   eps_local=eps_x
   err=1.d0
   DO WHILE( ABS(err) > 1e-6)
     bunch(bunch_number)%part(:)%cmp(4)=rnumber*eps_local/s_x
     err = calculate_emittance_x(bunch_number)-eps_x
     eps_local=eps_local-err
   ENDDO
 end subroutine initialise_emittance_x_weightedbunch

 subroutine initialise_emittance_y_weightedbunch(bunch_number,rnumber,nparticles,eps_y,s_y)
   integer, intent(in)    :: bunch_number,nparticles
   real(8), intent(in)    :: rnumber(nparticles),s_y,eps_y
   real(8) :: err,eps_local

   eps_local=eps_y
   err=1.d0
   DO WHILE( ABS(err) > 1e-6)
     bunch(bunch_number)%part(:)%cmp(5)=rnumber*eps_local/s_y
     err = calculate_emittance_y(bunch_number)-eps_y
     eps_local=eps_local-err
   ENDDO
 end subroutine initialise_emittance_y_weightedbunch

 subroutine initialise_energy_spread_weightedbunch(bunch_number,rnumber,nparticles,dgamma,gamma_m)
   integer, intent(in)    :: bunch_number,nparticles
   real(8), intent(in)    :: rnumber(nparticles),dgamma,gamma_m
   real(8) :: err,dgamma_local

   dgamma_local=dgamma
   err=1.d0
   DO WHILE( ABS(err) > 1d-8)
     bunch(bunch_number)%part(:)%cmp(6)=-(rnumber *0.01d0*dgamma_local*gamma_m + gamma_m)
     err = calculate_energy_spread(bunch_number)-(dgamma*0.01d0)
     dgamma_local=dgamma_local-err
   ENDDO
 end subroutine initialise_energy_spread_weightedbunch


 subroutine initialise_rmsdimension(bunch_number,rnumber,nparticles,sigma,direction)
   integer, intent(in)    :: bunch_number,nparticles
   real(8), intent(in)    :: rnumber(nparticles),sigma
   character(1), intent (in) :: direction
   real(8) :: err,sigma_local
   integer :: dir

   if(direction=='x') dir=1
   if(direction=='y') dir=2
   if(direction=='z') dir=3

   sigma_local=sigma
   err=1.d0
   DO WHILE( ABS(err) > 1e-6)
     bunch(bunch_number)%part(:)%cmp(dir)=rnumber*sigma_local
     err = sqrt(calculate_nth_moment(bunch_number,2,dir,'central'))-sigma
     sigma_local=sigma_local-err
   ENDDO
 end subroutine initialise_rmsdimension

 !--- --- ---!
 end module bunch_generation
 !--- --- ---!
