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

MODULE ion_background

 USE my_types
 USE use_my_types
 USE pstruct_data
 USE architect_class_structure
 USE utilities
 USE Compute_plasma_current

IMPLICIT NONE

CONTAINS

  subroutine initialise_ion_background
    integer tot_ions

    !--- *** ---!
    if(.not.ionisation%L_ionisation) then
      mesh(:,:)%ni_bck=mesh(:,:)%ne_bck
      return
    endif

    !--- *** ---!
    call from_ppc_to_npr_npz(ionisation%particle_per_cell,ionisation%np_per_r,ionisation%np_per_z)
    ionisation%tot_ions=(mesh_par%Nzm-2)*(mesh_par%Nxm-2)*ionisation%particle_per_cell
    write(*,'(A)')
    write(*,'(A)') 'Ionization module active :: static background ions'
    write(*,'(A,1I2)') 'np_per_z>',ionisation%np_per_z
    write(*,'(A,1I2)') 'np_per_r>',ionisation%np_per_r
    write(*,'(A,1I9)') 'total number of ions:',ionisation%tot_ions
    write(*,'(A)')

    allocate(static_ion(ionisation%tot_ions))
    static_ion(:)%cmp(7) = 0.

    call dispose_ions
    call from_particleZstart_to_nplasmai

  end subroutine initialise_ion_background

  !--------------------------------------!
  !--- dispose ions at t=0 --------------!
  !--------------------------------------!
  subroutine dispose_ions
    integer :: i,j,k,idx_z,idx_r,ion_counter=0
    real(8) :: hz,hr,z_coord,r_coord

    hz=mesh_par%dz/real(ionisation%np_per_z+1)
    hr=mesh_par%dr/real(ionisation%np_per_r+1)

    do i=1,mesh_par%Nzm-1
      do j=1,mesh_par%Nxm-1
        if( mesh_par%z_min_moving_um+(i-1)*mesh_par%dz/plasma%k_p<bck_plasma%z_coordinate_um(1) .and. &
            j<mesh_par%Nr_plasma ) then
            do idx_z=1,ionisation%np_per_z
              do idx_r=1,ionisation%np_per_r
                ion_counter=ion_counter+1
                z_coord=mesh_par%z_min_moving+(i-1)*mesh_par%dz+hz*idx_z
                r_coord=(j-1)*mesh_par%dr+hr*idx_r
                call initialise_ion(z_coord,r_coord,ionisation%mass_number,ionisation%atomic_number,1.d0,ion_counter)
              enddo
            enddo
        endif
      enddo
    enddo
  end subroutine dispose_ions


  !--------------------------------------!
  !---Initialise static background ion---!
  !--------------------------------------!
  subroutine initialise_ion(z_coord,r_coord,mass_number,atomic_number,Zstar,ion_id)
    real(8), intent(in) :: z_coord,r_coord,mass_number,atomic_number,Zstar
    integer, intent(in) :: ion_id
    integer :: i, j

    i=int((z_coord-mesh_par%z_min_moving)/mesh_par%dz)+2
    j=int(r_coord/mesh_par%dr)+2

    static_ion(ion_id)%cmp(1) = z_coord
    static_ion(ion_id)%cmp(2) = r_coord
    static_ion(ion_id)%cmp(3) = mass_number
    static_ion(ion_id)%cmp(4) = atomic_number
    static_ion(ion_id)%cmp(5) = Zstar
    static_ion(ion_id)%cmp(6) = background_density_value(i,j)
    static_ion(ion_id)%cmp(7) = 1.d0
  end subroutine initialise_ion




  !--->Subdivide ppc into np_per_radius np_per_longitudinal
  subroutine from_ppc_to_npr_npz(particle_per_cell,np_per_r,np_per_z)
  integer, intent(INOUT) :: particle_per_cell
  integer, intent(OUT) :: np_per_r,np_per_z
  integer i,number_of_factors
  integer, allocatable, dimension(:) :: factors

  !---case single particle---!
  if(particle_per_cell==1) then
    np_per_r=1
    np_per_z=1
    return
  endif

  !verify input 'ppc' are not prime numbers
  !do i=1,1
   !do while(ISPRIME(ppc(i)))
   do while(ISPRIME(particle_per_cell))
    !ppc(i)=ppc(i)+1
    particle_per_cell=particle_per_cell+1
   enddo
  !enddo

  !subdivide ppc into np_per_xc,yc,zc
  !do i=1,6
   !ALLOCATE(factors(ppc(i)/2))
   ALLOCATE(factors(particle_per_cell/2))

   CALL PRIMEFACTORS(particle_per_cell,factors,number_of_factors)
   np_per_r=PRODUCT(factors(1:number_of_factors/2))
   np_per_z=PRODUCT(factors(number_of_factors/2+1:number_of_factors))

   deallocate(factors)
  !enddo
  end subroutine from_ppc_to_npr_npz

  FUNCTION ISPRIME(num)
  INTEGER, INTENT(IN) :: num  !input number
  INTEGER :: i
  LOGICAL :: ISPRIME

  ISPRIME=.TRUE.

  Do i=2,num-1
   IF( MOD(num,i) == 0 ) then
    ISPRIME=.FALSE.
    EXIT
   ENDIF
  EndDo
  END FUNCTION ISPRIME

  SUBROUTINE PRIMEFACTORS(num, factors, number_factors)
  INTEGER, INTENT(IN) :: num  !input number
  INTEGER,INTENT(OUT), DIMENSION((num/2))::factors !Array to store factors
  INTEGER, INTENT(INOUT) :: number_factors
  INTEGER :: i, n
  i = 2  !Eligible factor
  number_factors = 1  !Number of factors
  n = num !store input number into a temporary variable
  DO
   IF (MOD(n,i) == 0) THEN !If i divides 2, it is a factor
    factors(number_factors) = i
    number_factors = number_factors+1
    n = n/i
   ELSE
    i = i+1     !Not a factor. Move to next number
   END IF
   IF (n == 1) THEN
    !Since f is incremented after a factor is found
    number_factors = number_factors-1  !its value will be one more than the number of factors
    !Hence the value of number_factors is decremented
    EXIT
   END IF
  END DO
  END SUBROUTINE PRIMEFACTORS

  !---------------------!
  !---rigid advection---!
  !---------------------!
  subroutine ion_advection
    call remove_outofboundaries_ions
  end subroutine ion_advection

  subroutine remove_outofboundaries_ions
		integer :: i
		do i=1,ionisation%tot_ions
			if(static_ion(i)%cmp(7)>0.d0 .and. static_ion(i)%cmp(1)>mesh_par%z_max_moving) static_ion(i)%cmp(7)=0.d0
		enddo
	end subroutine remove_outofboundaries_ions



	subroutine inject_ions
		integer :: i,j,k,idx_z,idx_r
    real(8) :: hz,hr,z_coord,r_coord

    hz=mesh_par%dz/real(ionisation%np_per_z+1)
    hr=mesh_par%dr/real(ionisation%np_per_r+1)

		i=2
    do j=2,mesh_par%Nr_plasma-1
	    do idx_z=1,ionisation%np_per_z
	      do idx_r=1,ionisation%np_per_r
					do k=1,ionisation%tot_ions
						if(static_ion(k)%cmp(7)<=0.) goto 111
					enddo
          111 continue
          z_coord=mesh_par%z_min_moving+(i-2)*mesh_par%dz+hz*idx_z
          r_coord=(j-2)*mesh_par%dr+hr*idx_r
          call initialise_ion(z_coord,r_coord,ionisation%mass_number,ionisation%atomic_number,1.d0,k)
          mesh(i,j)%ni_bck = mesh(i,j)%ni_bck + static_ion(k)%cmp(5)/real(ionisation%particle_per_cell) * static_ion(k)%cmp(6)
        enddo
      enddo
    enddo
    end subroutine inject_ions

    !---from particle Zstart to mesh%Zstar ---!
    subroutine from_particleZstart_to_meshZstar
      integer :: i,j,k

      if(.not.ionisation%L_ionisation) then
        mesh(:,:)%Zstar=1.d0
        return
      endif

      mesh(:,:)%Zstar=0.d0
      do k=1,ionisation%tot_ions
        if(static_ion(k)%cmp(7)>0.d0) then
          i=int((static_ion(k)%cmp(1)-mesh_par%z_min_moving)/mesh_par%dz)
          j=int( static_ion(k)%cmp(2)/mesh_par%dr)
          mesh(i,j)%Zstar = mesh(i,j)%Zstar + static_ion(k)%cmp(5)/real(ionisation%particle_per_cell)
        endif
      enddo
    end subroutine from_particleZstart_to_meshZstar


    !---from particle Zstart to mesh%ni_bck ---!
    subroutine from_particleZstart_to_nplasmai
      integer :: i,j,k

      if(.not.ionisation%L_ionisation) then
        mesh(:,:)%ni_bck=0.d0
        return
      endif

      mesh(:,:)%ni_bck=0.d0
      do k=1,ionisation%tot_ions
        if(static_ion(k)%cmp(7)>0.d0) then
          i=int((static_ion(k)%cmp(1)-mesh_par%z_min_moving)/mesh_par%dz)
          j=int( static_ion(k)%cmp(2)/mesh_par%dr)
          mesh(i,j)%ni_bck = mesh(i,j)%ni_bck + static_ion(k)%cmp(5)/real(ionisation%particle_per_cell) * static_ion(k)%cmp(6)
        endif
      enddo
    end subroutine from_particleZstart_to_nplasmai


END MODULE
