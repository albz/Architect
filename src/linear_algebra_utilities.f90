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

MODULE linear_algebra

USE my_types
USE use_my_types

	contains

	! LU Gauss Elimination
	subroutine GaussEliminationLU(AAA,bbb)

	real(8), dimension(2*bunchip%init_width_z*mesh_par%Nsample_z*(bunchip%init_width_r*mesh_par%Nsample_r/2), &
	2*bunchip%init_width_z*mesh_par%Nsample_z*(bunchip%init_width_r*mesh_par%Nsample_r/2)), intent(inout) :: AAA
	real(8), dimension((bunchip%init_width_r*mesh_par%Nsample_r/2)*2*bunchip%init_width_z*mesh_par%Nsample_z), intent(inout) :: bbb

  INTEGER n,i,j,k
	REAL(8) f,one

	n=2*bunchip%init_width_z*mesh_par%Nsample_z*(bunchip%init_width_r*mesh_par%Nsample_r/2)

    !---LU factorization---!
    DO i=1,n
        DO j=i+1,n
            f=-AAA(j,i)/AAA(i,i)
            AAA(j,i)=f
            DO k=i+1,n
                AAA(j,k)=f*AAA(i,k)+AAA(j,k)
            ENDDO
            bbb(j)=f*bbb(i)+bbb(j);
        ENDDO
    ENDDO

    !---backward substitution---!
    bbb(n)=bbb(n)/AAA(n,n)
    DO i=n-1,1,-1
        one=0.;
        DO j=i+1,n
        	one=one+AAA(i,j)*bbb(j)
		ENDDO
        bbb(i)=(bbb(i)-one)/AAA(i,i);
    ENDDO
	end subroutine GaussEliminationLU


	! Compute U matrix for Gauss Elimination
	subroutine Compute_U_matrix(AAA)
	real(8), dimension(2*bunchip%init_width_z*mesh_par%Nsample_z*(bunchip%init_width_r*mesh_par%Nsample_r/2), &
	2*bunchip%init_width_z*mesh_par%Nsample_z*(bunchip%init_width_r*mesh_par%Nsample_r/2)), intent(inout) :: AAA
    INTEGER n,i,j,k
	REAL(8) f,one

	n=2*bunchip%init_width_z*mesh_par%Nsample_z*(bunchip%init_width_r*mesh_par%Nsample_r/2)

    !---LU factorization---!
    !write(*,*)
    DO i=1,n
        DO j=i+1,n
            f=-AAA(j,i)/AAA(i,i)
            AAA(j,i)=f
            !write(*,*) 'f = ',f,j,i
            DO k=i+1,n
                AAA(j,k)=f*AAA(i,k)+AAA(j,k)
				!write(*,*) 'k=',k
            ENDDO
        ENDDO
    ENDDO

	end subroutine Compute_U_matrix

	! Backward substitution for Gauss Elimination
	subroutine Backward_substitution_Gauss_Elimination(AAA,BBB,ccc)
	real(8), dimension(2*bunchip%init_width_z*mesh_par%Nsample_z*(bunchip%init_width_r*mesh_par%Nsample_r/2), &
	2*bunchip%init_width_z*mesh_par%Nsample_z*(bunchip%init_width_r*mesh_par%Nsample_r/2)), intent(inout) :: AAA,BBB
	real(8), dimension((bunchip%init_width_r*mesh_par%Nsample_r/2)*2*bunchip%init_width_z*mesh_par%Nsample_z), intent(inout) :: ccc
    INTEGER n,i,j,k
	REAL(8) f,one

	! AAA is the coefficient matrix in the system AAA*x=ccc
	! BBB is the Upper triangular matrix from the LU decomposition of A

	n=2*bunchip%init_width_z*mesh_par%Nsample_z*(bunchip%init_width_r*mesh_par%Nsample_r/2)

	!---compute L^(-1) * ccc ---!
    DO i=1,n
        DO j=i+1,n
            f=-AAA(j,i)/AAA(i,i)
            AAA(j,i)=f
            DO k=i+1,n
                AAA(j,k)=f*AAA(i,k)+AAA(j,k)
            ENDDO
            ccc(j)=f*ccc(i)+ccc(j);
        ENDDO
    ENDDO

    !---backward substitution---!
    ccc(n)=ccc(n)/BBB(n,n)
    DO i=n-1,1,-1
        one=0.;
        DO j=i+1,n
        	one=one+BBB(i,j)*ccc(j)
		ENDDO
        ccc(i)=(ccc(i)-one)/BBB(i,i);
    ENDDO
	end subroutine Backward_substitution_Gauss_Elimination


!--- SOR Iterative Method ---!
	SUBROUTINE SOR(A2,B2,N,TOL,IterMax,W)
	REAL(8), intent(in)::A2(sim_parameters%dim_Laplacian,sim_parameters%dim_Laplacian)
	REAL(8), intent(inout)::B2(sim_parameters%dim_Laplacian)
	Real(8), intent(in)::tol
	Real(8), intent(in)::W
	real(8) :: X02(sim_parameters%dim_Laplacian),X(sim_parameters%dim_Laplacian),NORM,SUM1,SUM2
	Integer, intent(in)::N
	Integer, intent(in)::IterMax
	INTEGER::K=1,i,j,iter

		NORM=1.
		iter=1
		X=100.
		do while(NORM>TOL .and. iter<IterMax)

			DO I=1,N
				SUM1=0.0
				SUM2=0.0
				DO J=1,N
					IF (J.LT.I) SUM1=SUM1+A2(I,J)*X(J)
					IF (J.GT.I) SUM2=SUM2+A2(I,J)*X02(J)
				END DO
				X(I)=(1.0-W)*X02(I)+(W*(B2(I)-SUM1-SUM2))/A2(I,I)
			END DO

			NORM=MAXVAL(abs(X-X02))
			X02=X
			if(mod(iter,10)==0) write(*,'(A,1pe12.5,A,I4)') 'sor> norm:',norm, '  - iter: ',iter

			iter=iter+1
		end do !end-While
	B2=X02
	write(*,'(A,1pe12.5,A,I4)') 'SOR method ends. norm:',norm,'  - itermax: ',iter
	END SUBROUTINE


	SUBROUTINE SOR_sparse(A2,B2,N,TOL,IterMax,W)
	REAL(8), intent(in)::A2(sim_parameters%dim_Laplacian,sim_parameters%dim_Laplacian)
	REAL(8), intent(inout)::B2(sim_parameters%dim_Laplacian)
	Real(8), intent(in)::tol
	Real(8), intent(in)::W
	real(8) :: X02(sim_parameters%dim_Laplacian),X(sim_parameters%dim_Laplacian),NORM,SUM1,SUM2,diag
	Integer, intent(in)::N
	Integer, intent(in)::IterMax
	Integer, DIMENSION(:), ALLOCATABLE :: i_vector,j_vector
	Real(8), DIMENSION(:), ALLOCATABLE :: A_vector
	INTEGER::K=1,i,j,iter,count

	!---internal conversion from sparse to vector
	Allocate(A_vector(8*N))
	Allocate(i_vector(8*N))
	Allocate(j_vector(8*N))

	write(*,*) 'using the sparse matrix optimisation'

	count=0
	DO I=1,N
		DO J=1,N
			if(A2(I,J).ne.0.) then
				count=count+1
				A_vector(count)=A2(I,J)
				i_vector(count)=I
				j_vector(count)=J
			endIF
		endDO
	endDO
	i_vector(count+1)=i_vector(count)+1
	j_vector(count+1)=j_vector(count)+1
	write(*,*) 'Matrix vectorised. dim(',count,')'
	!---vectorised---!

		NORM=1.
		iter=1
		X=100.
		do while(NORM>TOL .and. iter<IterMax)
			SUM1=0.0
			SUM2=0.0
			DO I=1,count
				if(j_vector(I).lt.i_vector(I)) SUM1=SUM1+A_vector(I)*X(j_vector(I))
				if(j_vector(I).gt.i_vector(I)) SUM2=SUM2+A_vector(I)*X02(j_vector(I))
				if(j_vector(I).eq.i_vector(I)) diag=A_vector(I)
				if(i_vector(I+1).gt.i_vector(I)) then
					X(i_vector(I))=(1.0-W)*X02(i_vector(I))+(W*(B2(i_vector(I))-SUM1-SUM2))/diag
					SUM1=0.0
					SUM2=0.0
				endIF
			END DO

			NORM=MAXVAL(abs(X-X02))
			X02=X
			if(mod(iter,10)==0) write(*,'(A,1pe12.5,A,I4)') 'sor> norm:',norm, '  - iter: ',iter

			iter=iter+1
		end do !end-While
	B2=X02
	write(*,'(A,1pe12.5,A,I4)') 'SOR method ends. norm:',norm,'  - itermax: ',iter

	deallocate(A_vector)
	deallocate(i_vector)
	deallocate(j_vector)
	END SUBROUTINE SOR_sparse





	SUBROUTINE SOR_sparse_Matrix_sparse(A_vector,B,i_vector,j_vector,TOL,IterMax,W,dim_non_null_elements,dimB)
	Real(8), intent(in) :: tol,W
	Integer, intent(in) :: IterMax, dim_non_null_elements, dimB
	Integer, intent(in) :: i_vector(5*sim_parameters%dim_Laplacian), j_vector(5*sim_parameters%dim_Laplacian)
	Real(8), intent(in) :: A_vector(5*sim_parameters%dim_Laplacian)
	REAL(8), intent(inout) :: B(dimB)
	real(8) :: X02(dimB),X(dimB),NORM,SUM1,SUM2,diag
	INTEGER::K=1,i,j,iter

 		NORM=1.
 		iter=0
 		X = 0.
 		X02 = X
 		do while(NORM>TOL .and. iter<IterMax)
 			SUM1=0.0
 			SUM2=0.0
 			DO I=1,dim_non_null_elements
				!if(j_vector(I)==0) cycle
 				if(j_vector(I).lt.i_vector(I)) SUM1=SUM1+A_vector(I)*X(j_vector(I))
			  if(j_vector(I).gt.i_vector(I)) SUM2=SUM2+A_vector(I)*X02(j_vector(I))
 				if(j_vector(I).eq.i_vector(I)) diag=A_vector(I)
 				if(i_vector(I+1).gt.i_vector(I)) then
 					X(i_vector(I))=(1.0-W)*X02(i_vector(I))+(W*(B(i_vector(I))-SUM1-SUM2))/diag
 					SUM1=0.0
 					SUM2=0.0
 				endIF
 			END DO

 			NORM=MAXVAL(abs(X-X02))
 			X02=X
 			if(mod(iter,1000)==0) write(*,'(A,1pe12.5,A,I8)') 'sor> norm:',norm, '  - iter: ',iter

 			iter=iter+1
 		end do !end-While
 	B=X02
	write(*,'(A,1pe12.5,A,I8)') 'SOR method ends. norm:',norm,'  - itermax: ',iter

	END SUBROUTINE SOR_sparse_Matrix_sparse



	SUBROUTINE CG(A_vector,B,i_vector,j_vector,TOL,IterMax,W,dim_non_null_elements,dimB,dimL)
	Real(8), intent(in) :: tol,W
	Integer, intent(in) :: IterMax, dim_non_null_elements, dimB, dimL
	Integer, intent(in) :: i_vector(dim_non_null_elements), j_vector(dim_non_null_elements)
	Real(8), intent(in) :: A_vector(dim_non_null_elements)
	REAL(8), intent(inout) :: B(dimB)
	real(8) :: NORM,SUM1,SUM2,diag
	INTEGER :: kk,i,j,iter
	REAL(8) :: x(dimB), p(dimB), q(dimB), r(dimB), v(dimB)
	REAL(8) :: alpha, rho, rho0

	x = 0.D0
	v=0.d0
	do k=1,dim_non_null_elements
		i=i_vector(k)
		j=j_vector(k)
		v(j)=v(j)+A_vector(k)*x(j)
	enddo
	r = B - v
	p = r
	rho = dot_product(r,r)
	v=0.d0
	do k=1,dim_non_null_elements
		i=i_vector(k)
		j=j_vector(k)
		v(j)=v(j)+A_vector(k)*p(j)
	enddo
	q = v
	alpha = rho / dot_product(p,q)
	x = x + alpha*p
	r = r - alpha*q
	Do kk =1,IterMax
		rho0=rho
		rho=dot_product(r,r)
		p=r+(rho/rho0)*p
		v=0.d0
		do k=1,dim_non_null_elements
			i=i_vector(k)
			j=j_vector(k)
			v(j)=v(j)+A_vector(k)*p(j)
		enddo
		q=v
		alpha = rho / dot_product(p,q)
		x = x + alpha*p
		r = r - alpha*q
		if( MAXVAL(r)<TOL) EXIT
		write(*,'(A,1pe12.5,A,I8)') 'CG method norm:',MAXVAL(r),'  - iter: ',kk
	EndDO
	B=x
	END SUBROUTINE CG

	SUBROUTINE CG_sparse(A_vector,B,i_vector,j_vector,TOL,IterMax,WSOR,dim_non_null_elements,dimB)
	Real(8), intent(in) :: tol,WSOR
	Integer, intent(in) :: IterMax, dim_non_null_elements, dimB
	Integer, intent(in) :: i_vector(5*sim_parameters%dim_Laplacian), j_vector(5*sim_parameters%dim_Laplacian)
	Real(8), intent(inout) :: A_vector(5*sim_parameters%dim_Laplacian)
	REAL(8), intent(inout) :: B(dimB)
	real(8) :: x(dimB),p(dimB),q(dimB),r(dimB),w(dimB),utility(dimB),NORM,L(dim_non_null_elements)
	REAL(8) :: alpha, rho, rho0, beta, r_norm
	INTEGER :: k,i,j,iter,diag(dim_non_null_elements)

 		NORM=1.
 		iter=1

		x = 0.d0
		utility=0.d0
		Do k=1,dim_non_null_elements
			i=i_vector(k)
			j=j_vector(k)
			utility(i)=utility(i)+A_vector(k)*x(j)
		enddo
		r=B-utility
		r_norm=sqrt(dot_product(r,r))
		p=r
 		do while(NORM>TOL .and. iter<IterMax)
			utility=0.d0
			Do k=1,dim_non_null_elements
				i=i_vector(k)
				j=j_vector(k)
				utility(i)=utility(i)+A_vector(k)*p(j)
			enddo
			w=utility
			r_norm_inverse=1.d0/r_norm**2
			alpha =r_norm**2/dot_product(w,p)
			x=x+alpha*p
			r=r-alpha*w
			r_norm=sqrt(dot_product(r,r))
			beta=r_norm**2*r_norm_inverse
			p=r+beta*p
 			NORM=MAXVAL(abs(r))
 			if(mod(iter,100)==0) write(*,'(A,1pe12.5,A,I8)') 'CG> norm:',norm, '  - iter: ',iter
 			iter=iter+1
 		end do !end-While
 	B=x
	write(*,'(A,1pe12.5,A,I8)') 'CG method ends. norm:',norm,'  - itermax: ',iter
	END SUBROUTINE CG_sparse


	!***********************************************************************
	!**********************************************************
	!Modified (by Ema) from the sor subroutine taken from Numerical Recipes (f90)
	!Successive overrelaxation solution of equation (19.5.25) with Chebyshev acceleration. a , b ,
	!c , d , e , and f are input as the coefficients of the equation, each dimensioned to the grid
	!size J × J. u is input as the initial guess to the solution, usually zero, and returns with the
	!final value. rjac is input as the spectral radius of the Jacobi iteration, or an estimate of
	!it. Double precision is a good idea for J bigger than about 25.
	!
	!Plese note: differently from the algorithm presented in numerical recipes f90, this one
	! also works for a grid of dimension NxM when N and M can be even or odd (even-odd, even-even, odd-even, odd-odd).
	! The one presented in NRf90 only works for grids where N and M are both odd (when even, it skips the last lines/columns of the u matrix)
	!
	!**********************************************************
	!***********************************************************************
	SUBROUTINE sor_homo(a,b,c,d,e,f,u, omega_rule, omega_given, rjac, eps, maxits, &
	    minits, show_num_convergence, proceed_anyway_maxits ) !_homo stands for homogeneous (rhs of the equation is 0, except for what is given by Dirichlet nodes)
	  REAL(kind=8), DIMENSION(:,:), INTENT(INOUT) :: a,b,c,d,e,f
	  REAL(kind=8), DIMENSION(:,:), INTENT(INOUT) :: u
	  REAL(kind=8), INTENT(IN) :: rjac
	  integer(kind=4), intent(in) :: maxits, minits
	  real(kind=8):: u_norm, deltau_norm
	  real(kind=8), intent(in) :: eps, omega_given !eps: tolerance on the error; omega_given: given value for omega, it is used if omega_rule==2
	  integer, intent(in):: omega_rule ! if omega_rule==2 omega is imposed and equal to omega_given, otherwise if omega_rule==1 then it is computed from the rjac
	  logical, intent(in):: proceed_anyway_maxits
	  logical,intent(in):: show_num_convergence ! if true the program shows the numerical convergence data on the default output

	  REAL(kind=8), DIMENSION(size(a,1),size(a,2)) :: resid
	  INTEGER :: jmax,jm1,jm2,jm3, imax,im1,im2,im3, j_start, i_start
	  integer:: n,d1,d2
	  REAL(kind=8) :: omega
	  REAL(kind=8), DIMENSION(size(u,1),size(u,2)):: u_old

		!--- apply boundary conditions here ---!
		d1=size(a,1)
		d2=size(a,2)

		a(d1,:) =0.
		a(1,:)  =0.
		a(:,d2) =0.
		a(:,1)  =-1.
		f(:,1)  =0.

		b(d1,:) =0.
		b(1,:)  =0.
		b(:,d2) =-1.
		b(:,1)  =0.
		f(:,d2) =0.

		c(d1,:) =0.
		c(1,:)  =0.
		c(:,d2) =0.
		c(:,1)  =0.

		d(d1,:) =0.
		d(1,:)  =0.
		d(:,d2) =0.
		d(:,1)  =0.

		e(d1,:) =1.
		e(1,:)  =1.
		e(:,d2) =1.
		e(:,1)  =1.

	  ![Ema] Check whether all the sizes are the same and, if so, assigne the size value to jmax
	  jmax = size(a,2)
	  jm1 = jmax-1
	  jm2 = jmax-2
	  jm3 = jmax-3

	  imax = size(a,1)
	  im1 = imax-1
	  im2 = imax-2
	  im3 = imax-3

	  !since this subroutine is for a homogeneous problem, it is better to use as sopping criterium
	  ! |u-u_old|/|u| < eps, instead of using the criterium of the residual normalized to |rhs| lower than a threshold.
	  ! thus I comment out these lines twice, as to say that they stay here but they should probably never be used in this subroutine
	  !![Ema] This is a kind of norm of the matrix f also considering the contribution of the boundary nodes
	  !!anormf = sum(abs(f(1:imax,1:jmax))) ! this is with norm 1
	  !!anormf = norm2(f(1:imax,1:jmax)) ! this is with norm 2

	  if (omega_rule==1) then
	    omega = 1.0d0
	  else if (omega_rule==2) then
	    omega = omega_given
	  end if

	  ! I initialize the residual to 0, since all the boundary is currently Dririchlet boundary and there I know the solution, which means that there the residual is 0. On the other points the residual will be updated by the algorithm in the do-cycle
	  ! call init_matr_0(resid)
		resid=0.
	  !Compute initial norm of residual and terminate iteration when norm has been reduced by a
	  !factor EPS. This computation assumes initial u is zero (except for the dirichlet boundary).
	  !start the chessboard-like algorithm
	  do n=1,MAXITS
	  !First do the even-even and odd-odd squares of the grid, i.e., the red squares of the
	  !checkerboard:
	!------------------Compute the residual
	    !Non boundary nodes
	    !even-even:
	    resid(2:im1:2,2:jm1:2) = a(2:im1:2,2:jm1:2)*u(3:imax:2,2:jm1:2)+&
	      b(2:im1:2,2:jm1:2)*u(1:im2:2,2:jm1:2)+&
	      c(2:im1:2,2:jm1:2)*u(2:im1:2,3:jmax:2)+&
	      d(2:im1:2,2:jm1:2)*u(2:im1:2,1:jm2:2)+&
	      e(2:im1:2,2:jm1:2)*u(2:im1:2,2:jm1:2) - f(2:im1:2,2:jm1:2)
	    !odd-odd:
	    resid(3:im1:2,3:jm1:2) = a(3:im1:2,3:jm1:2)*u(4:imax:2,3:jm1:2)+&
	      b(3:im1:2,3:jm1:2)*u(2:im2:2,3:jm1:2)+&
	      c(3:im1:2,3:jm1:2)*u(3:im1:2,4:jmax:2)+&
	      d(3:im1:2,3:jm1:2)*u(3:im1:2,2:jm2:2)+&
	      e(3:im1:2,3:jm1:2)*u(3:im1:2,3:jm1:2) - f(3:im1:2,3:jm1:2)
	    ! For the boundary nodes (even-even and odd-odd all at once):
	    !top boundary
	    resid(1, 3:jm1:2) = a(1, 3:jm1:2)*u(2, 3:jm1:2)+&
	      c(1, 3:jm1:2)*u(1, 4:jmax:2)+&
	      d(1, 3:jm1:2)*u(1, 2:jm2:2)+&
	      e(1, 3:jm1:2)*u(1, 3:jm1:2) - f(1, 3:jm1:2)
	    !bottom boundary
	    j_start = mod(imax,2)+2
	    resid(imax, j_start:jm1:2) = b(imax, j_start:jm1:2)*u(im1, j_start:jm1:2)+&
	      c(imax, j_start:jm1:2)*u(imax, j_start+1:jmax:2)+&
	      d(imax, j_start:jm1:2)*u(imax, j_start-1:jm2:2)+&
	      e(imax, j_start:jm1:2)*u(imax, j_start:jm1:2) - f(imax, j_start:jm1:2)
	    !right boundary
	    i_start = mod(jmax,2)+2
	    resid(i_start:im1:2 ,jmax) = a(i_start:im1:2, jmax)*u(i_start+1:imax:2, jmax)+&
	      b(i_start:im1:2 ,jmax)*u(i_start-1:im2:2 ,jmax)+&
	      d(i_start:im1:2 ,jmax)*u(i_start:im1:2 ,jm1)+&
	      e(i_start:im1:2 ,jmax)*u(i_start:im1:2 ,jmax) - f(i_start:im1:2 ,jmax)
	    !left boundary
	    resid(3:im1:2 ,1) = a(3:im1:2 ,1)*u(4:imax:2 ,1)+&
	      b(3:im1:2, 1)*u(2:im2:2, 1)+&
	      c(3:im1:2, 1)*u(3:im1:2, 2)+&
	      e(3:im1:2, 1)*u(3:im1:2, 1) - f(3:im1:2, 1)
	    !corner top left
	    resid(1 ,1) = a(1 ,1)*u(2 ,1)+&
	      c(1, 1)*u(1, 2)+&
	      e(1, 1)*u(1, 1) - f(1, 1)
	    !corner bottom left
	    if (mod(imax,2)==1) then
	      resid(imax ,1) = b(imax ,1)*u(im1, 1)+&
	        c(imax ,1)*u(imax, 2)+&
	        e(imax ,1)*u(imax ,1) - f(imax ,1)
	    end if
	    !corner top right
	    if (mod(jmax,2)==1) then
	      resid(1, jmax) = a(1, jmax)*u(2, jmax)+&
	        d(1, jmax)*u(1, jm1)+&
	        e(1, jmax)*u(1, jmax) - f(1, jmax)
	    end if
	    !corner bottom right
	    if ( .not. xor(mod(jmax,2)==1,mod(imax,2)==1) ) then
	      resid(imax, jmax) = b(imax, jmax)*u(im1, jmax)+&
	        d(imax, jmax)*u(imax, jm1)+&
	        e(imax, jmax)*u(imax, jmax) - f(imax, jmax)
	    end if
	!---------------------------update u
	    !apply: u_new = u_old - omega*residual/e
	    u(2:imax:2,2:jmax:2) = u(2:imax:2,2:jmax:2) - omega*&
	      resid(2:imax:2,2:jmax:2)/e(2:imax:2,2:jmax:2)

	    u(1:imax:2,1:jmax:2) = u(1:imax:2,1:jmax:2)-omega*&
	      resid(1:imax:2,1:jmax:2)/e(1:imax:2,1:jmax:2)

	!---------------------------compute new omega
	    if (omega_rule==1) then
	      ![Ema] usage of function merge(tsource,fsource,mask): if mask is .true. it gives the elements in tsource, otherwise fsource.
	      !Updating omega according to Chebyshev acceleration (and not straight to the asymptotically optimal omega) (see page 859 of Numerical Recipes in F77).
	      omega = merge(1.0/(1.0-0.5*rjac**2), &
	        1.0/(1.0-0.25*rjac**2*omega), n == 1)
	    else if (omega_rule==2) then
	      ! do nothing, omega has just been set, I added this if-else to rememeber that there is also another option
	    end if

	    !Now do even-odd and odd-even squares of the grid, i.e., the black squares of the checker-
	    !board:
	!---------------------------compute the residuals
	    ! Non boundary nodes
	    ! odd-even
	    resid(3:im1:2,2:jm1:2) = a(3:im1:2,2:jm1:2)*u(4:imax:2,2:jm1:2)+&
	      b(3:im1:2,2:jm1:2)*u(2:im2:2,2:jm1:2)+&
	      c(3:im1:2,2:jm1:2)*u(3:im1:2,3:jmax:2)+&
	      d(3:im1:2,2:jm1:2)*u(3:im1:2,1:jm2:2)+&
	      e(3:im1:2,2:jm1:2)*u(3:im1:2,2:jm1:2)-f(3:im1:2,2:jm1:2)
	    ! even-odd
	    resid(2:im1:2,3:jm1:2) = a(2:im1:2,3:jm1:2)*u(3:imax:2,3:jm1:2)+&
	      b(2:im1:2,3:jm1:2)*u(1:im2:2,3:jm1:2)+&
	      c(2:im1:2,3:jm1:2)*u(2:im1:2,4:jmax:2)+&
	      d(2:im1:2,3:jm1:2)*u(2:im1:2,2:jm2:2)+&
	      e(2:im1:2,3:jm1:2)*u(2:im1:2,3:jm1:2)-f(2:im1:2,3:jm1:2)
	    !For the boundary nodes (even-odd and odd-even all at once):
	    !top boundary
	    resid(1, 2:jm1:2) = a(1, 2:jm1:2)*u(2, 2:jm1:2)+&
	      c(1, 2:jm1:2)*u(1, 3:jmax:2)+&
	      d(1, 2:jm1:2)*u(1, 1:jm2:2)+&
	      e(1, 2:jm1:2)*u(1, 2:jm1:2) - f(1, 2:jm1:2)
	    !bottom boundary
	    j_start = 3-mod(imax,2)
	    resid(imax, j_start:jm1:2) = b(imax, j_start:jm1:2)*u(im1, j_start:jm1:2)+&
	      c(imax, j_start:jm1:2)*u(imax, j_start+1:jmax:2)+&
	      d(imax, j_start:jm1:2)*u(imax, j_start-1:jm2:2)+&
	      e(imax, j_start:jm1:2)*u(imax, j_start:jm1:2) - f(imax, j_start:jm1:2)
	    !right boundary
	    i_start = 3-mod(jmax,2)
	    resid(i_start:im1:2 ,jmax) = a(i_start:im1:2, jmax)*u(i_start+1:imax:2, jmax)+&
	      b(i_start:im1:2 ,jmax)*u(i_start-1:im2:2 ,jmax)+&
	      d(i_start:im1:2 ,jmax)*u(i_start:im1:2 ,jm1)+&
	      e(i_start:im1:2 ,jmax)*u(i_start:im1:2 ,jmax) - f(i_start:im1:2 ,jmax)
	    !left boundary
	    resid(2:im1:2, 1) = a(2:im1:2, 1)*u(3:imax:2, 1)+&
	      b(2:im1:2, 1)*u(1:im2:2, 1)+&
	      c(2:im1:2, 1)*u(2:im1:2, 2)+&
	      e(2:im1:2, 1)*u(2:im1:2, 1) - f(2:im1:2, 1)
	    ! corner bottom left
	    if (mod(imax,2)==0) then
	      resid(imax ,1) = b(imax ,1)*u(im1, 1)+&
	        c(imax ,1)*u(imax, 2)+&
	        e(imax ,1)*u(imax ,1) - f(imax ,1)
	    end if
	    !corner top right
	    if (mod(jmax,2)==0) then
	      resid(1, jmax) = a(1, jmax)*u(2, jmax)+&
	        d(1, jmax)*u(1, jm1)+&
	        e(1, jmax)*u(1, jmax) - f(1, jmax)
	    end if
	    !corner bottom right
	    if ( xor(mod(jmax,2)==1,mod(imax,2)==1) ) then
	      resid(imax, jmax) = b(imax, jmax)*u(im1, jmax)+&
	        d(imax, jmax)*u(imax, jm1)+&
	        e(imax, jmax)*u(imax, jmax) - f(imax, jmax)
	    end if
	!----------------------------update u
	    u(1:imax:2,2:jmax:2) = u(1:imax:2,2:jmax:2)-omega*&
	      resid(1:imax:2,2:jmax:2)/e(1:imax:2,2:jmax:2)

	    u(2:imax:2,1:jmax:2) = u(2:imax:2,1:jmax:2)-omega*&
	      resid(2:imax:2,1:jmax:2)/e(2:imax:2,1:jmax:2)

	!----------------------------compute new omega
	    if (omega_rule==1) then
	      omega = 1.0/(1.0-0.25*rjac**2*omega)
	    else if (omega_rule==2) then
	      ! do nothing
	    end if

	!--------------- check the error (and whether the desired tollerance has been reached)
	    ! Compute the norm of the residual taking into account also the dirichlet nodes effect (i.e. start from the first node and end with the last, do not skip boundary nodes, "1:jmax" instead of "2:jm1")
	    u_norm = SQRT(SUM(u*u))!norm2(u) ! this is with norm 2
	    deltau_norm = SQRT(SUM((u-u_old)*(u-u_old)))
	    u_old = u

	    if (mod(n,1000)==0) write(*,'(A,1I6,A,1e14.5,A,1e14.5)') 'step :: ',n,' ---     err(norm2) ::',deltau_norm/u_norm,' ---      Absolute err(norm2) ::',deltau_norm

	    if ( (deltau_norm < eps*u_norm) .and. n>minits ) then
	      exit
	    end if

	  end do

	    !----- If the max number of iterations has been overcome, throw and error
	    if (n > MAXITS) then
	      if (.not. proceed_anyway_maxits) then
	        write(*,'(A)')'MAXITS exceeded in sor_omo'
					! stop
	      else
	        write(*,fmt="(A,1X,A)") "MAXITS exceeded, taking the last calculated value"
	      endif
	    else
	      write(*,fmt="(A,I8,A,1p1e14.5,A,1p1e14.5)") "SOR_homo converged at iteration ", n, '  with  err(norm2) ::',deltau_norm/u_norm,',      Absolute err(norm2) ::',deltau_norm
	    endif

	END SUBROUTINE sor_homo


END MODULE linear_algebra
