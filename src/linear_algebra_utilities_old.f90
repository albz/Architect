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

	real(8), dimension(2*bunch_initialization%init_width_z*mesh_par%Nsample_z*(bunch_initialization%init_width_r*mesh_par%Nsample_r/2), &
	2*bunch_initialization%init_width_z*mesh_par%Nsample_z*(bunch_initialization%init_width_r*mesh_par%Nsample_r/2)), intent(inout) :: AAA
	real(8), dimension((bunch_initialization%init_width_r*mesh_par%Nsample_r/2)*2*bunch_initialization%init_width_z*mesh_par%Nsample_z), intent(inout) :: bbb

  INTEGER n,i,j,k
	REAL(8) f,one

	n=2*bunch_initialization%init_width_z*mesh_par%Nsample_z*(bunch_initialization%init_width_r*mesh_par%Nsample_r/2)

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
	real(8), dimension(2*bunch_initialization%init_width_z*mesh_par%Nsample_z*(bunch_initialization%init_width_r*mesh_par%Nsample_r/2), &
	2*bunch_initialization%init_width_z*mesh_par%Nsample_z*(bunch_initialization%init_width_r*mesh_par%Nsample_r/2)), intent(inout) :: AAA
    INTEGER n,i,j,k
	REAL(8) f,one

	n=2*bunch_initialization%init_width_z*mesh_par%Nsample_z*(bunch_initialization%init_width_r*mesh_par%Nsample_r/2)

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
	real(8), dimension(2*bunch_initialization%init_width_z*mesh_par%Nsample_z*(bunch_initialization%init_width_r*mesh_par%Nsample_r/2), &
	2*bunch_initialization%init_width_z*mesh_par%Nsample_z*(bunch_initialization%init_width_r*mesh_par%Nsample_r/2)), intent(inout) :: AAA,BBB
	real(8), dimension((bunch_initialization%init_width_r*mesh_par%Nsample_r/2)*2*bunch_initialization%init_width_z*mesh_par%Nsample_z), intent(inout) :: ccc
    INTEGER n,i,j,k
	REAL(8) f,one

	! AAA is the coefficient matrix in the system AAA*x=ccc
	! BBB is the Upper triangular matrix from the LU decomposition of A

	n=2*bunch_initialization%init_width_z*mesh_par%Nsample_z*(bunch_initialization%init_width_r*mesh_par%Nsample_r/2)

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

END MODULE linear_algebra
