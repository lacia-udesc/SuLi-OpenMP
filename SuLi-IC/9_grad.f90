!Subrotina para calcular a pressão dinâmica pelo método do gradiente conjugado
!Referencia: Casulli (1992 e 1999)

!!! Implementação 15/04/2014
! Leonardo Romero Monteiro

!!! Modificações
! Leonardo Romero Monteiro - 13/01/2015

SUBROUTINE graddin()

	USE velpre
	USE parametros
    USE cond
	USE obst
	USE omp
	USE omp_lib

	IMPLICIT NONE

	!===================================================================================================================
	!DECLARADO SOMENTE NA SUBROTINA (ou não precisam de entrada)

	real(8) :: alfapr, alfamupr, alfadipr, betapr, betamupr
	real(8), dimension(nx1,ny,nz) :: matspri
	real(8), dimension(nx,ny1,nz) :: matsprj
	real(8), dimension(nx,ny,nz1) :: matsprk
	real(8), dimension(nx,ny,nz) :: matqpr, matapripos, mataprineg, mataprjpos, mataprjneg, mataprkpos, mataprkneg
	real(8), dimension(0:nx1,0:ny1,0:nz1) :: matdpr, matepr, erropr, erroppr, mppr

	!contadores
	integer :: i, j, k, cont

	!auxiliares
	real(8) :: aux1, aux2, aux3

	real(8), dimension(nx1,ny,nz) :: rhox
	real(8), dimension(nx,ny1,nz) :: rhoy
	real(8), dimension(nx,ny,nz1) :: rhoz

	!Inicialização de variáveis


	!===================================================================================================================
	!RESOLUÇÃO DO PROBLEMA
	!===================================================================================================================
	cont = 0

	!%%%!-- Método do Gradiente Conjugado - Para Pressão Dinâmica --!%%%!

	!Tempo da montagem inicial de matrizes p/ Fortran e OpenMP
	CALL cpu_time(start_outros2_f90)
	start_outros2_omp = omp_get_wtime()

	call interpx_cf(rho,nx,ny,nz,rhox) !(nx1,ny,nz)
	call interpy_cf(rho,nx,ny,nz,rhoy) !(nx,ny1,nz)
	call interpz_cf(rho,nx,ny,nz,rhoz) !(nx,ny,nz1)

	matspri = dt / (dx*dx*rhox)
	matsprj = dt / (dy*dy*rhoy)
	matsprk = dt / (dz*dz*rhoz)

	! Matrizes p e q
	do k = 1, nz
		do j = 1, ny
			do i = 1, nx
				matdpr(i,j,k) = matspri(i+1,j,k) + matspri(i,j,k) + matsprj(i,j+1,k) + matsprj(i,j,k) + matsprk(i,j,k+1) + matsprk(i,j,k)
				matqpr(i,j,k) = ( u(i,j,k) - u(i+1,j,k) )/dx + ( v(i,j,k)-v(i,j+1,k) )/dy  + (w(i,j,k) - w(i,j,k+1)) /dz
			enddo
		enddo
	enddo

	!dentro do obstáculo, divergência nula
	!do j = 1, ny
	!do i = 1, nx
		!matdpr(i,j,0:kw(i,j)) = 0.
	!	matqpr(i,j,0:kw(i,j)) = 0.
	!enddo
	!enddo

	! Condições de contorno, von Neumann
	if (ccx0.eq.0) then  ! Condição periódica
		matdpr(0,:,:)   = matdpr(nx,:,:)
		matdpr(nx1,:,:) = matdpr(1,:,:)
	else
		matdpr(0,:,:)   = matdpr(1,:,:)
		matdpr(nx1,:,:) = matdpr(nx,:,:)
	endif

	if (ccy0.eq.0) then  ! Condição periódica
		matdpr(:,0,:)   = matdpr(:,ny,:)
		matdpr(:,ny1,:) = matdpr(:,1,:)
	else
		matdpr(:,0,:)   = matdpr(:,1,:)
		matdpr(:,ny1,:) = matdpr(:,ny,:)
	endif

	matdpr(:,:,0)   = matdpr(:,:,1)
	matdpr(:,:,nz1) = matdpr(:,:,nz)
	
	! Normalização da pressão dinâmica
	do k = 0, nz1
		do j = 0, ny1
			do i = 0, nx1
				matepr(i,j,k) = prd1(i,j,k) * sqrt(matdpr(i,j,k))
			enddo
		enddo
	enddo

	CALL cpu_time(end_outros2_f90)
	end_outros2_omp = omp_get_wtime()

	write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
	soma_outros2_f90 = soma_outros2_f90 + (end_outros2_f90 - start_outros2_f90)
	soma_outros2_omp = soma_outros2_omp + (end_outros2_omp - start_outros2_omp)
	write(*,*) "Tempo individual de iniciar matrizes p/ Fortran:", end_outros2_f90 - start_outros2_f90
	write(*,*) "Tempo acumulado de iniciar matrizes p/ Fortran", soma_outros2_f90
	write(*,*) "Tempo individual de iniciar matrizes p/ OpenMP:", end_outros2_omp - start_outros2_omp
	write(*,*) "Tempo acumulado de iniciar matrizes p/ OpenMP", soma_outros2_omp

	!Tempo de norm. e 1_erro p/ Fortran e OpenMP
	CALL cpu_time(start_outros_f90)
	start_outros_omp = omp_get_wtime()

	!Normalização das matrizes e cálculo do primeiro erro
	erropr =   0.
	alfamupr = 0.
	do k = 1, nz
		do j = 1, ny
			do i = 1, nx

				matapripos(i,j,k) = matspri(i+1,j,k)/sqrt(matdpr(i,j,k)*matdpr(i+1,j,k))
				mataprineg(i,j,k) = matspri(i,j,k)/sqrt(matdpr(i,j,k)*matdpr(i-1,j,k))

				mataprjpos(i,j,k) = matsprj(i,j+1,k)/sqrt(matdpr(i,j,k)*matdpr(i,j+1,k))
				mataprjneg(i,j,k) = matsprj(i,j,k)/sqrt(matdpr(i,j,k)*matdpr(i,j-1,k))

				mataprkpos(i,j,k) = matsprk(i,j,k+1)/sqrt(matdpr(i,j,k)*matdpr(i,j,k+1))
				mataprkneg(i,j,k) = matsprk(i,j,k)/sqrt(matdpr(i,j,k)*matdpr(i,j,k-1))

				erropr(i,j,k) = matepr(i,j,k) - matapripos(i,j,k) * matepr(i+1,j,k) - mataprineg(i,j,k) * matepr(i-1,j,k) &
					- mataprjpos(i,j,k) * matepr(i,j+1,k) - mataprjneg(i,j,k) * matepr(i,j-1,k) &
					- mataprkpos(i,j,k) * matepr(i,j,k+1) - mataprkneg(i,j,k) * matepr(i,j,k-1) - matqpr(i,j,k)/sqrt(matdpr(i,j,k))

				alfamupr = alfamupr + erropr(i,j,k) * erropr(i,j,k)
			enddo
		enddo
	enddo

	if (ccx0.eq.0) then  ! Condição periódica
		erropr(0,:,:)   = erropr(nx,:,:)
		erropr(nx1,:,:) = erropr(1,:,:)
	else
		erropr(0,:,:)   = erropr(1,:,:)
		erropr(nx1,:,:) = erropr(nx,:,:)
	endif

	if (ccy0.eq.0) then  ! Condição periódica
		erropr(:,0,:)   = erropr(:,ny,:)
		erropr(:,ny1,:) = erropr(:,1,:)
	else
		erropr(:,0,:)   = erropr(:,1,:)
		erropr(:,ny1,:) = erropr(:,ny,:)
	endif

	erropr(:,:,0)   = erropr(:,:,1)
	erropr(:,:,nz1) = erropr(:,:,nz)

	erroppr = erropr

	CALL cpu_time(end_outros_f90)
	end_outros_omp = omp_get_wtime()

	write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
	soma_outros_f90 = soma_outros_f90 + (end_outros_f90 - start_outros_f90)
	soma_outros_omp = soma_outros_omp + (end_outros_omp - start_outros_omp)
	write(*,*) "Tempo individual de norm. e 1_erro p/ Fortran:", end_outros_f90 - start_outros_f90
	write(*,*) "Tempo acumulado de norm. e 1_erro p/ Fortran", soma_outros_f90
	write(*,*) "Tempo individual de norm. e 1_erro p/ OpenMP:", end_outros_omp - start_outros_omp
	write(*,*) "Tempo acumulado de norm. e 1_erro p/ OpenMP", soma_outros_omp

	!Tempo da redução do erro p/ Fortran e OpenMP
	CALL cpu_time(start_outros4_f90)
	start_outros4_omp = omp_get_wtime()

	!%%%%%%%%%%%%%   loop da redução do erro   %%%%%%%%%%%%%%!
	do while ((abs(alfamupr) > (0.0001/(nx*ny*nz))) .and. (cont < 10000) )

		if (cont == 9999) write(*,*) "pulou pressão; ", "erro =", abs(alfamupr)
			! OTIMIZAR CÓDIGO
			cont = cont +1

			!inicialização
			!alfapr   = 0.
			alfamupr = 0.
			alfadipr = 0.
			!betapr   = 0.
			betamupr = 0.
			!mppr     = 0.

			! Parâmetro mp e alfa

			do k = 1, nz
				do j = 1, ny
					do i = 1, nx
						mppr(i,j,k) = erroppr(i,j,k) - erroppr(i+1,j,k) * matapripos(i,j,k) - erroppr(i-1,j,k) * mataprineg(i,j,k) &
						- erroppr(i,j+1,k) * mataprjpos(i,j,k) - erroppr(i,j-1,k) * mataprjneg(i,j,k) & 
						- erroppr(i,j,k+1) * mataprkpos(i,j,k) - erroppr(i,j,k-1) * mataprkneg(i,j,k)
					enddo
				enddo
			enddo

			do k = 1, nz
				do j = 1, ny
					do i = 1, nx
						alfamupr = alfamupr + erropr(i,j,k) * erropr(i,j,k)
						!alfamupr = alfamupr + dot_product(erropr(i,j,k),erropr(i,j,k))				### DIMENSÃO INAPROPRIADA ###
					enddo
				enddo
			enddo

			do k = 1, nz
				do j = 1, ny
					do i = 1, nx
						alfadipr = alfadipr + erroppr(i,j,k) * mppr(i,j,k)
						!alfadipr = alfadipr + dot_product(erroppr(i,j,k),mppr(i,j,k))				### DIMENSÃO INAPROPRIADA ###
					enddo
				enddo
			enddo

			alfapr = alfamupr / alfadipr
	
			! Recálculo das matrizes e, erro e parâmetro beta

			do k = 1, nz
				do j = 1, ny
					do i = 1, nx
						matepr(i,j,k)  = matepr(i,j,k)  - alfapr * erroppr(i,j,k)
						!matepr(i,j,k) = matepr(i,j,k) - dot_product(alfapr,erroppr(i,j,k))				### DIMENSÃO INAPROPRIADA ###
						erropr(i,j,k) = erropr(i,j,k) - alfapr * mppr(i,j,k)
						!erropr(i,j,k) = erropr(i,j,k) - dot_product(alfapr,mppr(i,j,k))				### DIMENSÃO INAPROPRIADA ###
						betamupr = betamupr + erropr(i,j,k) * erropr(i,j,k)
						!betamupr = betamupr + dot_product(erropr(i,j,k),erropr(i,j,k))				### DIMENSÃO INAPROPRIADA ###
					enddo
				enddo
			enddo

			betapr = betamupr/alfamupr

			! Recálculo de errop
			do k = 1, nz
				do j = 1, ny
					do i = 1, nx
						erroppr(i,j,k) = erropr(i,j,k) + betapr * erroppr(i,j,k)
						!erroppr(i,j,k) = erropr(i,j,k) + dot_product(betapr,erroppr(i,j,k))				### DIMENSÃO INAPROPRIADA ###
					enddo
				enddo
			enddo

			! Condições de contorno

		if (ccx0.eq.0) then  ! Condição periódica
			matepr(0,:,:)   = matepr(nx,:,:)
			matepr(nx1,:,:) = matepr(1,:,:)
			erroppr(0,:,:)   = erroppr(nx,:,:)
			erroppr(nx1,:,:) = erroppr(1,:,:)
		else
			matepr(0,:,:)   = matepr(1,:,:)
			matepr(nx1,:,:) = matepr(nx,:,:)
			erroppr(0,:,:)   = erroppr(1,:,:)
			erroppr(nx1,:,:) = erroppr(nx,:,:)
		endif

		if (ccy0.eq.0) then  ! Condição periódica
			matepr(:,0,:)   = matepr(:,ny,:)
			matepr(:,ny1,:) = matepr(:,1,:)
			erroppr(:,0,:)   = erroppr(:,ny,:)
			erroppr(:,ny1,:) = erroppr(:,1,:)
		else
			matepr(:,0,:)   = matepr(:,1,:)
			matepr(:,ny1,:) = matepr(:,ny,:)
			erroppr(:,0,:)   = erroppr(:,1,:)
			erroppr(:,ny1,:) = erroppr(:,ny,:)
		endif

			matepr(:,:,0)   = matepr(:,:,1)
			matepr(:,:,nz1) = matepr(:,:,nz)

			erroppr(:,:,0)   = erroppr(:,:,1)
			erroppr(:,:,nz1) = erroppr(:,:,nz)
		
		write(*,*) "Contador: ", cont			!### PEDRO ###

	enddo		!De qual "do" É ESTE "enddo"? Atribuí como sendo do "do while".	### PEDRO ###
	! OTIMIZAR CÓDIGO

	CALL cpu_time(end_outros4_f90)
	end_outros4_omp = omp_get_wtime()

	write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
	soma_outros4_f90 = soma_outros4_f90 + (end_outros4_f90 - start_outros4_f90)
	soma_outros4_omp = soma_outros4_omp + (end_outros4_omp - start_outros4_omp)
	write(*,*) "Tempo individual da redução do erro p/ Fortran:", end_outros4_f90 - start_outros4_f90
	write(*,*) "Tempo acumulado da redução do erro p/ Fortran", soma_outros4_f90
	write(*,*) "Tempo individual da redução do erro p/ OpenMP:", end_outros4_omp - start_outros4_omp
	write(*,*) "Tempo acumulado da redução do erro p/ OpenMP", soma_outros4_omp

	! Desnormalização de matriz e para a pressão dinâmica
	do k = 0, nz+1
		do j = 0, ny+1
			do i = 0, nx+1
				prd1(i,j,k) = matepr(i,j,k) / sqrt(matdpr(i,j,k))
			enddo
		enddo
	enddo

END SUBROUTINE graddin

SUBROUTINE pressh()

	USE velpre
	USE parametros

	IMPLICIT NONE

	integer :: i, j, k
	real(8) :: intz

	!WORK Z-PENCILS
	do k = nz,1 ,-1
		do j = 1, ny
			do i = 1, nx
				if (k == nz) then !top boundary conditions
					prd1(i,j,k)   = dz*gz*rho(i,j,k)

				elseif (k == nz-1) then !top boundary conditions

					prd1(i,j,k)   = dz*gz*(rho(i,j,k)+ rho(i,j,k+1))*0.5 + prd1(i,j,k+1)

				elseif (k < nz-1) then

					intz = dz * (rho(i,j,k) + 4.*rho(i,j,k+1) + rho(i,j,k+2)) /3.
					prd1(i,j,k)   = gz * intz + prd1(i,j,k+2) 

				endif

			enddo
		enddo
	enddo

END SUBROUTINE pressh
