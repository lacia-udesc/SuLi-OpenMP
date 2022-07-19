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

	!integer, parameter :: Tipo1 = selected_real_kind(10,10)
	real(8) :: alfapr, alfamupr, alfadipr, betapr, betamupr
	real(8), dimension(nx1,ny,nz) :: matspri
	real(8), dimension(nx,ny1,nz) :: matsprj
	real(8), dimension(nx,ny,nz1) :: matsprk
	real(8), dimension(nx,ny,nz) :: matqpr, matapripos, mataprineg, mataprjpos, mataprjneg, mataprkpos, mataprkneg, mppr
	real(8), dimension(nx1+1,ny1+1,nz1+1) :: matdpr, matepr, erropr, erroppr

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
	!CALL cpu_time(start_outros2_f90)
	!start_outros2_omp = omp_get_wtime()

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
			matdpr(i+1,j+1,k+1) = matspri(i+1,j,k) + matspri(i,j,k) + matsprj(i,j+1,k) + matsprj(i,j,k) + matsprk(i,j,k+1) + matsprk(i,j,k)
				matqpr(i,j,k) = (u(i,j,k) - u(i+1,j,k))/dx + (v(i,j,k)-v(i,j+1,k))/dy  + (w(i,j,k) - w(i,j,k+1))/dz
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
			matdpr(1,:,:)   = matdpr(nx1,:,:)
			matdpr(nx1+1,:,:) = matdpr(2,:,:)
	else
			matdpr(1,:,:)   = matdpr(2,:,:)
			matdpr(nx1+1,:,:) = matdpr(nx1,:,:)
	endif

	if (ccy0.eq.0) then  ! Condição periódica
			matdpr(:,1,:)   = matdpr(:,ny1,:)
			matdpr(:,ny1+1,:) = matdpr(:,2,:)
	else
			matdpr(:,1,:)   = matdpr(:,2,:)
			matdpr(:,ny1+1,:) = matdpr(:,ny1,:)
	endif

		matdpr(:,:,1)   = matdpr(:,:,2)
		matdpr(:,:,nz1+1) = matdpr(:,:,nz1)
	
	! Normalização da pressão dinâmica
	do k = 0, nz1
		do j = 0, ny1
			do i = 0, nx1
			matepr(i+1,j+1,k+1) = prd1(i,j,k) * sqrt(matdpr(i+1,j+1,k+1))
			enddo
		enddo
	enddo

	!CALL cpu_time(end_outros2_f90)
	!end_outros2_omp = omp_get_wtime()

	!write(*,*) !"~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
	!soma_outros2_f90 = soma_outros2_f90 + (end_outros2_f90 - start_outros2_f90)
	!soma_outros2_omp = soma_outros2_omp + (end_outros2_omp - start_outros2_omp)
	!write(*,*) "Tempo individual de iniciar matrizes p/ Fortran:", end_outros2_f90 - start_outros2_f90
	!write(*,*) "Tempo acumulado de iniciar matrizes p/ Fortran", soma_outros2_f90
	!write(*,*) "Tempo individual de iniciar matrizes p/ OpenMP:", end_outros2_omp - start_outros2_omp
	!write(*,*) "Tempo acumulado de iniciar matrizes p/ OpenMP", soma_outros2_omp

	!Tempo de norm. e 1_erro p/ Fortran e OpenMP
	!CALL cpu_time(start_outros_f90)
	!start_outros_omp = omp_get_wtime()

	!Normalização das matrizes e cálculo do primeiro erro
	erropr =   0.
	alfamupr = 0.
	do k = 1, nz
		do j = 1, ny
			do i = 1, nx

			matapripos(i,j,k) = matspri(i+1,j,k)/sqrt(matdpr(i+1,j+1,k+1)*matdpr(i+2,j+1,k+1))
			mataprineg(i,j,k) = matspri(i,j,k)/sqrt(matdpr(i+1,j+1,k+1)*matdpr(i,j+1,k+1))

			mataprjpos(i,j,k) = matsprj(i,j+1,k)/sqrt(matdpr(i+1,j+1,k+1)*matdpr(i+1,j+2,k+1))
			mataprjneg(i,j,k) = matsprj(i,j,k)/sqrt(matdpr(i+1,j+1,k+1)*matdpr(i+1,j,k+1))

			mataprkpos(i,j,k) = matsprk(i,j,k+1)/sqrt(matdpr(i+1,j+1,k+1)*matdpr(i+1,j+1,k+2))
			mataprkneg(i,j,k) = matsprk(i,j,k)/sqrt(matdpr(i+1,j+1,k+1)*matdpr(i+1,j+1,k))

			erropr(i+1,j+1,k+1) = matepr(i+1,j+1,k+1) - matapripos(i,j,k) * matepr(i+2,j+1,k+1) - mataprineg(i,j,k) * matepr(i,j+1,k+1) &
				- mataprjpos(i,j,k) * matepr(i+1,j+2,k+1) - mataprjneg(i,j,k) * matepr(i+1,j,k+1) &
				- mataprkpos(i,j,k) * matepr(i+1,j+1,k+2) - mataprkneg(i,j,k) * matepr(i+1,j+1,k) - matqpr(i,j,k)/sqrt(matdpr(i+1,j+1,k+1))

			alfamupr = alfamupr + erropr(i+1,j+1,k+1) * erropr(i+1,j+1,k+1)
			enddo
		enddo
	enddo

	write(*,*) "Cálculo do primeiro erro", alfamupr

	if (ccx0.eq.0) then  ! Condição periódica
			erropr(1,:,:)   = erropr(nx+1,:,:)
			erropr(nx1+1,:,:) = erropr(2,:,:)
	else
			erropr(1,:,:)   = erropr(2,:,:)
			erropr(nx1+1,:,:) = erropr(nx1,:,:)
	endif

	if (ccy0.eq.0) then  ! Condição periódica
			erropr(:,1,:)   = erropr(:,ny+1,:)
			erropr(:,ny1+1,:) = erropr(:,2,:)
	else
			erropr(:,1,:)   = erropr(:,2,:)
			erropr(:,ny1+1,:) = erropr(:,ny1,:)
	endif

		erropr(:,:,1)   = erropr(:,:,2)
		erropr(:,:,nz1+1) = erropr(:,:,nz1)

	erroppr = erropr

	!CALL cpu_time(end_outros_f90)
	!end_outros_omp = omp_get_wtime()

	!write(*,*) !"~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
	!soma_outros_f90 = soma_outros_f90 + (end_outros_f90 - start_outros_f90)
	!soma_outros_omp = soma_outros_omp + (end_outros_omp - start_outros_omp)
	!write(*,*) "Tempo individual de norm. e 1_erro p/ Fortran:", end_outros_f90 - start_outros_f90
	!write(*,*) "Tempo acumulado de norm. e 1_erro p/ Fortran", soma_outros_f90
	!write(*,*) "Tempo individual de norm. e 1_erro p/ OpenMP:", end_outros_omp - start_outros_omp
	!write(*,*) "Tempo acumulado de norm. e 1_erro p/ OpenMP", soma_outros_omp

	!Tempo da redução do erro p/ Fortran e OpenMP
	CALL cpu_time(start_outros4_f90)
	start_outros4_omp = omp_get_wtime()

	!%%%%%%%%%%%%%   loop da redução do erro   %%%%%%%%%%%%%%!
	do while ((abs(alfamupr) > (0.0001/(nx*ny*nz))) .and. (cont < 10) )
	
		CALL cpu_time(fortran_start_grad_1)
		omp_start_grad_1 = omp_get_wtime()

		if (cont == 9999) write(*,*) "pulou pressão; ", "erro =", abs(alfamupr)
	
		! OTIMIZAR CÓDIGO
		cont = cont +1

		!inicialização
		alfapr   = 0.
		alfamupr = 0.
		alfadipr = 0.
		betapr   = 0.
		betamupr = 0.
		mppr = 0.

		! Parâmetro mp e alfa

		
		!!!	###########################################################################################################################################################
		
		do k = 1, nz
			do j = 1, ny
				do i = 1, nx
				mppr(i,j,k) = erroppr(i+1,j+1,k+1) &
				- erroppr(i+2,j+1,k+1) * matapripos(i,j,k) &
				- erroppr(i,j+1,k+1) * mataprineg(i,j,k) &
				- erroppr(i+1,j+2,k+1) * mataprjpos(i,j,k) & 
				- erroppr(i+1,j,k+1) * mataprjneg(i,j,k) & 
				- erroppr(i+1,j+1,k+2) * mataprkpos(i,j,k)
				- erroppr(i+1,j+1,k) * mataprkneg(i,j,k)
				enddo
			enddo
		enddo
		
		!!!	###########################################################################################################################################################
		
		do k = 1, nz
			do j = 1, ny
				do i = 1, nx
				alfamupr = alfamupr + erropr(i+1,j+1,k+1) * erropr(i+1,j+1,k+1)
				alfadipr = alfadipr + erroppr(i+1,j+1,k+1) * mppr(i,j,k)
				enddo
			enddo
		enddo
	
		!write(*,*) "Segundo Alfamupr", alfamupr

		alfapr = alfamupr / alfadipr

		!write(*,*) "Primeiro alfapr", alfapr
		
		CALL cpu_time(fortran_end_grad_1)
		omp_end_grad_1 = omp_get_wtime()

		soma_grad_1_f90 = soma_grad_1_f90 + (fortran_end_grad_1 - fortran_start_grad_1)
		soma_grad_1_omp = soma_grad_1_omp + (omp_end_grad_1 - omp_start_grad_1)

		CALL cpu_time(fortran_start_grad_2)
		omp_start_grad_2 = omp_get_wtime()

		! Recálculo das matrizes e, erro e parâmetro beta
		do k = 1, nz
			do j = 1, ny
				do i = 1, nx
				matepr(i+1,j+1,k+1) = matepr(i+1,j+1,k+1) - alfapr * erroppr(i+1,j+1,k+1)
				erropr(i+1,j+1,k+1) = erropr(i+1,j+1,k+1) - alfapr * mppr(i,j,k)
				betamupr = betamupr + erropr(i+1,j+1,k+1) * erropr(i+1,j+1,k+1)
				enddo
			enddo
		enddo
		
		betapr = betamupr/alfamupr

		! Recálculo de errop
		do k = 1, nz
			do j = 1, ny
				do i = 1, nx
				erroppr(i+1,j+1,k+1) = erropr(i+1,j+1,k+1) + betapr * erroppr(i+1,j+1,k+1)
				enddo
			enddo
		enddo

		CALL cpu_time(fortran_end_grad_2)
		omp_end_grad_2 = omp_get_wtime()

		soma_grad_2_f90 = soma_grad_2_f90 + (fortran_end_grad_2 - fortran_start_grad_2)
		soma_grad_2_omp = soma_grad_2_omp + (omp_end_grad_2 - omp_start_grad_2)

		CALL cpu_time(start_outros5_f90)
		start_outros5_omp = omp_get_wtime()

		! Condições de contorno

		if (ccx0.eq.0) then  ! Condição periódica
			do k = 1, nz2							! SERÁ QUE NÃO DEVERIA COMEÇAR O LOOP COM 0 AO INVÉS DE 1?
				do j = 1, ny2
					matepr(1,j,k) = matepr(nx+1,j,k)
					matepr(nx1+1,j,k) = matepr(2,j,k)
					erroppr(1,j,k) = erroppr(nx+1,j,k)
					erroppr(nx1+1,j,k) = erroppr(2,j,k)
				enddo
			enddo
		else
			do k = 1, nz2									! SERÁ QUE NÃO DEVERIA COMEÇAR O LOOP COM 0 AO INVÉS DE 1?
				do j = 1, ny2
					matepr(1,j,k) = matepr(2,j,k)
					matepr(nx1+1,j,k) = matepr(nx+1,j,k)
					erroppr(1,j,k) = erroppr(2,j,k)
					erroppr(nx1+1,j,k) = erroppr(nx+1,j,k)
				enddo
			enddo
		endif

		if (ccy0.eq.0) then  ! Condição periódica
			do k = 1, nz2									! SERÁ QUE NÃO DEVERIA COMEÇAR O LOOP COM 0 AO INVÉS DE 1?
				do i = 1, nx2
					matepr(i,1,k) = matepr(i,ny+1,k)
					matepr(i,ny1+1,k) = matepr(i,2,k)
					erroppr(i,1,k) = erroppr(i,ny+1,k)
					erroppr(i,ny1+1,k) = erroppr(i,2,k)
				enddo
			enddo			
		else
			do k = 1, nz2									! SERÁ QUE NÃO DEVERIA COMEÇAR O LOOP COM 0 AO INVÉS DE 1?
				do i = 1, nx2
					matepr(i,1,k) = matepr(i,2,k)
					matepr(i,ny1+1,k) = matepr(i,ny+1,k)
					erroppr(i,1,k) = erroppr(i,2,k)
					erroppr(i,ny1+1,k) = erroppr(i,ny+1,k)
				enddo
			enddo
		endif

		do j = 1, ny2
			do i = 1, nx2
					matepr(i,j,1) = matepr(i,j,2)
					matepr(i,j,nz1+1) = matepr(i,j,nz+1)
					erroppr(i,j,1) = erroppr(i,j,2)
					erroppr(i,j,nz1+1) = erroppr(i,j,nz+1)
			enddo
		enddo

		CALL cpu_time(end_outros5_f90)
		end_outros5_omp = omp_get_wtime()

		soma_outros5_f90 = soma_outros5_f90 + (end_outros5_f90 - start_outros5_f90)
		soma_outros5_omp = soma_outros5_omp + (end_outros5_omp - start_outros5_omp)
		
		write(*,*) "Contador		", "alfamupr		", "alfapr		", "betapr		", "betamupr		", "alfadipr		"			!### PEDRO ###
		write(*,*) cont, alfamupr, alfapr, betapr, betamupr, alfadipr						!### PEDRO ###

		write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"

	enddo
	! OTIMIZAR CÓDIGO

	CALL cpu_time(end_outros4_f90)
	end_outros4_omp = omp_get_wtime()

	write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
	write(*,*) "Tempo acumulado de 'Parâmetro mp e alfa' p/ Fortran", soma_grad_1_f90
	write(*,*) "Tempo acumulado de 'Parâmetro mp e alfa' p/ OpenMP", soma_grad_1_omp
	write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
	write(*,*) "Tempo acumulado de 'Recálculo de matrizes e erro' p/ Fortran", soma_grad_2_f90
	write(*,*) "Tempo acumulado de 'Recálculo de matrizes e erro' p/ OpenMP", soma_grad_2_omp
	write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
	write(*,*) "Tempo acumulado de 'Condições de contorno' p/ Fortran", soma_outros5_f90
	write(*,*) "Tempo acumulado de 'Condições de contorno' p/ OpenMP", soma_outros5_omp

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
			prd1(i,j,k) = matepr(i+1,j+1,k+1) / sqrt(matdpr(i+1,j+1,k+1))
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
