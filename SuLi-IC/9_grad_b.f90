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

	IMPLICIT NONE

	!===================================================================================================================
	!DECLARADO SOMENTE NA SUBROTINA (ou não precisam de entrada)

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

	!===================================================================================================================
	!RESOLUÇÃO DO PROBLEMA
	!===================================================================================================================
	cont = 0

	!%%%!-- Método do Gradiente Conjugado - Para Pressão Dinâmica --!%%%!

	matspri = dt / (dx*dx)
	matsprj = dt / (dy*dy)
	matsprk = dt / (dz*dz)


	call interpx_cf(rho,nx,ny,nz,rhox) !(nx1,ny,nz)
	call interpy_cf(rho,nx,ny,nz,rhoy) !(nx,ny1,nz)
	call interpz_cf(rho,nx,ny,nz,rhoz) !(nx,ny,nz1)

		do k = 1, nz
		do j = 1, ny
		do i = 1, nx1
			matspri(i,j,k) = matspri(i,j,k)/rhox(i,j,k)
		enddo
		do i = 1, nx
			matsprj(i,j,k) = matsprj(i,j,k)/rhoy(i,j,k)
			matsprk(i,j,k) = matsprk(i,j,k)/rhoz(i,j,k)
		enddo
		enddo
		j = ny1
		do i = 1, nx
			matsprj(i,j,k) = matsprj(i,j,k)/rhoy(i,j,k)
		enddo
		enddo

		k = nz1
		do j = 1, ny
		do i = 1, nx
			matsprk(i,j,k) = matsprk(i,j,k)/rhoz(i,j,k)
		enddo
		enddo


		! Matrizes p e q
		do k = 1, nz
		do j = 1, ny
		do i = 1, nx
			matdpr(i+1,j+1,k+1) = matspri(i+1,j,k) + matspri(i,j,k) + matsprj(i,j+1,k) + matsprj(i,j,k) + matsprk(i,j,k+1) + matsprk(i,j,k)
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

		!Normalização das matrizes s e cálculo do primeiro erro
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

		!%%%%%%%%%%%%%   loop da redução do erro   %%%%%%%%%%%%%%!
		do while ((abs(alfamupr) > (0.0001/(nx*ny*nz))) .and. (cont < 1000) )


		if (cont == 9999) write(*,*) "pulou pressão; ", "erro =", abs(alfamupr)

		cont = cont +1

			!inicialização
			alfapr   = 0.
			alfamupr = 0.
			alfadipr = 0.
			betapr   = 0.
			betamupr = 0.
			mppr     = 0.

			! Parâmetro mp e alfa

			do k = 1, nz
			do j = 1, ny
			do i = 1, nx
				mppr(i,j,k) = erroppr(i+1,j+1,k+1) - erroppr(i+2,j+1,k+1) * matapripos(i,j,k) &
				- erroppr(i+1,j+2,k+1) * mataprjpos(i,j,k) & 
				- erroppr(i+1,j+1,k+2) * mataprkpos(i,j,k)
			enddo
			enddo
			enddo		

			do k = 1, nz
			do j = 1, ny
			do i = 1, nx
				mppr(i,j,k) = mppr(i,j,k) - erroppr(i,j+1,k+1) * mataprineg(i,j,k) &
				- erroppr(i+1,j,k+1) * mataprjneg(i,j,k) & 
				- erroppr(i+1,j+1,k) * mataprkneg(i,j,k)
			enddo
			enddo
			enddo
			
			
			
			do k = 1, nz
			do j = 1, ny
			do i = 1, nx
				alfamupr = alfamupr + erropr(i+1,j+1,k+1) * erropr(i+1,j+1,k+1)
				alfadipr = alfadipr + erroppr(i+1,j+1,k+1) * mppr(i,j,k)
			enddo
			enddo
			enddo

			alfapr = alfamupr / alfadipr
			write(*,*) cont, alfapr
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

			! aqui fez mais sentido avançar todo o avanço, pois todos são deslocados
			do k = 2, nz+1
			do j = 2, ny+1
			do i = 2, nx+1
				erroppr(i,j,k) = erropr(i,j,k) + betapr * erroppr(i,j,k)
			enddo
			enddo
			enddo

			! Condições de contorno

		if (ccx0.eq.0) then  ! Condição periódica
			matepr(1,:,:)   = matepr(nx+1,:,:)
			matepr(nx1+1,:,:) = matepr(2,:,:)
			erroppr(1,:,:)   = erroppr(nx+1,:,:)
			erroppr(nx1+1,:,:) = erroppr(2,:,:)
		else
			matepr(1,:,:)   = matepr(2,:,:)
			matepr(nx1+1,:,:) = matepr(nx+1,:,:)
			erroppr(1,:,:)   = erroppr(2,:,:)
			erroppr(nx1+1,:,:) = erroppr(nx+1,:,:)
		endif

		if (ccy0.eq.0) then  ! Condição periódica
			matepr(:,1,:)   = matepr(:,ny+1,:)
			matepr(:,ny1+1,:) = matepr(:,2,:)
			erroppr(:,1,:)   = erroppr(:,ny+1,:)
			erroppr(:,ny1+1,:) = erroppr(:,2,:)
		else
			matepr(:,1,:)   = matepr(:,2,:)
			matepr(:,ny1+1,:) = matepr(:,ny+1,:)
			erroppr(:,1,:)   = erroppr(:,2,:)
			erroppr(:,ny1+1,:) = erroppr(:,ny+1,:)
		endif

			matepr(:,:,1)   = matepr(:,:,2)
			matepr(:,:,nz1+1) = matepr(:,:,nz+1)

			erroppr(:,:,1)   = erroppr(:,:,2)
			erroppr(:,:,nz1+1) = erroppr(:,:,nz+1)

		enddo
		!%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		! Desnormalização de matriz e para a pressão dinâmica
		do k = 0, nz+1
		do j = 0, ny+1
		do i = 0, nx+1
			prd1(i,j,k) = matepr(i+1,j+1,k+1) / sqrt(matdpr(i+1,j+1,k+1))
		enddo
		enddo
		enddo

	!===============================================================================================================

	END SUBROUTINE graddin


!######################################################################################
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





