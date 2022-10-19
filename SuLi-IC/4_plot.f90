!Subrotina definir as condições de contorno das velocidades e suas influências
!Referência: Gotoh, 2013

!Implementação em 15/04/2014
!Leonardo Romero Monteiro

!Modificações
!Leonardo Romero Monteiro em 12/01/2015

!Correções
!Leonardo Romero Monteiro & Mariana De Cesaro em 14/06/2022

!SUBROTINAS CONTIDAS NESTE MÓDULO: "plot_i", "plot_f", "plot_atrib" & "est"

SUBROUTINE plot_i()

	USE paodeplot_i

	IMPLICIT NONE

	integer :: i, j, k

	!Criar pastas, se já não existirem

	call system('mkdir -p arquivos')
	call system('mkdir -p dados')

	!Terminal inicial
	write(*,*) "Malha (nx x ny x nz):", nx, 'x', ny, 'x', nz
	!write(*,*) "nx = ", nx,"ny = , 'x',, n, 'x',, "nz = ", nz
	write(*,*) "Intervalo Espacial (dx x dy x dz):", dx, 'x', dy, 'x', dz
	!write(*,*) "dx = ", dx,"dy = ", dy, "dz = ", dz
	write(*,*) "Parâmetros temporais"
	write(*,*) "dt = ", dt,"dt_frame = ", dt_frame, "ts = ", ts
	write(*,*) "Parâmetros de modelo"
	write(*,*) "der = ", der,"advectivo = ", adv_type, "modelo de turb. = ", m_turb
	write(*,*) "Condições De Contorno:"
	write(*,*) "	ccxo		ccy0	ccz0"
	write(*,*) ccx0, ccy0, ccz0
	write(*,*) "	ccxf		ccyf	cczf"
	write(*,*) ccxf, ccyf, cczf
	!write(*,*) "	X0		Y0			Z0"
	!write(*,*) ccx0, ccy0, ccz0
	!write(*,*) "	Xf		Yf			Zf"
	!write(*,*) ccxf, ccyf, cczf
	!write(*) "ccx0 = ", ccx0 ,"ccxf = ", ccxf, "ccy0 = ", ccy0, "ccyf = ", ccyf, "ccz0 = ", ccz0, "cczf = ", cczf

	!Plotagens de condição inicial

	!Perfil longitudinal do desnível

	ii = 0
	do i = nx, 2, -1 !Tem que ser o último a variar
		do k = 1, nz
			do j = 1, ny
				if (ii == 0) then
					if (ls(i-1,j,k)*ls(i,j,k) <= 0.) then
						write(99998799,*) dt*it, dx*(i-0.5) + ls(i,j,k)
						ii = 1 ! fará com que saida de todos os loops
					endif
				endif
			enddo
		enddo
	enddo
	!close (unit=99998799)

	!Perfil longitudinal do desnível
	open (unit=99998800, action= 'write', file= 'dados//h2.txt', status= 'unknown') ! o status deveria ser para escrever em cima deste primeiro ..z
	ii = 0
	i = int(1.0/dx) 
	j = int(0.5/dy) 
	do k = nz, 2, -1
		if (ii == 0) then
			if (ls(i,j,k-1)*ls(i,j,k) <= 0.) then
				write(99998800,*) dt*it, dz*(k-0.5) + ls(i,j,k)
				ii = 1 !Fará com que saida de todos os loops
			endif
		endif
	enddo

	open (unit=99998801, action= 'write', file= 'dados//h4.txt', status= 'unknown') ! o status deveria ser para escrever em cima deste primeiro ..z
	ii = 0
	i = int(1.0/dx) 
	j = int(0.5/dy) 
	do k = nz, 2, -1
		if (ii == 0) then
			if (ls(i,j,k-1)*ls(i,j,k) <= 0.) then
				write(99998801,*) dt*it, dz*(k-0.5) + ls(i,j,k)
				ii = 1 !Fará com que saida de todos os loops
			endif
		endif
	enddo

	!close (unit=99998800)
	!close (unit=99998801)

	!Plotagens

	!Conservação de massa
	open (unit=100002, action= 'write', file= 'dados//conservacao_massa.txt', status= 'unknown')
	write(100002,*) "it*dt", " ", "vol_ini", " ", "vol_ins", " ","vol_ini-vol_ins", " ","divergencia"

	!Número de courant
	open (unit=9999991, action= 'write', file= 'dados//courant.txt', status= 'unknown')
	write(9999991,*) "it*dt", "  ", "ntal", "  ", "tal", "  ", "maxval(a)", " ", "maxval(d)", "  ", "loca", "  ", "locd"

	!Contagem temporal
	open (unit=200000, file='dados//contagem')
	write(200000,*) "it", " ", "it*dt", " ", "hoje", " ", "agora", " ", "duração da simulação (min)"

	call date_and_time(values = agora)
	agora1 = agora
	call cpu_time(t_i)
	!Contagem temporal
	write(200000,*) it, it*dt, agora, t_i - t_i

	!Perfil longitudinal do desnível
	open (unit=9999999, action= 'write', file= 'dados//resul.txt', status= 'unknown')
	j = floor(real(ny)*0.5)

	do i = 1, nx
		write(9999999,*) (real(i)-0.5)*dx, "arrumar"
	enddo
	
	close (unit=9999999)

	!Plotagens de planos (escolher um por simulação, fazer alterações dependendo do plano a ser extraído)
	cont = 10
	div = 0.

	x1(0,:,:) = -dx
	y1(:,0,:) = -dy
	z1(:,:,0) = -dz

	do k = 1, nz1
		do j = 1, ny1
			do i = 1, nx1
				uaux(i,j,k) = (u(i,j,k) + u(i,j-1,k) + u(i,j,k-1) + u(i,j-1,k-1)) * 0.25
				vaux(i,j,k) = (v(i,j,k) + v(i-1,j,k) + v(i,j,k-1) + v(i-1,j,k-1)) * 0.25
				waux(i,j,k) = (w(i,j,k) + w(i,j-1,k) + w(i-1,j,k) + w(i-1,j-1,k)) * 0.25

				x1(i,j,k) = x1(i-1,j,k) + dx
				y1(i,j,k) = y1(i,j-1,k) + dy
				z1(i,j,k) = z1(i,j,k-1) + dz
				x11(i,j,k) = x1(i,j,k)
				y11(i,j,k) = y1(i,j,k)
				z11(i,j,k) = z1(i,j,k)
			enddo
		enddo
	enddo

	if (obst_t == 0) then
		kaux = 0
	else
		kaux = 0.
		do j = 1, ny
			do i = 1, nx
				do k = 1, kw(i,j)-1 !-1 pois o kw está posicionado na aresta inferior e não no centro da célula
					!print *, ku(i,j)
					kaux(i,j,k)=1.
				enddo
			enddo
		enddo
	endif

	prdaux(1:nx,1:ny,1:nz) = prd1(1:nx,1:ny,1:nz)

	if (mms_t > 0) then
		CALL interpy_cf(sqrt(erro_u),nx1,ny,nz,auxz)
		CALL interpz_cf(auxz,nx1,ny1,nz,uaux)

		CALL interpx_cf(sqrt(erro_v),nx,ny1,nz,auxz)
		CALL interpz_cf(auxz,nx1,ny1,nz,vaux)

		CALL interpx_cf(sqrt(erro_w),nx,ny,nz1,auxy)
		CALL interpy_cf(auxy,nx1,ny,nz1,waux)

		prdaux = sqrt(erro_p)
	endif	

	if (t_plot == 0) then
		open (unit=cont, file='arquivos//campos_00010',form='unformatted',status='unknown')
		write(cont) real(x11,4),real(y11,4),real(z11,4),real(uaux,4), & 
		real(vaux,4),real(waux,4),real(ls,4),real(kaux,4) !,real(prdaux,4)
	elseif (t_plot == 1) then
		do k = 1, nz
			do j = 1, ny
				do i = 1, nx
					dudy(i,j,k) = ((u(i,j+1,k)-u(i,j,k))/dy)
					dudz(i,j,k) = ((u(i,j,k+1)-u(i,j,k))/dz) 
					dvdx(i,j,k) = ((v(i+1,j,k)-v(i,j,k))/dx)
					dvdz(i,j,k) = ((v(i,j,k+1)-v(i,j,k))/dz)
					dwdx(i,j,k) = ((w(i+1,j,k)-w(i,j,k))/dx)
					dwdy(i,j,k) = ((w(i,j+1,k)-w(i,j,k))/dy) 

					vorti(i,j,k) = dwdy(i,j,k) - dvdz(i,j,k)
					vortj(i,j,k) = dudz(i,j,k) - dwdx(i,j,k)
					vortk(i,j,k) = dvdx(i,j,k) - dudy(i,j,k)
				enddo
			enddo
		enddo

		call interpx_cf(nut,nx,ny,nz,xnuta) !(nx1,ny,nz)
		call interpy_cf(nut,nx,ny,nz,ynuta) !(nx,ny1,nz)
		call interpz_cf(nut,nx,ny,nz,znuta) !(nx,ny1,nz)

		open (unit=cont, file='arquivos//campos_00010',form='unformatted',status='unknown')
		write(cont) real(x11,4),real(y11,4),real(z11,4),real(uaux,4), & 
		real(vaux,4),real(waux,4),real(ls,4),real(kaux,4),real(prdaux,4), &
		real(vorti,4),real(vortj,4),real(vortk,4),real(xnuta,4),real(ynuta,4),real(znuta,4)
	endif

	close (unit=cont)

	!Chama paraview
	CALL visu ()

END SUBROUTINE plot_i

SUBROUTINE plot_f()

	USE paodeplot_f

	IMPLICIT NONE

	integer :: i, j, k

	!Cálculo para o a estimativa do tempo restante
	ciclo = (agora(5)-agora1(5)) * 60 * 60 + (agora(6)-agora1(6)) * 60 + (agora(7)-agora1(7)) + real(agora(8)-agora1(8))/1000
	prev = (prev*6 + (ts-it)*ciclo*1./(60.*60.))/(7.)
	agora1 = agora
	call date_and_time(values = agora)
	call cpu_time(t_a)
	nfil=it

	if(mod(it, ceiling(dt_frame/dt)).eq.0) then

		cont = cont + 1

		x1(0,:,:) = -dx
		y1(:,0,:) = -dy
		z1(:,:,0) = -dz

		do k = 1, nz1
			do j = 1, ny1
				do i = 1, nx1
					uaux(i,j,k) = (u(i,j,k) + u(i,j-1,k) + u(i,j,k-1) + u(i,j-1,k-1)) * 0.25
					vaux(i,j,k) = (v(i,j,k) + v(i-1,j,k) + v(i,j,k-1) + v(i-1,j,k-1)) * 0.25
					waux(i,j,k) = (w(i,j,k) + w(i,j-1,k) + w(i-1,j,k) + w(i-1,j-1,k)) * 0.25

					x1(i,j,k) = x1(i-1,j,k) + dx
					y1(i,j,k) = y1(i,j-1,k) + dy
					z1(i,j,k) = z1(i,j,k-1) + dz
					x11(i,j,k) = x1(i,j,k)
					y11(i,j,k) = y1(i,j,k)
					z11(i,j,k) = z1(i,j,k)
				enddo
			enddo
		enddo
	
		if (obst_t == 0) then
			kaux = 0
		else
			kaux=0.
			do j = 1, ny
				do i = 1, nx
					do k = 1, kw(i,j)-1!-1 prq o kw está posicionado na aresta inferior e não no centro da célula
						!print *, ku(i,j)
						kaux(i,j,k)=1.
					enddo
				enddo
			enddo
		endif

		do ifile = 1, cont
			dig1 =    ifile/10000 + 48
			dig2 = ( ifile - 10000*( ifile/10000 ) )/1000 + 48
			dig3 = ( ifile - 1000*( ifile/1000 ) )/100 + 48
			dig4 = ( ifile - 100*( ifile/100 ) )/10 + 48
			dig5 = ( ifile - 10*( ifile/10 ) )/1 + 48
			chits(1:5) = char(dig1)//char(dig2)//char(dig3)//char(dig4)//char(dig5)
		enddo

		prdaux(1:nx,1:ny,1:nz) = prd1(1:nx,1:ny,1:nz)

		if (mms_t > 0) then
			CALL interpy_cf(sqrt(erro_u),nx1,ny,nz,auxz)
			CALL interpz_cf(auxz,nx1,ny1,nz,uaux)

			CALL interpx_cf(sqrt(erro_v),nx,ny1,nz,auxz)
			CALL interpz_cf(auxz,nx1,ny1,nz,vaux)

			CALL interpx_cf(sqrt(erro_w),nx,ny,nz1,auxy)
			CALL interpy_cf(auxy,nx1,ny,nz1,waux)

			prdaux = sqrt(erro_p)
		endif

		if (t_plot == 0) then
			open (unit=cont, action= 'write', file= 'arquivos//campos_'//chits,form='unformatted',status='unknown')
			write(cont) real(x11,4),real(y11,4),real(z11,4),real(uaux,4), & 
			real(vaux,4),real(waux,4),real(ls,4),real(kaux,4)!,real(prdaux,4)
		elseif (t_plot == 1) then
			do k = 1, nz
				do j = 1, ny
					do i = 1, nx
						dudy(i,j,k) = ((u(i,j+1,k)-u(i,j,k))/dy)
						dudz(i,j,k) = ((u(i,j,k+1)-u(i,j,k))/dz) 
						dvdx(i,j,k) = ((v(i+1,j,k)-v(i,j,k))/dx)
						dvdz(i,j,k) = ((v(i,j,k+1)-v(i,j,k))/dz)
						dwdx(i,j,k) = ((w(i+1,j,k)-w(i,j,k))/dx)
						dwdy(i,j,k) = ((w(i,j+1,k)-w(i,j,k))/dy) 

						vorti(i,j,k) = dwdy(i,j,k) - dvdz(i,j,k)
						vortj(i,j,k) = dudz(i,j,k) - dwdx(i,j,k)
						vortk(i,j,k) = dvdx(i,j,k) - dudy(i,j,k)
					enddo
				enddo
			enddo

			call interpx_cf(nut,nx,ny,nz,xnuta) !(nx1,ny,nz)
			call interpy_cf(nut,nx,ny,nz,ynuta) !(nx,ny1,nz)
			call interpz_cf(nut,nx,ny,nz,znuta) !(nx,ny1,nz)

			!*********************************************************************
			open (unit=cont, action= 'write', file= 'arquivos//campos_'//chits,form='unformatted',status='unknown')
			write(cont) real(x11,4),real(y11,4),real(z11,4),real(uaux,4), & 
			real(vaux,4),real(waux,4),real(ls,4),real(kaux,4),real(prdaux,4), &
			real(vorti,4),real(vortj,4),real(vortk,4),real(xnuta,4),real(ynuta,4),real(znuta,4)
		endif

		close (unit=cont)

		!Contagem temporal
		!write(*,*) "it,", " ", "it*dt,", " ", "ciclo,", " ","tempo restante aproximado (horas),", " ","duração da simulação (min)"
		!write(*,*) it, it*dt, ciclo, prev, (t_a-t_i)/60.
		write(*,*)
		write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
		write(*,*) "it:", it
		write(*,*) "it*dt:", it*dt
		write(*,*) "Ciclo:", ciclo
		write(*,*) "Tempo restante aproximado (h):", 	prev
		write(*,*) "Tempo restante aproximado (min):", 	prev/60
		write(*,*) "Duração Da Simulação (s):", 	(t_a-t_i)
		write(*,*) "Duração Da Simulação (min):", 	(t_a-t_i)/60
	endif

	write(200000,*) it, it*dt, agora, (t_a-t_i)/60.

	!open (unit=99998799, action= 'write', file= 'dados//frente.txt', status= 'unknown')
	ii = 0

	do i = nx, 2, -1!Tem que ser o último a variar
		do k = 1, nz
			do j = 1, ny
				if (ii == 0) then
					if (ls(i-1,j,k)*ls(i,j,k) <= 0.) then
						write(99998799,*) dt*it, dx*(i-0.5) + ls(i,j,k)
						ii = 1 !Fará saída de todos os loops
					endif
				endif
			enddo
		enddo
	enddo

	ii = 0
	i = int(1.0/dx) 
	j = int(0.5/dy) 

	do k = nz, 2, -1
		if (ii == 0) then
			if (ls(i,j,k-1)*ls(i,j,k) <= 0.) then
				write(99998800,*) dt*it, dz*(k-0.5) + ls(i,j,k)
				ii = 1 !Fará saída de todos os loops
			endif
		endif
	enddo

	ii = 0
	i = int(2.66/dx) 
	j = int(0.5/dy) 

	do k = nz, 2, -1
		if (ii == 0) then
			if (ls(i,j,k-1)*ls(i,j,k) <= 0.) then
				write(99998801,*) dt*it, dz*(k-0.5) + ls(i,j,k)
				ii = 1 !Fará saída de todos os loops
			endif
		endif
	enddo

	!perfil velocidade horizontal seção 1 *************************************
	!	open (unit=100005, action= 'write', file= 'dados//perfil1.txt', status= 'unknown')

	!	i = int(0.8/dx)
	!	do k =1, nz
	!		write(100005,*) u(i,j,k), z1(i,j,k)-dz
	!	enddo	
	!	
	!	 close (unit=100005)
	!*********************************************************************

	!perfil velocidade horizontal seção 3 *************************************
	!	open (unit=100006, action= 'write', file= 'dados//perfil3.txt', status= 'unknown')

	!	i = int(1.2/dx)
	!	do k =1, nz
	!		write(100006,*) u(i,j,k), z1(i,j,k)-dz
	!	enddo	

	!	 close (unit=100006)
	!*********************************************************************

	!perfil velocidade horizontal seção 4 *************************************
	!	open (unit=100009, action= 'write', file= 'dados//perfil4.txt', status= 'unknown')

	!	i = int(1.4/dx)
	!	do k =1, nz
	!		write(100009,*) u(i,j,k), z1(i,j,k)-dz
	!	enddo	

	!	 close (unit=100009)
	!*********************************************************************

	!perfil velocidade horizontal seção 5 *************************************
	!	open (unit=100007, action= 'write', file= 'dados//perfil5.txt', status= 'unknown')

	!	i = int(1.9/dx)
	!	do k =1, nz
	!		write(100007,*) u(i,j,k), z1(i,j,k)-dz
	!	enddo	

	!	 close (unit=100007)
	!*********************************************************************

	!perfil velocidade horizontal seção 6 *************************************
	!	open (unit=100011, action= 'write', file= 'dados//perfil6.txt', status= 'unknown')

	!	i = int(2.4/dx)
	!	do k =1, nz
	!		write(100011,*) u(i,j,k), z1(i,j,k)-dz
	!	enddo	

	!	 close (unit=100011)
	!*********************************************************************

	!perfil velocidade horizontal seção 7 *************************************
	!	open (unit=100008, action= 'write', file= 'dados//perfil7.txt', status= 'unknown')

	!	i = int(2.6/dx)
	!	do k =1, nz
	!		write(100008,*) u(i,j,k), z1(i,j,k)-dz
	!	enddo	
	!	
	!	 close (unit=100008)
	!*********************************************************************

	!perfil velocidade horizontal seção 9 *************************************
	!	open (unit=100012, action= 'write', file= 'dados//perfil8.txt', status= 'unknown')

	!	i = int(2.8/dx)
	!	do k =1, nz
	!		write(100012,*) u(i,j,k), z1(i,j,k)-dz
	!	enddo	

	!	 close (unit=100012)

	CALL est(div)																		! ### PEDRO ###

END SUBROUTINE plot_f

SUBROUTINE plot_atrib()

	USE velpre
	USE obst

	open(unit=9, action= 'write', file= 'dados//atributos.txt', status= 'unknown')
	write(9,*) "u inicial = ", uinicial
	write(9,*) "dx = ", dx, ", dy = ", dy, ", dz = ", dz
	write(9,*) "nx = ", nx, "ny = ", ny, "nz = ", nz
	write(9,*) "dt = ", dt, "ts = ", ts, "duração da sim =", (t_a-t_i)/60.
	write(9,*) "amp = ", amp, "comp = ", comp
	write(9,*) "elev = ", elev
	close (unit=9)

	open(unit=9, action= 'write', file= 'dados//analise.txt', status= 'unknown') !ALTERAÇÃO BY PEDRO
	write(9,*) "Duração da simulação =", (t_a-t_i)/60
	!write(9,*) "Wtime =", suli_wtime
	close(unit=9)

END SUBROUTINE plot_atrib

SUBROUTINE est(div)

	USE disc, only: dt, dt_frame, it
	USE velpre
	USE ls_param

	IMPLICIT NONE
	real(8), dimension(nx,ny,nz) :: div
	integer :: i,j,k

	do k = 1, nz
		do j = 1, ny
			do i = 1, nx
				div(i,j,k) = (u(i+1,j,k) - u(i,j,k))/dx +   (v(i,j+1,k) - v(i,j,k))/dy +  (w(i,j,k+1) - w(i,j,k))/dz
			enddo
		enddo
	enddo

	!open (unit=100002, action= 'write', file= 'conservacao_massa.txt', status= 'unknown')
	write(100002,*) it*dt, vol_ini, vol_ins, vol_ini-vol_ins, maxval(abs(div))

	if(mod(it, ceiling(dt_frame/dt)).eq.0) then
		write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
		!write(*,*) "*vol. inicial", vol_ini, "vol. instant.", vol_ins, "erro (%)", (vol_ini-vol_ins)/vol_ini, "div.", maxval(abs(div))	
		write(*,*) "Vol. Inicial", vol_ini
		write(*,*) "Vol. Instant.", vol_ins
		write(*,*) "Erro (%)", (vol_ini-vol_ins)/vol_ini
		write(*,*) "Div.", maxval(abs(div))
	endif

END SUBROUTINE est
