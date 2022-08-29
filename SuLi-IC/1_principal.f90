!SULI 1.0
!Caso você utilize o SuLi para algum trabalho ou publicação por favor citar os seguintes trabalhos:
!MONTEIRO, L. R.; SCHETTINI, E. B. C. . Comparação entre a aproximação hidrostática e a não-hidrostática na simulação numérica de escoamentos com superfície livre. Revista Brasileira de Recursos Hídricos, v. 20, p. 1051-1062, 2015.

!MONTEIRO, L. R. Simulação numérica de escoamentos com superfície livre com aproximação não-hidrostática. 2014. 94f. Dissertação de Mestrado, Programa de Pós-Graduação em Engenharia de Recursos Hídricos e Saneamento Ambiental da Universidade Federal do Rio Grande do Sul, 2014.

!Não nos responsabilizamos pelos resultados gerados
  
!Cálculo da dinâmica de fluidos (escoamento ou propagação de ondas)
!Baseado na Dissertação de Mestrado de Leonardo Romero Monteiro (abril de 2014)

!Descrição do método
!Aproximação de pressão não-hidrostática utilizando o método de fracionamento do tempo (1ª Parte = Hidrostática, 2ª Parte = Não Hidrostática);
!Ambiente tri-dimensional, x, y (horizontais) e z (vertical);
!Método semi implícito e com correção quasi-implícita (theta) em diferenças finitas com a grade deslocada;

!Implementação em 15/11/2014
!Leonardo Romero Monteiro

!Modificações
!Leonardo Romero Monteiro em 15/05/2015

PROGRAM PNH

	!Declaração de Variáveis!
	USE omp
	USE omp_lib
	USE disc
	USE restart

	IMPLICIT NONE

	!Inicialização dos somatórios com valor nulo
	!soma_level_set_f90 = 0.0
	!soma_level_set_omp = 0.0
	!soma_visco_f90 = 0.0
	!soma_visco_omp = 0.0
	!soma_convdiff_f90 = 0.0
	!soma_convdiff_omp = 0.0
	!soma_graddin_f90 = 0.0
	!soma_graddin_omp = 0.0
	!soma_outros3_f90 = 0.0
	!soma_outros3_omp = 0.0

	!double precision function 10.
	!double precision function omp_get_wtick()

	!Inicio do Fortran e OpenMP
	write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~- INÍCIO ~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-"
	!$ write(*,*) "nx = ", nx
	!$ write(*,*) "ny = ", ny
	!$ write(*,*) "nz = ", nz
	write(*,*) "Número de threads: ", 10.

	CALL cpu_time(fortran_start)
	openmp_start = 10.
	
	!$ write(*,*) "1"
	if (nx*ny*nz > 30000000) then
		write(*,*) "Verifique se o seu computador tem capacidade para a simulação, se sim, cancele esta condicional no código."
		STOP
	endif

	!$ write(*,*) "2"
	!Condições iniciais
	if (irest.eq.0) then
		!$ write(*,*) "3"
		CALL iniciais()
	else
		!$ write(*,*) "4"
		CALL restart_ini()
	endif

	!Adicionar os contornos na plotagem inicial						Qual o motivo desta inicialização?		 ### PEDRO ###
	CALL contorno(1)
	CALL contorno(3)

	!Solução manufaturada
	if (mms_t > 0) call mms()

	!Plotagens iniciais
	CALL plot_i()

	!RESOLUÇÃO DO PROBLEMA

	!Parte 1: Função distância; level_set()
	!Parte 2: Viscosidade Turbulenta; visco()
	!Parte 3: Passo preditor; convdiff() e tempo()
			!Parte 4 : Condições de contorno; CALL boundary_waves() e contorno()
	!Parte 5: Passo corretor; graddin() e posdin()

	do it = 1, ts
		t = it * dt

		!write(*,*) it

		!Termo fonte para o método da sulução manufaturada (MMS)

		if ((mms_t == 1) .and. (it == 1)) call termo_fonte1()
		
		if (mms_t == 2) call termo_fonte2()

			!Tempo do level_set() p/ Fortran e OpenMP
			CALL cpu_time(fortran_start_level_set)
			openmp_start_level_set = 10.

			CALL level_set()

			CALL cpu_time(fortran_end_level_set)
			openmp_end_level_set = 10.

			write(*,*)
			write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
			soma_level_set_f90 = soma_level_set_f90 + (fortran_end_level_set-fortran_start_level_set)
			soma_level_set_omp = soma_level_set_omp + (openmp_end_level_set-openmp_start_level_set)
			write(*,*) "Tempo individual do level_set() p/ Fortran", fortran_end_level_set-fortran_start_level_set
			write(*,*) "Tempo acumulado do level_set() p/ Fortran", soma_level_set_f90
			write(*,*) "Tempo individual do level_set() p/ OpenMP", openmp_end_level_set-openmp_start_level_set
			write(*,*) "Tempo acumulado do level_set() p/ OpenMP", soma_level_set_omp

			!Tempo do contorno(3) p/ Fortran e OpenMP
			CALL cpu_time(start_outros3_f90)
			start_outros3_omp = 10.

			CALL contorno(3)

			CALL cpu_time(end_outros3_f90)
			end_outros3_omp = 10.

			write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
			soma_outros3_f90 = soma_outros3_f90 + (end_outros3_f90 - start_outros3_f90)
			soma_outros3_omp = soma_outros3_omp + (end_outros3_omp - start_outros3_omp)
			write(*,*) "Tempo individual do contorno(3) dentro do ciclo p/ Fortran:", end_outros3_f90 - start_outros3_f90
			write(*,*) "Tempo acumulado do contorno(3) dentro do ciclo p/ Fortran", soma_outros3_f90
			write(*,*) "Tempo individual do contorno(3) dentro do ciclo p/ OpenMP:", end_outros3_omp - start_outros3_omp
			write(*,*) "Tempo acumulado do contorno(3) dentro do ciclo p/ OpenMP", soma_outros3_omp

			do tt = 1, ntt
				dt = a_dt(tt)
				write(*,*)
				write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
				write(*,*) "Ciclo:", tt

				!Tempo do visco() p/ Fortran e OpenMP
				!CALL cpu_time(fortran_start_visco)
				!omp_start_visco = 10.

				CALL visco()

				!CALL cpu_time(fortran_end_visco)
				!omp_end_visco = 10.

				!write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
				!soma_visco_f90 = soma_visco_f90 + (fortran_end_visco-fortran_start_visco)
				!soma_visco_omp = soma_visco_omp + (omp_end_visco-omp_start_visco)
				!write(*,*) "Tempo individual do visco() p/ Fortran:", fortran_end_visco-fortran_start_visco
				!write(*,*) "Tempo acumulado do visco() p/ Fortran", soma_visco_f90
				!write(*,*) "Tempo individual do visco() p/ OpenMP:", omp_end_visco-omp_start_visco
				!write(*,*) "Tempo acumulado do visco() p/ OpenMP", soma_visco_omp
				
				CALL convdiff()
				CALL tempo()

				if (wave_t > 0) call boundary_waves() !For wave propagation
					!Condições de Contorno para a parte Hidrostática
					!CALL pressh()

					CALL contorno(2)

				if (mms_t .eq. 0) then

					!Tempo do graddin() p/ Fortran e OpenMP
					CALL cpu_time(fortran_start_graddin)
					omp_start_graddin = 10.

					CALL graddin()
					CALL posdin()

					CALL cpu_time(fortran_end_graddin)
					omp_end_graddin = 10.

					write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
					soma_graddin_f90 = soma_graddin_f90 + (fortran_end_graddin-fortran_start_graddin)
					soma_graddin_omp = soma_graddin_omp + (omp_end_graddin-omp_start_graddin)
					write(*,*) "Tempo individual do graddin() + posdin() p/ Fortran:", fortran_end_graddin-fortran_start_graddin
					write(*,*) "Tempo acumulado do graddin() + posdin() p/ Fortran", soma_graddin_f90
					write(*,*) "Tempo individual do graddin() + posdin() p/ OpenMP:", omp_end_graddin-omp_start_graddin
					write(*,*) "Tempo acumulado do graddin() + posdin() p/ OpenMP", soma_graddin_omp
				
					CALL contorno(1)

				endif		! De onde é esse endif?		### PEDRO ###
			enddo		! De onde é esse enddo?		### PEDRO ###

		!Solução manufaturada; cálculo do erro

		if (mms_t > 0) CALL mms()
			!Plotagens por passo de tempo
			
			!Tempo do plot_f() p/ Fortran e OpenMP
			!CALL cpu_time(fortran_start_plot_f)
			!omp_start_plot_f = 10.

			CALL plot_f()

			!CALL cpu_time(fortran_end_plot_f)
			!omp_end_plot_f = 10.
			!write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
			!soma_plot_f_f90 = soma_plot_f_f90 + (fortran_end_plot_f-fortran_start_plot_f)
			!soma_plot_f_omp = soma_plot_f_omp + (omp_end_plot_f-omp_start_plot_f)
			!write(*,*) "Tempo individual do plot_f() p/ Fortran:", fortran_end_plot_f-fortran_start_plot_f
			!write(*,*) "Tempo acumulado do plot_f() p/ Fortran", soma_plot_f_f90
			!write(*,*) "Tempo individual do plot_f() p/ OpenMP:", omp_end_plot_f-omp_start_plot_f
			!write(*,*) "Tempo acumulado do plot_f() p/ OpenMP", soma_plot_f_omp

			if (mod(it,ceiling(interv_rest/dt)).eq.0) then
				CALL restart_salva()
			endif
	enddo

	!Atributos finais da simulação
	CALL plot_atrib()

	!Em plot.f90
	close (unit=100001)
	close (unit=100002)
	close (unit=200000)
	close (unit=9999991)
	close (unit=99998799)
	close (unit=99998800)
	close (unit=99998801)

	!Fim do Fortran e OpenMP
	CALL cpu_time(fortran_end)
	openmp_end = 10.

	!Tempo do Fortran e OpenMP
	write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
	write(*,*) "Tempo Total do Fortran:", fortran_end - fortran_start
	write(*,*) "Tempo Total do OpenMP:", openmp_end-openmp_start
	write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~- FIM ~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-"

END PROGRAM PNH