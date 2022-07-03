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
	USE omp_lib
	USE disc
	USE restart

	IMPLICIT NONE

	!Declaração das variáveis para análise
	double precision :: fortran_start, fortran_end
	double precision :: openmp_start, openmp_end

	double precision :: fortran_start_level_set, fortran_end_level_set
	double precision :: openmp_start_level_set, openmp_end_level_set

	double precision :: fortran_start_visco, fortran_end_visco
	double precision :: omp_start_visco, omp_end_visco

	double precision :: fortran_start_convdiff, fortran_end_convdiff
	double precision :: omp_start_convdiff, omp_end_convdiff

	double precision :: fortran_start_graddin, fortran_end_graddin
	double precision :: omp_start_graddin, omp_end_graddin

	double precision :: fortran_start_plot_f, fortran_end_plot_f
	double precision :: omp_start_plot_f, omp_end_plot_f

	double precision :: fortran_start_plot_atri, fortran_end_plot_atri
	double precision :: omp_start_plot_atri, omp_end_plot_atri

	double precision :: start_outros_f90, end_outros_f90, start_outros2_f90, end_outros2_f90, start_outros3_f90, end_outros3_f90
	double precision :: start_outros_omp, end_outros_omp, start_outros2_omp, end_outros2_omp, start_outros3_omp, end_outros3_omp

	real :: soma_level_set_f90, soma_visco_f90, soma_convdiff_f90, soma_graddin_f90, soma_plot_f_f90
	real :: soma_level_set_omp, soma_visco_omp, soma_convdiff_omp, soma_graddin_omp, soma_plot_f_omp
	real :: soma_outros_f90, soma_outros2_f90, soma_outros3_f90
	real :: soma_outros_omp, soma_outros2_omp, soma_outros3_omp

	!Inicialização dos somatórios com valor nulo
	soma_level_set_f90 = 0.0
	soma_level_set_omp = 0.0
	soma_visco_f90 = 0.0
	soma_visco_omp = 0.0
	soma_convdiff_f90 = 0.0
	soma_convdiff_omp = 0.0
	soma_graddin_f90 = 0.0
	soma_graddin_omp = 0.0
	soma_outros_f90 = 0.0
	soma_outros_omp = 0.0
	soma_outros2_f90 = 0.0
	soma_outros2_omp = 0.0
	soma_outros3_f90 = 0.0
	soma_outros3_omp = 0.0

	!double precision function omp_get_wtime()
	!double precision function omp_get_wtick()

	!Inicio do Fortran e OpenMP
	write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~- INÍCIO ~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-"
	CALL cpu_time(fortran_start)
	openmp_start = omp_get_wtime()

	if (nx*ny*nz > 30000000) then
		write(*,*) "Verifique se o seu computador tem capacidade para a simulação, se sim, cancele esta condicional no código."
		STOP
	endif

	!Condições iniciais
	if (irest.eq.0) then
		CALL iniciais()
	else
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
			openmp_start_level_set = omp_get_wtime()

			CALL level_set()

			CALL cpu_time(fortran_end_level_set)
			openmp_end_level_set = omp_get_wtime()

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
			start_outros3_omp = omp_get_wtime()

			CALL contorno(3)

			CALL cpu_time(end_outros3_f90)
			end_outros3_omp = omp_get_wtime()

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
				CALL cpu_time(fortran_start_visco)
				omp_start_visco = omp_get_wtime()

				CALL visco()

				CALL cpu_time(fortran_end_visco)
				omp_end_visco = omp_get_wtime()

				write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
				soma_visco_f90 = soma_visco_f90 + (fortran_end_visco-fortran_start_visco)
				soma_visco_omp = soma_visco_omp + (omp_end_visco-omp_start_visco)
				write(*,*) "Tempo individual do visco() p/ Fortran:", fortran_end_visco-fortran_start_visco
				write(*,*) "Tempo acumulado do visco() p/ Fortran", soma_visco_f90
				write(*,*) "Tempo individual do visco() p/ OpenMP:", omp_end_visco-omp_start_visco
				write(*,*) "Tempo acumulado do visco() p/ OpenMP", soma_visco_omp

				!Tempo do convdiff() p/ Fortran e OpenMP
				CALL cpu_time(fortran_start_convdiff)
				omp_start_convdiff = omp_get_wtime()

				CALL convdiff()
				CALL tempo()

				CALL cpu_time(fortran_end_convdiff)
				omp_end_convdiff = omp_get_wtime()

				write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
				soma_convdiff_f90 = soma_convdiff_f90 + (fortran_end_convdiff-fortran_start_convdiff)
				soma_convdiff_omp = soma_convdiff_omp + (omp_end_convdiff-omp_start_convdiff)
				write(*,*) "Tempo individual do convdiff() p/ Fortran:", fortran_end_convdiff-fortran_start_convdiff
				write(*,*) "Tempo acumulado do convdiff() p/ Fortran", soma_convdiff_f90
				write(*,*) "Tempo individual do convdiff() p/ OpenMP:", omp_end_convdiff-omp_start_convdiff
				write(*,*) "Tempo acumulado do convdiff() p/ OpenMP", soma_convdiff_omp

				if (wave_t > 0) call boundary_waves() !For wave propagation
					!Condições de Contorno para a parte Hidrostática
					!CALL pressh()

					!Tempo do contorno(2) p/ Fortran e OpenMP
					CALL cpu_time(start_outros2_f90)
					start_outros2_omp = omp_get_wtime()

					CALL contorno(2)

					CALL cpu_time(end_outros2_f90)
					end_outros2_omp = omp_get_wtime()

					write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
					soma_outros2_f90 = soma_outros2_f90 + (end_outros2_f90 - start_outros2_f90)
					soma_outros2_omp = soma_outros2_omp + (end_outros2_omp - start_outros2_omp)
					write(*,*) "Tempo individual do contorno(2) dentro do ciclo p/ Fortran:", end_outros2_f90 - start_outros2_f90
					write(*,*) "Tempo acumulado do contorno(2) dentro do ciclo p/ Fortran", soma_outros2_f90
					write(*,*) "Tempo individual do contorno(2) dentro do ciclo p/ OpenMP:", end_outros2_omp - start_outros2_omp
					write(*,*) "Tempo acumulado do contorno(2) dentro do ciclo p/ OpenMP", soma_outros2_omp
				

				if (mms_t .eq. 0) then

					!Tempo do graddin() p/ Fortran e OpenMP
					CALL cpu_time(fortran_start_graddin)
					omp_start_graddin = omp_get_wtime()

					CALL graddin()
					CALL posdin()

					CALL cpu_time(fortran_end_graddin)
					omp_end_graddin = omp_get_wtime()

					write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
					soma_graddin_f90 = soma_graddin_f90 + (fortran_end_graddin-fortran_start_graddin)
					soma_graddin_omp = soma_graddin_omp + (omp_end_graddin-omp_start_graddin)
					write(*,*) "Tempo individual do graddin() p/ Fortran:", fortran_end_graddin-fortran_start_graddin
					write(*,*) "Tempo acumulado do graddin() p/ Fortran", soma_graddin_f90
					write(*,*) "Tempo individual do graddin() p/ OpenMP:", omp_end_graddin-omp_start_graddin
					write(*,*) "Tempo acumulado do graddin() p/ OpenMP", soma_graddin_omp
				
					!Tempo do contorno(1) p/ Fortran e OpenMP
					CALL cpu_time(start_outros_f90)
					start_outros_omp = omp_get_wtime()

					CALL contorno(1)

					CALL cpu_time(end_outros_f90)
					end_outros_omp = omp_get_wtime()

					write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
					soma_outros_f90 = soma_outros_f90 + (end_outros_f90 - start_outros_f90)
					soma_outros_omp = soma_outros_omp + (end_outros_omp - start_outros_omp)
					write(*,*) "Tempo individual do contorno(1) dentro do ciclo p/ Fortran:", end_outros_f90 - start_outros_f90
					write(*,*) "Tempo acumulado do contorno(1) dentro do ciclo p/ Fortran", soma_outros_f90
					write(*,*) "Tempo individual do contorno(1) dentro do ciclo p/ OpenMP:", end_outros_omp - start_outros_omp
					write(*,*) "Tempo acumulado do contorno(1) dentro do ciclo p/ OpenMP", soma_outros_omp
				endif		! De onde é esse endif?		### PEDRO ###
			enddo		! De onde é esse enddo?		### PEDRO ###

		!Solução manufaturada; cálculo do erro

		if (mms_t > 0) CALL mms()
			!Plotagens por passo de tempo
			!Tempo do plot_f() p/ Fortran e OpenMP
			CALL cpu_time(fortran_start_plot_f)
			omp_start_plot_f = omp_get_wtime()

			CALL plot_f()

			CALL cpu_time(fortran_end_plot_f)
			omp_end_plot_f = omp_get_wtime()
			write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
			soma_plot_f_f90 = soma_plot_f_f90 + (fortran_end_plot_f-fortran_start_plot_f)
			soma_plot_f_omp = soma_plot_f_omp + (omp_end_plot_f-omp_start_plot_f)
			write(*,*) "Tempo individual do plot_f() p/ Fortran:", fortran_end_plot_f-fortran_start_plot_f
			write(*,*) "Tempo acumulado do plot_f() p/ Fortran", soma_plot_f_f90
			write(*,*) "Tempo individual do plot_f() p/ OpenMP:", omp_end_plot_f-omp_start_plot_f
			write(*,*) "Tempo acumulado do plot_f() p/ OpenMP", soma_plot_f_omp

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
	openmp_end = omp_get_wtime()

	!Tempo do Fortran e OpenMP
	write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
	write(*,*) "Tempo Total do Fortran:", fortran_end - fortran_start
	write(*,*) "Tempo Total do OpenMP:", openmp_end-openmp_start
	write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~- FIM ~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-"

END PROGRAM PNH