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
	!soma_contorno3_f90 = 0.0
	!soma_contorno3_omp = 0.0

	!double precision function omp_get_wtime()
	!double precision function omp_get_wtick()

	!Inicio do Fortran e OpenMP
	write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~- INÍCIO ~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-"

	!$ CALL OMP_set_num_threads(6)		

	!$ CALL OMP_set_nested(.TRUE.)

	!# write(*,*) "Dynamic thread adjustment? ", OMP_GET_DYNAMIC()
	!$ write(*,*) "Nested? ", OMP_GET_NESTED()
	!$ write(*,*) "Número de threads início = ", OMP_GET_NUM_THREADS()
	!$ write(*,*) "Paralelo no início? ", OMP_in_parallel()

	write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~- INÍCIO ~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-"

	!$ CALL OMP_set_dynamic(.FALSE.)

	!$OMP PARALLEL

	!$ write(*,*) "Dynamic thread adjustment? ", OMP_GET_DYNAMIC()
	!$ write(*,*) "Nested? ", OMP_GET_NESTED()
	!$ write(*,*) "Número de threads início = ", OMP_GET_NUM_THREADS()
	!$ write(*,*) "Paralelo no início? ", OMP_in_parallel()
	write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~- INÍCIO ~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-"

	!$OMP END PARALLEL

	CALL cpu_time(fortran_start)
	!$ openmp_start = omp_get_wtime()
	
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

    CALL cpu_time(fortran_start_plot_i)
    !$ omp_start_plot_i = omp_get_wtime()

	CALL contorno(1)
	CALL contorno(3)

	!Solução manufaturada
	if (mms_t > 0) then
		write (*,*) "mms: ", mms_t
		call mms()
	else
		write (*,*) "mms: zero"
	endif

	!Plotagens iniciais

	CALL plot_i()

    CALL cpu_time(fortran_end_plot_i)
    !$ omp_end_plot_i = omp_get_wtime()

    !write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
    soma_plot_i_f90 = soma_plot_i_f90 + (fortran_end_plot_i-fortran_start_plot_i)
    soma_plot_i_omp = soma_plot_i_omp + (omp_end_plot_i-omp_start_plot_i)
    !write(*,*) "Tempo acumulado do plot_i() p/ Fortran", soma_plot_i_f90
    !write(*,*) "Tempo acumulado do plot_i() p/ OpenMP", soma_plot_i_omp

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

		if ((mms_t .eq. 1) .and. (it .eq. 1)) call termo_fonte1()
		
		if (mms_t .eq. 2) call termo_fonte2()

		!Tempo do level_set() p/ Fortran e OpenMP
		CALL cpu_time(fortran_start_level_set)
		!$ openmp_start_level_set = omp_get_wtime()

		CALL level_set()

		CALL cpu_time(fortran_end_level_set)
		!$ openmp_end_level_set = omp_get_wtime()

		!write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
		soma_level_set_f90 = soma_level_set_f90 + (fortran_end_level_set-fortran_start_level_set)
		soma_level_set_omp = soma_level_set_omp + (openmp_end_level_set-openmp_start_level_set)
		!write(*,*) "Tempo individual do level_set() p/ Fortran", fortran_end_level_set-fortran_start_level_set
		!write(*,*) "Tempo acumulado do level_set() p/ Fortran", soma_level_set_f90
		!write(*,*) "Tempo individual do level_set() p/ OpenMP", openmp_end_level_set-openmp_start_level_set
		!write(*,*) "Tempo acumulado do level_set() p/ OpenMP", soma_level_set_omp

		!Tempo do contorno(3) p/ Fortran e OpenMP
		CALL cpu_time(start_contorno3_f90)
		!$ start_contorno3_omp = omp_get_wtime()

		CALL contorno(3)

		CALL cpu_time(end_contorno3_f90)
		!$ end_contorno3_omp = omp_get_wtime()

		!write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
		soma_contorno3_f90 = soma_contorno3_f90 + (end_contorno3_f90 - start_contorno3_f90)
		soma_contorno3_omp = soma_contorno3_omp + (end_contorno3_omp - start_contorno3_omp)
		!write(*,*) "Tempo individual do contorno(3) dentro do ciclo p/ Fortran:", end_contorno3_f90 - start_contorno3_f90
		!write(*,*) "Tempo acumulado do contorno(3) dentro do ciclo p/ Fortran", soma_contorno3_f90
		!write(*,*) "Tempo individual do contorno(3) dentro do ciclo p/ OpenMP:", end_contorno3_omp - start_contorno3_omp
		!write(*,*) "Tempo acumulado do contorno(3) dentro do ciclo p/ OpenMP", soma_contorno3_omp

		do tt = 1, ntt
			dt = a_dt(tt)
			write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
			write(*,*) "Ciclo:", tt

			!Tempo do visco() p/ Fortran e OpenMP
			CALL cpu_time(fortran_start_visco)
			!$ omp_start_visco = omp_get_wtime()

			CALL visco()

			CALL cpu_time(fortran_end_visco)
			!$ omp_end_visco = omp_get_wtime()

			!write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
			soma_visco_f90 = soma_visco_f90 + (fortran_end_visco-fortran_start_visco)
			soma_visco_omp = soma_visco_omp + (omp_end_visco-omp_start_visco)
			!write(*,*) "Tempo individual do visco() p/ Fortran:", fortran_end_visco-fortran_start_visco
			!write(*,*) "Tempo acumulado do visco() p/ Fortran", soma_visco_f90
			!write(*,*) "Tempo individual do visco() p/ OpenMP:", omp_end_visco-omp_start_visco
			!write(*,*) "Tempo acumulado do visco() p/ OpenMP", soma_visco_omp

			CALL cpu_time(fortran_start_convdiff)
			!$ omp_start_convdiff = omp_get_wtime()

			CALL convdiff()

			CALL cpu_time(fortran_end_convdiff)
			!$ omp_end_convdiff = omp_get_wtime()

			!write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
			soma_convdiff_f90 = soma_convdiff_f90 + (fortran_end_convdiff-fortran_start_convdiff)
			soma_convdiff_omp = soma_convdiff_omp + (omp_end_convdiff-omp_start_convdiff)
			!write(*,*) "Tempo acumulado do convdiff() p/ Fortran", soma_convdiff_f90
			!write(*,*) "Tempo acumulado do convdiff() p/ OpenMP", soma_convdiff_omp				



			CALL cpu_time(fortran_start_tempo)
			!$ omp_start_tempo = omp_get_wtime()

			CALL tempo()

			CALL cpu_time(fortran_end_tempo)
			!$ omp_end_tempo = omp_get_wtime()

			!write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
			soma_tempo_f90 = soma_tempo_f90 + (fortran_end_tempo - fortran_start_tempo)
			soma_tempo_omp = soma_tempo_omp + (omp_end_tempo - omp_start_tempo)
			!write(*,*) "Tempo acumulado do tempo() p/ Fortran", soma_tempo_f90
			!write(*,*) "Tempo acumulado do tempo() p/ OpenMP", soma_tempo_omp



			if (wave_t > 0) call boundary_waves() !For wave propagation
				!Condições de Contorno para a parte Hidrostática
				!CALL pressh()


			CALL cpu_time(fortran_start_contorno2)
			!$ omp_start_contorno2 = omp_get_wtime()

			CALL contorno(2)

			CALL cpu_time(fortran_end_contorno2)
			!$ omp_end_contorno2 = omp_get_wtime()

			!write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
			soma_contorno2_f90 = soma_contorno2_f90 + (fortran_end_contorno2-fortran_start_contorno2)
			soma_contorno2_omp = soma_contorno2_omp + (omp_end_contorno2-omp_start_contorno2)
			!write(*,*) "Tempo acumulado do contorno(2) p/ Fortran", soma_contorno2_f90
			!write(*,*) "Tempo acumulado do contorno(2) p/ OpenMP", soma_contorno2_omp
			
			
			if (mms_t .eq. 0) then

				!Tempo do graddin() p/ Fortran e OpenMP
				CALL cpu_time(fortran_start_graddin)
				!$ omp_start_graddin = omp_get_wtime()
	
				CALL graddin()
	
				CALL cpu_time(fortran_end_graddin)
				!$ omp_end_graddin = omp_get_wtime()

				!write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
				soma_graddin_f90 = soma_graddin_f90 + (fortran_end_graddin-fortran_start_graddin)
				soma_graddin_omp = soma_graddin_omp + (omp_end_graddin-omp_start_graddin)
				!write(*,*) "Tempo individual do graddin() + posdin() p/ Fortran:", fortran_end_graddin-fortran_start_graddin
				!write(*,*) "Tempo acumulado do graddin() p/ Fortran", soma_graddin_f90
				!write(*,*) "Tempo individual do graddin() + posdin() p/ OpenMP:", omp_end_graddin-omp_start_graddin
				!write(*,*) "Tempo acumulado do graddin() p/ OpenMP", soma_graddin_omp
			


			
				CALL cpu_time(fortran_start_posdin)
				!$ omp_start_posdin = omp_get_wtime()

				CALL posdin()

				CALL cpu_time(fortran_end_posdin)
				!$ omp_end_posdin = omp_get_wtime()

				!write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
				soma_posdin_f90 = soma_posdin_f90 + (fortran_end_posdin - fortran_start_posdin)
				soma_posdin_omp = soma_posdin_omp + (omp_end_posdin - omp_start_posdin)
				!write(*,*) "Tempo acumulado do posdin() p/ Fortran", soma_posdin_f90
				!write(*,*) "Tempo acumulado do posdin() p/ OpenMP", soma_posdin_omp			





				CALL cpu_time(fortran_start_contorno1)
				!$ omp_start_contorno1 = omp_get_wtime()

				CALL contorno(1)

				CALL cpu_time(fortran_end_contorno1)
				!$ omp_end_contorno1 = omp_get_wtime()

				!write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
				soma_contorno1_f90 = soma_contorno1_f90 + (fortran_end_contorno1 - fortran_start_contorno1)
				soma_contorno1_omp = soma_contorno1_omp + (omp_end_contorno1 - omp_start_contorno1)
				!write(*,*) "Tempo acumulado do contorno(1) p/ Fortran", soma_contorno1_f90
				!write(*,*) "Tempo acumulado do contorno(1) p/ OpenMP", soma_contorno1_omp
				

			endif		! De onde é esse endif?		### PEDRO ###
		enddo		! De onde é esse enddo?		### PEDRO ###

		!Solução manufaturada; cálculo do erro

		if (mms_t > 0) CALL mms()
		!Plotagens por passo de tempo
		
		!Tempo do plot_f() p/ Fortran e OpenMP
		CALL cpu_time(fortran_start_plot_f)
		!$ omp_start_plot_f = omp_get_wtime()

		CALL plot_f()
	
		CALL cpu_time(fortran_end_plot_f)
		!$ omp_end_plot_f = omp_get_wtime()
		!write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
		soma_plot_f_f90 = soma_plot_f_f90 + (fortran_end_plot_f-fortran_start_plot_f)
		soma_plot_f_omp = soma_plot_f_omp + (omp_end_plot_f-omp_start_plot_f)
		!write(*,*) "Tempo individual do plot_f() p/ Fortran:", fortran_end_plot_f-fortran_start_plot_f
		!write(*,*) "Tempo acumulado do plot_f() p/ Fortran", soma_plot_f_f90
		!write(*,*) "Tempo individual do plot_f() p/ OpenMP:", omp_end_plot_f-omp_start_plot_f
		!write(*,*) "Tempo acumulado do plot_f() p/ OpenMP", soma_plot_f_omp

		if (mod(it,ceiling(interv_rest/dt)).eq.0) then
			CALL restart_salva()
		endif
	enddo

	!Atributos finais da simulação

    CALL cpu_time(fortran_start_plot_atri)
    !$ omp_start_plot_atri = omp_get_wtime()

    CALL plot_atrib()

    CALL cpu_time(fortran_end_plot_atri)
    !$ omp_end_plot_atri = omp_get_wtime()

    write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
    soma_plot_atri_f90 = soma_plot_atri_f90 + (fortran_end_plot_atri-fortran_start_plot_atri)
    soma_plot_atri_omp = soma_plot_atri_omp + (omp_end_plot_atri-omp_start_plot_atri)
    !write(*,*) "Tempo acumulado do plot_atri() p/ Fortran", soma_plot_atri_f90
    !write(*,*) "Tempo acumulado do plot_atri() p/ OpenMP", soma_plot_atri_omp

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
	!$ openmp_end = omp_get_wtime()

	!Tempo do Fortran e OpenMP
	write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~- FIM ~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-"
	write(*,*) soma_tempo_f90
	write(*,*) soma_tempo_omp
	write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
	write(*,*) soma_posdin_f90
	write(*,*) soma_posdin_omp	
	write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
	write(*,*) soma_visco_f90
	write(*,*) soma_visco_omp
	write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
	write(*,*) (soma_contorno1_f90 + soma_contorno2_f90 + soma_contorno3_f90)
	write(*,*) (soma_contorno1_omp + soma_contorno2_omp + soma_contorno3_omp)
	write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
	write(*,*) soma_convdiff_f90
	write(*,*) soma_convdiff_omp
	write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
	write(*,*) (soma_plot_i_f90 + soma_plot_f_f90 + soma_plot_atri_f90)
	write(*,*) (soma_plot_i_omp + soma_plot_f_omp + soma_plot_atri_omp)
	write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
	write(*,*) soma_level_set_f90
	write(*,*) soma_level_set_omp
	write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
	write(*,*) soma_graddin_f90
	write(*,*) soma_graddin_omp
	write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
	write(*,*) (fortran_end - fortran_start)
	write(*,*) (openmp_end - openmp_start)
	write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"

END PROGRAM PNH