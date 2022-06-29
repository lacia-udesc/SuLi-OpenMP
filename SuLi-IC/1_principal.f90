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

double precision :: fortran_start, fortran_end
double precision :: openmp_start, openmp_end

double precision :: fortran_start_level_set, fortran_end_level_set
double precision :: openmp_start_level_set, openmp_end_level_set

double precision :: fortran_start_convdiff, fortran_end_convdiff
double precision :: omp_start_convdiff, omp_end_convdiff

double precision :: fortran_start_visco, fortran_end_visco
double precision :: omp_start_visco, omp_end_visco

!double precision function omp_get_wtime()
!double precision function omp_get_wtick()

!Inicio do Fortran e OpenMP
write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~- INÍCIO ~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-"
call cpu_time(fortran_start)
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


!Adicionar os contornos na plotagem inicial
CALl contorno(1)
CALl contorno(3)


!Solução manufaturada
if (mms_t > 0) CALL mms()

!Plotagens iniciais
CALL plot_i()

!RESOLUÇÃO DO PROBLEMA

!Parte 1: Função distância; level_set()
!Parte 2: Viscosidade Turbulenta; visco()
!Parte 3: Passo preditor; convdiff() e tempo()
!Parte 4 : Condições de contorno; call boundary_waves() e contorno()
!Parte 5: Passo corretor; graddin() e posdin()

do it = 1, ts
t = it * dt

!write(*,*) it

!Termo fonte para o método da sulução manufaturada (MMS)
if ((mms_t == 1) .and. (it == 1)) call termo_fonte1()
if (mms_t == 2) call termo_fonte2()

!Tempo do level_set() p/ Fortran e OpenMP
call cpu_time(fortran_start_level_set)
openmp_start_level_set = omp_get_wtime()

CALL level_set()

call cpu_time(fortran_end_level_set)
openmp_end_level_set = omp_get_wtime()

write(*,*)
write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
write(*,*) "Tempo do level_set() p/ Fortran", fortran_end_level_set-fortran_start_level_set
write(*,*) "Tempo do level_set() p/ OpenMP", openmp_end_level_set-openmp_start_level_set

CALl contorno(3)

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

	!Tempo do visco() p/ Fortran e OpenMP
	write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
	write(*,*) "Tempo do convdiff() p/ Fortran:", fortran_end_visco-fortran_start_visco
	write(*,*) "Tempo do convdiff() p/ OpenMP:", omp_end_visco-omp_start_visco

	!Tempo do convdiff() p/ Fortran e OpenMP
	CALL cpu_time(fortran_start_convdiff)
	omp_start_convdiff = omp_get_wtime()

	CALL convdiff()

	CALL cpu_time(fortran_end_convdiff)
	omp_end_convdiff = omp_get_wtime()

	!Tempo do Fortran e OpenMP
	write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
	write(*,*) "Tempo do convdiff() p/ Fortran:", fortran_end_convdiff-fortran_start_convdiff
	write(*,*) "Tempo do convdiff() p/ OpenMP:", omp_end_convdiff-omp_start_convdiff

	CALL tempo()

	if (wave_t > 0) call boundary_waves() !For wave propagation
		!Condições de Contorno para a parte Hidrostática
		!CALL pressh()
	CALl contorno(2)

	if (mms_t .eq. 0) then
		CALL graddin()
		CALL posdin()
		CALl contorno(1)
	endif

 enddo


!Solução manufaturada; cálculo do erro
if (mms_t > 0) CALL mms()

!Plotagens por passo de tempo
CALL plot_f()

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


!Fim do Fortran
call cpu_time(fortran_end)

!Fim do OpenMP
openmp_end = omp_get_wtime()

!Tempo do Fortran
write(*,*) "Tempo Total do Fortran:", fortran_end-fortran_start 
!Tempo do OpenMP
write(*,*) "Tempo Total do OpenMP:", openmp_end-openmp_start
write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~- FIM ~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-"
End program PNH