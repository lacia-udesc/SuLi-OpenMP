!Características básicas do domínio, escoamento e método numérico SuLi-IC
module disc

    real(8),parameter ::  pi = acos(-1.) 

    !Discretizações espaciais em x e y (metros), discretização temporal (segundos)
    real(8),parameter :: dx = 0.02, dy = 0.02, dz = 0.02

    real(8) :: t, dt = 0.001, t_i, t_a
    !Número de células para x, y e z (-); número de pontos para x, y e z (-); tempo de simulação (segundos)

    !Número de tempo por arquivo plotado
    real(8),parameter :: dt_frame = 0.00001*10.

    integer,parameter :: nx=int(3./dx) , ny=int(1./dy), nz=int(1./dz)
    !integer,parameter :: nx=int(10) , ny=int(10), nz=int(10)
    !nz=int(10./dz1-0.1+0.5) porque a última célula é maior (0.5)
    integer,parameter :: nx1=nx+1, ny1=ny+1, nz1=nz+1, ts = ceiling(0.01/0.001)

    integer, parameter :: malha = nx*ny*nz

    integer,parameter :: nx2=nx1+1, ny2=ny1+1, nz2=nz1+1

    !Para fazer dz variável no espaço inicialmente criar uma função ...
    real(8),parameter :: uinicial = 0.2

    integer,parameter :: t_plot = 1 ! 0 = modo simples (velocidade, Level Set e IBM), 1 = modo completo (pressão, vorticidade, viscosidade)

    integer,parameter :: t_tempo = 2 ! 0 = Euler Explícito, 1 = RK 2, 2 = RK 3, 3 = AB2

    integer,parameter :: der = 3 ! 1 = upwind, 2 = centrado, 3 = upwind 2nd order (centrado só para advectivo clássico)
    integer,parameter :: adv_type = 1 ! 1 = advectivo clássico, 2 = rotacional, 3 = antissimétrico
        
    integer,parameter :: obst_t = 11 ! 0 = sem obst, 1 = dunas, 2 = dunas2, 3 = gaussiano3D, 4 = beji, 5 = delft degrau, 6 = delft 1_2, 7 = SBRH calombos e buracos, 8 = fennema1990, 9 = aureli2008, 10 = bd_koshizuka1995eKleefsman2005, 11= canal

    integer,parameter :: m_turb = 1 ! 0 = sem modelo, 1 = LES Smagorinsky-Lilly Clássico, 2 = LES Smagorinsky-Lilly Direcional

    integer,parameter :: esp_type = 1 ! 0 = sem camada esponja, 1 = leva em consideração a profundidade, 2 = não leva em consideração a profundidade, 3 = Método da Tangente Hiperbólica

    integer,parameter :: wave_t = 5 ! 0 = sem onda, 1 = Stokes I, 2 = Stokes II, 5 = Stokes V

    integer,parameter :: mms_t = 0  ! 0 = sem MMS, 1 = MMS permanente, 2 = MMS não permanente

    integer :: it, tt ,ntt

    real(8),dimension(3) :: a_dt

end module disc

module restart
    integer, parameter :: irest=0 !0=essa simulacao nao é um restart de outra 1 = o arquivo de restart vai ser lido e usado pra continuar a simulacao
    real :: interv_rest=100 !de quantas em quantas iteracoes salva o restart
end module restart

!Condicoes de contorno
module cond

    integer,parameter :: ccx0=3 !condicao de contorno parede x=0 --> 0 é periodico, 1 é free-slip, 2 é no-slip, 3 é prescrita, 4 é fluxo validacao
    integer,parameter :: ccxf=4 !condicao de contorno parede x=xf --> 0 é periodico, 1 é free-slip, 2 é no-slip, 3 é prescrita, 4 é saida livre
    !só pode usar condição periódica no final quando usar no começo e vice-versa
    integer,parameter :: ccy0=2 !condicao de contorno parede y=0 --> 0 é periodico, 1 é free-slip e 2 é no-slip, 3 é prescrita
    integer,parameter :: ccyf=2 !condicao de contorno parede y=yf --> 0 é periodico, 1 é free-slip e 2 é no-slip, 3 é prescrita
    integer,parameter :: ccz0=2 !condicao de contorno parede z=0 --> 1 é free-slip, 2 é no-slip, 3 é prescrita
    integer,parameter :: cczf=1 !condicao de contorno parede z=zf --> 1 é free-slip, 3 é prescrita

end module cond

module obst

    USE disc

    !Velocidade de fundo (m/s)
    real(8), dimension(0:nx1+1,0:ny+1,0:nz+1) :: ub
    real(8), dimension(0:nx+1,0:ny1+1,0:nz+1) :: vb
    real(8), dimension(0:nx+1,0:ny+1,0:nz1+1) :: wb

    !Obstáculo
    integer, dimension(0:nx1+1,0:ny+1) :: ku	
    !Indicam até que altura as velocidades tem que ser zeradas (até qual índice k)
    integer, dimension(0:nx+1,0:ny1+1) :: kv
    integer, dimension(0:nx+1,0:ny+1)  :: kw
    integer, parameter :: elev = 5. *dz                    
    !Comprimento a adicionar abaixo

    real(8), parameter :: amp = 0.25, comp = 1., fase = -pi/2. !amplitude, comprimento e fase da onda

end module obst

module velpre

    USE disc

    !Velocidades para x e z (m/s)
    real(8),dimension(0:nx1+1,0:ny+1,0:nz+1) :: u
    real(8),dimension(0:nx+1,0:ny1+1,0:nz+1) :: v
    real(8),dimension(0:nx+1,0:ny+1,0:nz1+1) :: w

    real(8),dimension(0:ny+1,0:nz+1)  :: bxx0, bxx1
    real(8),dimension(0:ny1+1,0:nz+1) :: bxy0
    real(8),dimension(0:ny+1,0:nz1+1) :: bxz0

    real(8),dimension(0:ny+1,0:nz+1)  :: bxxf, bxxf1
    real(8),dimension(0:ny1+1,0:nz+1) :: bxyf
    real(8),dimension(0:ny+1,0:nz1+1) :: bxzf

    real(8),dimension(0:nx1+1,0:nz+1)  :: byx0
    real(8),dimension(0:nx+1,0:nz+1) :: byy0, byy1
    real(8),dimension(0:nx+1,0:nz1+1) :: byz0

    real(8),dimension(0:nx1+1,0:nz+1)  :: byxf
    real(8),dimension(0:nx+1,0:nz+1) :: byyf, byyf1
    real(8),dimension(0:nx+1,0:nz1+1) :: byzf

    real(8),dimension(0:nx1+1,0:ny+1)  :: bzx0
    real(8),dimension(0:nx+1,0:ny1+1) :: bzy0
    real(8),dimension(0:nx+1,0:ny+1) :: bzz0, bzz1

    real(8),dimension(0:nx1+1,0:ny+1)  :: bzxf
    real(8),dimension(0:nx+1,0:ny1+1) :: bzyf
    real(8),dimension(0:nx+1,0:ny+1) :: bzzf, bzzf1

    !Pressão não-hidrostática (m²/s²)
    real(8),dimension(0:nx+1,0:ny+1,0:nz+1) :: prd1, prd0, prd
    real(8),dimension(nx,ny,nz) :: rho, ls_nu

    real(8) :: d_max, d_min, b_eta0, b_eta1

end module velpre

module vartempo

    USE disc
    real(8),dimension(nx1,ny,nz) :: Fu, u0, fu0, fu1
    real(8),dimension(nx,ny1,nz) :: Fv, v0, fv0, fv1
    real(8),dimension(nx,ny,nz1) :: Fw, w0, fw0, fw1

end module vartempo

module parametros

    !Parâmetros
    !Viscosidade cinemática (m²/s), coeficiente de chezy (m**(1/2)/s), aceleração da gravidade (m/s²) e implicitness parameter $Patnaik et al. 1987$ (-) 
    real(8), parameter :: tetah = 0.50, chezy = 75, decliv = 0.00001

    real(8), parameter :: gx = 9.80665 * sin(atan(decliv)) , gz = 9.80665 * cos(atan(decliv))

    !Real(8), parameter :: gx = 0. , gz = 0.

    !Wind stress coefficient (m/s) e velocidade do vento para x e y (m/s)
    !real(8) :: cwind, uwind, vwind

end module parametros

module tempo

    integer(8) :: cont

    ! hora e data
    integer :: agora(8), agora1(8)
    real(8) :: ciclo, prev

end module tempo

module wave_c

    USE disc

    real(8) :: p_w, n_w, c_w, l_w, a_w, f_w, h0_f, l0_w
    real(8) :: avel1, avel2, avel3, avel4, avel5
    real(8) :: aeta1, aeta2, aeta3, aeta4, aeta5
    real(8), dimension(0:nz+1) :: kp

end module wave_c

!Adicionado por Luísa Lucchese 07/16, variáveis para o LES
module smag

    use disc

    real(8), parameter :: csmag=0.13 ! 0.8 !
    real(8), dimension(nx,ny,nz)  :: nut 
    real(8), dimension(nx1,ny,nz) :: xnut 
    real(8), dimension(nx,ny1,nz) :: ynut 
    real(8), dimension(nx,ny,nz1) :: znut 

end module smag

module ls_param

    USE disc

    real(8), dimension(nx,ny,nz) :: ls, mod_ls, kurv, hs, ddlsdx, ddlsdy, ddlsdz, hsx, hsy, hsz
    real(8) :: dtau, alpha1, mi_f1, mi_f2, rho_f1, rho_f2 , vol_ini, vol_ins, ls_m, rho_m, sigma

    real(8), dimension(3) :: adtl, bdtl, gdtl
    real(8),parameter :: dx1 =  max(dx,dy,dz) !(dx+dy+dz)/3. !
    real(8) :: dt1 = 0.1 * dx1
    integer,parameter :: t_hs = 0 ! tipo de função heaviside

end module ls_param

!Variáveis para o LES
module mms_m

    use disc

    real(8), parameter :: a  = 0.!5
    real(8), parameter :: h0 = 2.

    ! dt / rho?
    real(8), parameter :: coef = 0. ! 1.
    real(8) :: erro_t ! erro médio acumulado
    real(8), dimension(nx,ny,nz)  :: tf_p, erro_p
    real(8), dimension(nx1,ny,nz) :: tf_u, erro_u, tf_px
    real(8), dimension(nx,ny1,nz) :: tf_v, erro_v, tf_py
    real(8), dimension(nx,ny,nz1) :: tf_w, erro_w, tf_pz

end module mms_m

module omp

	IMPLICIT NONE

	!Declaração das variáveis para análise
	double precision :: fortran_start, fortran_end
	double precision :: openmp_start, openmp_end

	double precision :: fortran_start_level_set, fortran_end_level_set
	double precision :: openmp_start_level_set, openmp_end_level_set

	double precision :: fortran_start_visco, fortran_end_visco
	double precision :: omp_start_visco, omp_end_visco

	double precision :: fortran_start_graddin, fortran_end_graddin
	double precision :: omp_start_graddin, omp_end_graddin

	double precision :: fortran_start_plot_i, fortran_end_plot_i
	double precision :: omp_start_plot_i, omp_end_plot_i

	double precision :: fortran_start_plot_f, fortran_end_plot_f
	double precision :: omp_start_plot_f, omp_end_plot_f

	double precision :: fortran_start_plot_atri, fortran_end_plot_atri
	double precision :: omp_start_plot_atri, omp_end_plot_atri

	double precision :: start_outros_f90, end_outros_f90, start_outros2_f90, end_outros2_f90, start_contorno3_f90, end_contorno3_f90
	double precision :: start_outros_omp, end_outros_omp, start_outros2_omp, end_outros2_omp, start_contorno3_omp, end_contorno3_omp
    
	double precision :: start_outros4_f90, end_outros4_f90, start_outros5_f90, end_outros5_f90
	double precision :: start_outros4_omp, end_outros4_omp, start_outros5_omp, end_outros5_omp

	double precision :: fortran_start_grad_1, fortran_end_grad_1, fortran_start_grad_2, fortran_end_grad_2 
	double precision :: omp_start_grad_1, omp_end_grad_1, omp_start_grad_2, omp_end_grad_2
    
	real :: soma_level_set_f90, soma_visco_f90, soma_graddin_f90, soma_plot_i_f90, soma_plot_f_f90, soma_plot_atri_f90
	real :: soma_level_set_omp, soma_visco_omp, soma_graddin_omp, soma_plot_i_omp, soma_plot_f_omp, soma_plot_atri_omp
	real :: soma_outros_f90, soma_outros2_f90, soma_contorno3_f90, soma_outros4_f90, soma_outros5_f90
	real :: soma_outros_omp, soma_outros2_omp, soma_contorno3_omp, soma_outros4_omp, soma_outros5_omp
	real :: soma_grad_1_f90, soma_grad_2_f90
	real :: soma_grad_1_omp, soma_grad_2_omp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    double precision :: fortran_start_convdiff, fortran_end_convdiff, fortran_start_posdin, fortran_end_posdin
	double precision :: omp_start_convdiff, omp_end_convdiff, omp_start_posdin, omp_end_posdin

    double precision :: fortran_start_contorno1, fortran_end_contorno1, fortran_start_contorno2, fortran_end_contorno2
	double precision :: omp_start_contorno1, omp_end_contorno1, omp_start_contorno2, omp_end_contorno2

    double precision :: fortran_start_tempo, fortran_end_tempo
	double precision :: omp_start_tempo, omp_end_tempo 

    real :: soma_convdiff_f90, soma_contorno2_f90, soma_posdin_f90, soma_contorno1_f90, soma_tempo_f90
	real :: soma_convdiff_omp, soma_contorno2_omp, soma_posdin_omp, soma_contorno1_omp, soma_tempo_omp
    
	!double precision function omp_get_wtime()
	!double precision function omp_get_wtick()

    !   Padrão de rastreio de processos   !

    !CALL cpu_time(fortran_start_processo)
    !!$ omp_start_processo = omp_get_wtime()

    !CALL processo()

    !CALL cpu_time(fortran_end_processo)
    !!$ omp_end_processo = omp_get_wtime()

    !write(*,*) "~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~"
    !soma_processo_f90 = soma_processo_f90 + (fortran_end_processo - fortran_start_processo)
    !soma_processo_omp = soma_processo_omp + (omp_end_processo - omp_start_processo)
    !write(*,*) "Tempo acumulado do processo() p/ Fortran", soma_processo_f90
    !write(*,*) "Tempo acumulado do processo() p/ OpenMP", soma_processo_omp
    !write(*,*) "Tempo individual do processo() p/ Fortran:", fortran_end_processo - fortran_start_processo
    !write(*,*) "Tempo individual do processo() p/ OpenMP:", omp_end_processo - omp_start_processo

end module omp

module paodemel

    USE disc

    IMPLICIT NONE

	real(8), dimension(nx1,ny,nz) :: rhox
	real(8), dimension(nx,ny1,nz) :: rhoy
	real(8), dimension(nx,ny,nz1) :: rhoz
    real(8) :: alfapr, alfamupr, alfadipr, betapr, betamupr
	real(8), dimension(nx1,ny,nz) :: matspri
	real(8), dimension(nx,ny1,nz) :: matsprj
	real(8), dimension(nx,ny,nz1) :: matsprk
	real(8), dimension(nx,ny,nz) :: matqpr, matapripos, mataprineg, mataprjpos, mataprjneg, mataprkpos, mataprkneg, mppr
	real(8), dimension(nx1+1,ny1+1,nz1+1) :: matdpr, matepr, erropr, erroppr

end module paodemel
