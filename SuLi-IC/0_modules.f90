!Características básicas do domínio, escoamento e método numérico SuLi-IC
module disc

    real(8),parameter ::  pi = acos(-1.) 

    !Discretizações espaciais em x e y (metros), discretização temporal (segundos)
    real(8),parameter :: dx = 0.02, dy = 0.02, dz = 0.02

    real(8) :: t, dt = 0.001, t_i, t_a
    !Número de células para x, y e z (-); número de pontos para x, y e z (-); tempo de simulação (segundos)

    !Número de tempo por arquivo plotado
    real(8),parameter :: dt_frame = 0.00001*10.

    integer,parameter :: nx=int(3.6/dx) , ny=int(1.2/dy), nz=int(1.2/dz)
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
    integer,parameter :: adv_type = 1 ! 1 = advectivo clássico, 2 = rotacional (apenas para der = 1), 3 = antissimétrico (apenas para der = 1)
        
    integer,parameter :: obst_t = 11 ! 0 = sem obst, 1 = dunas, 2 = dunas2, 3 = gaussiano3D, 4 = beji, 5 = delft degrau, 6 = delft 1_2, 7 = SBRH calombos e buracos, 8 = fennema1990, 9 = aureli2008, 10 = bd_koshizuka1995eKleefsman2005, 11= canal

    integer,parameter :: m_turb = 1 ! 0 = sem modelo, 1 = LES Smagorinsky-Lilly Clássico, 2 = LES Smagorinsky-Lilly Direcional

    integer,parameter :: esp_type = 1 ! 0 = sem camada esponja, 1 = leva em consideração a profundidade, 2 = não leva em consideração a profundidade, 3 = Método da Tangente Hiperbólica

    integer,parameter :: wave_t = 5 ! 0 = sem onda, 1 = Stokes I, 2 = Stokes II, 5 = Stokes V

    integer,parameter :: mms_t = 0  ! 0 = sem MMS, 1 = MMS permanente, 2 = MMS não permanente

    integer :: it, tt, ntt

    real(8),dimension(3) :: a_dt

end module disc

module restart
    integer, parameter :: irest=0 !0=essa simulacao nao é um restart de outra 1 = o arquivo de restart vai ser lido e usado pra continuar a simulacao
    real :: interv_rest=100 !de quantas em quantas iteracoes salva o restart
end module restart

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

    USE disc, only: nx, nx1, ny, ny1, nz, nz1, dz, pi

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

    USE disc, only: nx, nx1, nx2, ny, ny1, ny2, nz, nz1, nz2

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

    USE disc, only: nx, nx1, ny, ny1, nz, nz1

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

    USE disc, only: nz

    real(8) :: p_w, n_w, c_w, l_w, a_w, f_w, h0_f, l0_w
    real(8) :: avel1, avel2, avel3, avel4, avel5
    real(8) :: aeta1, aeta2, aeta3, aeta4, aeta5
    real(8), dimension(0:nz+1) :: kp

end module wave_c

!Adicionado por Luísa Lucchese 07/16, variáveis para o LES
module smag

   USE disc, only: nx, nx1, ny, ny1, nz, nz1

    real(8), parameter :: csmag=0.13 ! 0.8 !
    real(8), dimension(nx,ny,nz)  :: nut 
    real(8), dimension(nx1,ny,nz) :: xnut 
    real(8), dimension(nx,ny1,nz) :: ynut 
    real(8), dimension(nx,ny,nz1) :: znut 

end module smag

module ls_param
    
    USE disc, only: nx, ny, nz, dx, dy, dz

    real(8), dimension(nx,ny,nz) :: ls, mod_ls, kurv, hs, ddlsdx, ddlsdy, ddlsdz, hsx, hsy, hsz
    real(8) :: dtau, alpha1, mi_f1, mi_f2, rho_f1, rho_f2 , vol_ini, vol_ins, ls_m, rho_m, sigma

    real(8), dimension(3) :: adtl, bdtl, gdtl
    real(8),parameter :: dx1 =  max(dx,dy,dz) !(dx+dy+dz)/3. !
    real(8) :: dt1 = 0.1 * dx1
    integer,parameter :: t_hs = 0 ! tipo de função heaviside

end module ls_param

!Variáveis para o LES
module mms_m

    use disc, only: nx, nx1, ny, ny1, nz, nz1

    real(8), parameter :: a  = 0. !5
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
	real(8) :: fortran_start, fortran_end
	real(8) :: openmp_start, openmp_end

	real(8) :: fortran_start_level_set, fortran_end_level_set
	real(8) :: openmp_start_level_set, openmp_end_level_set

	real(8) :: fortran_start_visco, fortran_end_visco
	real(8) :: omp_start_visco, omp_end_visco

	real(8) :: fortran_start_graddin, fortran_end_graddin
	real(8) :: omp_start_graddin, omp_end_graddin

	real(8) :: fortran_start_plot_i, fortran_end_plot_i
	real(8) :: omp_start_plot_i, omp_end_plot_i

	real(8) :: fortran_start_plot_f, fortran_end_plot_f
	real(8) :: omp_start_plot_f, omp_end_plot_f

	real(8) :: fortran_start_plot_atri, fortran_end_plot_atri
	real(8) :: omp_start_plot_atri, omp_end_plot_atri

	real(8) :: start_outros_f90, end_outros_f90, start_outros2_f90, end_outros2_f90, start_contorno3_f90, end_contorno3_f90
	real(8) :: start_outros_omp, end_outros_omp, start_outros2_omp, end_outros2_omp, start_contorno3_omp, end_contorno3_omp

	real(8) :: start_outros4_f90, end_outros4_f90, start_outros5_f90, end_outros5_f90
	real(8) :: start_outros4_omp, end_outros4_omp, start_outros5_omp, end_outros5_omp

	real(8) :: fortran_start_grad_1, fortran_end_grad_1, fortran_start_grad_2, fortran_end_grad_2 
	real(8) :: omp_start_grad_1, omp_end_grad_1, omp_start_grad_2, omp_end_grad_2
    
	real(8) :: soma_level_set_f90, soma_visco_f90, soma_graddin_f90, soma_plot_i_f90, soma_plot_f_f90, soma_plot_atri_f90
	real(8) :: soma_level_set_omp, soma_visco_omp, soma_graddin_omp, soma_plot_i_omp, soma_plot_f_omp, soma_plot_atri_omp
	real(8) :: soma_outros_f90, soma_outros2_f90, soma_contorno3_f90, soma_outros4_f90, soma_outros5_f90
	real(8) :: soma_outros_omp, soma_outros2_omp, soma_contorno3_omp, soma_outros4_omp, soma_outros5_omp
	real(8) :: soma_grad_1_f90, soma_grad_2_f90
	real(8) :: soma_grad_1_omp, soma_grad_2_omp

    real(8) :: fortran_start_convdiff, fortran_end_convdiff, fortran_start_posdin, fortran_end_posdin
	real(8) :: omp_start_convdiff, omp_end_convdiff, omp_start_posdin, omp_end_posdin

    real(8) :: fortran_start_contorno1, fortran_end_contorno1, fortran_start_contorno2, fortran_end_contorno2
	real(8) :: omp_start_contorno1, omp_end_contorno1, omp_start_contorno2, omp_end_contorno2

    real(8) :: fortran_start_tempo, fortran_end_tempo
	real(8) :: omp_start_tempo, omp_end_tempo 

    real(8) :: soma_convdiff_f90, soma_contorno2_f90, soma_posdin_f90, soma_contorno1_f90, soma_tempo_f90
	real(8) :: soma_convdiff_omp, soma_contorno2_omp, soma_posdin_omp, soma_contorno1_omp, soma_tempo_omp
    
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

    USE disc, only: nx, nx1, ny, ny1, nz, nz1

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

	!contadores

	!auxiliares
	real(8) :: aux1, aux2, aux3

end module paodemel

module paodemod_ls11

	USE disc, only: nx, ny, nz, dx, dy, dz
	USE ls_param, only: ls, mod_ls, kurv, ddlsdx, ddlsdy, ddlsdz

	IMPLICIT NONE

	real(8) :: aux1, aux2
	real(8), dimension(nx,ny,nz) :: ta1,tb1,tc1,td1,te1,tf1,dlsdxa,dlsdya,dlsdza

end module paodemod_ls11

module paodecontorno

    USE disc, only: nx, nx1, ny, ny1, nz, nz1, mms_t, obst_t
    USE cond
    USE obst, only: ub, vb, wb, ku, kv, kw
    USE ls_param, only: ls
    USE velpre  	!prd0, prd, rho, ls_nu, d_max, d_min, b_eta0, b_eta1 não são usados

    real(8) :: zi, zj, zk
    integer :: niv, ii
    real(8),dimension(0:nx1+1,0:ny+1,0:nz+1) :: dpdx
    real(8),dimension(0:nx+1,0:ny1+1,0:nz+1) :: dpdy
    real(8),dimension(0:nx+1,0:ny+1,0:nz1+1) :: dpdz

end module paodecontorno

module paodeheaviside

	USE disc, only: nx, ny, nz, pi, dt, mms_t
	USE ls_param, only: ls, alpha1, dx1, mi_f1, mi_f2, hs, hsx, hsy, hsz, rho_f1, rho_f2, t_hs
	USE velpre, only: rho, ls_nu

	IMPLICIT NONE

	integer :: coefa1, ihs
	real(8) :: aux1, aux2, aux3, aux4
	real(8), dimension(nx,ny,nz) :: sy60, sy61,ta1,tb1,tc1,td1,te1,tf1

end module paodeheaviside

module paodeprd_corr

	USE disc, only: nx, nx1, ny, ny1, nz, nz1
	USE velpre, only: rho

	IMPLICIT NONE
	!Declarado também no programa

	real(8), dimension(nx1,ny,nz) :: rhox
	real(8), dimension(nx,ny1,nz) :: rhoy
	real(8), dimension(nx,ny,nz1) :: rhoz

end module paodeprd_corr

module paodemms

	USE disc, only: nx, ny, nz, nx1, ny1, nz1, dx, dy,dz, it, t, dt, mms_t, ts, dt_frame, pi
	USE velpre, only: u, v, w, prd1
	USE mms_m
	USE ls_param, only: ls

	IMPLICIT NONE
	
	integer :: nxc, nyc, nzc, npc

	real(8), dimension(nx1,ny,nz) :: u_m
	real(8), dimension(nx,ny1,nz) :: v_m
	real(8), dimension(nx,ny,nz1) :: w_m
	real(8), dimension(nx,ny,nz)  :: p_m
	real(8), dimension(nx,ny)     :: h_m
	real(8), dimension(nx1,ny)    :: h_mx
	real(8), dimension(nx,ny1)    :: h_my
	real(8), dimension(nx1,ny,nz) :: lsx
	real(8), dimension(nx,ny1,nz) :: lsy
	real(8), dimension(nx,ny,nz1) :: lsz

	real(8) :: x,y,z,h, hpi
	real(8) :: erro_u1, erro_v1, erro_w1, erro_p1

end module paodemms

module paodeplot_i

    USE disc, only: nx, nx1, ny, ny1, nz, nz1, dx, dy, dz, dt, dt_frame, it, ts, der, adv_type, m_turb, mms_t, obst_t, t_i, t_plot
	USE ls_param, only: ls, vol_ini, vol_ins
	USE velpre, only: u, v, w, prd1
	USE tempo, only: cont, agora, agora1
	USE smag, only: nut
	USE obst, only: ku, kw
	USE mms_m, only: a, erro_p, erro_u, erro_v, erro_w
	USE cond

	IMPLICIT NONE

	integer :: ii

	real(8), dimension(nx1,ny1,nz1) :: uaux, vaux, waux, x11, y11, z11
	real(8), dimension(nx,ny,nz) :: dudy, dudz, dvdx, dvdz, dwdx, dwdy
	real(8), dimension(nx,ny,nz)    ::nutaux, prdaux, div, kaux, vorti, vortj, vortk, lasux
	real(8), dimension(nx1,ny,nz) :: xnuta
	real(8), dimension(nx,ny1,nz) :: ynuta	
	real(8), dimension(nx,ny,nz1) :: znuta
	real(8), dimension(nx1,ny,nz1) :: auxy
	real(8), dimension(nx1,ny1,nz) :: auxz
	real(8), dimension(0:nx1,0:ny1,0:nz1) :: x1, y1, z1

end module paodeplot_i

module paodelevel_set

	USE disc, only: nx, ny, nz, dx, dy, dz, dt, mms_t
	USE ls_param, only: ls, mod_ls, hs, vol_ini, vol_ins, adtl, bdtl, gdtl, dt1

	IMPLICIT NONE
	
	real(8),dimension(nx,ny,nz) :: sy7_ls,gx_ls,ta1_ls,sy7_ls1,gx_ls1,ta1_ls1
	real(8) :: aux1, aux2, dtaux

end module paodelevel_set

module paodeconv_weno

	USE disc, only: nx, ny, nz
	USE ls_param, only: ls
	USE velpre, only: u, v, w

	IMPLICIT NONE

	real(8),dimension(nx,ny,nz) :: ta1, tb1, tc1, td1, te1, tf1
	real(8) :: apos, aneg, bpos, bneg, cpos, cneg

end module paodeconv_weno

module paodereinic_weno

	USE disc, only: nx, ny, nz
	USE ls_param, only: ls, dx1, dt1, alpha1

	IMPLICIT NONE

	real(8),dimension(nx,ny,nz) :: sy1, sy4, func_s, ddd, ta1, tb1, tc1, td1, te1, tf1
	real(8),dimension(nx,ny,nz) :: lsaux,ls0
	real(8) :: error
	integer :: l, il, nr,ihs,itrl
	real(8) :: mod_ls1, aux1, aux2

end module paodereinic_weno

module paodevisco

	USE disc, only: nx, ny, nz, nx1, ny1, nz1, dx, dy, dz, m_turb
	USE smag
	USE velpre, only: u, v, w, rho, ls_nu

	IMPLICIT NONE
  
	real(8) :: p1, p2, p3
							
	real(8), dimension(nx,ny,nz) :: aux1

	real(8), dimension(nx1,ny,nz) :: dudx_i, dudy_i, dudz_i
	real(8), dimension(nx,ny1,nz) :: dvdx_i, dvdy_i, dvdz_i
	real(8), dimension(nx,ny,nz1) :: dwdx_i, dwdy_i, dwdz_i

	real(8), dimension(nx,ny,nz) :: dudx, dudy, dudz
	real(8), dimension(nx,ny,nz) :: dvdx, dvdy, dvdz
	real(8), dimension(nx,ny,nz) :: dwdx, dwdy, dwdz

	real(8), dimension(nx1,ny,nz) :: dudx_x, dudy_x, dudz_x, dwdx_x, dvdx_x
	real(8), dimension(nx,ny1,nz) :: dudy_y, dvdx_y, dwdy_y, dvdz_y, dvdy_y
	real(8), dimension(nx,ny,nz1) :: dwdy_z, dvdz_z, dwdz_z, dwdx_z, dudz_z

	real(8), dimension(nx1,ny1,nz) :: dvdx_a, dudy_a
	real(8), dimension(nx1,ny,nz1) :: dwdx_a, dudz_a 
	real(8), dimension(nx,ny1,nz1) :: dwdy_a, dvdz_a

	real(8) :: aux2,aux3,aux4,deltag,nuor

	real(8), dimension(nx,ny1,nz) :: ynut_a
	real(8), dimension(nx1,ny,nz) :: xnut_a
	real(8), dimension(nx,ny,nz1) :: znut_a

	real(8), dimension(nx1,ny,nz) :: ls_nux
	real(8), dimension(nx,ny1,nz) :: ls_nuy
	real(8), dimension(nx,ny,nz1) :: ls_nuz

end module paodevisco

module paodeconvdiff

	USE disc, only: nx, ny, nz, nx1, ny1, nz1, dx, dy, dz, adv_type
	USE velpre, only: u, v, w, rho, prd0, prd1
	USE parametros, only: chezy, gx, gz
	USE smag, only: nut, xnut, ynut, znut
	USE ls_param, only: kurv, hsx, hsy, hsz, sigma
	USE vartempo, only: Fu, Fv, Fw
	USE mms_m, only: a, tf_u, tf_v, tf_w
	USE obst, only: ub, vb, wb

	IMPLICIT NONE

	! auxiliares de velocidades: velocidades lagrangianas
	real(8), dimension(nx1,ny,nz) :: uint
	real(8), dimension(nx,ny1,nz) :: vint
	real(8), dimension(nx,ny,nz1) :: wint

	real(8), dimension(0:nx1,ny,nz) :: ama
	real(8), dimension(nx,0:ny1,nz) :: bmb
	real(8), dimension(nx,ny,0:nz1) :: dmd
	real(8), dimension(nx1,ny1,nz)  :: amb, bma
	real(8), dimension(nx1,ny,nz1)  :: amd, dma
	real(8), dimension(nx,ny1,nz1)  :: bmd, dmb

	!real(8), dimension(0:nx1,ny1,nz)  :: vx
	!real(8), dimension(0:nx1,ny,nz1)  :: wx

	!real(8), dimension(nx1,0:ny1,nz)  :: uy
	!real(8), dimension(nx,0:ny1,nz1)  :: wy

	!real(8), dimension(nx1,ny,0:nz1)  :: uz
	!real(8), dimension(nx,ny1,0:nz1)  :: vz

	real(8), dimension(nx1,ny,nz) :: rhox, dhsdx
	real(8), dimension(nx,ny1,nz) :: rhoy, dhsdy
	real(8), dimension(nx,ny,nz1) :: rhoz, dhsdz, epis_z

	!contadores
	integer :: ntal
	real(8)    :: tal

	!auxiliares
	real(8) :: aux1, aux2

end module paodeconvdiff

module paodeclassico

	USE disc, only: nx, ny, nz, nx1, ny1, nz1, dx, dy,dz, der
	USE velpre, only: u, v, w

	IMPLICIT NONE

	real(8), dimension(nx1,ny,nz) :: dudx, dudy, dudz, bma, dma, dudxa, dudya, dudza
	real(8), dimension(nx,ny1,nz) :: dvdx, dvdy, dvdz, amb, dmb, dvdxa, dvdya, dvdza
	real(8), dimension(nx,ny,nz1) :: dwdx, dwdy, dwdz, amd, bmd, dwdxa, dwdya, dwdza
	real(8), dimension(nx,ny,nz)  :: aux

	!contadores

	!auxiliares
	real(8) :: aux1, aux2

end module paodeclassico

module paoderotacional

	USE disc, only: nx, ny, nz, nx1, ny1, nz1, dx, dy,dz, der
	USE velpre, only: u, v, w

	IMPLICIT NONE

	real(8), dimension(0:nx1+1,0:ny1+1,0:nz1+1) :: ap, an, bma, dma
	real(8), dimension(0:nx1+1,0:ny1+1,0:nz1+1) :: bp, bn, amb, dmb
	real(8), dimension(0:nx1+1,0:ny1+1,0:nz1+1) :: dp, dn, amd, bmd

	!
	real(8), dimension(nx1,ny1,nz1) :: dudx, dvdx, dwdx
	real(8), dimension(nx1,ny1,nz1) :: dudy, dvdy, dwdy
	real(8), dimension(nx1,ny1,nz1) :: dudz, dvdz, dwdz

	!
	real(8) :: aa, bb, dd

	!contadores
	integer :: ai, bi, di

	!plotagem
	real(8) :: acont, bcont, dcont 
	integer :: loca(3), locb(3), locd(3)

	!auxiliares
	real(8) :: aux1, aux2

end module paoderotacional

module paodeantissim

	USE disc, only: nx, ny, nz, nx1, ny1, nz1, dx, dy,dz, der
	USE velpre, only: u, v, w

	IMPLICIT NONE

	real(8), dimension(0:nx1+1,0:ny1+1,0:nz1+1) :: ap, an
	real(8), dimension(0:nx1+1,0:ny1+1,0:nz1+1) :: bp, bn
	real(8), dimension(0:nx1+1,0:ny1+1,0:nz1+1) :: dp, dn

	!
	real(8), dimension(nx1,ny1,nz1) :: dudx, dvdx, dwdx
	real(8), dimension(nx1,ny1,nz1) :: dudy, dvdy, dwdy 
	real(8), dimension(nx1,ny1,nz1) :: dudz, dvdz, dwdz
	real(8), dimension(nx1,ny1,nz1) :: bma, dma,amb, dmb, amd, bmd

	!
	real(8) :: aa, bb, dd

	!contadores
	integer ::  ai, bi, di

	!plotagem
	real(8) :: acont, bcont, dcont 
	integer :: loca(3), locb(3), locd(3)

	!auxiliares
	real(8) :: aux1, aux2

end module paodeantissim

module paodeboundary_waves

	USE disc, only: nx, ny, nz, nx1, ny1, nz1, dx, dz,t, wave_t
	USE wave_c, only: n_w, a_w, f_w, h0_f, avel1, avel2, avel3, avel4, avel5, aeta1, aeta2, aeta3, aeta4, aeta5, kp
	USE velpre, only: u, v, w, bxx0, bxy0, bxz0, bxx1
	USE ls_param, only: ls

	IMPLICIT NONE

	real(8) :: aux1, aux2, aux3, aux4, aux5, h_fa, l_wa

	real(8), dimension(0:nx) :: h_f

	real(8),dimension(0:nx1+1,0:ny+1,0:nz+1) :: u1
	real(8),dimension(0:nx+1,0:ny1+1,0:nz+1) :: v1
	real(8),dimension(0:nx+1,0:ny+1,0:nz1+1) :: w1

	real(4),dimension(0:nx+1,0:ny+1,0:nz+1) :: ls1

endmodule paodeboundary_waves

module paodeplot_f

	USE disc, only: nx, ny, nz, nx1, ny1, nz1, dx, dy,dz, ts, it, dt_frame, dt, mms_t, obst_t, t_a, t_i, t_plot
	USE ls_param, only: ls
	USE velpre, only: u, v, w, prd1
	USE tempo
	USE smag, only: nut
	USE obst, only: ku, kw
	USE mms_m, only: a, erro_p, erro_u, erro_v, erro_w

	IMPLICIT NONE
	!Declarado também no programa

	real(8), dimension(nx1,ny1,nz1) :: uaux, vaux, waux, x11, y11, z11
	real(8), dimension(nx,ny,nz) :: dudy, dudz, dvdx, dvdz, dwdx, dwdy
	real(8), dimension(nx,ny,nz)    ::nutaux, prdaux, div, kaux, vorti, vortj, vortk
	real(8), dimension(nx1,ny,nz1) :: auxy
	real(8), dimension(nx1,ny1,nz) :: auxz
	real(8), dimension(0:nx1,0:ny1,0:nz1) :: x1, y1, z1
	integer :: ifile, nfil, ii

	!Número do arquivo de saída
	integer :: dig1, dig2, dig3, dig4, dig5

	!Nome do arquivo de saída
	character(5) chits

	real(8), dimension(nx1,ny,nz) :: xnuta
	real(8), dimension(nx,ny1,nz) :: ynuta	
	real(8), dimension(nx,ny,nz1) :: znuta

end module paodeplot_f

module paodeposdin

	USE disc, only: nx, ny, nz, nx1, ny1, nz1, dx, dy,dz, dt
	USE velpre, only: u, v, w, rho, prd, prd0, prd1
	USE mms_m, only: tf_px, tf_py, tf_pz
	
	IMPLICIT NONE

	real(8) :: aux1, aux2
	real(8), dimension(nx1,ny,nz) :: rhox, hs_x
	real(8), dimension(nx,ny1,nz) :: rhoy, hs_y
	real(8), dimension(nx,ny,nz1) :: rhoz, hs_z

end module paodeposdin
