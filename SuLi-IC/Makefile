#======================================================================
# Makefile for SuLi
#=======================================================================

# GNU
FC = gfortran  # Fortran Compiler
FFLAGS = -O3 -fdefault-real-8 -fdefault-double-8 -mcmodel=medium -fcheck=all -fopenmp # -g -Wimplicit # #  #large #-floop-optimize -funroll-loops -ftree-fre -finline-functions -fpeel-loops -funsafe-math-optimizations -fprefetch-loop-arrays # -Wall -Wextra -std=f2008
EXE = rodar_3d

#0_modules.o: FFLAGS+= -fopenmp
#1_principal.o: FFLAGS+= -fopenmp
#2_iniciais.o: FFLAGS+= -fopenmp
#3_contorno.o: FFLAGS+= -fopenmp
#4_plot.o: FFLAGS+= -fopenmp
#5_viscosidade.o: FFLAGS+= -fopenmp
#6_convdiff.o: FFLAGS+= -fopenmp
#7_advectivo.o: FFLAGS+= -fopenmp
#8_derivadas.o: FFLAGS+= -fopenmp
#9_grad.o: FFLAGS+= -fopenmp
#10_pos.o: FFLAGS+= -fopenmp
#11_diversos.o: FFLAGS+= -fopenmp
#12_paraview_3D.o: FFLAGS+= -fopenmp
#13_mms.o: FFLAGS+= -fopenmp
#14_level_set.o: FFLAGS+= -fopenmp

# List of source files
#SRC = 	$(wildcard *.f90)		ALteração meio sus		### PEDRO ###
SRC = 0_modules.f90 1_principal.f90 2_iniciais.f90 3_contorno.f90 4_plot.f90 5_viscosidade.f90 6_convdiff.f90 7_advectivo.f90 8_derivadas.f90 9_grad.f90 10_pos.f90 11_diversos.f90 12_paraview_3D.f90 13_mms.f90 14_level_set.f90

#OBJ = ${SRC:.f90=.o}	#Substitui .f90 por .o
OBJ = $(patsubst %.f90,%.o,$(SRC))

rodar_3d:

$(EXE): $(OBJ)
	$(FC) -o $@ $^ $(FFLAGS)

%.o: %.f90				#Regra: criação de *.o depende da existência de *.f90
	$(FC) -c -o $@ $^ $(FFLAGS)

clean: 
	rm -f *.o *.mod

allclean:
	rm -f *.o *.mod rodar_3d
	rm -rf arquivos
	rm -rf dados
	rm -f saida
	rm -f simulacao_
	rm -f 0_modules.o
	rm -f fort.99998799
