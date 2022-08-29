#======================================================================
# Makefile for SuLi
#=======================================================================

# GNU
FC = gfortran  # Fortran Compiler
FFLAGS = -O3 -fdefault-real-8 -fdefault-double-8 -mcmodel=medium -fcheck=all # #  #large #-floop-optimize -funroll-loops -ftree-fre -finline-functions -fpeel-loops -funsafe-math-optimizations -fprefetch-loop-arrays # -Wall -Wextra -std=f2008
EXE = rodar_3d

9_grad.o: FFLAGS+= -fopenmp

# List of source files
#SRC = 	$(wildcard *.f90)		ALteração meio sus		### PEDRO ###
SRC = 0_modules.f90 1_principal.f90 2_iniciais.f90 3_contorno.f90 4_plot.f90 5_viscosidade.f90 6_convdiff.f90 7_advectivo.f90 8_derivadas.f90 9_grad.f90 10_pos.f90 11_diversos.f90 12_paraview_3D.f90 13_mms.f90 14_level_set.f90

#OBJ = ${SRC:.f90=.o}	#Substitui .f90 por .o
OBJ = $(patsubst %.f90,%.o,$(SRC))

rodar_3d: $(EXE)

$(EXE): $(OBJ)
	$(FC) -o $@ $^ $(FFLAGS) -fopenmp

%.o: %.f90				#Regra: criação de *.o depende da existência de *.f90
	$(FC) -c -o $@ $^ $(FFLAGS)

clean: 
	rm -f *.o *.mod

allclean:
	rm -f *.o *.mod rodar_3d
	rm -rf arquivos
	rm -rf dados
	rm -f saida
	rm -f 0_modules.o
	rm -f fort.99998799
