#======================================================================
# Makefile for SuLi
#=======================================================================

# GNU
FC = gfortran
OPTFC = -O3 -fdefault-real-8 -fdefault-double-8 -mcmodel=medium # -fcheck=all  #  #large #-floop-optimize -funroll-loops -ftree-fre -finline-functions -fpeel-loops -funsafe-math-optimizations -fprefetch-loop-arrays 
CC = gcc

# Intel
#FC = ifort
#OPTFC =  -fpp -O3 -xHost -IPO -autodouble
#CC = icc

FLAGS = -O3
LIBS = -fopenmp	#Pedro

# List of source files
SRC = 	$(wildcard *.f90)

#-----------------------------------------------------------------------
# Normally no need to change anything below

OBJ =	$(SRC:.f90=.o)

rodar_3d : $(OBJ)
	$(FC) -o $@ $(OBJ) $(LIBS)

%.o : %.f90
	$(FC) $(OPTFC) $(OPTIONS) -c $<

.PHONY: clean

clean:
	rm -f *.o *.mod 

allclean:
	rm -f *.o *.mod rodar_3d
	rm -rf arquivos
	rm -rf dados
	rm -f saida
	rm -f fort.99998799