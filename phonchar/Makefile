FC = gfortran

FLAGS = -I. -ffree-line-length-none -fopenmp

# debug
#FLAGS+= -Wall -Wextra -Wconversion -pedantic -O -fcheck=all -g -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid -ffree-line-length-none -fopenmp

# opt
FLAGS+= -O3

EXE = phonchar


OBJECTS = modules.o main.o show_logo.o init.o read_poscar.o read_eigen.o \
          print_frac.o rotvec.o calc_dispchar.o calc_eigchar.o calc_dirchar.o \
          deallocate_all.o credits.o

all: $(OBJECTS)
	$(FC) -o $(EXE) $(FLAGS) $(OBJECTS)

clean:
	-rm -f *__genmod.f90 *.o *.mod $(EXE)

%.o: %.f90
	$(FC) $(FLAGS) -c $<
