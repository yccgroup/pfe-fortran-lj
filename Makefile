# compiler
FC = gfortran

# compile flags
FCFLAGS = -g -c -O3 -Wall -march=native
#FCFLAGS = -g -c -O -Wall -fcheck=all -ffpe-trap=invalid,zero,overflow

# link flags
FLFLAGS =

# source files and objects
SRCS = const.f90 utilities.f90 LJparticles.f90 sampling.f90 pfe.f90 debug.f90 Main.f90
OBJS = $(patsubst %.f90, %.o, $(SRCS))

# program name
PROGRAM = pfe

all: pfe

$(PROGRAM): $(OBJS)
	$(FC) $(FLFLAGS) -o $@ $+

%.o: %.f90
	$(FC) $(FCFLAGS) -o $@ $<

clean:
	rm -f *.o *.mod pfe

