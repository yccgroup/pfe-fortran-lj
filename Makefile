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

TestPot_SRCS = const.f90 LJparticles.f90 modopenmm.f90 TestPot.f90
TestPot_OBJS = $(patsubst %.f90, %.o, $(TestPot_SRCS))

TestSort_SRCS = utilities.f90 sort.f90 TestSort.f90
TestSort_OBJS = $(patsubst %.f90, %.o, $(TestSort_SRCS))

# program name
PROGRAM = pfe
TestPot_PROGRAM = TestPot
TestSort_PROGRAM = TestSort

all: $(PROGRAM) $(TestPot_PROGRAM) $(TestSort_PROGRAM)

$(PROGRAM): $(OBJS)
	$(FC) $(FLFLAGS) -o $@ $+

$(TestPot_PROGRAM): $(TestPot_OBJS)
	$(FC) $(FLFLAGS) -o $@ $+

$(TestSort_PROGRAM): $(TestSort_OBJS)
	$(FC) $(FLFLAGS) -o $@ $+

%.o: %.f90
	$(FC) $(FCFLAGS) -o $@ $<

clean:
	rm -f *.o *.mod $(PROGRAM) $(TestPot_PROGRAM) $(TestSort_PROGRAM)

