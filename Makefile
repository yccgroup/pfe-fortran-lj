# compiler
FC = gfortran

# compile flags
FCFLAGS = -g -c -O2

# link flags
FLFLAGS =

# source files and objects
SRCS = const.f90 utilities.f90 LJparticles.f90 sampling.f90 nsrafep.f90 debug.f90 Main.f90
OBJS = $(patsubst %.f90, %.o, $(SRCS))

# program name
PROGRAM = rafep

all:
	make rafep

$(PROGRAM): $(OBJS)
	$(FC) $(FLFLAGS) -o $@ $+

%.o: %.f90
	$(FC) $(FCFLAGS) -o $@ $<

clean:
	rm -f *.o *.mod rafep

