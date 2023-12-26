# PBC or FBC
BC := PBC

# compiler
FC = gfortran

# compile flags
FCFLAGS = -g -c -O2

# link flags
FLFLAGS =

# source files and objects
SRCS = tempar.f90 utilities.f90 $(BC)LJparticles.f90 sampling.f90 nsrafep.f90 debug.f90 Main.f90
OBJS = $(patsubst %.f90, %.o, $(SRCS))

# program name
PROGRAM = rafep-$(BC)

all:
	make rafep BC=PBC
	make rafep BC=FBC

.PHONY: rafep
rafep: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(FC) $(FLFLAGS) -o $@ $+

%.o: %.f90
	$(FC) $(FCFLAGS) -o $@ $<

clean:
	rm -f *.o *.mod rafep-PBC rafep-FBC

