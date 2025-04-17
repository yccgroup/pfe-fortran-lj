# PFE-Fortran-LJ

This code is a Fortran implementation of the [PFE](https://dx.doi.org/10.1063/5.0237340)
method for a system of Lennard-Jones particles in a fixed volume. The code can perform
the following:

- Collect a number of samples from the canonical distribution of
  the Lennard-Jones potential, via Metropolis Monte-Carlo.
  Alternatively, energy samples can be read in that were produced
  by an external program.
- Estimate the value of the partition function Z (at the given
  temperature) from the samples, via the PFE method.
- Find the optimal value for the PFE energy cutoff `E*`.
- Calculate a reference value for Z by computing the complete
  density of states via nested sampling.


## Prerequisites

The code is written in modern Fortran, using some features from
F2003 and F2008. Therefore, a reasonably recent Fortran compiler
is required (e.g. Gfortran 8 or later).  The compilation is
controlled through a Makefile, please use GNU make to compile.

The code is not parallelized and only requires a moderate amount
of RAM, so it can be run on any standard desktop or laptop machine.


## Compilation

The supplied Makefile works fine when using Gfortran as compiler.
You can edit the Makefile and adjust the following variables:

* `FC`: the Fortran compiler
* `FCFLAGS`: flags for compiling the object files
* `FLFLAGS`: flags for linking the executable

Run `make pfe` to compile just the `pfe` main program. Running `make`
will additionally produce some test programs that were used for
internal testing during code development (these are undocumented).


## Running the Code

To run the code, an input file is needed (see below) and must be named
"input.dat". We recommend to place this and the "pfe" executable into
a dedicated directory. The code can then be run via

    ./pfe

Some results are printed to standard output, therefore it is recommended
to redirect it so that you have a record, e.g. via

    ./pfe | tee pfe.log

The run will produce some additional output files, documented below.


## Input File

The input file "input.dat" specifies all parameters for the computation.
The format is simple, with one parameter per line. A commented sample
input file is provided. Please note that the comments are ignored and
carry no semantic information. Therefore, you must not reorder the lines,
delete any lines, or introduce any extra lines!

The full list of parameters is as follows, where FP, INT, STR
denote the type of the parameter as floating-point (real), integer,
or string (no spaces allowed).

* The system temperature, in Kelvin (FP).
* The number of Lennard-Jones (LJ) particles (INT).
* The mass of one LJ particle, in amu (FP).
* The LJ epsilon parameter, in kcal/mol (FP).
* The LJ sigma parameter, in Å (FP).
* The side length of the cubic box that the particles are confined to, in Å (FP).
  Periodic boundary conditions are applied throughout.
* The LJ potential cutoff distance, in Å (FP).
  This must be less than half the side length of the box.
  The inter-particle potential is set to zero above this distance. Note that
  the LJ potential is vertically shifted such that it goes continuously to
  zero at the cutoff distance.
* The random seed (INT).
  If set to zero, the seed is initialized with random data retrieved from
  the OS environment.
* The number of steps for the Monte Carlo (MC) equilibration phase (INT).
* The maximum step size for the MC equilibration phase, in Å (FP).
  In each MC step, a single particle is randomly moved by displacing each
  coordinate by an amount uniformly chosen between minus and plus this
  step size.
* The output frequency for the MC equilibration phase (INT).
  (Actually, this parameter has no significant effect.)
* The number of steps for the MC production run (INT).
* The maximum step size for the MC production run, in Å (FP).
* The output frequency for the MC production run (INT).
  This is the number of MC steps (both rejected and accepts steps count)
  that need to have passed before a new configuration / data point is output.
* The filename for the MC data (STR).
  This file contains, in binary format, a record of the essential system parameters,
  all recorded MC configurations, and their energies. For details, see subroutine
  `Write` in code file `sampling.f90`.
* The filename for the energy data (STR).
  This file contains, in text format, only the energies of the recorded MC configurations,
  one per line, in kJ/mol.
* The `bothend` flag which determines whether only the highest energies should be cut off
  (`.FALSE.`) or also the lowest energies (`.TRUE.`). See below for an extended description.
  It is recommended to leave this set to `.FALSE`.
* The number of bins for the energy distribution (INT).
  Only used if `bothend` is `.TRUE.`.
* The percentage of highest energies to discard if the automated search for the
  optimal `E*` fails (FP).
  Only used if `bothend` is `.FALSE.`.
* The number of walkers for the nested sampling (NS) run (INT).
* The number of extra NS relaxation steps (INT). See below for a discussion of
  the details of our NS implementation.
* The number of NS equlibration steps (INT).
* The maximum step size for the NS steps, in Å (FP).
* The fraction `p` for determining the new NS energy level (FP).
  Must be positive and smaller than one.
* The initial energy ceiling for NS, in kcal/mol (FP).
  This should be large enough so that a randomly chosen configuration is likely
  to have an energy lower than this.
* The job to run (STR). These will be described in detail below. One of:
  * `MC`: Run Metropolis Monte Carlo to gather samples.
  * `PFE`: Read in energy samples, and run nested sampling to compute the volume
    correction term, and estimate Z via our PFE method.
  * `Exact`: Run nested sampling to compute the complete density of states (DoS),
    and compute Z by integrating the DoS with the Boltzmann factor.
  * `POT`: Compute and output the LJ dimer potential.


## Job Types

### `MC`

This job runs Metropolis Monte Carlo to gather samples from the canonical
distribution of the system according to the parameters from the input file.

Output files:

* `Energy.dat` (default name, but can be changed in input file):
  the energy data of the samples (see above)
* `MC.bin` (default name, but can be changed in input file):
  the MC data (see above)
* `Traj.dat`: for each sample...
  * the current step (`Time = ...`)
  * for each atom, its x/y/z coordinates in Å (one atom per line)

Some information is written to standard output, including:

* `lambda_th`: the thermal deBroglie wavelength (in Å)
* `Acceptance rate`: the acceptance rate of the Metropolis sampling
  (between zero and one)
* `Average energy`: the average energy of the whole system (in kJ/mol)


### `PFE`

This job reads in a number of energy samples (from `Energy.dat` by default),
determines a suitable value for the energy cutoff `E*`, runs nested sampling to
compute the volume correction term, and estimates Z via our PFE method.

Output files:

* `Statistics.dat`: for each nested sampling iteration...
  * the iteration counter
  * the number of outliers encountered in this iteration
  * the total number of relaxation steps, aggregated over the outliers
  * the total number of propagation steps, aggregated over the outliers
  * the average MC acceptance rate for relaxing the outliers
  * the average MC acceptance rate for propagating the outliers
  * the "move distance" average, across the outliers; this measures how
    far an outlier moves in configuration space during the propagation
    phase (in Å)
  * the "move distance" standard deviation, across the outliers
* `Levels.dat`: for each nested sampling level...
  * the energy (in kJ/mol)
  * the log of the relative configuration space volume (dimensionless)
  * the relative configuration space volume (dimensionless);
    this is relative to the full configuration space `L^3N`

Some progress and information are written to standard output, including:

* `PartFunc: Emax / Emin`: the maximum/minimum values among the energy samples (in kcal/mol)
* `PartFunc: Estar`: value of `E*` as determined by the algorithm from Section IIB in our paper
   (note that this might fail if the energy data is too noisy, in which case `E*` is chosen
    by cutting off a given fraction of the highest energy samples)
* `initial conformation generation success rate`: the fraction of initially
   generated random configurations whose system energy falls below the initial
   energy ceiling
* `logVolume`: the log of the relative configuration space volume after the final NS iteration
* `volume`: the relative configuration space volume after the final NS iteration
* `err`: estimate of the relative error of `volume` (and for the absolute error of `logVolume`)
* `PFE lnQ`: the log of the configurational part of the partition function (as computed via
   PFE, cf. Eq.(10) in our paper)
* `PFE lnP`: the log of the momentum part of the partition function (analytical formula,
   cf. Eq.(3) in our paper)
* `PFE lnZ`: the log of the full partition function (sum of `lnP` and `lnQ`)
* `PFE Err VErr Err+VErr`: estimates for relative errors of various terms:
  * `Err`: for the `<f>` term (cf. Eq.(12))
  * `VErr`: for the volume term
  * `Err+VErr`: for their fraction; this is also the absolute error of `PFE lnQ` and thus
    of `PFE lnZ` (as the `lnP` term is analytic)
* `Estar`: the actual value of `E*` used
* `Cutoff %`: the percentage of energy samples that lie above `E*`


### `Exact`

This job is for calculating reference data, by computing the complete density
of states (DoS) via nested sampling down to very low energies, and then
integrating the DoS times the Boltzmann factor to obtain Z.

The DoS gets computed down to a low-cutoff energy `E_low`, which is obtained by
first performing an energy minimization of the whole system, and then
adding kT/1000.  The energy minimization is performed via Monte Carlo,
running for the number of steps specified for the production run,
where moves are only accepted if they lower the system energy.

Output files:

* `Statistics.dat`: as for `PFE`
* `Levels.dat`: as for `PFE`

Some progress and information are written to standard output, including:

* `Emin`: the system energy after energy minimization (in kcal/mol)
* `initial conformation generation success rate`: as for `PFE`
* `logVolume`: the log of the relative configuration space volume after the final NS iteration
* `volume`: the relative configuration space volume after the final NS iteration
* `err`: estimate of the relative error of `volume` (and for the absolute error of `logVolume`)
* `NS lnP`: the log of the momentum part of the partition function (analytical formula,
   cf. Eq.(3) in our paper)
* `NS lnQ`: the log of the configurational part of the partition function (obtained by
   integrating over the density of states)
* `NS lnZ`: the log of the full partition function (sum of `lnP` and `lnQ`)

### `POT`

This job was used to verify the implementation of the LJ potential function.
The number of atoms in the system must be exactly two for this job.

Output Files:

* `pot.dat`:
  * column 1: distance between the atoms (in nm)
  * column 2: inter-particle potential energy, directly based on distance (in kcal/mol)
  * column 3: system energy, taking into account the PBC (in kcal/mol);
    can differ slightly from column 2 due to floating-point rounding errors
    during the calculation of the inter-particle distance

