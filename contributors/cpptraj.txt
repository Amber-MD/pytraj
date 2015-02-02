About CPPTRAJ
=============

CPPTRAJ is a program designed to load and analyze molecular dynamics
trajectories and relevant data sets derived from their analysis.

CPPTRAJ is a rewrite of the PTRAJ trajectory analysis code in C++ that is part
of AmberTools. The intent is to make the code more readable, leak-free, and
thread-safe. The biggest functional change from PTRAJ is the ability to load and
process trajectories with different topology files in the same run.

Disclaimer and Copyright
========================

CPPTRAJ is Copyright (c) 2010-2014 Daniel R. Roe. 
The terms for using, copying, modifying, and distributing CPPTRAJ are 
specified in the file LICENSE.

Contributors to CPPTRAJ
=======================
Lead Author: Daniel R. Roe (daniel.r.roe@gmail.com)
             Department of Medicinal Chemistry
             University of Utah, Salt Lake City, UT.

  CPPTRAJ is based on PTRAJ by Thomas E. Cheatham, III (University of Utah,
Salt Lake City, UT, USA) and many routines from PTRAJ have been adapted for 
use by CPPTRAJ, including (but not limited to) code used in the following 
classes:
  - Analysis_CrankShaft
  - Analysis_Statistics
  - Action_DNAionTracker
  - Action_RandomizeIons
  - Action_Principal
  - Action_Grid
  - Grid
  - Action_Image
  - ImageRoutines

Contributors
============
  - James Maier (Stony Brook University, Stony Brook, NY, USA)
    Code for calculating J-couplings (used in Action_Jcoupling).

  - Jason M. Swails (University of Florida, Gainesville, FL, USA)
    Action_LIE, Analysis_RunningAvg, Action_Volmap, Grid OpenDX output.

  - Jason M. Swails (University of Florida, Gainesville, FL, USA)
  - Guanglei Cui (GlaxoSmithKline, Upper Providence, PA, USA)
    Action_SPAM

  - Mark J. Williamson (Unilever Centre for Molecular Informatics, 
      Department of Chemistry, Cambridge, UK).
    Action_GridFreeEnergy

  - Hannes H. Loeffler (STFC Daresbury, Scientific Computing Department,
                        Warrington, WA4 4AD, UK)
    Action_Density, Action_OrderParameter, Action_PairDist

  - Crystal N. Nguyen (University of California, San Diego).
  - Romelia F. Salomon (University of California, San Diego).
    Action_Gist

  The following people wrote code in PTRAJ which was eventually adapted for
use in CPPTRAJ:
  - Holger Gohlke (Heinrich-Heine-University, Düsseldorf, Germany)
  - Alrun N. Koller (Heinrich-Heine-University, Düsseldorf, Germany) 
    Original implementation of matrix/vector functionality in PTRAJ, including
    matrix diagonalization, IRED analysis, eigenmode analysis, and vector time 
    correlations.

  - Holger Gohlke (Heinrich-Heine-University, Düsseldorf, Germany)
    Original code for DSSP (secstruct).

  - Michael Crowley (University of Southern California, Los Angeles, CA, USA)
    Original code for dealing with truncated octahedral unit cells.

  - Viktor Hornak (Merck, NJ, USA)
    Original code for mask expression parser.

  - John Mongan (UCSD, San Diego, CA, USA)
    Original implementation of the Amber NetCDF trajectory format.

  - Hannes H. Loeffler (STFC Daresbury, Scientific Computing Department,
                        Warrington, WA4 4AD, UK)
    Diffusion calculation code adapted for use in Action_STFC_Diffusion.

Documentation
=============
  The main documentation for CPPTRAJ usage is in the AmberTools user manual,
available with any AmberTools distribution in `$AMBERHOME/doc`. There is also
limited help for commands in interactive mode:
```
help [<command>]
```
`help` with no arguments lists all known commands.

  Code documentation can be generated via Doxygen by typing `make docs`. This
will install HTML and Latex documentation at doc/html/index.html and in 
the `doc/latex` respectively. A limited developers guide is available in
Lyx format in `doc/CpptrajDevlopmentGuide.lyx`.


Installation & Testing
======================
Run `./configure --help` for the complete list of configure options.

  `./configure gnu` should be adequate to set up compilation for most systems.
For systems without BLAS/LAPACK/ARPACK and/or NETCDF libraries installed,
the `-amberlib` flag can be specified to use the ones already compiled in
an AmberTools installation (`$AMBERHOME` must be set), e.g.
'./configure -amberlib gnu'. For multicore systems, the `-openmp` flag can
be specified to enable OpenMP parallelization, e.g. `./configure -openmp gnu`.

The configure script by default sets everything up to link dynamically. The
`-static` flag can be used to force static linking. If linking errors are
encountered you may need to specify library locations using the `--with-LIB=`
options. For example, to use NetCDF compiled in /opt/netcdf use the option 
`--with-netcdf=/opt/netcdf`. Alternatively, individual libraries can be 
disabled with the '-noLIB' options.

After config.h has been successfully generated, `make install` will
compile and place the cpptraj binary in the `bin/` subdirectory.
