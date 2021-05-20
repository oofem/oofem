# OOFEM.org
![Build workflow](https://github.com/oofem/oofem/actions/workflows/build.yml/badge.svg)
[![License: LGPL v2.1](https://img.shields.io/badge/License-LGPL%20v2.1-blue.svg)](https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4339630.svg)](https://doi.org/10.5281/zenodo.4339630)


OOFEM is parallel, object-oriented finite element code for solving mechanical, transport and fluid mechanics problems. 

OOFEM is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.

Copyright (C) 1993 - 2021   Borek Patzak
    
## Getting Started
### What is here
The source directory tree holds source code to the OOFEM package.  
```
  OOFEM_TOP_DIR
  |
  |-- doc - contains the "User's guide", sources to generate "Reference manual", 
  |         documents describing the input file specifications, element and
  |         material libraries, and other useful documents. 
  |
  |-- src - source files of all oofem modules
  |   |
  |   |-- oofemlib - sources of the core part of OOFEM, the OOFEMlib module.
  |   |
  |   |-- sm       - sources of structural analysis module.
  |   |
  |   |-- tm       - sources of transport problem module.
  |   |
  |   |-- fm       - sources of fluid mechanics module.
  |   |
  |   |-- dss      - includes the sources for Direct Sparse Solver (DSS),
  |   |              contributed by R. Vondracek)
  |   |
  |   |-- main     - contains the sources of main() and implementation of some 
  |                  global functions for oofem, oofeg.
  |
  |-- tools   - sources for several utility programs.
  |
  |-- tests   - contains several tests, which are useful to verify
  |             the program functionality.
  |
  |-- bindings - holds sources to generate OOFEM bindings to Python programming language.
```


### Pre-requisites

* The oofem requires the CMake cross-platform build system and C++ compiler with STL support (Standard Template Library).

* The oofem contains interface to IML++ library. It is the C++ templated library of modern iterative methods for solving both symmetric and 
non-symmetric linear systems of equations, written by Roldan Pozo. It can be downloaded from http://math.nist.gov/iml++/. 

* Parallel support for distributed memory requires MPI library to be installed. If you do not have any, we recommend to use Open MPI. 
This is a freely available, high-performance, and portable implementation of MPI (http://www.open-mpi.org/).

* For high performance linear solvers, OOFEM can use [PETSc](https://www.mcs.anl.gov/petsc/), [SuperLU](http://crd-legacy.lbl.gov/~xiaoye/SuperLU/), MKL [Pardiso](https://software.intel.com/en-us/mkl-developer-reference-fortran-intel-mkl-pardiso-parallel-direct-sparse-solver-interface) or [Pardiso-project.org]( http://www.pardiso-project.org/) solvers. 

* For high performance eigenvalue solvers, oofem uses SLEPc. The toolkit can be downloaded from SLEPc project home page 
(http://www.grycap.upv.es/slepc). 

## Installation

* unpack sources 
* create an out-of-tree build directory
```
   mkdir -p ~/build/debug
```
* configure the target
```
   cd ~/build/debug; cmake PATH_TO_OOFEM_SOURCES
```
   where PATH_TO_OOFEM_SOURCES is the path to OOFEM source directory,
   created in step 0 (~/oofem-2.2, for example). 
* compile OOFEM
```
  make
```

### Running the tests
To run the tests, go to your build directory and run ctest
```
ctest
```


### Additional instructions 
* Instead of cmake you can use ccmake which uses an ncurses interface, 
  or cmake-gui for a GUI. Use the command make help for a list of all targets. 

* You can find detailed installation (including installation on Windows)
  instruction on OOFEM wiki (http://www.oofem.org/wiki/doku.php?id=installation)

* To get support check out oofem wiki (www.oofem.org/wiki) and
  oofem forum (www.oofem.org/forum) as well.


### Running oofem 
The oofem program prints out the available options when run without
any option. To run a specific job, you can enter
```
oofem -f input_file_name
```
## Documentation
You may find OOFEM documentation [here](http://oofem.org/en/manual).



## Further information
Please consult oofem home page (http://www.oofem.org) for 
more documentation, manuals and new versions.

## Authors

See the list of [contributors](http://www.oofem.org/doku.php?id=en:credits) who participated in this project.

## Acknowledgments
http://www.oofem.org/doku.php?id=en:funding




