Introduction
============

OOFEM is a general purpose, object-oriented Finite Element code for solving mechancal. transport and fluid mechanics problems, primarily focused on academic/research use. It has been continuously developped since 1997. At present it consist of more than 230 thousands lines of source C++ code. It has acive development team and large user base around the world. The code has been sucessfully applied to solution of several industrial problems.

OOFEM is trying to be one of the best FEM solvers. To make this a reality, we focus on this and do not try to provide extensive pre and post processing capabilities. Here we rely on third party external tools.

OOFEM as a FEM solver solves the problems described by a set of partial differential equations. On input it expects the discretized domain and corresponding problem parameters, initial and boundary conditions. On output, it provides the primary unknown fields as well as other secondary fields of interest. The results can be exported in many ways to facilitate post-processing.

OOFEM is free software, distributed under GNU LGPL license.



General Features
----------------
* Structural, heat & mass transfer, and CFD modules
* Support for parallel processing on shared and massively parallelmachines, computer clusters, Dynamic load balancing
* Direct & iterative solvers (IML, Spooles, PETSc, SuperLU, Pardisso)
* Full restart capability, support for adaptive and staggered analyses
* Postprocessing: X-Windows post-processing, export to VTK
