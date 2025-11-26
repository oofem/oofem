
First Steps
###########

This section covers the basic usage of oofem python bindings. 
It covers the installation, consisting of generating the python bindings and demonstrates the basic usage.

Installation
============

Prerequisities
--------------
The pybind11 library is included as submodule inside oofem repository. To build oofem with python bindings requires **python-dev** or **python3-dev** packages to be installed.

Generate binding code
---------------------
Configure oofem target to build python bindings

.. code-block:: bash

   cd oofem.build
   cmake path_to_oofem_git_repository -D USE_PYTHON_BINDINGS=ON
   make
   
Building the above project will produce a binary module file that can be imported to Python.

Basic usage
===========
In the previous section we have built the oofem python module. Assuming that the compiled module is located in the
current directory, the following interactive Python session shows how to
load and execute the example:

.. code-block:: pycon

    $ python3
    Python 2.7.10 (default, Aug 22 2015, 20:33:39)
    Python 3.6.8 (default, Oct  7 2019, 12:59:55) 
    [GCC 8.3.0] on linux
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import oofempy
    >>> dr=oofempy.OOFEMTXTDataReader("patch010.in")
    >>> problem=oofempy.InstanciateProblem(dr, oofempy.problemMode.processor, False, None, False)
    >>> problem.init()
    >>> problem.solveYourself()
    Computing initial guess
    StaticStructural :: solveYourselfAt - Solving step 1, metastep 1, (neq = 3)
    NRSolver: Iteration ForceError
    ----------------------------------------------------------------------------
    NRSolver: 0      D_u:  0.000e+00
    Checking rules...
    EngngModel info: user time consumed by solution step 1: 0.00s
    >>> problem.terminateAnalysis()
    ANALYSIS FINISHED
    Real time consumed: 000h:00m:45s
    User time consumed: 000h:00m:00s




