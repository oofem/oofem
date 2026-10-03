Getting started
===============

To get OOFEM installed on your system, please follow the :ref:`installation`
instructions first.

OOFEM is a console application; it should be executed from the command line
with arguments specifying the path to the input file. When no arguments are
provided, the help is printed on the output.

.. code-block:: bash

    $ ./oofem

.. code-block:: text

    ____________________________________________________
      ____  ____  __________  ___
     / __ \/ __ \/ __/ __/  |/  /___  _______
    / /_/ / /_/ / _// _// /|_/ // _ \/ __/ _ \
    \____/\____/_/ /___/_/  /_(_)___/_/  \_, /
    Copyright (C) 1994-2025 Borek Patzak /__/
    ____________________________________________________

    Options:

      -v  prints oofem version
      -f  (string) input file name
      -r  (int) restarts analysis from given step
      -ar (int) restarts adaptive analysis from given step
      -l  (int) sets treshold for log messages (Errors=0, Warnings=1,
                Relevant=2, Info=3, Debug=4)
      -rn turns on renumbering
      -qo (string) redirects the standard output stream to given file
      -qe (string) redirects the standard error stream to given file
      -c  creates context file for each solution step
      -p  runs in parallel mode (MPI)
      -t  (int) sets number of threads for OpenMP parallelization
      -m  shows solution status monitor output,
          redirecting standard output to files,
          oofem.stdout, oofem.stderr by default, use -qo -qe to override
      -Werror  promote warnings to errors

    Copyright (C) 1994-2025 Borek Patzak
    This is free software; see the source for copying conditions.  There is NO
    warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

The most useful options in daily work are

``-f <file>``
    The input file to solve. This is the only required option.

``-l <level>``
    How verbose the console log should be. ``-l 2`` keeps only the relevant
    messages, ``-l 4`` adds debugging detail.

``-qo <file>``, ``-qe <file>``
    Redirect the console log (standard output) and the error stream to files.
    Note that these control the *log*; the result output file is a separate
    thing, named in the input file itself (see :ref:`understanding_output`).

``-c`` and ``-r <step>``
    ``-c`` saves a context (restart) file after each solution step; ``-r``
    restarts a previously interrupted analysis from the given step. Both are
    needed for long running analyses.

``-t <n>``
    Number of OpenMP threads. Requires a build configured with
    ``USE_OPENMP_PARALLEL``.

``-p``
    Run as an MPI process of a distributed memory parallel analysis. Requires a
    build configured with ``USE_MPI_PARALLEL`` and a partitioned input.

``-Werror``
    Abort on the first warning instead of continuing. Useful when preparing a
    new model, because warnings about ignored or misspelled input fields then
    cannot pass unnoticed.

Running ``oofem -v`` prints the version together with the list of compiled-in
modules and the exact commit the binary was built from, which is what you
should quote when reporting a problem on the forum:

.. code-block:: bash

    $ ./oofem -v

.. code-block:: text

    OOFEM version 3.0 (x86_64-Linux, fm tm sm sm_elements sm_materials mpm IML++)
    Git RepoURL: https://github.com/oofem/oofem.git
        Branch: devel, Hash: v3.0

Running the first analysis
--------------------------
The problem to be solved is fully described in the input file. The structure of
the input is explained in the :ref:`understanding_input` section and fully
documented in the
`OOFEM Input manual <http://www.oofem.org/resources/doc/oofemInput/html/index.html>`_.

To run oofem with a specific input, use the ``-f`` option followed by the path
to the input file. For illustration, we demonstrate the execution using the
`beam2d_1.in <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/beam2d_1.in>`_
test, which is part of the OOFEM test suite and is located in the
``tests/regression/sm`` directory.

.. code-block:: bash

    $ ./oofem -f /home/user/oofem.src/tests/regression/sm/beam2d_1.in

.. code-block:: text

    ____________________________________________________
      ____  ____  __________  ___
     / __ \/ __ \/ __/ __/  |/  /___  _______
    / /_/ / /_/ / _// _// /|_/ // _ \/ __/ _ \
    \____/\____/_/ /___/_/  /_(_)___/_/  \_, /
    Copyright (C) 1994-2025 Borek Patzak /__/
    ____________________________________________________
    Job: /home/user/oofem.src/tests/regression/sm/beam2d_1.in
    Computing initial guess

    StaticStructural :: solveYourselfAt - Solving Metastep 1,  Step 1, Starting Time 0.000000e+00, Time Increment 1.000000e+00, Final Time 1.000000e+00,  (neq = 15)
    NRSolver: Iteration ForceError
    ----------------------------------------------------------------------------
    NRSolver: 0      D_u:  1.000e+00  D_w:  1.000e+00  R_v:  1.000e+00
    NRSolver: 1      D_u:  4.005e-06  D_w:  7.393e-06  R_v:  3.385e-16
    EngngModel info: user time consumed by solution step 1: 0.00s

    ...

    ANALYSIS FINISHED

    Real time consumed: 000h:00m:00s
    User time consumed: 000h:00m:00s
    Total 0 error(s) and 0 warning(s) reported

The lines starting with ``NRSolver:`` report the convergence of the equilibrium
iterations, one line per iteration, with the relative error of each unknown
group. The last line of the run is a summary you should always check: an
analysis that finishes with reported errors, or with warnings you did not
expect, has probably not solved the problem you intended.

By default, a text output file with the results is created (as specified in the
input file). In this case, the ``beam2d_1.out`` file is created in the current
working directory; its content is described in :ref:`understanding_output`.
Note that the paths in the input file, the output file name included, are
interpreted relative to the current working directory, not to the location of
the input file.

Where to go next
----------------
* :ref:`understanding_input` — how the input file is structured
* :ref:`understanding_output` — how to read the results
* :ref:`postprocessing` — how to get the results into Paraview
* :ref:`python` — how to drive the same analysis from Python
* The ``tests/regression`` directory of the source distribution contains
  several hundred small, working input files covering nearly every element,
  material and analysis type. Finding a test close to your problem and
  adapting it is usually the fastest way to build a new model.
