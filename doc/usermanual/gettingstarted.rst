Getting started
===============

To get OOFEM installed on yout system, please follow :ref:`installation` instructions first.

OOFEM is console application, it should be executed from command line with arguments specifying the path to input file. When no arguments are privided, the help is printed on output.

.. code-block:: bash

    $ ./oofem

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

    Copyright (C) 1994-2017 Borek Patzak
    This is free software; see the source for copying conditions.  There is NO
    warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

The problem to be solved is fully described in input file. The structure of input is explained in :ref:`understanding_input` sectiona and fully documented in `OOFEM Input manual <http://www.oofem.org/resources/doc/oofemInput/html/oofemInput.html>`_.

To run oofem with specific input, use ``-f`` option followed by a system path to input file. For illustration, we demonstrate the execution using ``beam2d_1.in`` test which is a part of OOFEM test suite and is located in ``tests/sm`` directory.

.. code-block:: bash

    $ ./oofem -f /home/user/oofem.git/tests/sm/beam2d_1.in
        ____________________________________________________
                    OOFEM - Finite Element Solver
                Copyright (C) 1994-2017 Borek Patzak
        ____________________________________________________
        Computing initial guess
        
        StaticStructural :: solveYourselfAt - Solving step 1, metastep 1, (neq = 15)

        ...

        ANALYSIS FINISHED

        Real time consumed: 000h:00m:00s
        User time consumed: 000h:00m:00s
        Total 0 error(s) and 0 warning(s) reported

By default, the text output file with results is created (as specified in input file). In this case, the ``beam2d_1.out`` file is created in the current working directory.



