.. _understanding_output:

Understanding output file
=========================

By default, oofem produces output in the form of readable text file, called output file.
Its content can be controlled to some extend in ``OutputManager`` record of input file. One can filter output for selected solution steps, and specific sets of nodes and elements, for example. 

An Example
^^^^^^^^^^
Consider the same linear elastic analysis of beam structure, as in previous section.
We first run the solver with the example 

.. code-block:: bash

    $ ./oofem -f /home/user/oofem.git/tests/sm/beam2d_1.in 

Upon the successful execution, the text output file containing simulation results is created (the name of output file is specified in input file). In our case, the ``beam2d_1.out`` file is created in the current working directory.
The output file can be inspected by any text editor.

Header section
--------------
Every output file starts with Header containing information on solver version, job name and 
starting date and time of the analysis

::

    ############################################################## www.oofem.org ###
           ####    ####   ######  ######  #    #
         #    #  #    #  #       #       ##  ##
        #    #  #    #  #####   #####   # ## #
       #    #  #    #  #       #       #    #
      #    #  #    #  #       #       #    #   OOFEM ver. 2.6
      ####    ####   #       ######  #    #    Copyright (C) 1994-2020 Borek Patzak
    ################################################################################


    Starting analysis on: Wed Jun  9 09:32:17 2021

    Homework www sm40 no. 1

Solution step section(s)
------------------------
The output file continues with output for each solution step of the analysis. 
This consists of simple header indicating the solution step time

::

    ==============================================================
     Output for time 1.00000000e+00
    ==============================================================

The output for solution step consists of output of problem domain(s)

::

    Output for domain   1
    
consisting of output for all nodes and elements. We start with nodal output. 
The output for each node starts with ``Node`` keyword, followed by node label and number (in parenthesis).
This is followed by indented block containing the output for every degree of freedom of the node.
DOF meaning is identified by a integer code after ``dof`` keyword. The codes are defined by DofIDItem enum. 
For example, the mechanical unknowns have following codes: ``1`` for displacement in x direction, ``2`` for displacement in y direction, ``3`` for displacement in z direction, ``4`` for rotation around x axis, ``5`` for rotation around y axis and ``6`` for rotation around z axis (see src/oofemlib/dofiditem.h for full definition).
Our 2D beam structure is in xz plane, so the relevant DOFs are displacement in x, displacement in z and rotation around y, corresponding to DOF codes 1,3,5. 
The ``d`` code dof output means that actual value is printed, for some analyses, also velocities and accelerations of the corresponding unknown can be printed.

::

    DofManager output:
    ------------------
        Node           1 (       1):
            dof 1   d -1.37172495e-03
            dof 3   d  0.00000000e+00
            dof 5   d -2.38787802e-05
        Node           2 (       2):
            dof 1   d -1.37172495e-03
            dof 3   d  2.03123137e-14
            dof 5   d -1.01549259e-06
        Node           3 (       3):
            dof 1   d -1.37172495e-03
            dof 3   d  4.20823351e-05
            dof 5   d  0.00000000e+00
        Node           4 (       4):
            dof 1   d -1.75286359e-03
            dof 3   d  5.50267187e-04
            dof 5   d  1.05205838e-05
        Node           5 (       5):
            dof 1   d -1.34016320e-03
            dof 3   d  0.00000000e+00
            dof 5   d  4.07440098e-04
        Node           6 (       6):
            dof 1   d  0.00000000e+00
            dof 3   d  0.00000000e+00
            dof 5   d -0.00000000e+00

The nodal output is followed by element output. The actual format depends on particular
element, but generally internal variables at each integration point are
reported. In case of beam element, the local displacements and end forces are printed as well.

::

    Element output:
    ---------------
    beam element 1 (       1) :
    local displacements  -1.3717e-03 0.0000e+00 -2.3879e-05 -1.3717e-03 2.0312e-14 -1.0155e-06
    local end forces     0.0000e+00 -8.9375e+00 0.0000e+00 0.0000e+00 -1.5062e+01 -7.3499e+00
    GP  1.1  :  strains  0.0000e+00 0.0000e+00 0.0000e+00 3.6323e-05 0.0000e+00 0.0000e+00 -2.4500e-33 0.0000e+00
                stresses 0.0000e+00 0.0000e+00 0.0000e+00 4.2897e+00 0.0000e+00 0.0000e+00 -3.0625e+00 0.0000e+00
    GP  1.2  :  strains  0.0000e+00 0.0000e+00 0.0000e+00 2.0106e-05 0.0000e+00 0.0000e+00 -2.4500e-33 0.0000e+00
                stresses 0.0000e+00 0.0000e+00 0.0000e+00 2.3745e+00 0.0000e+00 0.0000e+00 -3.0625e+00 0.0000e+00
    GP  1.3  :  strains  0.0000e+00 0.0000e+00 0.0000e+00 -1.0531e-06 0.0000e+00 0.0000e+00 -2.4500e-33 0.0000e+00
                stresses 0.0000e+00 0.0000e+00 0.0000e+00 -1.2437e-01 0.0000e+00 0.0000e+00 -3.0625e+00 0.0000e+00
    GP  1.4  :  strains  0.0000e+00 0.0000e+00 0.0000e+00 -1.7270e-05 0.0000e+00 0.0000e+00 -2.4500e-33 0.0000e+00
                stresses 0.0000e+00 0.0000e+00 0.0000e+00 -2.0396e+00 0.0000e+00 0.0000e+00 -3.0625e+00 0.0000e+00
    ...

For structural analyses the reaction table is reported:

::

    R E A C T I O N S  O U T P U T:
    _______________________________

    Node        1 iDof  3 reaction -8.9375e+00    [bc-id: 1]
    Node        3 iDof  5 reaction  0.0000e+00    [bc-id: 2]
    Node        5 iDof  3 reaction -1.8750e+01    [bc-id: 1]
    Node        6 iDof  1 reaction  1.8000e+01    [bc-id: 3]
    Node        6 iDof  3 reaction -2.0312e+01    [bc-id: 3]
    Node        6 iDof  5 reaction -5.3999e+01    [bc-id: 3]

Finally, the solution time for every time step is reported.

::

  User time consumed by solution step 1: 0.004 [s]

The output is then repeated for each solution step of the problem.
Finally, the accumulated total solution time is reported.

::

  Finishing analysis on: Wed Jun  9 09:32:17 2021

  Real time consumed: 000h:00m:00s
  User time consumed: 000h:00m:00s

