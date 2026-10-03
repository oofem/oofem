.. _understanding_output:

Understanding output file
=========================

By default, oofem produces output in the form of a readable text file, called
the output file. Its content can be controlled to some extent in the
``OutputManager`` record of the input file. One can filter the output for
selected solution steps and for specific sets of nodes and elements, for
example.

Do not confuse the output file with the console log. The console log reports
the progress of the solution (which step is being solved, how the equilibrium
iterations converge) and can be redirected with the ``-qo`` and ``-qe``
options; the output file contains the results and is named in the first line of
the input file.

Controlling the content
-----------------------
The ``OutputManager`` record accepts the following selectors:

``tstep_all``
    Output for all solution steps.

``tstep_step <n>``
    Output for every *n*-th solution step.

``tsteps_out <size> <list>``
    Output only for the explicitly listed solution steps.

``dofman_all`` / ``element_all``
    Output for all dof managers / all elements.

``dofman_output <size> <list>`` / ``element_output <size> <list>``
    Output only for the listed dof managers / elements. Ranges may be given in
    the ``{(first last)}`` form.

``dofman_except <size> <list>`` / ``element_except <size> <list>``
    Suppress the output for the listed dof managers / elements.

Restricting the output to the steps and components of interest is worth doing
for larger models — for a mesh of any size the complete output of every
integration point in every step runs into hundreds of megabytes and is not
readable anyway. In that situation an export module (see
:ref:`postprocessing`) is the better tool.

An Example
^^^^^^^^^^
Consider the same linear elastic analysis of the beam structure as in the
previous section. We first run the solver with the example

.. code-block:: bash

    $ ./oofem -f /home/user/oofem.src/tests/regression/sm/beam2d_1.in

Upon successful execution, the text output file containing the simulation
results is created (the name of the output file is specified in the input
file). In our case, the ``beam2d_1.out`` file is created in the current working
directory. The output file can be inspected in any text editor.

Header section
--------------
Every output file starts with a header containing information on the solver
version, the repository and commit the solver was built from, the job
description and the starting date and time of the analysis.

::

    ############################################################## www.oofem.org ###
      ____  ____  __________  ___
     / __ \/ __ \/ __/ __/  |/  /___  _______
    / /_/ / /_/ / _// _// /|_/ // _ \/ __/ _ \
    \____/\____/_/ /___/_/  /_(_)___/_/  \_, /  OOFEM ver. 3.0
                                        /___/   Copyright (C) 1994-2025 Borek Patzak
    ################################################################################
    Git URL:https://github.com/oofem/oofem.git,
        branch:devel, hash:v3.0


    Starting analysis on: Wed Jun  9 09:32:17 2021

    Homework www sm40 no. 1

The ``Git URL`` block identifies the exact sources the results were produced
with, which makes an old output file traceable; quote it when reporting a
problem. The last line is the job description record copied verbatim from the
input file.

Solution step section(s)
------------------------
The output file continues with the output for each solution step of the
analysis. This consists of a simple header indicating the solution step time

::

    ==============================================================
    Output for time 1.00000000e+00
    ==============================================================

The output for a solution step consists of the output of the problem domain(s)

::

    Output for domain   1

consisting of the output for all nodes and elements. We start with the nodal
output. The output for each node starts with the ``Node`` keyword, followed by
the node label and the internal number (in parentheses). This is followed by an
indented block containing the output for every degree of freedom of the node.
The DOF meaning is identified by an integer code after the ``dof`` keyword. The
codes are defined by the ``DofIDItem`` enum. For example, the mechanical
unknowns have the following codes: ``1`` for displacement in the x direction,
``2`` for displacement in the y direction, ``3`` for displacement in the z
direction, ``4`` for rotation around the x axis, ``5`` for rotation around the
y axis and ``6`` for rotation around the z axis (see
`src/core/dofiditem.h <https://github.com/oofem/oofem/blob/devel/src/core/dofiditem.h>`_
for the full definition). Our 2D beam structure is in the xz plane, so the
relevant DOFs are displacement in x, displacement in z and rotation around y,
corresponding to DOF codes 1, 3 and 5. The ``d`` code in the dof output means
that the actual value is printed; for some analyses the velocities (``v``) and
accelerations (``a``) of the corresponding unknown can be printed as well.

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

The nodal output is followed by the element output. The actual format depends
on the particular element, but generally the internal variables at each
integration point are reported, prefixed by ``GP`` and the integration rule and
point numbers. In the case of a beam element, the local displacements and end
forces are printed as well.

::

    Element output:
    ---------------
    beam element 1 (       1) :
      local displacements  -1.3717e-03 0.0000e+00 -2.3879e-05 -1.3717e-03 2.0312e-14 -1.0155e-06
      local end forces     0.0000e+00 -8.9375e+00 0.0000e+00 0.0000e+00 -1.5062e+01 -7.3499e+00
      GP  1.1  :  strains  +0.0000e+00 +0.0000e+00 +0.0000e+00 +3.6323e-05 +0.0000e+00 +0.0000e+00 -2.4500e-33 +0.0000e+00
                  stresses +0.0000e+00 +0.0000e+00 +0.0000e+00 +4.2897e+00 +0.0000e+00 +0.0000e+00 -3.0625e+00 +0.0000e+00
      GP  1.2  :  strains  +0.0000e+00 +0.0000e+00 +0.0000e+00 +2.0106e-05 +0.0000e+00 +0.0000e+00 -2.4500e-33 +0.0000e+00
                  stresses +0.0000e+00 +0.0000e+00 +0.0000e+00 +2.3745e+00 +0.0000e+00 +0.0000e+00 -3.0625e+00 +0.0000e+00
      GP  1.3  :  strains  +0.0000e+00 +0.0000e+00 +0.0000e+00 -1.0531e-06 +0.0000e+00 +0.0000e+00 -2.4500e-33 +0.0000e+00
                  stresses +0.0000e+00 +0.0000e+00 +0.0000e+00 -1.2437e-01 +0.0000e+00 +0.0000e+00 -3.0625e+00 +0.0000e+00
      GP  1.4  :  strains  +0.0000e+00 +0.0000e+00 +0.0000e+00 -1.7270e-05 +0.0000e+00 +0.0000e+00 -2.4500e-33 +0.0000e+00
                  stresses +0.0000e+00 +0.0000e+00 +0.0000e+00 -2.0396e+00 +0.0000e+00 +0.0000e+00 -3.0625e+00 +0.0000e+00
    ...

The stress and strain components are reported in the element's own local
coordinate system; for a beam element these are generalized quantities (forces
and moments rather than stresses). Their number and ordering follow from the
stress-strain mode of the element and are documented per element type in the
`OOFEM Element Library Manual <http://www.oofem.org/resources/doc/elementlibmanual/html/elementlibmanual.html>`_.
The number of ``GP`` blocks follows from the element type and its integration
rule.

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

Only the DOFs constrained by a Dirichlet boundary condition appear here, and
``[bc-id: n]`` identifies which boundary condition of the input file is
responsible. Checking that the reactions sum up to the applied load is the
quickest sanity check on a new model.

Finally, the solution time for every time step is reported.

::

  User time consumed by solution step 1: 0.004 [s]

The output is then repeated for each solution step of the problem. Finally, the
accumulated total solution time is reported.

::

  Finishing analysis on: Wed Jun  9 09:32:17 2021

  Real time consumed: 000h:00m:00s
  User time consumed: 000h:00m:00s

If the output file ends without this closing block, the analysis did not finish
— it either diverged or was interrupted — and the results of the last step
present in the file should not be trusted.
