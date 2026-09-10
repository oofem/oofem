MPM Symbolic Elements (MPM Module)
==================================


In this section we describe the elements supporting the Symbolic mode of Multiphysics
Module (MPM). The symbolic mode is used to define the weak form of the problem in terms
of variables, terms and integrals. The symbolic MPM elements support different types of
approximations. The additional nodes are generated authomatically and this process
ensures the reauired continuity of interpolations among element boundaries. Also, the
needed DOFs are generated automatically, depending on the configured variables.

.. toctree::
   :maxdepth: 2

   q1-element
   l1-element
