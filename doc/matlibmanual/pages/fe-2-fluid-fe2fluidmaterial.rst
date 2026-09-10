FE 2 fluid - FE2FluidMaterial
=============================


.. _FE2FluidMaterial:

FE² fluid - FE2FluidMaterial
----------------------------

Constitutive model of multiscale fluids.
The macroscale fluid behavior is determined by a Representative Volume Element (RVE) which is solved for in each integration point.
Some requirements are put on the RVE, such as the it must use the MixedGradientPressure boundary condition along its outer boundary. This boundary condition is equivalent to that of Dirichlet boundary condition in classical homogenization, only adjusted for a mixed control; deviatoric gradient + pressure, instead of the full gradient only.
Worth noting is that the macroscale (and microscale) behavior can still be compressible,  in particular is the RVE contains pores.
The triangular and tetrahedral Stokes' flow elements both support compressible behavior, but it perhaps be applicable other types of flow, as long as the extended mixed control option is used when calling the material routine (the function with gradient and pressure as input). If the RVE doesn't contain any pores (such that the response turns out to be incompressible) then material should work for all problem classes.
This could for example be the case of a oil-water mixture.
The model parameters are summarized in :numref:`FE2FluidMaterial_table`.

.. table:: FE² fluid material - summary.
   :name: FE2FluidMaterial_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | FE\ :sup:`2` Fluid material                                                                  |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`FE2FluidMaterial` :elemparam:`num{in}` :elemparam:`d{rn}`                         |
   |                    | :elemparam:`inputfile{s}`                                                                    |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`num` material model number                                                         |
   |                    | - :param:`d` unused                                                                          |
   |                    | - :param:`inputfile` input file for RVE problem                                              |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 2d, 3d flow                                                                                  |
   +--------------------+----------------------------------------------------------------------------------------------+

.. clearpage
