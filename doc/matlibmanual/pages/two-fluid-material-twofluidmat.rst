
.. _TwoFluidMaterial:

Two-fluid material - TwoFluidMat
================================

Material coupling the behaviour of two particular materials based on
rule of mixture. The weighting factor is VOF fraction.
The model parameters are summarized in :numref:`TwoFluidMaterial_table`.

.. table:: Two-Fluid material - summary.
   :name: TwoFluidMaterial_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Two-Fluid material                                                                           |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`TwoFluidMat` :elemparam:`num{in}` :elemparam:`mat{ia}`                            |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`num` material model number                                                         |
   |                    | - :param:`mat` integer array containing two numbers representing numbers of material models  |
   |                    |   of which the receiver is composed. Material with index 0 is a material, that is fully      |
   |                    |   active in a cell with VOF=0, material with index 1 is a material fully active in a cell    |
   |                    |   with VOF=1.                                                                                |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 2d flow                                                                                      |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `benchmark/fm/bdam7.oofem.in                                                                 |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/benchmark/fm/bdam7.oofem.in>`_,  |
   |                    | `fm/scctest01.in                                                                             |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/fm/scctest01.in>`_               |
   +--------------------+----------------------------------------------------------------------------------------------+
