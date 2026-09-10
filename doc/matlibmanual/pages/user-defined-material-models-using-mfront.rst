
.. _user-defined-material-models-using-mfront:

User-defined material models using MFront
=========================================

A user-defined material created with MFront can be used if MGIS is installed
(https://github.com/thelfer/MFrontGenericInterfaceSupport) and linked to OOFEM with cmake before the compilation.
Run cmake with following flags and specify the location of the MGIS cmake files:

.. code-block:: bash

   -DUSE_MFRONT=ON -DMFrontGenericInterface_DIR=/path/to/mgis/cmake/

Once OOFEM is compiled with this setting, the file /tests/smmfront/mfront01.in will be included in the tests.
It requires the file

/tests/smmfront/plasticityIsotropicLinearHardeningPlasticity.mfront to be processed by MFront before the tests are run.

See /tests/smmfront/readme_install.txt for more detailed instructions.

.. table:: MFrontUserMaterial material model -- summary.
   :name: MFrontUserMaterial_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | MFrontUserMaterial allows use of user-defined material models based on MFront                |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`MFrontUserMaterial` :elemparam:`d{rn}` :elemparam:`modelname{s}`                  |
   |                    | :elemparam:`libpath{s}` :elemparam:`properties{dc}`                                          |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`num` material model number                                                         |
   |                    | - :param:`d` specific weight                                                                 |
   |                    | - :param:`modelname` name of the material model from the shared library to be used           |
   |                    | - :param:`libpath` path to the shared library of the material created with MFront            |
   |                    | - :param:`properties` are not used now, which yields in "properties 0" on input for now      |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat                                                                                        |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `smmfront/mfront01.in                                                                        |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/smmfront/mfront01.in>`_          |
   +--------------------+----------------------------------------------------------------------------------------------+
