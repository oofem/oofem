.. _InitModulesSec:

Initialization modules
======================

Initialization modules initialise the state variables from data computed
previously by external software.  The number of initialization module records is
declared by :param:`ninitmodules` in the analysis record; see
:ref:`AnalysisRecord`.

.. record:: Record syntax common to all initialization modules

   :descitem:`EntType` :elemparam:`initfile{s}`

**Parameters**

:param:`initfile{s}`
    Path to the file containing the initialization data.

.. note::

   Unlike most string attributes, this file name is given **without** quotes.

Gauss point initialization module
---------------------------------

Currently the only initialization module.

.. record::

   :descitem:`GPInitModule` :elemparam:`initfile{s}`

Format of the initialization file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Each Gauss point is represented by one line.
* The Gauss points must be given in a specific order, by element number and
  then Gauss point number, matching the mesh specified in the later sections of
  the input file.
* Each line referring to a Gauss point contains:

  .. code-block:: none

     elnum gpnum coords ng var_1_id values_1 ... var_ng_id values_ng

``elnum`` (in)
    Element number.

``gpnum`` (in)
    Gauss point number.

``coords`` (ra)
    Coordinates of the Gauss point.

``ng`` (in)
    Number of groups of variables that follow.

``var_i_id`` (in)
    Identification number of variable group *i*, according to the
    ``InternalStateType`` enumeration in ``src/core/internalstatetype.h``.

``values_i`` (ra)
    Values of the variables in group *i*.

Example
~~~~~~~

.. code-block:: none

   37 4 3 0.02 0.04 0.05 3 52 1 0.23 62 1 0.049 1 6 0 -2.08e+07 0 0 0 0

This says that Gauss point number 4 of element number 37 has coordinates
:math:`x=0.02`, :math:`y=0.04`, :math:`z=0.05`, and that initial values are
given for 3 groups of variables:

* the first group, variable ID 52, is of type ``IST_DamageScalar`` and holds 1
  variable — it is a scalar — with the value 0.23;
* the second group, ID 62, is of type ``IST_CumPlasticStrain`` and holds 1
  variable with the value 0.049;
* the third group, ID 1, is of type ``IST_StressTensor`` and holds 6 variables —
  the stress components :math:`\sigma_x`, :math:`\sigma_y` and so on — with the
  values 0, -2.08e+07, 0, 0, 0, 0.
