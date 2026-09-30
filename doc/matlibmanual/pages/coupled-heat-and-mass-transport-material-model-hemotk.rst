
.. _hemotk:

Coupled heat and mass transport material model - HeMotk
=======================================================

Coupled heat and mass transfer material model.
Source: T. Krejci doctoral thesis; Bazant and Najjar, 1972;
Pedersen, 1990. Assumptions: water vapor is the only driving
mechanism; relative humidity is from range 0.2 - 0.98 (I and II
regions). The model parameters are summarized
in :numref:`hemotk_table`.

.. table:: Coupled heat and mass transfer material model - summary.
   :name: hemotk_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Coupled heat and mass transfer material model                                                |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`HeMotk` :elemparam:`num{in}` :elemparam:`d{rn}` :elemparam:`a_0{rn}`              |
   |                    | :elemparam:`nn{rn}` :elemparam:`phi_c{rn}` :elemparam:`delta_wet{rn}` :elemparam:`w_h{rn}`   |
   |                    | :elemparam:`n{rn}` :elemparam:`a{rn}` :elemparam:`latent{rn}` :elemparam:`c{rn}`             |
   |                    | :elemparam:`rho{rn}` :elemparam:`chi_eff{rn}` :elemparam:`por{rn}` :elemparam:`rho_gws{rn}`  |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`num` material model number                                                         |
   |                    | - :param:`d`, :param:`rho` material density                                                  |
   |                    | - :param:`a_0` constant (obtained from experiments) a_0 [Bazant and Najjar, 1972]            |
   |                    | - :param:`nn` constant-exponent (obtained from experiments) n [Bazant and Najjar, 1972]      |
   |                    | - :param:`phi_c` constant-relative humidity (obtained from experiments) phi_c [Bazant and    |
   |                    |   Najjar, 1972]                                                                              |
   |                    | - :param:`delta_wet` constant-water vapor permeability (obtained from experiments) delta_wet |
   |                    |   [Bazant and Najjar, 1972]                                                                  |
   |                    | - :param:`w_h` constant water content (obtained from experiments) w_h [Pedersen, 1990]       |
   |                    | - :param:`n` constant-exponent (obtained from experiments) n [Pedersen, 1990]                |
   |                    | - :param:`a` constant (obtained from experiments) A [Pedersen, 1990]                         |
   |                    | - :param:`latent` latent heat of evaporation                                                 |
   |                    | - :param:`c` thermal capacity                                                                |
   |                    | - :param:`chi_eff` effective thermal conductivity                                            |
   |                    | - :param:`por` porosity                                                                      |
   |                    | - :param:`rho_gws` saturation volume density                                                 |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | _2dHeMo                                                                                      |
   +--------------------+----------------------------------------------------------------------------------------------+

.. clearpage
