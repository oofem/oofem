
.. _sec-coupled-heat-mass-transport:

Coupled heat and mass transport material model - HeMoKunzel
===========================================================

The presented formulation is based on the work of Kuenzel [Kuenzel]_. 
The model is suitable for problems with dominating water diffusion
and negligible water convection. The governing equations for temperature and humidity reads

.. math::
   :label: eq-governing-heat

   \frac{\partial Q}{\partial t} &=& \frac{\partial Q}{\partial T} \frac{\partial T}{\partial t} = C_v \frac{\partial T}{\partial t} = -\nabla q_T = \nabla \left( \lambda \nabla T \right ) + h_v \nabla \left( \delta_p \nabla (H p_{sat}) \right )\\
   \frac{\partial w}{\partial t} &=& \frac{\partial w}{\partial H} \frac{\partial H}{\partial t} = -\nabla q_H = \nabla \left( D_H \nabla H + \delta_p \nabla (H p_{sat}) \right )

.. table:: Parameters from Kunzel's model.
   :name: tab-parameters-kunzel

   +--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | Symbol                                                                                                                                                                                                       | Unit                                                                                                                               | Description                                                                                                                                                                                                                                                                                                 |
   +==============================================================================================================================================================================================================+====================================================================================================================================+=============================================================================================================================================================================================================================================================================================================+
   | :math:`T` :math:`H`                                                                                                                                                                                          | (K) (-)                                                                                                                            | Temperature Relative humidity 0-1                                                                                                                                                                                                                                                                           |
   +--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | :math:`\frac{\partial Q}{\partial T} \approx C_v` :math:`\frac{\partial w}{\partial H}` :math:`Q` :math:`q_T` :math:`\lambda` :math:`h_v` :math:`\delta_p` :math:`p_{sat}` :math:`w` :math:`D_H` :math:`D_w` | J/K/m\ :math:`^3` kg/m\ :math:`^3` J/m\ :math:`^3` W/m\ :math:`^2` W/m/K J/kg kg/m/s/Pa Pa kg/m\ :math:`^3` kg/m/s m\ :math:`^2`/s | Heat storage capacity per volume Moisture storage capacity - sorption isotherm Total amount of heat in unit volume Heat flux Thermal conductivity Evaporation enthalpy of water Water vapour permeability Water vapour saturation pressure Moisture content Liquid conduction coefficient Water diffusivity |
   +--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

Numerical solution leads to the system of equations

.. math::
   :label: final_eq

   \left[ \begin{array}{cc}
   \tenss{K}_{TT} & \tenss{K}_{TH} \\
   \tenss{K}_{HT} & \tenss{K}_{HH}
   \end{array} \right]
   \left\{ \begin{array}{c}
   \tenss{r}_T \\
   \tenss{r}_H
   \end{array} \right\} + 
   \left[ \begin{array}{cc}
   \tenss{C}_{TT} & \tenss{C}_{TH} \\
   \tenss{C}_{HT} & \tenss{C}_{HH}
   \end{array} \right]
   \left\{ \begin{array}{c}
   \dot{\tenss{r}}_T \\
   \dot{\tenss{r}}_H
   \end{array} \right\} = 
   \left\{ \begin{array}{c}
   \tenss{q}_{T} \\
   \tenss{q}_{H}
   \end{array} \right\},

where

.. math::

   \tenss{K}_{TT} = \int_{\Omega} \tenss{B}^T k_{TT} \tenss{B}{\rm d}\Omega,\qquad
   \tenss{K}_{TH} = \int_{\Omega} \tenss{B}^T k_{Th} \tenss{B}{\rm d}\Omega,\qquad
   \tenss{K}_{HT} = \int_{\Omega} \tenss{B}^T k_{hT} \tenss{B}{\rm d}\Omega,\qquad
   \tenss{K}_{HH} = \int_{\Omega} \tenss{B}^T k_{hh} \tenss{B}{\rm d}\Omega,\qquad

.. math::

   \tenss{C}_{TT} = \int_{\Omega} \tenss{N}^T c_{TT}  \tenss{N} {\rm d}\Omega,\qquad
   \tenss{C}_{TH} = \int_{\Omega} \tenss{N}^T c_{Th}  \tenss{N} {\rm d}\Omega,\qquad
   \tenss{C}_{HT} = \int_{\Omega} \tenss{N}^T c_{hT}  \tenss{N} {\rm d}\Omega,\qquad
   \tenss{C}_{HH} = \int_{\Omega} \tenss{N}^T c_{hh}  \tenss{N} {\rm d}\Omega,\qquad

.. math::

   \tenss{q}_T = \int_{\Gamma_2} \tenss{N}^T  \overline{q}_{T}{\rm d}\Gamma,\qquad
   \tenss{q}_H = \int_{\Gamma_2} \tenss{N}^T  \overline{q}^{h}{\rm d}\Gamma,\qquad

where

.. math::
   :label: k_TT

   k_{TT} &=& \lambda(w) + h_v \cdot \delta_p(T) \cdot H \cdot  \frac{\Delta p_{sat}}{\Delta T}(T),\\
   k_{TH} &=& h_v \cdot \delta_p(T) \cdot p_{sat}(T),\\
   k_{HT} &=& \delta_p(T) \cdot H \cdot \frac{\Delta p_{sat}}{\Delta T}(T),\\
   k_{HH} &=& D_w(H) \cdot \frac{\Delta w}{\Delta H}(H) + \delta_p(T) \cdot p_{sat}(T),\\
   c_{TT} &=& C_s \cdot \rho + C_w \cdot w,\\
   c_{TH} &=& 0,\\
   c_{HT} &=& 0,\\
   c_{HH} &=& \frac{\Delta w}{\Delta H}(H).

Note, that conductivity matrix \tenss{K} is unsymmetric hence unsymmetric matrix storage needs to be used (smtype).

The model parameters are summarized in :numref:`hemokunzel_table`.

.. table:: Coupled heat and mass transfer material model Kunzel - summary.
   :name: hemokunzel_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Coupled heat and mass transfer material model                                                |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`HeMoKunzel` :elemparam:`num{in}` :elemparam:`d{rn}` :elemparam:`iso_type{in}`     |
   |                    | :elemparam:`iso_wh{rn}` :elemparam:`mu{rn}` :elemparam:`permeability_type{in}`               |
   |                    | :elemparam:`A{rn}` :elemparam:`lambda0{rn}` :elemparam:`b{rn}` :elemparam:`cs{rn}`           |
   |                    | :optelemparam:`pl{rn}` :optelemparam:`rhoH2O{rn}` :optelemparam:`cw{rn}`                     |
   |                    | :optelemparam:`hv{rn}`                                                                       |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`num` material model number                                                         |
   |                    | - :param:`d` bulk density of dry building material [kg/m\ :math:`^3`]                        |
   |                    | - :param:`iso_type`\ =0 is isotherm from Hansen needing :param:`iso_n`, :param:`iso_a`, =1   |
   |                    |   is Kunzel which needs :param:`iso_b`                                                       |
   |                    | - :param:`iso_wh` maximum adsorbed water content [kg/m\ :math:`^3`]                          |
   |                    | - :param:`mu` water vapor diffusion resistance [-]                                           |
   |                    | - :param:`permeability_type` =0 is Multilin_h needing :param:`perm_h`, :param:`perm_Dw(h)`,  |
   |                    |   =1 is Multilin_wV needs :param:`perm_wV`, :param:`perm_DwwV`, =2 is Kunzelperm needs       |
   |                    |   :param:`A` as water absorption coefficient [kg/m/s\ :math:`^0.5`]                          |
   |                    | - :param:`lambda0`\ rn thermal conductivity [W/m/K]                                          |
   |                    | - :param:`b` thermal conductivity supplement [-]                                             |
   |                    | - :param:`cs` specific heat capacity of the building material [J/kg/K]                       |
   |                    | - :optparam:`pl` ambient air pressure [Pa], default = 101325                                 |
   |                    | - :optparam:`rhoH2O` water density [kg/m3], default = 1000                                   |
   |                    | - :optparam:`cw` specific heat capacity of liquid water, default = 4183                      |
   |                    | - :optparam:`hv` latent heat of water phase change [J/kg], default = 2.5e+6                  |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | _2dHeMo, _3dHeMo                                                                             |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `tm/HeMoKunzel_1.in                                                                          |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/HeMoKunzel_1.in>`_            |
   +--------------------+----------------------------------------------------------------------------------------------+
