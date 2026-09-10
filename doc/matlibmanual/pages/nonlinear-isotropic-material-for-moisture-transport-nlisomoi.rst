Nonlinear isotropic material for moisture transport - NlIsoMoisture
===================================================================



.. _sec:NlIsoMoistureMat:

Nonlinear isotropic material for moisture transport -- NlIsoMoisture
--------------------------------------------------------------------

This is a more general model for nonlinear moisture transport in isotropic 
porous materials, based on a nonlinear sorption isotherm (relation between 
the pore relative humidity :math:`h` and the water content :math:`w`) 
and on a humidity-dependent
moisture permeability.
The governing differential equation solves water mass balance in a unit volume [kg/m\ :sup:`3`/s] and reads

.. math::
   :label: nlisomoisture:governing

   k(h) \frac{\partial h}{\partial t} = \nabla \cdot \left[ c(h) \nabla h \right] + w_n \frac{\partial \alpha}{\partial t}

where 
:math:`k(h)` [kg/m\ :sup:`3`] is the humidity-dependent moisture capacity (:math:`k(h)=\frac{\partial w}{\partial h}` which is derivative of the moisture content :math:`w(h)` [kg/m\ :sup:`3`] with respect to the relative humidity), :math:`c(h)`
[kg/m/s] is the moisture permeability and the sink term 
:math:`w_n \frac{\partial \alpha}{\partial t}` corresponds to non-evaporable water loss due to hydration. :math:`w_n` is
non-evaporable water content for complete hydration per m\ :sup:`3` [kg/m\ :sup:`3`] and :math:`\alpha` is degree of hydration [-]. 
For the majority of cements, 1 kg of cement consumes approximately 0.23 kg of non-evaporable water at complete hydration.

So far, six different functions for the **sorption
isotherm** have been implemented  (in fact, what matters for the model
is not the isotherm itself but its derivative---the moisture capacity):

.. _nlisomoisture:

Isotherms
~~~~~~~~~

- Isotherm proposed by **Kuenzel** [Kuenzel]_ (``isothermType=3``) 
  in the form

.. math::
   :label: nlisomoisture:kuenzel

   w(h) = w_f \frac{(b-1)h}{b-h}

where :math:`w_f` [kg/m\ :sup:`3`] is the moisture content at free saturation and
:math:`b` is a dimensionless approximation factor greater than 1.

- Isotherm proposed by **Hansen** [Hansen]_ (``isothermType=4``) in the form

.. math::
   :label: nlisomoisture:hansen

   u(h) = u_h \left(1- \frac{\ln h}{A} \right)^{-1/n}

characterizes the amount of adsorbed water 
by the moisture ratio :math:`u` [kg/kg]. To obtain the
moisture content :math:`w`, it is necessary to multiply the moisture ratio
by the density of the solid phase. In (:eq:`nlisomoisture:hansen`),
:math:`u_h` is the maximum hygroscopically
bound water by adsorption, and :math:`A` and :math:`n` are constants
obtained by fitting of experimental data.

- The **BSB** isotherm [BSB]_ (``isothermType=5``) is an improved
  version of the famous BET isotherm. It is expressed in terms
  of the moisture ratio

.. math::
   :label: nlisomoisture:BSB

   u(h) = \frac{C k V_m h}{(1-k h)(1+(C-1)k h)}

where :math:`V_m` is the monolayer capacity, and :math:`C` depends on the absolute
temperature :math:`T` and on the difference between the heat of adsorption and
condensation. Empirical formulae for estimation of the parameters can
be found in [Xi]_. Note that these formulae hold quite accurately
for cement paste only; a reduction of the moisture ratio is necessary 
if the isotherm
should be applied for concrete.

- **Bilinear** isotherm (``isothermType=6``) is defined by its
  moisture capacity :emph:`capa` for relative humidity less than
  :emph:`hx`. Parameter :emph:`wf` defines moisture content at full
  saturation, :math:`h = 1`. Convergence is substantially improved by a smooth transition on
  interval :math:`(hx-dx, hx+dx)`.
  Similarly to the linear isotherm the moisture content can be
  adjusted by parameter :emph:`isooffset`.

The present implementation covers three functions for 
**moisture permeability**:


Moisture Permeability
~~~~~~~~~~~~~~~~~~~~~

- **Piecewise linear** permeability  (``permeabilityType=0``) 
  is defined by two arrays
  with the values of pore relative humidity :math:`perm\_h` and the
  corresponding values of moisture content :math:`perm\_c(h)`. The arrays must be of
  the same size.

The **Bažant-Najjar** permeability function (:math:`permeabilityType=1`) is given by the same formula (:eq:`BNmodel:diffusivity`) as the diffusivity in 
Section :ref:`sec:BazantNajjarMoistureMat`. All parameters have a similar
meaning as in (:eq:`BNmodel:diffusivity`) but :math:`c1` is now the
moisture permeability at full saturation [kg/m\ :math:`\cdot`\ s].

The permeability function proposed by **Xi et al.** [Xi]_ (:math:`permeabilityType=2`) reads

.. math::
   :label: nlisomoisture:Xi

   c(h) = \alpha_h + \beta_h \left[ 1 - 2^{-10^{\gamma_h(h-1)}} \right]

where :math:`\alpha_h`, :math:`\beta_h` and :math:`\gamma_h` are parameters that can be
evaluated using empirical mixture-based formulae presented in [Xi]_. 
However, if those formulae are used outside the range of water-cement 
ratios for which they were calibrated, the
permeability may become negative. Also the physical units are unclear.

Note that the Bajant-Najjar 
model from Section :ref:`sec:BazantNajjarMoistureMat` 
can be obtained as a special
case of the present model if :math:`permeabilityType` is set to 1
and :math:`isothermType` is set to 0. 
The ratio :math:`c1/moistureCapacity` then corresponds to the diffusivity
parameter :math:`C_1` from :eq:`BNmodel:governing`.

The model parameters are summarized in :numref:`NlIsoMoistureMat`.

.. table:: Nonlinear isotropic material for moisture transport
   :name: NlIsoMoistureMat

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Nonlinear isotropic material for moisture transport                                          |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`NlIsoMoistureMat` :elemparam:`num{in}` :elemparam:`d{rn}`                         |
   |                    | :elemparam:`isothermType{in}` :elemparam:`permeabilityType{in}` :optelemparam:`rhodry{rn}`   |
   |                    | :optelemparam:`capa{rn}` :optelemparam:`iso_h{ra}` :optelemparam:`iso_w(h){ra}`              |
   |                    | :optelemparam:`dd{rn}` :optelemparam:`wf{rn}` :optelemparam:`b{rn}` :optelemparam:`uh{rn}`   |
   |                    | :optelemparam:`A{rn}` :optelemparam:`nn{rn}` :optelemparam:`c{rn}` :optelemparam:`k{rn}`     |
   |                    | :optelemparam:`Vm{rn}` :optelemparam:`hx{rn}` :optelemparam:`dx{rn}`                         |
   |                    | :optelemparam:`perm_h{ra}` :optelemparam:`perm_c(h){ra}` :optelemparam:`c1{rn}`              |
   |                    | :optelemparam:`n{rn}` :optelemparam:`alpha0{rn}` :optelemparam:`hc{rn}`                      |
   |                    | :optelemparam:`alphah{rn}` :optelemparam:`betah{rn}` :optelemparam:`gammah{rn}`              |
   |                    | :optelemparam:`wn{rn}` :optelemparam:`alpha{expr}`                                           |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`num` material model number                                                         |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`isothermType` isotherm function as listed above (0, 1, ...6)                       |
   |                    | - :param:`permeabilityType` moisture permeability function as listed above (0, 1, 2)         |
   |                    | - :param:`rhodry` [kg/m\ :math:`^3`] density of dry material (for :math:`isothermType=4` and |
   |                    |   5)                                                                                         |
   |                    | - :param:`capa` [kg/m\ :math:`^3`] moisture capacity (for :math:`isothermType=0` and 6)      |
   |                    | - :param:`isooffset` [kg/m\ :math:`^3`] moisture capacity at zero relative humidity (for     |
   |                    |   :math:`isothermType=0` and 6)                                                              |
   |                    | - :param:`iso_h` [-] humidity array (for :math:`isothermType=1`)                             |
   |                    | - :param:`iso_w(h)` [kg/m\ :math:`^3`] moisture content array (for :math:`isothermType=1`)   |
   |                    | - :param:`dd` [-] parameter (for :math:`isothermType=2`)                                     |
   |                    | - :param:`wf` [kg/m\ :math:`^3`] is the moisture content at free saturation (for             |
   |                    |   :math:`isothermType=3` and 6)                                                              |
   |                    | - :param:`b` [-] parameter (for :math:`isothermType=3`)                                      |
   |                    | - :param:`uh` [kg/kg] maximum hygroscopically bound water by adsorption (for                 |
   |                    |   :math:`isothermType=4`)                                                                    |
   |                    | - :param:`A` [-] parameter (for :math:`isothermType=4`)                                      |
   |                    | - :param:`n` [-] parameter (for :math:`isothermType=4`)                                      |
   |                    | - :param:`Vm` (for :math:`isothermType=5`)                                                   |
   |                    | - :param:`k` (for :math:`isothermType=5`)                                                    |
   |                    | - :param:`C` (for :math:`isothermType=5`)                                                    |
   |                    | - :param:`hx` [-] transition relative humidity (for :math:`isothermType=6`)                  |
   |                    | - :param:`dx` [-] length of the smooth transition (for :math:`isothermType=6`)               |
   |                    | - :param:`perm_c(h)` [kg m\ :math:`^{-1}` s\ :math:`^{-1}`] moisture permeability array (for |
   |                    |   :math:`permeabilityType=0`)                                                                |
   |                    | - :param:`c1` [kg m\ :math:`^{-1}` s\ :math:`^{-1}`] moisture permeability at full           |
   |                    |   saturation (for :math:`permeabilityType=1`)                                                |
   |                    | - :param:`n` [-] exponent (for :math:`permeabilityType=1`)                                   |
   |                    | - :param:`alpha0` [-] ratio between minimum and maximum diffusivity (for                     |
   |                    |   :math:`permeabilityType=1`)                                                                |
   |                    | - :param:`hc` [-] relative humidity at which the diffusivity is exactly between its minimum  |
   |                    |   and maximum value (for :math:`permeabilityType=1`)                                         |
   |                    | - :param:`alphah` [kg m\ :math:`^{-1}` s\ :math:`^{-1}`] (for :math:`permeabilityType=2`)    |
   |                    | - :param:`betah` [kg m\ :math:`^{-1}` s\ :math:`^{-1}`] (for :math:`permeabilityType=2`)     |
   |                    | - :param:`gammah` [-] (for :math:`permeabilityType=2`)                                       |
   |                    | - :param:`wn` [kg m\ :math:`^{-3}`] nonevaporable water content per m\ :math:`^3` of         |
   |                    |   concrete, default 0.23 kg/kg of cement                                                     |
   |                    | - :param:`alpha` [-] function of degree of hydration                                         |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | _2dHeat                                                                                      |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `tm/nlisomoisture01.in                                                                       |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/nlisomoisture01.in>`_,        |
   |                    | `tm/nlisomoisture02.in                                                                       |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/nlisomoisture02.in>`_,        |
   |                    | `tm/nlisomoisture03.in                                                                       |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/nlisomoisture03.in>`_         |
   +--------------------+----------------------------------------------------------------------------------------------+


.. _Cemhyd:

CemhydMat
~~~~~~~~~

CemhydMat represents a hydrating material based on CEMHYD3D model version 3.0,
developed at NIST [NISTIR7232]_. The model represents a digital hydrating
microstructure, driven with cellular automata rules and combined with cement
chemistry. Ordinary Portland cement is treated without any difficulties, blended
cements are usually decomposed into hydrating Portland contribution and intert
secondary cementitious material. The microstructure size can be from
:math:`10\times10\times10` to over :math:`200\times200\times200` :math:`\mu\text{m}`. For standard
computations the size :math:`50\times50\times50` suffices.

Each material instance creates an independent microstructure. It is also
possible to enforce having different microstructures in each integration point.
The hydrating model is coupled with temperature and averaging over shared
elements within one material instance occurs during the solution. Such approach
allows domain partitioning to many CemhydMat instances, depending on expected
accuracy or computational speed. A more detailed description with engineering
examples was published [Smilauer:09]_. :numref:`Cemhydmat_table` summarizes input parameters.

.. table:: Cemhydmat - summary
   :name: Cemhydmat_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Cemhyd - hydrating material                                                                  |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`CemhydMat` :elemparam:`num{in}` :elemparam:`d{rn}` :elemparam:`k{rn}`             |
   |                    | :elemparam:`c{rn}` :elemparam:`file{s}` [:elemparam:`eachGP{in}`]                            |
   |                    | [:elemparam:`densityType{in}`] [:elemparam:`conductivityType{in}`]                           |
   |                    | [:elemparam:`capacityType{in}`] [:elemparam:`castingtime{rn}`] [:elemparam:`nowarnings{ia}`] |
   |                    | [:elemparam:`scaling{ra}`] [:elemparam:`reinforcementDegree{rn}`]                            |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`num` material model number                                                         |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`k` Conductivity                                                                    |
   |                    | - :param:`c` Specific heat capacity                                                          |
   |                    | - :param:`file` XML input file for cement microstructure and concrete composition            |
   |                    | - :param:`eachGP` 0 (default) no separate microstructures in each GP, 1 assign separate      |
   |                    |   microstructures to each GP                                                                 |
   |                    | - :param:`densityType` 0 (default) get density from OOFEM input file, 1 get it from XML      |
   |                    |   input file                                                                                 |
   |                    | - :param:`conductivityType` 0 (default) get constant conductivity from OOFEM input file, 1   |
   |                    |   compute as :math:`\lambda = \textrm{k} (1.33-0.33\alpha)` [Ruiz:01]_                       |
   |                    | - :param:`capacityType` 0 (default) get capacity, 1 according to Bentz, 2 according to XML   |
   |                    |   and CEMHYD3D routines                                                                      |
   |                    | - :param:`castingtime` optional casting time of concrete, from which hydration takes place.  |
   |                    |   Absolute time is used.                                                                     |
   |                    | - :param:`nowarnings` supresses warnings when material data are out of standard ranges. The  |
   |                    |   array of size 4 represent entries for density, conductivity, capacity, temperature.        |
   |                    |   Nonzero values mean supression.                                                            |
   |                    | - :param:`scaling` components in the array scale density, conductivity, capacity in this     |
   |                    |   order. :param:`nowarnings` are checked before scaling.                                     |
   |                    | - :param:`reinforcementDegree` specifies the area fraction of reinforcement. Typical values  |
   |                    |   is 0.015. Steel reinforcement slightly increases concrete conductivity and slightly        |
   |                    |   decreases its capacity. Thermal properties of steel are considered 20 W/m/K and 500        |
   |                    |   J/kg/K.                                                                                    |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | _2dHeat, _3dHeat                                                                             |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `tmcemhyd/cemhyd01.in                                                                        |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/tmcemhyd/cemhyd01.in>`_,         |
   |                    | `tmcemhyd/cemhyd02.in                                                                        |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/tmcemhyd/cemhyd02.in>`_          |
   +--------------------+----------------------------------------------------------------------------------------------+

The input XML file specifies the details about cement and concrete composition.
It is possible to start all simulations from the scratch, i.e. with the
reconstruction of digital microstructure. Alternatively, the digital
microstructure can be provided directly in two files; one for chemical phases,
the second for particle's IDs. The XML input file can be created with the CemPy
package, obtainable from
http://mech.fsv.cvut.cz/\~smilauer/index.php?id=software. The CemPy package
alleviates tedious preparation of particle size distribution etc.

The linear solver (specified as NonStationaryProblem) performs well when the
time integration step is small enough (order of minutes) and heat capacity,
conductivity and density remain constant. If not so, use of nonlinear solver is
strongly suggested (specified as NlTransientTransportProblem).

.. _Affinity1:

HydratingConcreteMat
~~~~~~~~~~~~~~~~~~~~

Simple hydration models based on chemical affinity are implemented. The models calculate the degree of hydration of cement, :math:`\alpha`, which can be scaled to the level of concrete when providing the corresponding amount of cement in concrete. Blended cements can be considered as well, either by separating supplementary cementitious materials from pure Portland clinker or by providing parameters for the evolution of hydration degree and potential heat. The released heat from the cement paste is obtained from

.. math::

   Q(t) = \alpha Q_{pot},

where the potential hydration heat, :math:`Q_{pot}`, is expressed in kJ/kg of cement and for pure Portland cement is around 500 kJ/kg.

Evolution of hydration degree under isothermal curing conditions is approximated by several models. Scaling from a reference temperature to arbitrary temperature is based on the Arrhenius equation, which coincides with the maturity method approach. The equivalent time, :math:`t_e`, is defined as the time under constant reference (isothermal) temperature

.. math::

   t_e(T_0) &= t(T) k_{rate},\\
   k_{rate} &= \exp\left[\frac{E_a}{R}\left(\frac{1}{T_0} - \frac{1}{T}\right)\right],

where :math:`t` is real time, :math:`T` is the arbitrary constant temperature of hydration, :math:`T_0` is a reference temperature, :math:`R` is the universal gas constant (8.314 Jmol\ :math:`^{-1}`\ K\ :math:`^{-1}`), and :math:`E_a` is the apparent activation energy. Due to the varying history of temperature, an incremental solution is adopted. Linear and nonlinear nonstationary solvers are supported for all :param:`hydrationmodeltype`'s. Hydration models are evaluated at intrinsic time in each time step. Usually, intrinsic time is in the middle of the time step.

The :param:`hydrationmodeltype` = 1 is based on an exponential approximation of hydration degree [Schindler:2005]_. An equivalent time increment is added in each time step. Thus, all the thermal history is stored in the equivalent time

.. math::

   \alpha(t_e) = \alpha_\infty \exp\left(-\left[\frac{\tau}{t_e}\right]^\beta\right)

where three parameters :math:`\tau`, :math:`\beta`, and :math:`\alpha_\infty` are needed. Some meaningful parameters are provided in [Schindler:2005]_, e.g., :math:`\tau=26\cdot3600=93600` s, :math:`\beta=0.75`, and :math:`\alpha_\infty=0.90`.

The :param:`hydrationmodeltype` = 2 is inspired by Cervera et al. [Cervera:99]_, who proposed an analytical form of the normalized affinity which was refined in [Gawin:06a]_. A slightly modified formulation is proposed here. The affinity model is formulated for a reference temperature 25°C

.. math::
   :label: eq:affinity2

   \frac{\partial \alpha}{\partial t} &=& \tilde{A}_{25}(\alpha) k_{rate}= B_1 \left( \frac{B_2}{\alpha_\infty} + \alpha \right ) \left( \alpha_\infty - \alpha \right) f_s \exp\left(-\bar{\eta}\frac{\alpha}{\alpha_\infty}\right) k_{rate}\\
   \alpha &>& DoH_1 \Rightarrow f_s = 1+P_1(\alpha - DoH_1)~\mathrm{else}~f_s = 1 

where :math:`B_1, B_2` are coefficients to be calibrated, :math:`\alpha_\infty` is the ultimate hydration degree and :math:`\bar{\eta}` represents microdiffusion of free water through formed hydrates. The function :math:`f_s` adds additional peak which may occur in slag-rich blended cements with two parameters :math:`DoH_1,P_1`. The solution proceeds incrementally, where :math:`\alpha` is the unknown. During one macroscopic time step, :eq:`eq:affinity2` needs to be integrated in finer inner steps. This is controlled with two optional variables; :param:`maxmodelintegrationtime` specifies maximum integration time in the loop while :param:`minmodeltimestepintegrations` specifies minimum number of integration steps.

.. table:: HydratingConcreteMat - summary of affinity hydration models.
   :name: Affinity1_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | HydratingConcreteMat                                                                         |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`HydratingConcreteMat` :elemparam:`num{in}` :elemparam:`d{rn}` :elemparam:`k{rn}`  |
   |                    | :elemparam:`c{rn}` :elemparam:`hydrationmodeltype{in}` :elemparam:`Qpot{rn}`                 |
   |                    | :elemparam:`masscement{rn}` [:elemparam:`activationenergy{rn}`]                              |
   |                    | :elemparam:`reinforcementdegree{rn}`] [:elemparam:`densitytype{in}`]                         |
   |                    | [:elemparam:`conductivitytype{in}`] [:elemparam:`capacityType{in}`]                          |
   |                    | [:elemparam:`minModelTimeStepIntegrations{in}`] [:elemparam:`maxmodelintegrationtime{rn}`]   |
   |                    | [:elemparam:`castingTime{rn}`]                                                               |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`num` material model number                                                         |
   |                    | - :param:`d` material density about 2300 kg/m\ :math:`^3`                                    |
   |                    | - :param:`k` Conductivity about 1.7 W/m/K                                                    |
   |                    | - :param:`c` Specific heat capacity about 870 J/kg/K                                         |
   |                    | - :param:`hydrationmodeltype` 1 is exponential model from \refeqeq:affinity1, 2 is affinity  |
   |                    |   model from \refeqeq:affinity2                                                              |
   |                    | - :param:`Qpot` Potential heat of hydration, about 500 kJ/kg of cement                       |
   |                    | - :param:`masscement` Cement mass per 1m\ :math:`^3` of concrete, about 200-450              |
   |                    | - :param:`activationenergy` Arrhenius activation energy, 38400. (default)                    |
   |                    | - :param:`DoHinf` degree of hydration at infinite time                                       |
   |                    | - :param:`B1,B2,eta,DoH1,P1` parameters from \refeqeq:affinity2-\refeqeq:affinity3           |
   |                    | - :param:`reinforcementDegree` specifies the area fraction of reinforcement. Typical values  |
   |                    |   is 0.015. Steel reinforcement slightly increases concrete conductivity and slightly        |
   |                    |   decreases its capacity. Thermal properties of steel are considered 20 W/m/K and 500        |
   |                    |   J/kg/K.                                                                                    |
   |                    | - :param:`densityType` 0 (default)                                                           |
   |                    | - :param:`conductivityType` 0 (default), 1 compute as :math:`\lambda = \textrm{k}            |
   |                    |   (1.33-0.33\alpha)` [Ruiz:01]_                                                              |
   |                    | - :param:`capacityType` 0 (default)                                                          |
   |                    | - :param:`minModelTimeStepIntegrations` Minimum integrations per time step in affinity model |
   |                    |   30 (default)                                                                               |
   |                    | - :param:`maxmodelintegrationtime` Maximum integration time step in affinity model 36000 s   |
   |                    |   (default)                                                                                  |
   |                    | - :param:`castingtime` optional casting time of concrete, from which hydration takes place,  |
   |                    |   in s                                                                                       |
   |                    | - :param:`relmatage` material age at casting time, 0 s (default)                             |
   |                    | - :param:`scaling` components in the array scale density, conductivity, capacity in this     |
   |                    |   order. :param:`nowarnings` are checked before scaling.                                     |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | _2dHeat, _3dHeat                                                                             |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `tm/TwoStepCasting_01.in                                                                     |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/TwoStepCasting_01.in>`_,      |
   |                    | `tm/hydratingConcreteMat01.in                                                                |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/hydratingConcreteMat01.in>`_, |
   |                    | `tm/hydratingConcreteMat02.in                                                                |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/hydratingConcreteMat02.in>`_, |
   |                    | `tm/hydratingConcreteMat03.in                                                                |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/hydratingConcreteMat03.in>`_, |
   |                    | `tm/hydratingConcreteMat04.in                                                                |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/hydratingConcreteMat04.in>`_, |
   |                    | `tm/hydratingConcreteMat05.in                                                                |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/hydratingConcreteMat05.in>`_, |
   |                    | `tm/hydratingConcreteMat06.in                                                                |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/hydratingConcreteMat06.in>`_, |
   |                    | `tm/hydratingConcreteMat07.in                                                                |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/hydratingConcreteMat07.in>`_  |
   +--------------------+----------------------------------------------------------------------------------------------+

.. _hydration_comparison:

.. figure:: /figures/Mokra_OOFEM_affinity_time.png
   :width: 70%
   :alt: Performance of implemented hydration models

   Performance of implemented hydration models: exponential model from ``eq:affinity1``, affinity model from :eq:`eq:affinity2`, CEMHYD3D model from Subsection :ref:`Cemhyd`.

   Parameters for exponential model according to ``eq:affinity1`` are :math:`\tau=26\cdot3600=93600` s, :math:`\beta=0.75`, :math:`\alpha_\infty=0.90`. Parameters for affinity model according to :eq:`eq:affinity2` are :math:`B_1=3.519e-4` s\ :sup:`-1`, :math:`B_2=8.0e-7`, :math:`\eta=7.4`, :math:`\alpha_\infty=0.85`.
