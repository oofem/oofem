Material models specific to concrete
====================================


.. _mazars-model:

Mazars damage model for concrete - MazarsModel
----------------------------------------------

This isotropic damage model assumes that the stiffness degradation is isotropic, i.e., stiffness moduli corresponding to different directions decrease proportionally and independently of direction of loading.

It introduces two damage parameters :math:`\omega_t` and :math:`\omega_c` that are computed from the same equivalent strain using two different damage functions :math:`g_t` and :math:`g_c`. The :math:`g_t` is identified from the uniaxial tension tests, while :math:`g_c` from compressive test. The damage parameter for general stress states :math:`\omega` is obtained as a linear combination of :math:`\omega_t` and :math:`\omega_c`:
:math:`\omega=\alpha_t g_t + \alpha_c g_c`, where the coefficients :math:`\alpha_t` and :math:`\alpha_c` take into account the character of the stress state.

The damaged stiffness tensor is expressed as
:math:`\mbf{D}=(1-\omega)\mbf{D}^e`.

Damage evolution law is postulated in an explicit form, relating damage parameter and scalar measure of largest reached strain level in material, taking into account the principle of preserving of fracture energy :math:`G_f`. The equivalent strain, i.e., a scalar measure of the strain level is defined as norm from positive principal strains.

The model description and parameters are summarized in :numref:`maz_table`.

.. table:: Mazars damage model -- summary
   :name: maz_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Mazars damage model for concrete                                                             |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`MazarsModel` :elemparam:`d{rn}` :elemparam:`E{rn}` :elemparam:`n{rn}`             |
   |                    | :elemparam:`e0{rn}` :elemparam:`ac{rn}` [:elemparam:`bc{rn}`] [:elemparam:`beta{rn}`]        |
   |                    | :elemparam:`at{rn}` :optelemparam:`bt{rn}` [:elemparam:`hreft{rn}`] [:elemparam:`hrefc{rn}`] |
   |                    | [:elemparam:`version{in}`] [:elemparam:`tAlpha{rn}`] [:elemparam:`equivstraintype{in}`]      |
   |                    | [:elemparam:`maxOmega{rn}`]                                                                  |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`num` material model number                                                         |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`E` Young modulus                                                                   |
   |                    | - :param:`n` Poisson ratio                                                                   |
   |                    | - :param:`e0` max effective strain at peak                                                   |
   |                    | - :param:`ac`,:param:`bc` material parameters related to the shape of uniaxial compression   |
   |                    |   curve (A sample set used by Saouridis is :math:`A_c = 1.34, B_c = 2537`                    |
   |                    | - :param:`beta` coefficient reducing the effect of damage under response under shear.        |
   |                    |   Default value set to 1.06                                                                  |
   |                    | - :param:`at`, :optparam:`bt` material parameters related to the shape of uniaxial tension   |
   |                    |   curve. Meaning dependent on :param:`version` parameter.                                    |
   |                    | - :param:`hreft`, :param:`hrefc` reference characteristic lengths for tension and            |
   |                    |   compression. The material parameters are specified for element with these characteristic   |
   |                    |   lengths. The current element then will have the same COD (Crack Opening Displacement) as   |
   |                    |   reference one.                                                                             |
   |                    | - :param:`version` Model variant. if 0 specified, the original form :math:`g_t=              |
   |                    |   1.0-(1.0-A_t)*\varepsilon_0/\kappa - A_t*\exp(-B_t*\ (\kappa-\varepsilon_0));` of tension  |
   |                    |   damage evolution law is used, if equal 1, the modified law used which asymptotically tends |
   |                    |   to zero :math:`g_t =                                                                       |
   |                    |   1.0-(\varepsilon_0/\kappa)*\exp((\varepsilon_0-\kappa)/\varepsilon_f)`                     |
   |                    | - :param:`tAlpha` thermal dilatation coefficient                                             |
   |                    | - :param:`equivstraintype` see Tab. :numref:`id_table`                                       |
   |                    | - :param:`maxOmega` limit maximum damage, use for convergency improvement                    |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat, PlaneStress, PlaneStrain, 1dMat                                                       |
   +--------------------+----------------------------------------------------------------------------------------------+


Nonlocal Mazars damage model for concrete - MazarsModelnl
---------------------------------------------------------


The nonlocal variant of Mazars damage model for concrete.
Model based on nonlocal averaging of equivalent strain.
The nonlocal averaging acts as a powerful localization limiter. The bell-shaped nonlocal averaging function is used.
The model description and parameters are summarized in :numref:`maznl_table`.

.. table:: Nonlocal Mazars damage model -- summary
   :name: maznl_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Nonlocal Mazars damage model for concrete                                                    |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`MazarsModelnl` :elemparam:`r{rn}` :elemparam:`E{rn}` :elemparam:`n{rn}`           |
   |                    | :elemparam:`e0{rn}` :elemparam:`ac{rn}` :elemparam:`bc{rn}` :elemparam:`beta{rn}`            |
   |                    | :elemparam:`version{in}` :elemparam:`at{rn}` :optelemparam:`bt{rn}` :elemparam:`r{rn}`       |
   |                    | :elemparam:`tAlpha{rn}`                                                                      |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`num` material model number                                                         |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`E` Young modulus                                                                   |
   |                    | - :param:`n` Poisson ratio                                                                   |
   |                    | - :param:`maxOmega` limit maximum damage, use for convergency improvement                    |
   |                    | - :param:`tAlpha` thermal dilatation coefficient                                             |
   |                    | - :param:`version` Model variant. if 0 specified, the original form :math:`g_t=              |
   |                    |   1.0-(1.0-A_t)*\varepsilon_0/\kappa - A_t*\exp(-B_t*\ (\kappa-\varepsilon_0));` of tension  |
   |                    |   damage evolution law is used, if equal 1, the modified law used which asymptotically tends |
   |                    |   to zero :math:`g_t = 1.0-(\varepsilon_0/\kappa)*\exp((\varepsilon_0-\kappa)/A_t)`          |
   |                    | - :param:`ac`,:param:`bc` material parameters related to the shape of uniaxial compression   |
   |                    |   curve (A sample set used by Saouridis is :math:`A_c = 1.34, B_c = 2537`                    |
   |                    | - :param:`at`, :optparam:`bt` material parameters related to the shape of uniaxial tension   |
   |                    |   curve. Meaning dependent on :param:`version` parameter.                                    |
   |                    | - :param:`beta` coefficient reducing the effect of damage under response under shear.        |
   |                    |   Default value set to 1.06                                                                  |
   |                    | - :param:`r` parameter specifying the width of nonlocal averaging zone                       |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat, PlaneStress, PlaneStrain, 1dMat                                                       |
   +--------------------+----------------------------------------------------------------------------------------------+


CebFip78 model for concrete creep with aging - CebFip78
-------------------------------------------------------


Implementation of aging viscoelastic model for concrete creep
according to the CEB-FIP Model Code.
The model parameters are summarized
in :numref:`cebfip_table`.


CebFip78 material model -- summary.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. table:: CebFip78 material model -- summary.
   :name: cebfip_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | CebFip78 model for concrete creep with aging                                                 |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`CebFip78` :elemparam:`n{rn}` :elemparam:`relMatAge{rn}` :elemparam:`E28{rn}`      |
   |                    | :elemparam:`fibf{rn}` :elemparam:`kap_a{rn}` :elemparam:`kap_c{rn}` :elemparam:`kap_tt{rn}`  |
   |                    | :elemparam:`u{rn}`                                                                           |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`num` material model number                                                         |
   |                    | - :param:`E28` Young modulus at age of 28 days [MPa]                                         |
   |                    | - :param:`n` Poisson ratio                                                                   |
   |                    | - :param:`fibf` basic creep coefficient                                                      |
   |                    | - :param:`kap_a` coefficient of hydrometric conditions                                       |
   |                    | - :param:`kap_c` coefficient of type of cement                                               |
   |                    | - :param:`kap_tt` coeficient of temperature effects                                          |
   |                    | - :param:`u` surface imposed to environment [:math:`mm^2`]; temporary here; should be in     |
   |                    |   crosssection level                                                                         |
   |                    | - :param:`relmatage` relative material age                                                   |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat, PlaneStress, PlaneStrain, 1dMat, 2dPlateLayer,2dBeamLayer, 3dShellLayer               |
   +--------------------+----------------------------------------------------------------------------------------------+


Double-power law model for concrete creep with aging - DoublePowerLaw
---------------------------------------------------------------------

Implementation of aging viscoelastic model for concrete creep
with compliance function given by the double-power law.
The model parameters are summarized
in :numref:`doublepowerlaw_table`.

.. table:: Double-power law model -- summary.
   :name: doublepowerlaw_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Double-power law model for concrete creep with aging                                         |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`DoublePowerLaw` :elemparam:`n{rn}` :elemparam:`relMatAge{rn}`                     |
   |                    | :elemparam:`E28{rn}` :elemparam:`fi1{rn}` :elemparam:`m{rn}` :elemparam:`n{rn}`              |
   |                    | :elemparam:`alpha{rn}`                                                                       |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`num` material model number                                                         |
   |                    | - :param:`E28` Young modulus at age of 28 days [MPa]                                         |
   |                    | - :param:`n` Poisson ratio                                                                   |
   |                    | - :param:`fibf` basic creep coefficient                                                      |
   |                    | - :param:`m` coefficient                                                                     |
   |                    | - :param:`n` coefficient                                                                     |
   |                    | - :param:`alpha` coeficient                                                                  |
   |                    | - :param:`relmatage` relative material age                                                   |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat, PlaneStress, PlaneStrain, 1dMat, 2dPlateLayer,2dBeamLayer, 3dShellLayer               |
   +--------------------+----------------------------------------------------------------------------------------------+


Eurocode 2 model for concrete creep and shrinkage - EC2CreepMat
---------------------------------------------------------------

Implementation of aging viscoelastic model for concrete creep
according to Eurocode 2 for concrete structures
The model parameters are summarized
in :numref:`ec2creep_table`.

According to EC2, the compliance function is defined using the creep coefficient :math:`\varphi` as

.. math::

   J(t,t') = \frac{1}{E(t')} + \frac{\varphi(t,t')}{ 1.05 E_{cm}}

where :math:`E_{cm}` is the mean elastic modulus at the age of 28 days and :math:`E(t')` is the elastic modulus at the age of loading, :math:`t'`.

Current implementation supports only linear creep which is valid only for stresses below 0.45 of the characteristic compressive strength at the time of loading.

The elastic modulus at age :math:`t` (in days) is defined as

.. math::

   E(t) = \left[ \exp \left( s \left(1-\sqrt{28/t} \right) \right)
   \right]^{0.3} E_{cm}

where :math:`s` is a cement-type dependent constant (0.2 for class R, 0.25 for class N and finally 0.38 for type S), and the mean secant elastic modulus at 28 days can be estimated from the mean compressive strength

.. math::

   E_{cm} = 22 \: \left(0.1 f_{cm,28} \right)^{0.3} 

where :math:`f_{cm,28}` is in MPa and the resulting modulus is in GPa.

The creep coefficient is given by the expression from Annex B of the standard.

.. math::

   \varphi(t,t') = \varphi_{RH} \times \frac{16.8}{f_{cm}} \times \frac{1}{0.1 + t'^{0.2}}  \left( \frac{t-t'}{\beta_H + t - t'_T} \right)^{0.3}

with 

.. math::

   \varphi_{RH} = \left( 1 + \alpha_1 \frac{1-h_{env}}{0.1 \: (h_0)^{1/3}} 
   \right) \times \alpha_2\\
   \beta_H = 1.5 \left(1 + \left( 1.2 h_{env} \right)^{18} \right) h_0 +
   250 \: \alpha_3 \leq 1500 \: \alpha_3

and 

.. math::

   \alpha_1 = \alpha_2 = \alpha_3  = 1 \quad \mathrm{for} \: f_{cm} \leq
   35 \: \mathrm{MPa}

else 

.. math::

   \alpha_1 = \left( \frac{35}{f_{cm}} \right)^{0.7}\\
   \alpha_2 = \left( \frac{35}{f_{cm}} \right)^{0.2}\\
   \alpha_3 = \left( \frac{35}{f_{cm}} \right)^{0.5}

and :math:`t'_T` is the temperature-adusted age according to B.9 in the code. It is computed automatically if string :param:`temperatureDependent` appears in the input record.

The shrinkage deformation is additively split into two parts: drying shrinkage :math:`\varepsilon_{sh,d}` and autogenous shrinkage :math:`\varepsilon_{sh,a}`. Drying shrinkage strain at time :math:`t` is computed from

.. math::

   \varepsilon_{sh,d} = \frac{t-t_0}{t-t_0 + 0.04 \: h^{3/2}} \: k_h \: \varepsilon_{sh,d,0}

where :math:`t_0` is the duration of curing, :math:`k_h` is a thickness-dependent parameter and

.. math::

   \varepsilon_{sh,d,0} = 1.3175 \times 10^{-6} \left[ (220 + 110 \alpha_{ds1} ) \exp \left( -0.1 \alpha_{ds2} \: f_{cm} \right) \right] \left(1- h_{env}^3\right)

with 
:math:`\alpha_{ds1} = 3` for cement class *S*, 4 for class *\ N*, and 6 for class *\ R*, and :math:`\alpha_{ds2} = 0.13` for cement class *\ S*, 0.12 for class *\ N*, and 0.11 for class *\ R*.

Autogenous shrinkage strain can be computed as

.. math::

   \varepsilon_{sh,a} = 2.5 \: (f_{cm} - 18) \left[ 1- \exp \left( -0.2 \sqrt{t} \right)  \right] \times 10^{-6}

.. table:: EC2Creep material model -- summary.
   :name: ec2creep_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | EC2CreepMat model for concrete creep and shrinkage                                           |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`EC2CreepMat` :elemparam:`n{rn}` :optelemparam:`begOfTimeOfInterest{rn}`           |
   |                    | :optelemparam:`endOfTimeOfInterest{rn}` :elemparam:`relMatAge{rn}`                           |
   |                    | :optelemparam:`timeFactor{rn}` :elemparam:`stiffnessFactor{rn}` :optelemparam:`tAlpha{rn}`   |
   |                    | :elemparam:`fcm28{rn}` :elemparam:`t0{rn}` :elemparam:`cemType{in}` :optelemparam:`henv{rn}` |
   |                    | :elemparam:`h0{rn}` :elemparam:`shType{in}` :optelemstring:`spectrum`                        |
   |                    | :optelemstring:`temperatureDependent`                                                        |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`num` material model number                                                         |
   |                    | - :param:`n` Poisson's ratio                                                                 |
   |                    | - :param:`begOfTimeOfInterest` determines the shortest time which is reasonably well         |
   |                    |   captured by the approximated compliance function (default value is 0.1); the units are the |
   |                    |   time units of the analysis                                                                 |
   |                    | - :param:`endOfTimeOfInterest` determines the longest time which is reasonably well captured |
   |                    |   by the approximated compliance function (if not provided it is read from the engineering   |
   |                    |   model); the units are the time units of the analysis                                       |
   |                    | - :param:`relMatAge` time shift used to specify the age of material on the begging of the    |
   |                    |   analysis, the meaning is the material age at time t = 0;                                   |
   |                    | - :param:`timeFactor` scaling factor transforming the actual time into appropriate units     |
   |                    |   needed by the formulae of the eurocode. For analysis in days timeFactor = 1, for analysis  |
   |                    |   in seconds timeFactor = 86,400.                                                            |
   |                    | - :param:`stiffnessFactor` scaling factor transforming predicted stiffness into appropriate  |
   |                    |   units of the analysis, for analysis in MPa stiffnessFactor = 1.e6 (default), for Pa        |
   |                    |   stiffnessFactor = 1                                                                        |
   |                    | - :param:`fcm28` mean compressive strength measured on cylinders at the age of 28 days in    |
   |                    |   MPa                                                                                        |
   |                    | - :param:`t0` duration of curing [day] (this is relevant only for drying shrinkage, not for  |
   |                    |   creep)                                                                                     |
   |                    | - :param:`cemType` type of cement, 1 = class *R*, 2 = class *\ N*, 3 = class *\ S*           |
   |                    | - :param:`henv` ambient relative humidity expressed as decimal                               |
   |                    | - :param:`h0` effective member thickness in [mm] calculated according to EC2 as 2\           |
   |                    |   :math:`\times`\ A/u where A is the cross-section are and u is the cross-section perimeter  |
   |                    |   exposed to drying                                                                          |
   |                    | - :param:`shType` shrinkage type; 0 = no shrinkage, 1 = both drying and autogenous           |
   |                    |   shrinkage, 2 = drying shrinkage only, 3 = autogenous shrinkage only                        |
   |                    | - :param:`spectrum` this string switches on evaluation of the moduli of the aging Kelvin     |
   |                    |   chain using the retardation spectrum of the compliance function, otherwise (default        |
   |                    |   option) the least-squares method is used                                                   |
   |                    | - :param:`temperatureDependent` turns on the influence of temperature on concrete maturity   |
   |                    |   (equivalent age concept) by default this option is not activated.                          |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat, PlaneStress, PlaneStrain, 1dMat, 2dPlateLayer,2dBeamLayer, 3dShellLayer               |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `sm/EC2creep.in                                                                              |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/EC2creep.in>`_,               |
   |                    | `sm/EC2creep_casting.in                                                                      |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/EC2creep_casting.in>`_,       |
   |                    | `sm/EC2shrinkage.in                                                                          |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/EC2shrinkage.in>`_            |
   +--------------------+----------------------------------------------------------------------------------------------+

.. _B3-and-MPS-models-for-concrete-creep-with-aging:

B3 and MPS models for concrete creep with aging
-----------------------------------------------

Model B3 is an aging viscoelastic model for concrete creep and shrinkage, developed by Prof. Bažant and coworkers. In OOFEM it is implemented in three different ways.

The first version, “B3mat”, is kept in OOFEM for compatibility. It is based on an aging Maxwell chain. The moduli of individual units in the chain are evaluated in each step using the least-squares method.

The second, more recent version, is referred to as “B3solidmat”. Depending on the specified input it exploits either a non-aging Kelvin chain combined with the solidification theory, or an aging Kelvin chain. It is extended to the so-called Microprestress-Solidification theory (MPS), which in this implementation  takes into account only the effects of variable humidity on creep; the effects of temperature on creep are not considered. 
The underlying rheological chain consists of four serially coupled components. The solidifying Kelvin chain represents short-term creep; it is serially coupled with a non-aging elastic spring that reflects instantaneous deformation. Long-term creep is captured by an aging dashpot with viscosity dependent on the microprestress, the evolution of which is affected by changes of humidity. The last  unit describes volumetric deformations (shrinkage and thermal strains).
Drying creep is incorporated either by the 
“averaged cross-sectional approach”, or by the “point approach”.

MPS
~~~

The latest version is denoted as “MPS” and is based on the microprestress-solidification theory [Bazant-89-I]_ [Bazant-97-I]_ [JirHav14a]_. The rheological model consists of the same four components as in “B3solidmat”, but now the implemented exponential algorithm is designed especially for the solidifying Kelvin chain, which is a special case of an aging Kelvin chain. This model takes into account both humidity and temperature effects on creep. Drying creep is incorporated exclusively by the so-called “point approach”. The model can operate in four modes, controlled by the keyword :math:`CoupledAnalysisType`. The first mode (:math:`CoupledAnalysisType = 0`) solves only the basic creep and runs as a single problem, while the remaining three modes need to be run as a staggered problem with humidity and/or temperature analysis preceding the mechanical problem; both humidity and temperature fields are read when :math:`CoupledAnalysisType = 1`, only the field of relative humidity is taken into account when :math:`CoupledAnalysisType = 2` and finally, only temperature when :math:`CoupledAnalysisType = 3`.

The basic creep is in the microprestress-solidification theory influenced by the same four parameters :math:`q_1` - :math:`q_4` as in the model B3. Values of these parameters can be estimated from the composition of concrete mixture and its compressive strength using the following empirical formulae:

.. math::

   q_1 &= 126.77 \bar{f_c}^{-0.5} \hspace{5 mm} [10^{-6}/\mbox {MPa}]\\
   q_2 &= 185.4 c^{0.5} \bar{f_c}^{-0.9} \hspace{5 mm} [10^{-6}/\mbox {MPa}]\\
   q_3 &= 0.29 \left(w/c\right)^4 q_2 \hspace{5 mm} [10^{-6}/\mbox {MPa}]\\
   q_4 &= 20.3 \left(a/c\right)^{-0.7} \hspace{5 mm} [10^{-6}/\mbox {MPa}]  

Here, :math:`\bar{f_c}` is the average compressive cylinder strength at age of 28 days [MPa], :math:`a`, :math:`w` and :math:`c` is the weight of aggregates, water and cement per unit volume of concrete [kg/m^3].

.. _spring:

Spring
~~~~~~

The non-aging spring stiffness represents the asymptotic modulus of the material; it is equal to :math:`1/q_1`.

.. _kelvin-chain:

Kelvin Chain
~~~~~~~~~~~~

The solidifying Kelvin chain is composed of :math:`M` Kelvin units with fixed retardation times :math:`\tau_{\mu}`, :math:`\mu = 1, 2, \dots, M`, which form a geometric progression with quotient 10. The lowest retardation time :math:`\tau_1` is equal to :math:`0.3\; \text{begoftimeofinterest}`, the highest retardation time :math:`\tau_M` is bigger than :math:`0.5\; \text{endoftimeofinterest}`. The chain also contains a spring with stiffness :math:`E_0^\infty` (a special case of Kelvin unit with zero retardation time).

Moduli :math:`E_{\mu}^\infty` of individual Kelvin units are determined such that the chain provides a good approximation of the non-aging micro-compliance function of the solidifying constituent, :math:`\Phi(t-t') = q_2 \ln \left( 1 + \left( \left(t-t' \right)/\lambda_0 \right)^n \right)`, where :math:`\lambda_0 = 1` day and :math:`n = 0.1`.

The technique based on the continuous retardation spectrum leads to the following formulae:

.. math::

   \frac{1}{E_0^\infty} &= q_2 \ln\left(1+\tilde{\tau_0}\right) - \frac{q_2 \tilde{\tau_0} }{10 \left( 1 + \tilde{\tau_0} \right)} \quad \text{where} \quad \tilde{\tau_0} = \left( \frac{2 \tau_1}{ \sqrt{10}}  \right)^{0.1} \\
   \frac{1}{E_\mu^\infty} &=  (\ln 10) \frac{q_2 \tilde{\tau}_\mu \left(0.9 + \tilde{\tau}_\mu \right)}{10 \left( 1 + \tilde{\tau}_\mu \right)^2} \quad \text{where} \quad \tilde{\tau}_\mu = \left(2 \tau_\mu \right)^{0.1}, \quad \mu = 1,2,\dots M

Viscosities :math:`\eta_{\mu}^\infty` of individual Kelvin units are obtained from the simple relation :math:`\eta_{\mu}^\infty=\tau_{\mu}/E_{\mu}^\infty`. A higher accuracy is reached if all retardation times are multiplied by the factor 1.35 and the last modulus :math:`E_M` is divided by 1.2.

.. _solidification:

Solidification
~~~~~~~~~~~~~~

The actual viscosities :math:`\eta_{\mu}` and stiffnesses :math:`E_{\mu}` of the solidifying chain change in time according to :math:`\eta_{\mu}(t) = v(t) \eta_{\mu}^\infty` and :math:`E_{\mu}(t) = v(t) E_{\mu}^\infty`, where

.. math::

   v(t)= \frac{1}{\frac{q_3}{q_2} + \left( \frac{\lambda_0}{ t} \right)^m}

is the volume growth function, and exponent :math:`m = 0.5`. In the case of variable temperature or humidity, the actual age of concrete :math:`t` is replaced by the equivalent time :math:`t_e`, which is obtained by integrating :math:`(``\ dtedt``)`.

.. _microprestress:

Microprestress
~~~~~~~~~~~~~~

Evolution of viscosity of the aging dashpot is governed by the differential
equation

.. math::
   :label: eq:mps_viscosity

   \dot{\eta} + \frac{1}{\mu_S T_0} \left | T \frac{\dot{h}}{h} - \kappa_T k_T \dot{T} \right | \left( \mu_S \eta \right)^{p/(p-1)} = \frac{\psi_S}{q_4}

where :math:`h` is the relative pore humidity, :math:`T` is the absolute temperature [K], :math:`T_0 = 298` K is the room temperature, and parameter :math:`p = 2`. Parameter :math:`k_T` is different for monotonically increasing and for cyclic temperature, and is defined as

.. math::
   :label: eq:mps_temperature

   \kappa_T = \begin{cases}
   k_{Tm} & \text{if } T = T_{\max} \text{ and } \dot{T} > 0 \\
   k_{Tc} & \text{if } T < T_{\max} \text{ or } \dot{T} \leq 0 
   \end{cases}

in which :math:`k_{Tm}` [-] and :math:`k_{Tc}` [-] are new parameters and :math:`T_{\max}` is the maximum temperature attained in the previous history of the material point.

Equation (:eq:`eq:mps_viscosity`) differs from the one presented in the original work; it replaces the differential equation for microprestress, which is not used here. The evolution of viscosity can be captured directly, without the need for microprestress. What matters is only the relative humidity and temperature and their rates. Parameters :math:`c_0` and :math:`k_1` of the original MPS theory are replaced by :math:`\mu_S = c_0 T_0^{p-1} k_1^{p-1} q_4 (p-1)^p`. The initial value of viscosity is defined as :math:`\eta(t_0) = t_0/q_4`, where :math:`t_0` is age of concrete at the onset drying or when the temperature starts changing, in the present implementation it is set :math:`relMatAge` which corresponds to the material age when the material is cast.

As mentioned above, under variable humidity and temperature conditions the physical time :math:`t` in function :math:`v(t)` describing evolution of the solidified volume is replaced by the equivalent time :math:`t_e`. In a similar spirit, :math:`t` is replaced by the solidification time :math:`t_s` in the equation describing creep of the solidifying constituent, and by the reduced time :math:`t_r` in equation :math:`\frac{d \varepsilon_f}{d t_r} = \sigma / \eta(t)` relating the flow strain rate to the stress.

Factors transforming the physical time :math:`t` into :math:`t_e`, :math:`t_r` and :math:`t_s` are defined as follows:

.. math::
   :label: dtedt

   \frac{dt_e}{dt} &=& \psi_e(t) = \beta_{eT}(T(t))\, \beta_{eh}(h(t))\\
   \frac{dt_r}{dt} &=& \psi_r(t) = \beta_{rT}(T(t))\, \beta_{rh}(h(t))\\
   \frac{dt_s}{dt} &=& \psi_s(t) = \beta_{sT}(T(t))\, \beta_{sh}(h(t))

Functions describing the influence of temperature have the form 

.. math::

   \beta_{eT}(T) &=& \exp \left[ \frac{Q_e}{R}\left( \frac{1}{T_0} - \frac{1}{T} \right) \right]\\
   \beta_{rT}(T) &=& \exp \left[ \frac{Q_r}{R}\left( \frac{1}{T_0} - \frac{1}{T} \right) \right]\\
   \beta_{sT}(T) &=& \exp \left[ \frac{Q_s}{R}\left( \frac{1}{T_0} - \frac{1}{T} \right) \right]

motivated by the rate process theory. :math:`R` is the universal gas constant and :math:`Q_e`, :math:`Q_r`, :math:`Q_s` are activation energies for hydration, viscous processes and microprestress relaxation, respectively. Only the ratios :math:`Q_e/R`, :math:`Q_r/R` and :math:`Q_s/R` have to be specified.

Functions describing the influence of humidity have the form 

.. math::

   \beta_{eh}(h) &=& \frac{1}{1+\left[\alpha_e \left( 1-h\right) \right]^4}\\
   \beta_{rh}(h) &=& \alpha_r + \left( 1 - \alpha_r \right) h^2\\
   \beta_{sh}(h) &=& \alpha_s + \left( 1 - \alpha_s \right) h^2

where :math:`\alpha_e`, :math:`\alpha_r` and :math:`\alpha_s` are parameters. 

At sealed conditions (or :math:`CoupledAnalysisType = 2`) the auxiliary coefficients :math:`\beta_{eT} = \beta_{rT} = \beta_{sT} = 1` while at room temperature (or with :math:`CoupledAnalysisType = 3`) factors :math:`\beta_{e,h} = \beta_{r,h}  = \beta_{s,h} = 1`.

Both the size effect on drying creep as well as its delay behind drying shrinkage can be addressed through parameter :math:`p` in the governing equation ``mps_viscosity``. In the experiments, the average (cross-sectional) drying creep is decreasing with specimen size. Unfortunately, for the standard value :math:`p=2`, the MPS model exhibits the opposite trend. For :math:`p = \infty` the size effect disappears, and for :math:`p < 0` it corresponds to the experiments. It should be noted that for negative or infinite values of :math:`p` the underlying theory loses its original physical background. 

If the experimental data of drying creep measured on different sizes are missing, the exponent :math:`p = \infty` can be taken as a realistic estimate.

When parameter :math:`p` is changed from its recommended value :math:`p = 2` it is advantageous to rewrite the governing differential equation to the following form

.. math::
   :label: eq:mps_reformulate

   \dot \eta_f + \frac{k_3}{T_0} \left| T \frac{\dot h}{h}  - \kappa_T  \dot T \right|  \eta_f^{\tilde{p}} = \frac{\psi_S}{q_4}

with newly introduced parameters

.. math::
   :label: eq:mps_p

   \tilde{p} = p / (p-1)

.. math::
   :label: eq:mps_k3

   k_3 = \mu_S^{\frac{1}{p-1}}

With the “standard value” :math:`p = 2` (reverted size effect) the new parameter :math:`\tilde{p}` is also 2, for :math:`p < 0` (correct size effect) :math:`\tilde{p} < 1` and finally with :math:`p = \infty` the first parameter :math:`\tilde{p} = 1` and the second parameter :math:`k_3` becomes dimensionless. The other advantage of :math:`\tilde{p} = 1` is that the governing differential equation becomes linear and can be solved directly.


The rate of thermal strain is expressed as :math:`\dot{\varepsilon}_T = \alpha_T \dot{T}` and the rate of drying shrinkage strain as :math:`\dot{\varepsilon}_{sh} = k_{sh} \dot{h}`, where both :math:`\alpha_T` and :math:`k_{sh}` are assumed to be constant in time and independent of temperature and humidity.

There are two options to simulate the autogenous shrinkage in the MPS material model. The first one is according to the B4 model

.. math::
   :label: eq:auto_shr_b4

   \varepsilon_{sh,au,B4}(t_e) = \varepsilon_{sh,au,B4}^\infty \left[ 1 + \left( \frac{ \tau_{au} }{t_e} \right) ^ {w/c/0.38} \right]^{-4.5}

and the second one is proposed in the Model Code 2010

.. math::
   :label: eq:auto_shr_fib

   \varepsilon_{sh,au,fib}(t_e) = \varepsilon_{sh,au,fib}^\infty \left( 1- \exp \left( -0.2 \: \sqrt{t_e} \right) \right)

For the normally hardening cement, the ultimate value of the autogenous shrinkage can be estimated from the composition using the empirical formula of the B4 model

.. math::
   :label: eq:auto_ultimate_B4

   \varepsilon_{sh,au,B4}^\infty = -210 \times 10^{-6}  \left ( \frac{a/c}{6} \right ) ^{-0.75} \left ( \frac{w/c}{0.38} \right ) ^ {-3.5}

Similarly, in the Model Code, the ultimate shrinkage strain can be estimated from the mean concrete strength at the age of 28 days and from the cement grade 

.. math::
   :label: eq:auto_ultimate_fib

   \varepsilon_{sh,au,fib}^\infty = -\alpha_{as} \left( \frac{ 0.1 f_{cm} } { 6 + 0.1 f_{cm} }\right) ^{2.5} \times 10^{-6}

where :math:`\alpha_{as}` is 600 for cement grades 42.5 R and 52.5 R and N, 700 for 32.5 R and 42.5 N and 800 for 32.5 N.

The model description and parameters are summarized
in :numref:`b3_table` for “B3mat”, in :numref:`b3solid_table` for “B3solidmat”, and in :numref:`mps_table` for “MPS”. 
Since some model parameters are determined from the composition
and strength using empirical formulae, it is necessary to use the
specified units (e.g.\ compressive strength always in MPa, irrespectively
of the units used in the simulation for stress). 
For “B3mat” and “B3solidmat” it is strictly 
required to use the specified units in the material input record (stress always in MPa, time in days etc.). 
The “MPS” model is almost unit-independent, except for 
:math:`\bar{f}_c` in MPa and :math:`c` in kg/m\ :sup:`3`, which are used in empirical formulae.

.. table:: B3 creep and shrinkage model -- summary
   :name: b3_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | B3 material model for concrete aging                                                         |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`B3mat` :elemparam:`d{rn}` :elemparam:`n{rn}` :elemparam:`talpha{rn}`              |
   |                    | :optelemparam:`begoftimeofinterest{rn}` :optelemparam:`endoftimeofinterest{rn}`              |
   |                    | :elemparam:`timefactor{rn}` :elemparam:`relMatAge{rn}` :optelemparam:`mode{in}`              |
   |                    | :elemparam:`fc{rn}` :elemparam:`cc{rn}` :elemparam:`w/c{rn}` :elemparam:`a/c{rn}`            |
   |                    | :elemparam:`t0{rn}` :elemparam:`q1{rn}` :elemparam:`q2{rn}` :elemparam:`q3{rn}`              |
   |                    | :elemparam:`q4{rn}` :elemparam:`shmode{in}` :elemparam:`ks{rn}` :elemparam:`vs{rn}`          |
   |                    | :elemparam:`hum{rn}` :optelemparam:`alpha1{rn}` :optelemparam:`alpha2{rn}`                   |
   |                    | :elemparam:`kt{rn}` :elemparam:`EpsSinf{rn}` :elemparam:`q5{rn}` :elemparam:`es0{rn}`        |
   |                    | :elemparam:`r{rn}` :elemparam:`rprime{rn}` :elemparam:`at{rn}` :elemparam:`w_h{rn}`          |
   |                    | :elemparam:`ncoeff{rn}` :elemparam:`a{rn}`                                                   |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`num` material model number                                                         |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`n` Poisson ratio                                                                   |
   |                    | - :param:`talpha` coefficient of thermal expansion                                           |
   |                    | - :param:`begoftimeofinterest` optional parameter; lower boundary of time interval with good |
   |                    |   approximation of the compliance function [day]; default 0.1 day                            |
   |                    | - :param:`endoftimeofinterest` optional parameter; upper boundary of time interval with good |
   |                    |   approximation of the compliance function [day]                                             |
   |                    | - :param:`timefactor` scaling factor transforming the simulation time units into days        |
   |                    | - :param:`relMatAge` relative material age [day]                                             |
   |                    | - :param:`mode` if :math:`mode = 0` (default value) creep and shrinkage parameters are       |
   |                    |   predicted from composition; for :math:`mode = 1` parameters must be user-specified.        |
   |                    | - :param:`fc` 28-day mean cylinder compression strength [MPa]                                |
   |                    | - :param:`cc` cement content of concrete [kg/m\ :math:`^{3}`]                                |
   |                    | - :param:`w/c` ratio (by weight) of water to cementitious material                           |
   |                    | - :param:`a/c` ratio (by weight) of aggregate to cement                                      |
   |                    | - :param:`t0` age when drying begins [day]                                                   |
   |                    | - :param:`q1-q4` parameters of B3 model for basic creep [1/TPa]                              |
   |                    | - :param:`shmode` shrinkage mode; :math:`0=` no shrinkage; :math:`1=` average shrinkage (the |
   |                    |   following parameters must be specified: :param:`ks`, :param:`vs`, :param:`hum` and         |
   |                    |   additionally :param:`alpha1` :param:`alpha2` for :math:`mode = 0` and :param:`kt`          |
   |                    |   :param:`EpsSinf` :param:`q5` :param:`t0` for :math:`mode = 1`; :math:`2=` point shrinkage  |
   |                    |   (needed: :param:`es0`, :param:`r`, :param:`rprime`, :param:`at`, :param:`w_h`,             |
   |                    |   :param:`ncoeff`, :param:`a`)                                                               |
   |                    | - :param:`ks` cross-section shape factor [-]                                                 |
   |                    | - :param:`vs` volume to surface ratio [m]                                                    |
   |                    | - :param:`hum` relative humidity of the environment [-]                                      |
   |                    | - :param:`alpha1` shrinkage parameter -- influence of cement type [-]                        |
   |                    | - :param:`alpha2` shrinkage parameter -- influence of curing type [-]                        |
   |                    | - :param:`kt` shrinkage parameter [day/m\ :math:`^2`]                                        |
   |                    | - :param:`EpsSinf` shrinkage parameter [10\ :math:`^{-6}`]                                   |
   |                    | - :param:`q5` drying creep parameter [1/TPa]                                                 |
   |                    | - :param:`es0` final shrinkage at material point                                             |
   |                    | - :param:`r`, :param:`rprime` coefficients                                                   |
   |                    | - :param:`at` oefficient relating stress-induced thermal strain and shrinkage                |
   |                    | - :param:`w_h`, :param:`ncoeff`, :param:`a` sorption isotherm parameters obtained from       |
   |                    |   experiments [Pedersen, 1990]                                                               |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat, PlaneStress, PlaneStrain, 1dMat, 2dPlateLayer,2dBeamLayer, 3dShellLayer               |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `sm/trussb3_creep.in                                                                         |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/trussb3_creep.in>`_,          |
   |                    | `sm/trussb3_relax.in                                                                         |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/trussb3_relax.in>`_           |
   +--------------------+----------------------------------------------------------------------------------------------+

For illustration, sample input records for the material considered in Example 3.1 of the creep book by Bažant and Jirásek is presented. The concrete mix is composed of 170 kg/m\ :sup:`3`\ of water, 450 kg/m\ :sup:`3`\ of type-I cement and 1800 kg/m\ :sup:`3`\ of aggregates, which corresponds to ratios :math:`w/c=0.3778` and :math:`a/c=4`. The compressive strength is :math:`\bar{f}_c=45.4` MPa. The concrete slab of thickness 200 mm is cured in air with initial protection against drying until the age of 7 days. Subsequently, the slab is exposed to an environment with relative humidity of 70\%. The following input record can be used for the first version of the model (B3mat):

.. code-block:: text

   B3mat 1 n 0.2 d 0. talpha 1.2e-5 relMatAge 28. fc 45.4 cc 450. w/c 0.3778 a/c 4. t0 7. timefactor 1. alpha1 1. alpha2 1.2 ks 1. hum 0.7 vs 0.1 shmode 1

Parameter :math:`\alpha_1=1` corresponds to type-I cement, parameter :math:`\alpha_2=1.2` to curing in air, parameter :math:`k_s=1` to an infinite slab. The volume-to-surface ratio is in this case equal to one half of the slab thickness and must be specified in meters, independently of the length units that are used in the finite element analysis (e.g., for nodal coordinates). The value of relMatAge must be specified in days. Parameter :tt:`relMatAge 28.` means that time 0 of the analysis corresponds to concrete age 28 days. If material B3mat is used, the finite element analysis must use days as the units of time (not only for relMatAge, but in general, e.g.\ for the time increments).

If only the basic creep (without shrinkage) should be computed, then the material input record reduces to following:

.. code-block:: text

   B3mat 1 n 0.2 d 0. talpha 1.2e-5 relMatAge 28. fc 45.4 cc 450. w/c 0.3778 a/c 4. t0 7. timefactor 1. shmode 0


B3 SOLID MAT


Parameters

- :param:`num` material model number
- :param:`d` material density
- :param:`n` Poisson ratio
- :param:`talpha` coefficient of thermal expansion

- :param:`mode` optional parameter; if :math:`mode = 0` (default), parameters :math:`q1-q4` are predicted from composition of the concrete mixture (parameters :math:`fc`, :math:`cc`, :math:`w/c`, :math:`a/c` and :math:`t0` need to be specified). Otherwise values of parameters :math:`q1-q4` are expected.
- :param:`EmoduliMode` optional parameter; analysis of retardation spectrum (:math:`=0`, default value) or least-squares method (:math:`=1`) is used for evaluation of Kelvin units moduli
- :param:`Microprestress` :math:`0=` basic creep; :math:`1=` drying creep (must be run as a staggered problem with preceding analysis of humidity diffusion. Parameter :param:`shm` must be equal to 3. The following parameters must be specified: :param:`c0`, :param:`c1`, :param:`tS0`, :param:`w_h`, :param:`ncoeff`, :param:`a`)

:param:`shmode` shrinkage mode;
:math:`0`: no shrinkage;
:math:`1`: average shrinkage (the following parameters must be specified:
:param:`ks`, :param:`vs`, :param:`hum` and additionally :param:`alpha1`, :param:`alpha2` for :math:`mode = 0`
and :param:`kt`, :param:`EpsSinf`, :param:`q5` for :math:`mode = 1`;
:math:`2`: point shrinkage (needed: :param:`es0`, :param:`r`, :param:`rprime`, :param:`at`), :param:`w_h`, :param:`ncoeff`, :param:`a`;
:math:`3`: point shrinkage based on MPS theory (needed: parameter :param:`kSh` or value of :param:`kSh` can be approximately determined if following parameters are given: :param:`inithum`, :param:`finalhum`, :param:`alpha1`, and :param:`alpha2`)

:param:`begoftimeofinterest` optional parameter; lower boundary of time interval with good approximation of the compliance function [day]; default 0.1 day
:param:`endoftimeofinterest` optional parameter; upper boundary of time interval with good approximation of the compliance function [day]

:param:`timefactor` scaling factor transforming the simulation time units into days
:param:`relMatAge` relative material age [day]

:param:`fc` 28-day mean cylinder compression strength [MPa]
:param:`cc` cement content of concrete mixture [kg/m\ :sup:`3`\ ]
:param:`w/c` water to cement ratio (by weight)
:param:`a/c` aggregate to cement ratio (by weight)
:param:`t0` age of concrete when drying begins [day]
:param:`q1`, :param:`q2`, :param:`q3`, :param:`q4` parameters (compliances) of B3 model for basic creep [1/TPa]

:param:`c0` MPS theory parameter [MPa\ :sup:`-1`\ day\ :sup:`-1`\ ]
:param:`c1` MPS theory parameter [MPa]
:param:`tS0` MPS theory parameter - time when drying begins [day]
:param:`w_h`, :param:`ncoeff`, :param:`a` sorption isotherm parameters obtained from experiments [Pedersen, 1990]

:param:`ks` cross section shape factor [-]
:param:`alpha1` optional shrinkage parameter - influence of cement type (optional parameter, default value is 1.0)
:param:`alpha2` optional shrinkage parameter - influence of curing type (optional parameter, default value is 1.0)
:param:`hum` relative humidity of the environment [-]
:param:`vs` volume to surface ratio [m]
:param:`q5` drying creep parameter [1/TPa]
:param:`kt` shrinkage parameter [day/m\ :sup:`2`\ ]
:param:`EpsSinf` shrinkage parameter [10\ :sup:`-6`\ ]

:param es0: final shrinkage at material point
:param at: coefficient relating stress-induced thermal strain and shrinkage
:param rprime, r: coefficients

:param kSh: influences magnitude of shrinkage in MPS theory [-]
:param inithum, finalhum: [-], if provided, approximate value of :math:`kSh` can be computed

Supported modes: 3dMat, PlaneStress, PlaneStrain, 1dMat, 2dPlateLayer, 2dBeamLayer, 3dShellLayer

B3solid creep and shrinkage model -- summary.

.. table:: B3solid creep and shrinkage model -- summary.
   :name: b3solid_table

   +--------------------+----------------------------------------------------------------------------------------------+
   |                    | - :param:`fc` 28-day mean cylinder compression strength [MPa]                                |
   |                    | - :param:`cc` cement content of concrete mixture [kg/m\ :math:`^{3}`]                        |
   |                    | - :param:`w/c` water to cement ratio (by weight)                                             |
   |                    | - :param:`a/c` aggregate to cement ratio (by weight)                                         |
   |                    | - :param:`t0` age of concrete when drying begins [day]                                       |
   |                    | - :param:`q1`, :param:`q2`, :param:`q3`, :param:`q4` parameters (compliances) of B3 model    |
   |                    |   for basic creep [1/TPa]                                                                    |
   |                    | - :param:`c0` MPS theory parameter [MPa\ :math:`^{-1}` day\ :math:`^{-1}`]                   |
   |                    | - :param:`c1` MPS theory parameter [MPa]                                                     |
   |                    | - :param:`tS0` MPS theory parameter - time when drying begins [day]                          |
   |                    | - :param:`w_h`, :param:`ncoeff`, :param:`a` sorption isotherm parameters obtained from       |
   |                    |   experiments [Pedersen, 1990]                                                               |
   |                    | - :param:`ks` cross section shape factor [-]                                                 |
   |                    | - :param:`alpha1` optional shrinkage parameter - influence of cement type (optional          |
   |                    |   parameter, default value is 1.0)                                                           |
   |                    | - :param:`alpha2` optional shrinkage parameter - influence of curing type (optional          |
   |                    |   parameter, default value is 1.0)                                                           |
   |                    | - :param:`hum` relative humidity of the environment [-]                                      |
   |                    | - :param:`vs` volume to surface ratio [m]                                                    |
   |                    | - :param:`q5` drying creep parameter [1/TPa]                                                 |
   |                    | - :param:`kt` shrinkage parameter [day/m\ :math:`^2`]                                        |
   |                    | - :param:`EpsSinf` shrinkage parameter [10\ :math:`^{-6}`]                                   |
   |                    | - :param:`es0` final shrinkage at material point                                             |
   |                    | - :param:`at` coefficient relating stress-induced thermal strain and shrinkage               |
   |                    | - :param:`rprime`, :param:`r` coefficients                                                   |
   |                    | - :param:`kSh` influences magnitude of shrinkage in MPS theory [-]                           |
   |                    | - :param:`inithum` [-], :param:`finalhum` [-] if provided, approximate value of :param:`kSh` |
   |                    |   can be computed                                                                            |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat, PlaneStress, PlaneStrain, 1dMat, 2dPlateLayer, 2dBeamLayer, 3dShellLayer              |
   +--------------------+----------------------------------------------------------------------------------------------+

Additional parameters
~~~~~~~~~~~~~~~~~~~~~

depend on the specific type of analysis:

1. Computing basic creep only, shrinkage not considered,
parameters :math:`q_i` estimated from composition.

.. code-block:: text

   mode 0 fc 45.4 cc 450.  w/c 0.3778 a/c 4.\\ t0 7. microprestress 0 shmode 0

2. Computing basic creep only, shrinkage not considered,
parameters :math:`q_i` specified by the user.

.. code-block:: text

   mode 1 q1 18.81 q2 126.9 q3 0.7494 q4 7.692\\ microprestress 0 shmode 0

3. Computing basic creep only, shrinkage handled using the sectional approach,
parameters estimated from composition.

.. code-block:: text

   mode 0 fc 45.4 cc 450. w/c 0.3778 a/c 4. t0 7.\\ microprestress 0 shmode 1 ks 1. alpha1 1. alpha2 1.2 hum 0.7 vs 0.1

4. Computing basic creep only, shrinkage handled using the sectional approach,
parameters specified by the user.

.. code-block:: text

   mode 1 q1 18.81 q2 126.9 q3 0.7494 q4 7.692\\
   microprestress 0 shmode 1 ks 1. q5 326.7 kt  28025 EpsSinf 702.4 t0 7. hum 0.7 vs 0.1

5. Computing basic creep only, shrinkage handled using the point approach (B3),
parameters specified by the user.

.. code-block:: text

   mode 1 q1 18.81 q2 126.9 q3 0.7494 q4 7.692\\ 
   microprestress 0 shmode 2 es0 ... r ... rprime ... at ...

6. Computing drying creep, shrinkage handled using the point approach (MPS),
parameters :math:`q_i` estimated from composition.

.. code-block:: text

   mode 0 fc 45.4 w/c 0.3778 a/c 4. t0 7. microprestress 1 \\
   shmode 3 c0 1. c1 0.2 tS0 7. w\_h 0.0476 ncoeff 0.182 a 4.867  kSh 1.27258e-003


:param mode: optional parameter; if :math:`mode = 0` (default), parameters :math:`q1-q4` are predicted from the composition of the concrete mixture (parameters fc, cc, w/c, a/c, and stiffnessfactor need to be specified). Otherwise, values of parameters :math:`q1-q4` are expected.
:param CoupledAnalysisType: :math:`0 =` basic creep; :math:`1 =` (default) drying creep, shrinkage, temperature transient creep, and creep at elevated temperature; :math:`2 =` drying creep, shrinkage; :math:`3 =` temperature transient creep and creep at elevated temperature; for choice \#1, 2, 3, the problem must be run as a staggered problem with preceding analysis of humidity and/or temperature distribution. Following parameters must be specified: :param:`mus` or :param:`k3` (according to exponent :param:`p`), :param:`kTm` (compulsory for choice \#3 otherwise optional).
:param lambda0: scaling factor equal to 1.0 day in time units of analysis (e.g. 86400 if the analysis runs in seconds).
:param begoftimeofinterest: lower boundary of time interval with good approximation of the compliance function; default value = 0.01 :math:`\lambda0`.
:param endoftimeofinterest: upper boundary of time interval with good approximation of the compliance function; default value = 10000. :math:`\lambda0`.
:param timefactor: scaling factor, for mps material must be equal to 1.0.
:param relMatAge: relative material age = age at time when the material is cast in the structure.



.. table:: MPS theory---summary
   :name: mps_table

   +--------------------+----------------------------------------------------------------------------------------------+
   |                    | - :param:`fc` 28-day standard cylinder compression strength [MPa]                            |
   |                    | - :param:`cc` cement content of concrete mixture [kg m\ :math:`^{-3}`]                       |
   |                    | - :param:`w/c` water to cement weight ratio                                                  |
   |                    | - :param:`a/c` aggregate to cement weight ratio                                              |
   |                    | - :param:`stiffnessfactor` scaling factor converting “predicted” parameters :math:`q_1` -    |
   |                    |   :math:`q_4` into proper units (e.g. 1.0 if stiffness is measured in Pa, 1.e6 for MPa)      |
   |                    | - :param:`q1`, :param:`q2`, :param:`q3`, :param:`q4` parameters of B3 model for basic creep  |
   |                    | - :param:`p` and :param:`p_tilde` replaceable parameters in the governing equation for       |
   |                    |   viscosity, default value is :math:`p = 2`                                                  |
   |                    | - :param:`mus` parameter governing to the evolution of viscosity; for exponent :math:`p=2`,  |
   |                    |   :math:`\mu_S = c_0 c_1 q_4` [Pa\ :math:`^{-1}` s\ :math:`^{-1}`]                           |
   |                    | - :param:`k3` dimensionless parameter governing to the evolution of viscosity replacing      |
   |                    |   :param:`mus` in the special case when :math:`p > 100` (then :math:`p` is automatically set |
   |                    |   to :math:`\infty` which is equivalent to :math:`p\_tilde = 1`)                             |
   |                    | - :param:`ksh` parameter relating rate of shrinkage to rate of humidity [-], default value   |
   |                    |   is 0.0, i.e. no shrinkage                                                                  |
   |                    | - :param:`t0` time of the first temperature or humidity change                               |
   |                    | - :param:`alphaE` constant, default value 10.                                                |
   |                    | - :param:`alphaR` constant, default value 0.1                                                |
   |                    | - :param:`alphaS` constant, default value 0.1                                                |
   |                    | - :param:`QEtoR` activation energy ratio, default value 2700. K                              |
   |                    | - :param:`QRtoR` activation energy ratio, default value 5000. K                              |
   |                    | - :param:`QStoR` activation energy ratio, default value 3000. K                              |
   |                    | - :param:`kTm` replaces :math:`\ln{h}` in the governing equation for viscosity               |
   |                    | - :param:`kTc` controls creep at cyclic temperature                                          |
   |                    | - :param:`alpha_as` and :param:`eps_cas0` control the ultimate value of autogenous shrinkage |
   |                    |   according to *fib* Model Code 2010; this ultimate value can be provided either directly    |
   |                    |   via :param:`eps_cas0` (negative value for contraction) or using :param:`alpha_as`, see     |
   |                    |   equation (:eq:`eq:auto_ultimate_fib`)                                                      |
   |                    | - :param:`b4_eps_au_infty`, :param:`b4_tau_au`, :param:`b4_alpha`, :param:`b4_r_t`,          |
   |                    |   :param:`b4_cem_type` control the evolution and the ultimate value of autogenous shrinkage  |
   |                    |   according to model B4; all parameters are predicted from composition if the                |
   |                    |   :param:`b4_cem_type` is provided; however, this prediction can be manually overridden      |
   |                    | - :param:`temperInCelsius` this string enables to run the supplementary transport problem    |
   |                    |   with temperature in Celsius instead of Kelvin                                              |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat, PlaneStress, PlaneStrain, 1dMat, 2dPlateLayer, 2dBeamLayer, 3dShellLayer              |
   +--------------------+----------------------------------------------------------------------------------------------+

Finally consider the same conditions for “MPS material”.
In all the examples below,
the input record with the material description can start by
``mps 1 d 2420. n 0.2 talpha 12.e-6 referencetemperature 296.``

Additional parameters
depend on the specific type of analysis:
1. Computing basic creep only, shrinkage not considered,
parameters :math:`q_i` estimated from composition and simulation time in
days and stiffnesses in MPa.

.. code-block:: text

   mode 0 fc 45.4 cc 450. w/c 0.3778 a/c 4. stiffnessFactor 1.e6
   timefactor 1. lambda0 1. begoftimeofinterest 1.e-2
   endoftimeofinterest 3.e4  relMatAge 28.  CoupledAnalysisType 0.

2. Computing basic creep only, shrinkage not considered,
parameters :math:`q_i` specified by user, simulation time in seconds and stiffnesses in Pa.

.. code-block:: text

   mode 1 q1 18.81e-12 q2 126.9e-12 q3 0.7494e-12 q4 7.6926e-12
   timefactor 1. lambda0 86400. begoftimeofinterest 864.
   endoftimeofinterest 2.592e9 relMatAge 2419200. CoupledAnalysisType 0.

3. Computing both basic and drying creep,
parameters :math:`q_i` specified by user, simulation time in seconds and stiffnesses in MPa.

.. code-block:: text

   mode 1 q1 18.81e-6 q2 126.9e-6 q3 0.7494e-6 q4 7.6926e-6
   timefactor 1. lambda0 86400. begoftimeofinterest 864.
   endoftimeofinterest 2.592e9  relMatAge 2419200. CoupledAnalysisType
   1.
   ksh 0.0004921875. t0 2419200. kappaT 0.005051 mus 4.0509259e-8

Final recommendations:

- to simulate basic creep without shrinkage it is possible to use
  all three models

  - B3mat with :math:`shmode = 0`
  - B3Solidmat with :math:`shmode = 0` and :math:`microprestress = 0`
  - MPS with :math:`CoupledAnalysisType = 0`
- to simulate drying creep with shrinkage using “sectional approach”, only the
  first material model (B3mat) is suitable can be used (with :math:`shmode = 1`)
- to simulate drying creep without shrinkage using “sectional approach”, only the
  first material model (B3mat) is suitable (with :math:`shmode = 1`, :math:`mode = 1`
  and :math:`EpsSinf = 0.0`)
- to simulate drying creep with shrinkage using “point approach”
  according to B3 model there are two options:

  - B3mat (with :math:`shmode = 2`)
  - B3Solidmat (with :math:`shmode = 2`)
    In order to suppress shrinkage set :math:`es0 = 0.0`

to simulate drying creep with shrinkage using “point approach”
according to MPS model there are two options:

   - B3Solidmat (with :math:`shmode = 3`)
   - MPS with :math:`CoupledAnalysisType = 1`

In order to suppress shrinkage set :math:`kSh = 0.0`


MPS damage model
----------------


This model extends the model based on the Microprestress-Solidification theory, described in previous section and summarized in :numref:`mps_table` for “MPS”, by tensile cracking.

This extension uses the isotropic damage model with Rankine definition of the equivalent deformation
defined as the biggest principal effective stress divided by the
elastic modulus. The softening law cannot deal directly with the
strain because due to creep strain increases and this would lead to
further softening without additional loading.

Two different approaches are implemented. The first one, default, (\#1) reduces
the stiffness only in the directions of tension (in case the tensile
strength is exceeded). A full stiffness is restored in compression and
after unloading from tension.  

The other approach (\#2) is the standard isotropic damage model which
reduces the stiffness equally in all directions independently of
loading. This approach leads to faster convergence because the secant
stiffness can be used instead of the incremental viscoelastic
stiffness which must be used in the first approach. The second
approach becomes useful when the loading is monotonic or when the
benefit of the accelerated computation prevails over the consequences
of the reduced/underestimated stiffness in compression. 

Proper energy dissipation is guaranteed by the crack-band approach.

The following algorithm is used to compute the stress vector (in each time step and until the iteration criteria are met):


.. _compute-effective-stress:

Compute effective stress
~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

   \sigma_{\mathrm{eff},k+1} = \sigma_{\mathrm{eff},k} + \bar{E}_k \boldsymbol{D}_\nu (\Delta \varepsilon_k - \Delta \varepsilon_k'' - \Delta \varepsilon_{sh,k} - \Delta \varepsilon_{T,k})

where :math:`\sigma_{\mathrm{eff},k}` is the effective stress vector in the preceding step, :math:`\bar{E}_k` is the incremental stiffness, :math:`\boldsymbol{D}_\nu` is a unit stiffness matrix and :math:`\Delta \varepsilon_k`, :math:`\Delta \varepsilon_k''`, :math:`\Delta \varepsilon_{sh,k}`, :math:`\Delta \varepsilon_{T,k}` are the increments of the total strain, strain due to creep, shrinkage strain (sum of drying and autogenous shrinkage) and thermal strain, respectively.

.. _compute-principal-effective-stresses:

Compute principal effective stresses :math:`\sigma_{\mathrm{eff,1}}`, :math:`\sigma_{\mathrm{eff,2}}`, :math:`\sigma_{\mathrm{eff,3}}`

.. _evaluate-equivalent-deformation:

Evaluate equivalent deformation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

   \tilde \varepsilon = \mathrm{max}(\sigma_{\mathrm{eff,1}}, \sigma_{\mathrm{eff,2}}, \sigma_{\mathrm{eff,3}})/E

.. _evaluate-corresponding-damage:

Evaluate corresponding damage
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if :math:`\tilde \varepsilon > \varepsilon_0`

.. math::

   \omega = 1 - \frac{\varepsilon_0}{\tilde \varepsilon} \exp \left(-\frac{\tilde \varepsilon - \varepsilon_0}{\varepsilon_{c,f} - \varepsilon_0} \right)

where :math:`\varepsilon_{c,f} = \frac{G_f}{f_t \: h}`

and :math:`h` is the characteristic length of the finite element in the direction of the biggest principal stress.

.. _compute-principal-nominal-stresses:

Compute principal nominal stresses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   :depth: 1

   - approach \#1:
     for :math:`i = 1, 2, 3`
     if :math:`\sigma_{\mathrm{eff},i} > 0`, :math:`\sigma_{i} = (1-\omega) \sigma_{\mathrm{eff},i}`
     else :math:`\sigma_{i} = \sigma_{\mathrm{eff},i}`

   - approach \#2:
     for :math:`i = 1, 2, 3`
     :math:`\sigma_{i} = (1-\omega) \sigma_{\mathrm{eff},i}`

.. _construct-stress-vector:

Construct the stress vector in the original configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

   \boldsymbol {\sigma} = \boldsymbol{T} \boldsymbol {\sigma}_{princ}

where :math:`\boldsymbol{T}` is the stress transformation matrix.

The tensile strength and fracture energy can be defined either as fixed, time-independent values (parameters :math:`ft`, :math:`gf`) or as variable, depending on the equivalent hydration time :math:`t_e`. The latter case is activated by keyword :math:`timeDepFracturing`; tensile strength and fracture energy are then linked to the current value of the mean compressive strength using (slightly modified) formulae from the :sl:`fib` Model Code 2010.

For :math:`0 \: \mathrm{MPa} \leq f_{cm}(t) \leq 20` MPa the tensile strength is computed from a linear function 

.. math::

   f_{tm}(t) = 0.015 \cdot 12^{2/3} \cdot f_{cm}(t) = 0.07862 f_{cm}(t)

which starts at the origin and at :math:`f_{cm}(t) = 20` MPa connects to expression

.. math::

   f_{tm}(t) = 0.3 (f_{cm}(t) - 8 \: \mathrm{MPa})^{2/3} 

which is valid for :math:`20 \: \mathrm{MPa} \leq f_{cm}(t) \leq 58` MPa.
Finally, for :math:`f_{cm}(t) > 58` MPa 

.. math::

   f_{tm}(t) = 2.12 \: \ln \left( 1 + 0.1 f_{cm}(t) \right ) 

Following the guidelines from the Model Code 2010 (section 5.1.5.2),
the fracture energy :math:`G_f` can be estimated from the mean compressive
strength at 28 days using an empirical formula

.. math::
   :label: gf

   G_{f,28} = 73 \times \left( f_{cm,28} \right)^{0.18}

In this function, the mean compressive strength must be provided in MPa and
the resulting fracture energy is in N/m.

The growth in fracture
energy is approximately proportional to the tensile strength. 
The current value of fracture energy is very simply computed as 

.. math::
   :label: gf_time

   G_f(t) = G_{f,28} \times f_{tm}(t) / f_{tm,28}(28)

The development of the mean compressive strength in time is adopted
entirely from the Model Code 2010 (see section 5.1.9.1)

.. math::
   :label: fcm

   f_{cm}(t) = f_{cm,28} \exp{\left( fib\_s \: \left( 1 - \sqrt{28/t_e} \right) \right) }

where :math:`fib\_s` is a cement-type dependent coefficient (0.38 for cement 32.5 N; 0.25 for 32.5 R and 42.5 N; 0.2 for 42.5 R and 52.5 N and R), :math:`t_e` is the equivalent age and :math:`f_{cm,28}` is the mean compressive strength at the age of 28 days. 

For some concretes the prediction formulae do not provide accurate predictions, therefore it is possible to scale the evolution of tensile strength and fracture energy by providing their 28-day values :param:`ft28` and :param:`gf28`. 

The model description and parameters are summarized in :numref:`mps_dam_table`.

.. table:: MPS damage--summary.
   :name: mps_dam_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | MPS damage model for concrete creep with cracking                                            |
   |                    | **(additionaly all parameters from Table :numref:`mps_table`)**                              |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`MPSDamMat` :optelemparam:`ft{rn}` :optelemparam:`gf{rn}`                          |
   |                    | :optelemstring:`timeDepFracturing` :optelemparam:`fib_s{rn}` :optelemparam:`ft28{rn}`        |
   |                    | :optelemparam:`gf28{rn}` :optelemparam:`damlaw{in}` :optelemstring:`isotropic`               |
   |                    | :optelemparam:`checksnapback{in}`                                                            |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`ft` tensile strength (constant in time)                                            |
   |                    | - :param:`gf` fracture energy (constant in time)                                             |
   |                    | - :param:`timeDepFracturing` string activating *fib* MC 2010 prediction for tensile strength |
   |                    |   and fracture energy + their time evolution                                                 |
   |                    | - :param:`ft28` manual override for the *fib* MC 2010 prediction of hydration-dependent      |
   |                    |   tensile strength                                                                           |
   |                    | - :param:`gf28` manual override for the *fib* MC 2010 prediction of hydration-dependent      |
   |                    |   fracture energy                                                                            |
   |                    | - :param:`fib_s` cement-type dependent coefficient (compulsory with                          |
   |                    |   :param:`timeDepFracturing`)                                                                |
   |                    | - :param:`damageLaw` traction-separation law: 0 = exponential softening (default), 1 =       |
   |                    |   linear softening, 6 = damage disabled                                                      |
   |                    | - :param:`isotropic` string activating same reduction of stiffness after the onset of        |
   |                    |   cracking both in in compression and tension (enables to use secant stiffness matrix for    |
   |                    |   faster convergence)                                                                        |
   |                    | - :param:`checkSnapBack` switch for snap-back checking; 1 = activated (default), 0 =         |
   |                    |   deactiavated                                                                               |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat, PlaneStress, PlaneStrain                                                              |
   +--------------------+----------------------------------------------------------------------------------------------+


Microplane model M4 - Microplane\_M4
------------------------------------


Model M4 covers inelastic behavior of concrete under complex triaxial stress states. It is based on the microplane concept and can describe softening. However, objectivity with respect to element size is not ensured -- the parameters need to be manually adjusted to the element size. Since the tangent stiffness matrix is not available, elastic stiffness is used. This can lead to a very slow convergence when used within an implicit approach. The model parameters are summarized in :numref:`m4_table`.

.. table:: Microplane model M4 -- summary.
   :name: m4_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | M4 material model                                                                            |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`Microplane_M4` :elemparam:`nmp{in}` :elemparam:`c3{rn}` :elemparam:`c20{rn}`      |
   |                    | :elemparam:`k1{rn}` :elemparam:`k2{rn}` :elemparam:`k3{rn}` :elemparam:`k4{rn}`              |
   |                    | :elemparam:`E{rn}` :elemparam:`n{rn}`                                                        |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`nmp` number of microplanes, supported values are 21, 28 and 61                     |
   |                    | - :param:`n` Poisson ratio                                                                   |
   |                    | - :param:`E` Young modulus                                                                   |
   |                    | - :param:`c3`,:param:`c20`, :param:`k1`, :param:`k2`, :param:`k3`, :param:`k4` model         |
   |                    |   parameters                                                                                 |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat                                                                                        |
   +--------------------+----------------------------------------------------------------------------------------------+

.. _sec:cdpm:

Damage-plastic model for concrete - ConcreteDPM
-----------------------------------------------

This model, developed by Grassl and Jirásek for failure of concrete under general triaxial stress, is described in detail in [GraJir]_. It belongs to the class of damage-plastic models with yield condition formulated in terms of the effective stress :math:`\bar{\vsig} = \mD_e : (\veps - \veps_p)`. The stress-strain law is postulated in the form

.. math::
   :label: cdpm1-1

   \vsig = (1 - \omega) \bar{\vsig} = (1 - \omega) \mD_e : (\veps - \veps_p)

where :math:`\mD_e` is the elastic stiffness tensor and :math:`\omega` is a scalar damage parameter. The plastic part of the model consists of a three-invariant yield condition, nonassociated flow rule, and pressure-dependent hardening law. For simplicity, damage is assumed to be isotropic. In contrast to pure damage models with damage driven by the total strain, here the damage is linked to the evolution of plastic strain.

The *yield surface* is described in terms of the cylindrical coordinates in the principal effective stress space (Haigh-Westergaard coordinates), which are the volumetric effective stress :math:`\bar{\sigma}_{\rm V} = {I_1(\bar{\vsig})}/{3}`, the norm of the deviatoric effective stress :math:`\bar{\rho} = \sqrt{2 J_2(\bar{\vsig})}`, and the Lode angle :math:`\theta` defined by the relation

.. math::
   :label: lodeangle

   \cos 3 \theta = \frac{3 \sqrt{3}}{2}\frac{J_3}{J_2^{3/2}}

where :math:`J_2` and :math:`J_3` are the second and third deviatoric invariants.

The yield function

.. math::
   :label: eq:yieldSurface

   f_{\rm p}(\bar\sigma_{\rm V},\bar{\rho},\bar{\theta};\kappa_{\rm p})&=&\left(\left[1-q_{\rm{h}}(\kappa_{\rm p})\right]\left( \frac{\bar{\rho}} {\sqrt{6}\fc} + \frac{\bar{\sigma}_{\rm V}} {\fc} \right)^2 + \sqrt{\frac{3}{2}} \frac {\bar{\rho}}{\fc} \right)^2 +\\
   && 
   +m_0 q_{\rm{h}}^2(\kappa_{\rm p}) \left(\frac{\bar{\rho}r(\bar{\theta}) }{\sqrt{6}\fc} + \frac{\bar{\sigma}_{\rm V}}{\fc} \right) - q_{\rm{h}}^2(\kappa_{\rm p})

depends on the effective stress (which enters in the form of cylindrical coordinates) and on the hardening variable :math:`\kappa_{\rm p}` (which enters through a dimensionless variable :math:`q_{\rm h}`). Parameter :math:`\fc` is the uniaxial compressive strength. 

Note that, under uniaxial compression characterized by axial stress :math:`\bar{\sigma}<0`, we have :math:`\bar{\sigma}_{\rm V}=\bar{\sigma}/3`, :math:`\bar{\rho}=-\sqrt{2/3}\,\bar{\sigma}` and :math:`\bar{\theta}=60^o`. The yield function then reduces to :math:`f_{\rm p}=(\bar{\sigma}/\fc)^2-q_{\rm{h}}^2`. This means that function :math:`q_{\rm{h}}` describes the evolution of the uniaxial compressive yield stress normalized by its maximum value, :math:`\fc`.

The evolution of the yield surface during hardening is presented in :numref:`fig:surfaceMerDev`. The parabolic shape of the meridians (:numref:`fig:surfaceMerDev`\ a) is controlled by the hardening variable :math:`q_{\rm h}` and the friction parameter :math:`m_0`. The initial yield surface is closed, which allows modeling of compaction under highly confined compression. The initial and intermediate yield surfaces have two vertices on the hydrostatic axis but the ultimate yield surface has only one vertex on the tensile part of the hydrostatic axis and opens up along the compressive part of the hydrostatic axis. The deviatoric sections evolve as shown in :numref:`fig:surfaceMerDev`\ b, and their final shape at full hardening is a rounded triangle at low confinement and almost circular at high confinement. The shape of the deviatoric section is controlled by the Willam-Warnke function

.. math::
   :label: bf103

   r(\theta)=\frac{4(1-e^2)\cos^2\theta+(2e-1)^2}{2(1-e^2)\cos\theta+(2e-1)\sqrt{4(1-e^2)\cos^2\theta+5e^2-4e}}

The eccentricity parameter :math:`e` that appears in this function, as well as the friction parameter :math:`m_0`, are calibrated from the values of uniaxial and equibiaxial compressive strengths and uniaxial tensile strength.

.. figure:: /figures/meridians.png
   :align: center
   :width: 100%
   :alt: Evolution of the yield surface during hardening
   :name: fig:surfaceMerDev

   Evolution of the yield surface during hardening:
   a) meridional section, b) deviatoric section for a constant volumetric effective stress of :math:`\bar{\sigma}_{\rm V} = -\fc/3`

The maximum size of the elastic domain is attained when the variable :math:`q_{\rm h}` is equal to one (which is its maximum value, as follows from the hardening law, to be specified in :eq:`eq:hardeningLaw`). The yield surface is then described by the equation

.. math::
   :label: eq:failureSurface

   f_{\rm p}\left(\bar{\sigma}_{\rm V},\bar{\rho},\bar{\theta};1\right) \equiv \frac{3}{2} \frac{\bar{\rho}^2}{\fc^2} + m_0 \left(\frac{\bar{\rho} }{\sqrt{6}\fc}r(\bar{\theta}) + \frac{\bar{\sigma}_{\rm V}}{\fc} \right) - 1 = 0

The **flow rule** 

.. math::
   :label: eq:flowRuleDPM

   \dot{\mbf{\eps}}_{\rm p} = \dot{\lambda} \frac{\partial g_{\rm p}}{\partial\bar{\sigma}} 

is non-associated, which means that the yield function :math:`f_{\rm p}` and the plastic potential 

.. math::
   :label: eq:plasticPotential

   g_{\rm p}(\bar{\sigma}_{\rm V},\bar{\rho};\kappa_{\rm p}) &= \left(\left[1-q_{\rm{h}}(\kappa_{\rm p})\right] \left( \frac{\bar{\rho}} {\sqrt{6}\fc} + \frac{\bar{\sigma}_{\rm V}}{\fc} \right)^2 + \sqrt{\frac{3}{2}} \frac {\bar{\rho}}{\fc} \right)^2 + \\
   &+q_{\rm{h}}^2(\kappa_{\rm p}) \left( \frac{m_0 \bar{\rho}}{\sqrt{6}\fc} + \frac{m_{\rm g}(\bar{\sigma}_{\rm V})}{\fc} \right)

do not coincide and, therefore, the direction of the plastic flow :math:`\partial g_{\rm p}/\partial\bar{\vsig}` is not normal to the yield surface. The ratio of the volumetric and the deviatoric parts of the flow direction is controlled by function :math:`m_{\rm g}`, which depends on the volumetric stress and is defined as 

.. math::
   :label: eq:mg

   m_{\rm g}(\bar{\sigma}_{\rm V})=A_{\rm g}B_{\rm g}\fc \exp{\frac{\bar{\sigma}_{\rm V} -\ft/3}{B_{\rm g}\fc}}

where :math:`A_{\rm g}` and :math:`B_{\rm g}` are model parameters that are determined from certain assumptions on the plastic flow in uniaxial tension and compression.

The dimensionless variable :math:`q_{\rm h}` that appears in the yield function (:eq:`eq:yieldSurface`) is a function of the hardening variable :math:`\kappa_{\rm p}`. It controls the size and shape of the yield surface and, thereby, of the elastic domain. The **hardening law** is given by

.. math::
   :label: eq:hardeningLaw

   \qh(\kappa_{\rm p}) = \left\{ \begin{array}{ll}
   {\qh}_0 + (1-{\qh}_0)\kappa_{\rm p}({\kappa_{\rm p}}^2 - 3\kappa_{\rm p}+3) & \mbox{if $\kappa_{\rm p} < 1$} \\
   1 & \mbox{if $\kappa_{\rm p} \ge 1$}
   \end{array}
   \right.

The initial inclination of the hardening curve (at :math:`\kappa_{\rm p}=0`) is positive and finite, and the inclination at peak (i.e., at :math:`\kappa_{\rm p} = 1`) is zero.

The evolution law for the hardening variable, (In the original paper [GraJir]_, equation (:numref:)eq:dotkappap`) was written with :math:`\cos^2\bar{\theta}` instead of :math:`(2\cos\bar{\theta})^2`, but all the results presented in that paper were computed with OOFEM using an implementation based on (:eq:`eq:dotkappap`).`, is given by:

.. math::
   :label: eq:dotkappap

   \dot{\kappa}_{\rm {p}} = \frac {\| \dot{\veps}_{\rm{p}}\|}{x_{\rm{h}}\left(\bar{\sigma}_{\rm V} \right)}(2\cos\bar{\theta})^2 

This sets the rate of the hardening variable equal to the norm of the plastic strain rate scaled by a hardening ductility measure:

.. math::
   :label: eq:hardeningMeasure

   x_{\rm{h}}\left(\bar{\sigma}_{\rm V}\right) = \left\{ \begin{array}{ll}
   A_{\rm{h}} - \left(A_{\rm{h}} - B_{\rm{h}} \right)\exp{\left(-R_{\rm h}(\bar{\sigma}_{\rm V})/C_{\rm h}\right)} & \mbox{if $R_{\rm h}(\bar{\sigma}_{\rm V}) \geq 0 $} \\[5mm]
   E_{\rm h}\exp({R_{\rm h}(\bar{\sigma}_{\rm V})/F_{\rm h}})+D_{\rm h} & \mbox{if $R_{\rm h}(\bar{\sigma}_{\rm V}) < 0$}
   \end{array}
   \right.

The dependence of the scaling factor :math:`x_{\rm{h}}` on the volumetric effective stress :math:`\bar{\sigma}_{\rm V}` is constructed such that the model response is more ductile under compression.

The variable

.. math::

   R_{\rm h}(\bar{\sigma}_{\rm V}) = -\frac{\bar{\sigma}_{\rm V}}{\fc}-\frac{1}{3}

is a linear function of the volumetric effective stress.

Model parameters :math:`A_{\rm h}, B_{\rm h}, C_{\rm h}` and :math:`D_{\rm h}` are calibrated from the values of strain at peak stress under uniaxial tension, uniaxial compression and triaxial compression, whereas the parameters

.. math::

   E_{\rm h} &= B_{\rm h} - D_{\rm h}
   \\
   F_{\rm h} &= \frac{\left(B_{\rm h} -D_{\rm h}\right)C_{\rm h}}{B_{\rm h} - A_{\rm h}}

are determined from the conditions of a smooth transition between the two parts of equation (:eq:`eq:hardeningMeasure`) at :math:`R_{\rm h} = 0`.

For the present model, the *evolution of damage* starts after full saturation of plastic hardening, i.e., at :math:`\kappa_{\rm p}=1`. This greatly facilitates calibration of model parameters, because the strength envelope is fully controlled by the plastic part of the model and damage affects only the softening behavior.

In contrast to pure damage models, damage is assumed to be driven by the plastic strain, more specifically by its volumetric part, which is closely related to cracking. To slow down the evolution of damage under compressive stress states, the damage-driving variable :math:`\kappa_{\rm d}` is not set equal to the volumetric plastic strain, but it is defined incrementally by the rate equation

.. math::
   :label: eq:equivStrain

   \dot{\kappa}_{\rm d} = \left\{ \begin{array}{ll}
   0 & \mbox{if $\kappa_{\rm p}< 1$} \\
   {\rm Tr}(\dot{\veps}_{\rm p})/x_{\rm s}\left(\bar{\sigma}_{\rm V}\right) & \mbox{if $\kappa_{\rm p}\ge 1$}
   \end{array}
    \right.

where

.. math::
   :label: eq:damageDuctility

   x_{\rm s}(\bar{\sigma}_{\rm V}) = \left\{\begin{array}{ll}
   1 + A_{\rm s}R_{\rm s}^2(\bar{\sigma}_{\rm V}) & \mbox{if $R_{\rm s}(\bar{\sigma}_{\rm V}) < 1$}\\[5mm]
   1 - 3 A_{\rm s} + 4 A_{\rm s} \sqrt{R_{\rm s}(\bar{\sigma}_{\rm V})} & \mbox{if $R_{\rm s}(\bar{\sigma}_{\rm V}) \geq 1$}
   \end{array}
   \right .

is a softening ductility measure. Parameter :math:`A_{\rm s}` is determined from the softening response in uniaxial compression. The dimensionless variable :math:`R_{\rm s}=\dot{\varepsilon}_{\rm pV}^{-}/\dot{\varepsilon}_{\rm {pV}}` is defined as the ratio between the “negative” volumetric plastic strain rate

.. math::
   :label: eq:xs

   \dot{\varepsilon}_{\rm pV}^{-} = \displaystyle{\sum_{I = 1}^3 \langle-\dot{\varepsilon}_{\rm{p}I}\rangle}

and the total volumetric plastic strain rate :math:`\dot{\varepsilon}_{\rm {pV}}`. This ratio depends only on the flow direction :math:`\partial g_{\rm p}/\partial\bar{\vsig}`, and thus :math:`R_{\rm s}` can be shown to be a unique function of the volumetric effective stress. In (:eq:`eq:xs`), :math:`\dot{\varepsilon}_{\rm{p}I}` are the principal components of the rate of plastic strains and :math:`\left<\cdot\right>` denotes the McAuley brackets (positive-part operator).

For uniaxial tension, for instance, all three principal plastic strain rates are nonnegative, and so :math:`\dot{\varepsilon}_{\rm pV}^{-} =0`, :math:`R_{\rm s}=0` and :math:`x_{\rm s}=1`. This means that under uniaxial tensile loading we have :math:`\kappa_{\rm d}=\kappa_{\rm p}-1`. On the other hand, under compressive stress states the negative principal plastic strain rates lead to a ductility measure :math:`x_{\rm s}` greater than one and the evolution of damage is slowed down. It should be emphasized that the flow rule for this specific model is constructed such that the volumetric part of plastic strain rate at the ultimate yield surface cannot be negative.

.. table:: Damage-plastic model for concrete -- summary
   :name: dpm_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Damage-plastic model for concrete                                                            |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`ConcreteDPM` :elemparam:`d{rn}` :elemparam:`E{rn}` :elemparam:`n{rn}`             |
   |                    | :elemparam:`tAlpha{rn}` :elemparam:`ft{rn}` :elemparam:`fc{rn}` :elemparam:`wf{rn}`          |
   |                    | :elemparam:`Gf{rn}` :elemparam:`ecc{rn}` :elemparam:`kinit{rn}` :elemparam:`Ahard{rn}`       |
   |                    | :elemparam:`Bhard{rn}` :elemparam:`Chard{rn}` :elemparam:`Dhard{rn}` :elemparam:`Asoft{rn}`  |
   |                    | :elemparam:`helem{rn}` :elemparam:`href{rn}` :elemparam:`dilation{rn}`                       |
   |                    | :elemparam:`yieldtol{rn}` :elemparam:`newtoniter{in}`                                        |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`d` material density                                                                |
   |                    | - :param:`E` Young modulus                                                                   |
   |                    | - :param:`n` Poisson ratio                                                                   |
   |                    | - :param:`tAlpha` thermal dilatation coefficient                                             |
   |                    | - :param:`ft` uniaxial tensile strength :math:`f_t`                                          |
   |                    | - :param:`fc` uniaxial compressive strength                                                  |
   |                    | - :param:`wf` parameter :math:`w_f` that controls the slope of the softening branch (serves  |
   |                    |   for the evaluation of :math:`\eps_f=w_f/h` to be used in (:eq:`eq:damageLaw`))             |
   |                    | - :param:`Gf` fracture energy, can be specified instead of wf, it is converted to            |
   |                    |   :math:`w_f=G_f/f_t`                                                                        |
   |                    | - :param:`ecc` eccentricity parameter :math:`e` from (:eq:`bf103`), optional, default        |
   |                    |   value 0.525                                                                                |
   |                    | - :param:`kinit` parameter :math:`{\qh}_0` from (:eq:`eq:hardeningLaw`), optional,           |
   |                    |   default value 0.1                                                                          |
   |                    | - :param:`Ahard` parameter :math:`A_{\rm h}` from (:eq:`eq:hardeningMeasure`), optional,     |
   |                    |   default value 0.08                                                                         |
   |                    | - :param:`Bhard` parameter :math:`B_{\rm h}` from (:eq:`eq:hardeningMeasure`), optional,     |
   |                    |   default value 0.003                                                                        |
   |                    | - :param:`Chard` parameter :math:`C_{\rm h}` from (:eq:`eq:hardeningMeasure`), optional,     |
   |                    |   default value 2                                                                            |
   |                    | - :param:`Dhard` parameter :math:`D_{\rm h}` from (:eq:`eq:hardeningMeasure`), optional,     |
   |                    |   default value :math:`10^{-6}`                                                              |
   |                    | - :param:`Asoft` parameter :math:`A_{\rm s}` from (:eq:`eq:damageDuctility`), optional,      |
   |                    |   default value 15                                                                           |
   |                    | - :param:`helem` element size :math:`h`, optional (if not specified, the actual element size |
   |                    |   is used)                                                                                   |
   |                    | - :param:`href` reference element size :math:`h_{ref}`, optional (if not specified, the      |
   |                    |   standard adjustment of the damage law is used)                                             |
   |                    | - :param:`dilation` dilation factor (ratio between lateral and axial plastic strain rates in |
   |                    |   the softening regime under uniaxial compression), optional, default value -0.85            |
   |                    | - :param:`yieldtol` tolerance for the implicit stress return algorithm, optional, default    |
   |                    |   value :math:`10^{-10}`                                                                     |
   |                    | - :param:`newtoniter` maximum number of iterations in the implicit stress return algorithm,  |
   |                    |   optional, default value 100                                                                |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat                                                                                        |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `sm/con1dpm1.in                                                                              |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/con1dpm1.in>`_,               |
   |                    | `sm/con1dpm2.in                                                                              |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/con1dpm2.in>`_,               |
   |                    | `sm/con1dpm3.in                                                                              |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/con1dpm3.in>`_                |
   +--------------------+----------------------------------------------------------------------------------------------+

The relation between the damage variable :math:`\omega` and the internal variable :math:`\kappa_{\rm d}` (maximum level of equivalent strain) is assumed to have the exponential form

.. math::
   :label: eq:damageLaw

   \omega = 
   1-\exp\left(-\kappa_{\rm d}/\varepsilon_f\right) 

where :math:`\varepsilon_f` is a parameter that controls the slope of the softening curve. In fact, equation (:eq:`eq:damageLaw`) is used by the nonlocal version
of the damage-plastic model, with :math:`\kappa_{\rm d}` replaced by its weighted
spatial average (not yet available in the public version of OOFEM).
For the local model, it is necessary to adjust softening according
to the element size, otherwise the results would suffer by pathological 
mesh sensitivity. It is assumed that localization takes place at the
peak of the stress-strain diagram, i.e., at the onset of damage. 
After that, the strain is decomposed into the distributed
part, which corresponds to unloading from peak, and the localized part,
which is added if the material is softening.
The localized part of strain 
is transformed into an equivalent crack opening, :math:`w`, which is 
under uniaxial tension linked
to the stress by the exponential law

.. math::
   :label: eq:damageLaw2

   \sigma = (1-\omega) f_t =
   f_t \exp\left(-w/w_f\right) 

Here, :math:`f_t` is the uniaxial tensile strength and :math:`w_f` is the characteristic
crack opening, playing a similar role to :math:`\varepsilon_f`.
Under uniaxial tension, the localized strain can be expressed as the sum
of the post-peak plastic strain (equal to variable :math:`\kappa_{\rm d}`) 
and the unloaded part of elastic strain (equal to :math:`\omega f_t/E`).
Denoting the effective element size as :math:`h`, we can write

.. math::

   w = h(\kappa_{\rm d}+\omega f_t/E)

and substituting this into (:eq:`eq:damageLaw2`), we obtain a nonlinear
equation

.. math::
   :label: eq:damageLaw4

   1-\omega = \exp\left(-\frac{h}{w_f}\left(\kappa_{\rm d}+\omega f_t/E\right)\right)

from which the damage variable :math:`\omega` corresponding to the given internal
variable :math:`\kappa_{\rm d}` can be computed by Newton iteration.
The effective element size :math:`h` is obtained by projecting the element onto
the direction of the maximum principal strain at the onset of cracking,
and afterwards it is held fixed. The evaluation of :math:`\omega` from
:math:`\kappa_{\rm d}` is no longer explicit, but the resulting load-displacement curve
of a bar under uniaxial tension is totally independent of the mesh size.
A simpler approach would be to use (:eq:`eq:damageLaw`) with 
:math:`\varepsilon_f = w_f/h`, but then the scaling would not be perfect and
the shape of the load-displacement curve (and also the dissipated energy)
would slightly depend on the mesh size. With the present approach,
the energy per unit sectional area dissipated under uniaxial tension 
is exactly :math:`G_f = w_f f_t`. The input parameter controling the damage law
can be either the characteristic crack opening :math:`w_f`, 
or the fracture energy :math:`G_f`. If both are specified, :math:`w_f` is used and
:math:`G_f` is ignored. If only :math:`G_f` is specified, :math:`w_f` is set to :math:`G_f/f_t`.

The onset of damage corresponds to the peak of the stress-strain diagram under proportional loading, when the ratios of the stress components are fixed. This is the case, e.g., for uniaxial tension, uniaxial compression, or shear under free expansion of the material (with zero normal stresses). However, for shear under confinement, the shear stress can rise even after the onset of damage, due to increasing hydrostatic pressure, which increases the mobilized friction. It has been observed that the standard approach leads to strong sensitivity of the peak shear stress to the element size. To reduce this pathological effect, a modified approach has been implemented. The second-order work (product of stress increment and strain increment) is checked after each step, and the element-size dependent adjustment of the damage law is applied only after the second-order work becomes negative. Up to this stage, the damage law corresponds to a fixed reference element size, which is independent of the actual size of the element. This size is set by the optional parameter :param:`href`. If this parameter is not specified, the standard approach is used. For testing purposes, one can also specify the actual element size, :param:`helem`, as a “material property”. If this parameter is not specified, the element size is computed for each element separately and represents its actual size.

The damage-plastic model contains 15 parameters, but only 6 of them need to be actually calibrated for different concrete types, namely Young's modulus :math:`E`, Poisson's ratio :math:`\nu`, tensile strength :math:`f_{\rm t}`, compressive strength :math:`f_{\rm c}`, parameter :math:`w_f` (or fracture energy :math:`G_f`), and parameter :math:`A_{\rm s}` in the ductility measure (:eq:`eq:damageDuctility`) of the damage model. The remaining parameters can be set to their default values specified in [GraJir]_.

The model parameters are summarized in :numref:`dpm_table`. Note that it is possible to specify the “size” of finite element, :math:`h`, which (if specified) replaces the actual element size in (:eq:`eq:damageLaw4`). The usual approach is to consider :math:`h` as the actual element size (evaluated automatically by OOFEM), in which case the optional parameter :math:`h` is missing (or is set to 0., which has the same effect in the code). However, for various studies of mesh sensitivity, it is useful to have the option of specifying :math:`h` as an input “material” value.

If the element is too large, it may become too brittle and local snap-back
occurs in the stress-strain diagram, which is not acceptable. 
In such a case, an error message is issued and the program execution
is terminated. The maximum admis
sible element size 

.. math::

   h_{\max} = \frac{EG_f}{f_t^2} = \frac{Ew_f}{f_t}

happens to be equal to Hillerborg's characteristic material length. 
For typical concretes it is in the order of a few hundred mm. 
If the condition :math:`h < h_{\max}` is violated, the mesh needs to be refined.
Note that the effective element size :math:`h` is obtained by projecting the element.
For instance, if the element is a cube of edge length 100 mm, its effective
size in the direction of the body diagonal can be 173 mm. 

.. _sec:cdpm2:

Damage-plastic model for concrete - CDPM2
-----------------------------------------

This model is an extension of the ConcreteDPM presented in :ref:`sec:cdpm`.
CDPM2 has been developed by Grassl and coworkers for modelling the failure of concrete for both static and dynamic loading. It is is described in detail in [GraXenNys13]_.
The main differences between CDPM2 and ConcreteDPM are that in CDPM2

- the plasticity part exhibits hardening once damage is active,
- two independent damage parameters, which describe tensile and compressive damage, are introduced.
  These extensions provide more flexibility, and the shapes
  of stress-strain curves under various loading scenarios, including unloading, can be better adapted. The price to pay
  is an increased number of parameters and thus more difficult calibration.

CDPM2 deals with the same definition of effective stress as ConcreteDPM, namely

.. math::

   \bar{\boldsymbol{\sigma}} = \boldsymbol{D}_e:(\boldsymbol{\varepsilon}-\boldsymbol{\varepsilon}_p)

which means that the effective stress is computed using the plastic part of the model. The effect of damage is incorporated in the transformation from effective to nominal stress. While ConcreteDPM uses a simple reduction factor :math:`1-\omega` (see :eq:`cdpm1-1`), CDPM2 splits the effective stress into the positive part, :math:`\bar{\boldsymbol{\sigma}}_{\rm t}`, and the negative part, :math:`\bar{\boldsymbol{\sigma}}_{\rm c}`, and then uses two damage variables, :math:`\omega_t` and :math:`\omega_c` (for tension and for compression). The nominal stress is thus evaluated as

.. math::
   :label: cdpm2:efftonominal

   \boldsymbol{\sigma} = \left(1-\omega_{\rm t}\right) \bar{\boldsymbol{\sigma}}_{\rm t} + \left(1-\omega_{\rm c}\right) \bar{\boldsymbol{\sigma}}_{\rm c}

The values of :math:`\omega_{\rm t}` and :math:`\omega_{\rm c}` range from 0 (undamaged) to 1 (fully damaged). The decomposition of the effective stress tensor into the positive and negative parts is based on principal values. If :math:`\bar{\sigma}_I`, :math:`I=1,2,3`, are the principal effective stresses and :math:`\boldsymbol{n}_I`, :math:`I=1,2,3`, are the corresponding principal directions, the effective stress tensor

.. math::

   \bar{\boldsymbol{\sigma}} = \sum_{I=1}^3 \bar{\sigma}_I\,\boldsymbol{n}_I\otimes\boldsymbol{n}_I

is decomposed into

.. math::

   \bar{\boldsymbol{\sigma}}_t &=& \sum_{I=1}^3 \langle\bar{\sigma}_I\,\rangle\boldsymbol{n}_I\otimes\boldsymbol{n}_I \\
   \bar{\boldsymbol{\sigma}}_c &=& -\sum_{I=1}^3 \langle-\bar{\sigma}_I\,\rangle\boldsymbol{n}_I\otimes\boldsymbol{n}_I

where :math:`\langle\ldots\rangle` are Macaulay brackets (positive-part operator).

Let us now describe in more detail the **plastic part** of the model. 

The **yield condition** is formulated in terms of the effective stress :math:`\bar{\boldsymbol{\sigma}}` and therefore is unaffected by damage. The yield function is given by

.. math::
   :label: eqcdpm2:yieldSurface

   \begin{split}
   & f_{\rm p}(\bar{\sigma}_{\rm V},\bar\rho,\bar\theta;\kappa_{\rm p})=  \left\{\left[1-q_{\rm{h1}}(\kappa_{\rm p})\right]\left( \frac{\bar{\rho}} {\sqrt{6}f_{\rm c}} + \frac{\bar{\sigma}_{\rm V}} {f_{\rm c}} \right)^2 + \sqrt{\frac{3}{2}} \frac {\bar{\rho}}{f_{\rm c}} \right\}^2 \\
   & +m_0 q^2_{\rm{h1}}(\kappa_{\rm p})q_{\rm{h2}}(\kappa_{\rm p}) \left[\frac{\bar{\rho} }{\sqrt{6}f_{\rm c}}r(\cos{\bar{\theta}}) + \frac{\bar{\sigma}_{\rm V}}{f_{\rm c}} \right] - q_{\rm{h1}}^2(\kappa_{\rm p}) q_{{\rm h2}}^2(\kappa_{\rm p})
   \end{split}

where :math:`\bar\sigma_{\rm V}` is the volumetric effective stress, :math:`\bar\rho` is the norm of the deviatoric effective stress and :math:`\bar\theta` is the Lode angle (another invariant of the effective stress tensor). Parameter :math:`f_{\rm c}` is the uniaxial compressive strength. Evolution of the yield surface is controlled by two dimensionless stress-like hardening variables, :math:`q_{h1}` and :math:`q_{h2}`, which are considered as functions of one dimensionless strain-like hardening variable, :math:`\kappa_p`.

The **hardening laws** read

.. math::
   :label: eq:hardeningLawOne

   q_{\rm h1}(\kappa_{\rm p}) &=& 
   \left\{ \begin{array}{ll} q_{\rm h0} + \left(1-q_{\rm h0}\right) \left( \kappa_{\rm p}^3 - 3 \kappa_{\rm p}^2 + 3 \kappa_{\rm p} \right) - H_{\rm p} \left(\kappa_{\rm p}^3 - 3 \kappa_{\rm p}^2 + 2 \kappa_{\rm p}\right) & \mbox{if $\kappa_{\rm p} < 1$} \\
   1 & \mbox{if $\kappa_{\rm p} \ge 1$}
   \end{array}
   \right\}

.. math::
   :label: eq:hardeningLawTwo

   q_{\rm h2}(\kappa_{\rm p}) &=& \left\{ \begin{array}{ll} 1 & \mbox{if $\kappa_{\rm p} < 1$} \\
   1 + H_{\rm p} (\kappa_{\rm p} - 1)  & \mbox{if $\kappa_{\rm p} \ge 1$}
   \end{array}
   \right\}

During the initial stage of hardening, :math:`\kappa_p` grows from 0 to 1, :math:`q_{h1}` grows from its initial value :math:`q_{h0}` to 1 while :math:`q_{h2}` remains equal to 1. No damage is induced and the yield surface is expanding. A change of behavior occurs when :math:`\kappa_p` becomes equal to 1. For the original ConcreteDPM model, there would not be any further evolution of the yield surface. In CDPM2, linear hardening with constant plastic modulus :math:`H_p` is assumed. Variable :math:`q_{h1}` remains equal to 1 while :math:`q_{h2}` grows proportionally to the increments of :math:`\kappa_p`. If :math:`H_p` is set to zero (i.e., :math:`q_{\rm h2}` remains equal to 1), the same behavior as for the ConcreteDPM is obtained.

The shape of the deviatoric section is controlled by the Willam-Warnke function

.. math::
   :label: eq:rFunction

   r(\cos{\bar{\theta}}) = \frac{4(1-e^2)\cos^2{\bar{\theta}} + (2e-1)^2}{2(1-e^2)\cos{\bar{\theta}} + (2e-1)\sqrt{4(1-e^2)\cos^2{\bar{\theta}}+5e^2 -4e}} 

Here, :math:`e` is the eccentricity parameter.

The friction parameter :math:`m_0` is given by 

.. math::
   :label: eq:frictionM

   m_0 = \dfrac{3 \left(f_{\rm c}^2 - f_{\rm t}^2\right)}{f_{\rm c}f_{\rm t}} \dfrac{e}{e+1}

where :math:`f_{\rm t}` is the tensile strength.

Hardening variable :math:`\kappa_p` is a generalized form of cumulative plastic strain. Its rate corresponds to the rate of the norm of plastic strain rate scaled by a factor that depends on the stress state. This is described by the evolution equation

.. math::
   :label: eq:kappap

   \dot{\kappa}_p = \frac{\Vert\dot{\boldsymbol{\varepsilon}}_p\Vert}{x_h(\bar{\sigma}_V)}\left(2\cos\bar{\theta}\right)^2

in which :math:`\bar{\sigma}_V` is the hydrostatic part of effective stress and :math:`\bar{\theta}` is the Lode angle. Function :math:`x_h` has a complicated form, designed such that :math:`\kappa_p` grows faster for positive hydrostatic stress and more slowly under negative hydrostatic stress. Consequently, the material response is brittle in tension and less brittle or even ductile in compression, especially under confinement. The form of function :math:`x_h` suggested in [GraXenNys13]_ is

.. math::
   :label: eq:xh

   x_h(\bar{\sigma}_V) = \left\{\begin{array}{ll}
   A_h-(A_h-B_h)\exp(-R_h(\bar{\sigma}_V)/C_h) & \hskip 10mm \mbox{if } R_h(\bar{\sigma}_V) \ge 0
   \\[3mm]
   D_h + E_h\exp(R_h(\bar{\sigma}_V)/F_h) & \hskip 10mm \mbox{if } R_h(\bar{\sigma}_V) < 0
   \end{array}\right.

where

.. math::

   R_h(\bar{\sigma}_V) = -\frac{\bar{\sigma}_V}{f_c}-\frac{1}{3}

At the peak of the uniaxial compressive stress-strain curve, we have :math:`\bar{\sigma}_V=-f_c/3` and thus :math:`R_h=0`. Positive values of :math:`R_h` correspond to confined compression. Under compression with no confinement, :math:`x_h` is equal to :math:`B_h`, and under high confinement it is close to :math:`A_h`. Parameter :math:`C_h` controls the transition from :math:`B_h` to :math:`A_h` and if it is high, the transition occurs at higher compressive stresses. In this way, parameters :math:`A_h`, :math:`B_h` and :math:`C_h` can be used to tune up the sensitivity of the model to confining stresses. Parameter :math:`D_h` affects the behavior under combined tension and compression, or under uniaxial or multiaxial tension.

Parameters :math:`E_h` and :math:`F_h` are not considered as independent.

To ensure a smooth transition between cases with positive and negative :math:`R_h`, they need to be set to :math:`E_h=B_h-D_h` and :math:`F_h=(B_h-D_h)C_h/(A_h-B_h)`. Note that all these parameters are dimensionless. In [GraXenNys13]_ it was suggested to calibrate parameters :math:`A_h`, :math:`B_h`, :math:`C_h` and :math:`D_h` from the strains at peak stress under uniaxial tension, uniaxial compression and triaxial compression. Of course, four parameters cannot be uniquely determined from just three experimental values. It is desirable to use experimental data for triaxial compression at several levels of confinement. From the foregoing discussion of the role of individual parameters it follows that :math:`A_h` should be greater than :math:`B_h`, to obtain a positive effect of confinement on ductility, and parameter :math:`D_h` should be smaller than :math:`B_h`, to get a more brittle behavior in tension than in compression.

The **flow rule** is postulated in the non-associated format (Equation :eq:`eq:flowRuleDPM`)
with a plastic potential given by the same expression (Equation :eq:`eq:plasticPotential`)
as for the ConcreteDPM model, just with hardening variable :math:`q_h` formally replaced by :math:`q_{h1}`.
The direction of plastic flow is given by the gradient of the plastic potential,
and it can be split into the volumetric and deviatoric parts:

.. math::
   :label: eq:generalFlow

   \boldsymbol{m} = \dfrac{\partial g}{\partial \bar{\boldsymbol{\sigma}} } = \frac{\partial g}{\partial \bar{\sigma}_{\rm V}} \frac{\partial \bar{\sigma}_{\rm V}}{\partial \bar{\boldsymbol{\sigma}}} + \frac{\partial g}{\partial \bar{\rho}} \frac{\partial \bar{\rho} }{\partial \bar{\boldsymbol{\sigma}}}

Taking into account that :math:`\partial\bar{\sigma}_{\rm V}/\partial\bar{\boldsymbol{\sigma}}=\boldsymbol{\delta}/3` and :math:`\partial\bar{\rho}/\partial\bar{\boldsymbol{\sigma}}=\bar{\mbf{s}}/\bar\rho`, restricting attention to the onset of the post-peak regime (in which :math:`q_{\rm h1} = q_{\rm h2} = 1`) and differentiating the plastic potential (Equation :eq:`eq:plasticPotential`), we rewrite equation (Equation :eq:`eq:generalFlow`) as

.. math::
   :label: eq:detailedFlow

   \mbf{m} = \dfrac{\partial g}{\partial \bar{\boldsymbol{\sigma}}} =  \frac{{\partial}{m_{\rm g}}}{{\partial}{\bar{\sigma}_{\rm V}}} \frac{\boldsymbol{\delta}}{3f_{\rm c}} 
   + \left(\frac{3}{f_{\rm c}} + \frac{m_0}{\sqrt{6}\bar\rho} \right)\frac{\bar{\mbf{s}}}{f_{\rm c}}

Since the flow rule is non-associative, the direction of the plastic flow is not normal to the yield surface. This is important for concrete since an associative flow rule would give an overestimated maximum stress for passive confinement.

**Damage** is initiated when the equivalent strain reaches the threshold :math:`\varepsilon_{0} = f_{\rm t}/E`,
which represents the limit elastic strain under uniaxial tension. 
The evolution of damage variables depends
on certain internal variables linked to the elastic and plastic strains
computed using the plastic part of the model.

The **tensile damage** variable :math:`\omega_t` is a function of three auxiliary internal variables, :math:`\kappa_{dt}`, :math:`\kappa_{dt1}`, and :math:`\kappa_{dt2}`.

Variable :math:`\kappa_{dt}` is the maximum previously reached value of equivalent strain

.. math::
   :label: cdpm2:eps

   \tilde\varepsilon = \frac{\varepsilon_0m_0}{2}\left(\frac{\bar\rho r(\cos\bar\theta)}{\sqrt{6}f_c}+\frac{\bar\sigma_V}{f_c}\right) + \sqrt{\frac{\varepsilon_0^2m_0^2}{4}\left(\frac{\bar\rho r(\cos\bar\theta)}{\sqrt{6}f_c}+\frac{\bar\sigma_V}{f_c}\right)^2+\frac{3\varepsilon_0^2\bar\rho^2}{2f_c^2}}    

where :math:`r` is the function of the Lode angle specified in :eq:`eq:rFunction` and :math:`m_0` is the friction parameter specified in :eq:`eq:frictionM`.

Expression :eq:`cdpm2:eps` has been derived from the yield condition :math:`f_{\rm p} = 0` by setting :math:`q_{\rm {h1}}=1` and :math:`q_{\rm{h2}} = \tilde{\varepsilon}/\varepsilon_{0}`, which leads to a quadratic equation.

Variables :math:`\kappa_{dt1}` and :math:`\kappa_{dt2}` are defined by the rate equations

.. math::

   \dot\kappa_{dt1} &= \left\{\begin{array}{ll}
   \displaystyle\frac{\Vert\dot{\boldsymbol\varepsilon}_p\Vert}{x_s(\bar\rho,\bar\sigma_V)} & \hskip 10mm \mbox{if } \dot\kappa_{dt}>0\;\mbox{ and }\; \kappa_{dt}\ge \varepsilon_0
   \\[3mm]
   0 & \hskip 10mm \mbox{if } \dot\kappa_{dt}=0\;\mbox{ or }\; \kappa_{dt}< \varepsilon_0
   \end{array}\right.
   \\
   \dot\kappa_{dt2} &=  \frac{\dot\kappa_{dt}}{x_s(\bar\rho,\bar\sigma_V)}
   

where 

.. math::
   :label: eq:xs2

   x_s(\bar\rho,\bar\sigma_V) =
   1+(A_s-1)\frac{\sqrt{6}\langle-\bar\sigma_V\rangle}{\bar\rho}

is a ductility measure. For positive hydrostatic stress :math:`\bar\sigma_V`, we have :math:`\langle-\bar\sigma_V\rangle=0` and :math:`x_s=1`. Under compression, :math:`x_s` increases (provided that parameter :math:`A_s` is greater than 1), which slows down the evolution of variables :math:`\kappa_{dt1}` and :math:`\kappa_{dt2}` and consequently of damage. For uniaxial compression, we have :math:`\bar\sigma_V/\bar\rho=-1/\sqrt{6}` and formula :eq:`eq:xs2` gives :math:`x_s=A_s`.

The specific form of the dependence of damage variable :math:`\omega_t` on :math:`\kappa_{dt}`, :math:`\kappa_{dt1}`, and :math:`\kappa_{dt2}` is dictated by the desired shape of the softening part of the uniaxial tensile stress-strain diagram. For instance, bilinear softening is obtained with

.. math::
   :label: cdpm2:omegat

   \omega_t = \left\{\begin{array}{ll}
   \dfrac{(E\kappa_{dt}-f_t)\varepsilon_{f1}-(\sigma_1-f_t)\kappa_{dt1}}{E\kappa_{dt}\varepsilon_{f1}+(\sigma_1-f_t)\kappa_{dt2}} & \hskip 10mm \mbox{ if } 0\le\kappa_{dt1}+\omega_t\kappa_{dt2}\le\varepsilon_{f1}
   \\[3mm]
     \dfrac{E\kappa_{dt}(\varepsilon_f-\varepsilon_{f1})+\sigma_1(\kappa_{dt1}-\varepsilon_f)}{E\kappa_{dt}(\varepsilon_f-\varepsilon_{f1})-\sigma_1\kappa_{dt2}} & \hskip 10mm \mbox{ if } \varepsilon_{f1}<\kappa_{dt1}+\omega_t\kappa_{dt2}\le\varepsilon_{f}
     \\[3mm]
     0 & \hskip 10mm \mbox{ if } \varepsilon_{f}<\kappa_{dt1}+\omega_t\kappa_{dt2}
   \end{array}\right\}

where :math:`\varepsilon_{f1}` and :math:`\sigma_1` is the strain and stress at the knee point of the softening diagram (where the two linear segments meet), and :math:`\varepsilon_{f}` is the strain at which the stress vanishes. The choice of a bilinear softening diagram for tension makes it possible to obtain a steep initial descent followed by a long tail, which would not be that easy with the exponential diagram used in the ConcreteDPM model. In the local version of the model, objectivity with respect to mesh size must be enforced by adjustment of the softening diagram. Numerical parameters :math:`\varepsilon_{f1}` and :math:`\varepsilon_{f}` are then derived from material parameters :math:`w_{f1}` and :math:`w_f` that define the cohesive law in terms of the dependence of stress on crack opening, :math:`w`. Transformation of crack opening into the smeared cracking strain then leads to parameters :math:`\varepsilon_{f1}=w_{f1}/h` and :math:`\varepsilon_{f}=w_f/h` where :math:`h` is the estimated thickness of the computationally resolved localized band.

The **compressive damage** variable :math:`\omega_c` is a function of three auxiliary internal variables, :math:`\kappa_{dc}`, :math:`\kappa_{dc1}`, and :math:`\kappa_{dc2}`.

Variable :math:`\kappa_{dc}` is defined by the rate equation

.. math::

   \dot\kappa_{dc} = \alpha_c\dot{\tilde\varepsilon}

where :math:`\tilde\varepsilon` is the already defined equivalent strain, and the scaling factor :math:`\alpha_c` depends on the type of stress state and is given by

.. math::

   \alpha_c = \frac{\sum_{I=1}^3 \langle-\bar\sigma_I\rangle^2}{\sum_{I=1}^3 \bar\sigma_I^2}

If all principal effective stresses :math:`\sigma_I` are positive (or at least non-negative), we have :math:`\alpha_c=0`, and if all principal effective stresses are negative (or at least non-positive), we have :math:`\alpha_c=1`. For the zero stress state, :math:`\alpha_c` is undefined, but for that particular state it is not needed.

Two other auxiliary internal variables that affect compressive damage are defined by the rate equations

.. math::
   :label: cdpm2:kappa_dc1

   \dot\kappa_{dc1} = \left\{
   \begin{array}{ll}
   \displaystyle\frac{\alpha_c\beta_c\vert\dot{\boldsymbol\varepsilon}_p\vert}{x_s} & \hskip 10mm \mbox{if } \dot\kappa_{dt}>0\;\mbox{ and }\; \kappa_{dt}\ge \varepsilon_0
   \\[3mm]
   0 & \hskip 10mm \mbox{if } \dot\kappa_{dt}=0\;\mbox{ or }\; \kappa_{dt}< \varepsilon_0
   \end{array}\right\}

.. math::
   :label: cdpm2:kappa_dc2

   \dot\kappa_{dc2} =  \frac{\dot\kappa_{dc}}{x_s}

where factor

.. math::
   :label: cdpm2:betac

   \beta_c =\frac{f_tq_{h2}\sqrt{2/3}}{\bar\rho\sqrt{1+2D_f^2}}

depends on the plastic hardening variable :math:`q_{h2}` and effective stress invariant :math:`\bar\rho` and provides a smooth transition from pure damage to damage-plastic softening, which can occur during cyclic loading. The dependence of the compressive damage :math:`\omega_c` on the internal variables is defined such that the relation between the stress and the inelastic strain, :math:`\varepsilon_{\rm i}`, under uniaxial compression is exponential and is given by

.. math::
   :label: cdpm2:sigma

   \begin{array}{ll}
   \sigma = f_{\rm c} \exp\left(-\dfrac{\varepsilon_{\rm i}}{\varepsilon_{\rm fc}}\right) & \mbox{if \( 0 < \varepsilon_{\rm i} \)}
   \end{array}

where :math:`\varepsilon_{\rm fc}` is a parameter controls the initial inclination of the softening curve. The inelastic strain is defined as the difference between the total strain, :math:`\varepsilon`, and the elastic strain, :math:`\sigma/E`. Replacing :math:`\sigma` by the right-hand side of (:eq:`cdpm2:sigma`), we obtain a nonlinear equation

.. math::
   :label: cdpm2:epsi

   \varepsilon_{\rm i} + \dfrac{f_{\rm c}}{E} \exp\left(-\dfrac{\varepsilon_{\rm i}}{\varepsilon_{\rm fc}}\right) = \varepsilon

from which :math:`\varepsilon_{\rm i}` can be evaluated iteratively if the total strain :math:`\varepsilon` is given. The use of different damage evolution rules for tension and compression is one important improvement over ConcreteDPM. In contrast to tension, the evolution law for compressive damage is considered to be independent of the element size, since compressive failure is often characterized by mesh-independent zones of damage growth.

.. table:: Summary of parameters of CDPM2
   :name: CDPM2-params
   :align: center

   +--------------------------+------------------+------+-----------------+----------------------------------------------------------------------------------------+
   | parameter                | OOFEM identifier | unit | default         | meaning                                                                                |
   +==========================+==================+======+=================+========================================================================================+
   | :math:`E`                | E                | Pa   |                 | Young's modulus                                                                        |
   +--------------------------+------------------+------+-----------------+----------------------------------------------------------------------------------------+
   | :math:`\nu`              | n                | -    |                 | Poisson's ratio                                                                        |
   +--------------------------+------------------+------+-----------------+----------------------------------------------------------------------------------------+
   | :math:`f_t`              | ft               | Pa   |                 | uniaxial tensile strength                                                              |
   +--------------------------+------------------+------+-----------------+----------------------------------------------------------------------------------------+
   | :math:`f_c`              | fc               | Pa   |                 | uniaxial compression strength                                                          |
   +--------------------------+------------------+------+-----------------+----------------------------------------------------------------------------------------+
   | :math:`w_f`              | wf               | m    |                 | critical crack opening                                                                 |
   +--------------------------+------------------+------+-----------------+----------------------------------------------------------------------------------------+
   | :math:`e`                | ecc              | -    | 0.525           | eccentricity, eq.~(:eq:`eq:rFunction`)--(:eq:`eq:frictionM`)                           |
   +--------------------------+------------------+------+-----------------+----------------------------------------------------------------------------------------+
   | :math:`q_{h0}`           | kinit            | -    | 0.3             | initial value of hardening variable :math:`q_{h1}`, eq.~(:eq:`eq:hardeningLawOne`)     |
   +--------------------------+------------------+------+-----------------+----------------------------------------------------------------------------------------+
   | :math:`H_p`              | hp               | -    | 0.5             | hardening modulus (for the last stage, eq.~(:eq:`eq:hardeningLawTwo`))                 |
   +--------------------------+------------------+------+-----------------+----------------------------------------------------------------------------------------+
   | :math:`D_f`              | dilation         | -    | 0.85            | dilation factor, eq.~(:eq:`cdpm2:betac`)                                               |
   +--------------------------+------------------+------+-----------------+----------------------------------------------------------------------------------------+
   | :math:`A_h`              | Ahard            | -    | 0.08            | hardening parameter, eq.~(:eq:`eq:xh`)                                                 |
   +--------------------------+------------------+------+-----------------+----------------------------------------------------------------------------------------+
   | :math:`B_h`              | Bhard            | -    | 0.003           | hardening parameter, eq.~(:eq:`eq:xh`)                                                 |
   +--------------------------+------------------+------+-----------------+----------------------------------------------------------------------------------------+
   | :math:`C_h`              | Chard            | -    | 2               | hardening parameter, eq.~(:eq:`eq:xh`)                                                 |
   +--------------------------+------------------+------+-----------------+----------------------------------------------------------------------------------------+
   | :math:`D_h`              | Dhard            | -    | :math:`10^{-6}` | hardening parameter, eq.~(:eq:`eq:xh`)                                                 |
   +--------------------------+------------------+------+-----------------+----------------------------------------------------------------------------------------+
   | :math:`A_s`              | Asoft            | -    | 15              | softening parameter, eq.~(:eq:`eq:xs`)                                                 |
   +--------------------------+------------------+------+-----------------+----------------------------------------------------------------------------------------+
   | :math:`\varepsilon_{fc}` | efc              | -    | :math:`10^{-4}` | softening parameter for compression, eq.~(:eq:`cdpm2:sigma`)                           |
   +--------------------------+------------------+------+-----------------+----------------------------------------------------------------------------------------+
   | stype                    |                  | -    | 1               | type of softening (1=bilinear)                                                         |
   +--------------------------+------------------+------+-----------------+----------------------------------------------------------------------------------------+
   | :math:`w_{f1}/w_f`       | wf1              | -    | 0.15            | softening parameter for tension                                                        |
   +--------------------------+------------------+------+-----------------+----------------------------------------------------------------------------------------+
   | :math:`\sigma_1/f_t`     | ft1              | -    | 0.3             | softening parameter for tension, eq.~(:eq:`cdpm2:omegat`)                              |
   +--------------------------+------------------+------+-----------------+----------------------------------------------------------------------------------------+

From the presented description of the model equations, it is clear that CDPM2 uses a large number of parameters. In [GraXenNys13]_ it was recommended to adjust only a few basic parameters, most of which have a certain physical meaning, and to set all the other parameters to their default values. The physical parameters that can be adjusted depending on the specific type of concrete are Young's modulus, Poisson's ratio, uniaxial tensile and compression strengths, and the critical crack opening, which controls the tensile fracture energy. A summary of all model parameters and of their recommended default values is provided in Table :numref:`CDPM2-params`. The table shows the symbol used for the parameter in equations describing the theoretical background and also the corresponding identifier used in OOFEM input files.

Parameters listed on the first five lines of Table :numref:`CDPM2-params` must be specified by the user. Parameter :math:`w_f` corresponds to the opening of a cohesive crack at which the stress transmitted by the crack vanishes. It is closely related to the tensile fracture energy, :math:`G_{Ft}`. The shape of the cohesive curve for tensile failure is affected by the choice of the type of softening. The default choice is bilinear softening, which requires, in addition to :math:`f_t` and :math:`w_f`, two more parameters that control the position of the point at which the slope of the diagram changes. This point is characterized by crack opening :math:`w_{f1}` and stress :math:`\sigma_1`, but the input parameters for OOFEM are dimensionless ratios :math:`w_{f1}/w_f` and :math:`\sigma_1/f_t` with default values 0.15 and 0.3. In general, the fracture energy is evaluated as :math:`G_{Ft}=(f_tw_{f1} +\sigma_1 w_f)/2`, and for the default values this gives :math:`G_{Ft}=f_tw_f/4.444`, which means that parameter :math:`w_f` can be evaluated from the tensile strength and fracture energy using formula :math:`w_f=4.444\,G_{Ft}/f_t`.

The eccentricity parameter :math:`e` could be evaluated from :math:`f_t` and :math:`f_c` combined with the biaxial compression strength, :math:`f_{bc}`, using formula

.. math::

   e=    \frac{1+\epsilon}{2-\epsilon}

where 

.. math::

   \epsilon=\frac{f_t}{f_{bc}}\frac{f_{bc}^2-f_c^2}{f_c^2-f_t^2}

The default value is not mentioned in  [GraXenNys13]_; in OOFEM it is taken as :math:`e=0.525`. 

Parameter :math:`H_p` is a dimensionless hardening modulus that sets the ratio between the increments of variables :math:`q_{h2}` and :math:`\kappa_p` during the last stage of the plastic hardening process. The value recommended in  [GraXenNys13]_ is :math:`H_p=0.01`, but the default value used in OOFEM is 0.5.

The dilation parameter :math:`D_f` controls the ratio between the transversal and axial plastic strain rates in uniaxial compression and thus can be used to control the volumetric expansion in the softening regime. The default value is set to 0.85, as suggested in [GraJir]_, which means that the transversal plastic strain rate will be :math:`-0.85` times the (negative) axial plastic strain rate. 

The default values of :math:`q_{h0}=0.3`, :math:`A_h=0.08`, :math:`B_h=0.003`, :math:`C_h=2`, :math:`D_h=10^{-6}` and :math:`\varepsilon_{fc}=10^{-4}` are taken from [GraXenNys13]_. Parameter :math:`A_s` was varied in [GraXenNys13]_ in the range from 1.5 to 15. In OOFEM, the default value of :math:`A_s` is taken as 15, which was recommended in [GraJir]_.

Concrete is strongly :math:`\bf{rate-dependent}`. If the loading rate is increased, the tensile and compressive strength increase and are more prominent in tension than in compression. The dependency of the stress-strain response on the dynamic increase factor is taken into account by an additional variable :math:`\alpha_{\rm r}`. The rate dependency is included by scaling both the equivalent strain rate and the inelastic strain. The rate parameter is defined by

.. math::

   \alpha_{\rm r} = (1-X) \alpha_{\rm rt} + X \alpha_{\rm rc}

where :math:`X` is the continuous compression measure (= 1 means only compression, = 0 means only tension).

The dynamic increase factor (DIF) for strength is computed according to Model Code 2010.

The :math:`\bf{parameters of CDPM2}` for OOFEM input file are summarised in :numref:`cdpm2_table`.
By default, two damage variables are used, i.e., the effective stress is transformed into nominal stress using equation (:eq:`cdpm2:efftonominal`). If parameter :em:`damflag` is set to 3, only :math:`\omega_t` is used and the effective stress is reduced by :math:`1-\omega_t`.

.. table:: CDPM2 -- summary
   :name: cdpm2_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | CDPM2                                                                                        |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`con2dpm` :elemparam:`d{rn}` :elemparam:`E{rn}` :elemparam:`n{rn}`                 |
   |                    | :elemparam:`tAlpha{rn}` :elemparam:`ft{rn}` :elemparam:`fc{rn}` :elemparam:`wf{rn}`          |
   |                    | [:elemparam:`stype{in}`] [:elemparam:`ft1{rn}`] [:elemparam:`wf1{rn}`]                       |
   |                    | [:elemparam:`efc{rn}`] [:elemparam:`ecc{rn}`] [:elemparam:`kinit{rn}`]                       |
   |                    | [:elemparam:`Ahard{rn}`] [:elemparam:`Bhard{rn}`] [:elemparam:`Chard{rn}`]                   |
   |                    | [:elemparam:`Dhard{rn}`] [:elemparam:`Asoft{rn}`] [:elemparam:`helem{rn}`]                   |
   |                    | [:elemparam:`dilation{rn}`] [:elemparam:`hp{rn}`] [:elemparam:`damflag{in}`]                 |
   |                    | [:elemparam:`sratetype{in}`] :elemparam:`eratetype{in}` [:elemparam:`yieldtol{rn}`]          |
   |                    | [:elemparam:`newtoniter{in}`]                                                                |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`d` material density                                                                |
   |                    | - :param:`E` Young modulus                                                                   |
   |                    | - :param:`n` Poisson ratio                                                                   |
   |                    | - :param:`tAlpha` thermal dilatation coefficient                                             |
   |                    | - :param:`ft` uniaxial tensile strength                                                      |
   |                    | - :param:`fc` uniaxial compressive strength                                                  |
   |                    | - :param:`wf` parameter that controls the slope of the softening branch to be used           |
   |                    | - :param:`stype` allows to choose different types of softening laws, default value 1: ; 0 -  |
   |                    |   linear softening ; 1 - bilinear softening ; 2 - exponential softening                      |
   |                    | - :param:`ft1` parameter for the bilinear softening law (stype = 1) defining the ratio       |
   |                    |   between intermediate stress and tensile strength, optional, default value 0.3              |
   |                    | - :param:`wf1` parameter for the bilinear softening law (stype = 1) defining the             |
   |                    |   intermediate crack opening, optional, default value 0.15                                   |
   |                    | - :param:`efc` parameter for exponential softening law in compression, optional, default     |
   |                    |   value :math:`100 \times 10^{-6}`                                                           |
   |                    | - :param:`ecc` eccentricity parameter, optional, default value 0.525                         |
   |                    | - :param:`kinit` initial value of hardening law, optional, default value 0.3                 |
   |                    | - :param:`Ahard` hardening parameter, optional, default value 0.08                           |
   |                    | - :param:`Bhard` hardening parameter, optional, default value 0.003                          |
   |                    | - :param:`Chard` hardening parameter, optional, default value 2                              |
   |                    | - :param:`Dhard` hardening parameter, optional, default value :math:`1\times 10^{-6}`        |
   |                    | - :param:`Asoft` softening parameter, optional, default value 15                             |
   |                    | - :param:`helem` element size :math:`h`, optional (if not specified, the actual element size |
   |                    |   is used)                                                                                   |
   |                    | - :param:`dilation` dilation factor, optional, default value 0.85                            |
   |                    | - :param:`hp` hardening modulus, optional, default value 0.5                                 |
   |                    | - :param:`damflag` flag which allows for the use of different damage formulations. Default   |
   |                    |   value is 1. If set to zero, no damage is used. If set to 3, only the tensile damage        |
   |                    |   varialbe is used.                                                                          |
   |                    | - :param:`sratetype` Dynamic increase factor (DIF) for strength ; 0 - DIF = 1 (no rate       |
   |                    |   effect) ; 1 - First branch of Model Code 2010 DIF for strength. ; 2 - First and second     |
   |                    |   branch of Model Code 2010 DIF for strength                                                 |
   |                    | - :param:`eratetype` Dynamic increase factor (DIF) for fracture energy. Only used if         |
   |                    |   sratetype is not equal to 0. ; 0 - DIF = 1 ; 1 - DIF from sratetype ; 2 - DIF squared from |
   |                    |   sratetype                                                                                  |
   |                    | - :param:`yieldtol` tolerance for the implicit stress return algorithm, optional, default    |
   |                    |   value :math:`1 \times 10^{-6}`                                                             |
   |                    | - :param:`newtoniter` maximum number of iterations in the implicit stress return algorithm,  |
   |                    |   optional, default value 100                                                                |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat, PlaneStrain, 1dMat                                                                    |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `sm/con2dpm1.in                                                                              |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/con2dpm1.in>`_,               |
   |                    | `sm/con2dpm10.in                                                                             |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/con2dpm10.in>`_,              |
   |                    | `sm/con2dpm2.in                                                                              |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/con2dpm2.in>`_,               |
   |                    | `sm/con2dpm3.in                                                                              |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/con2dpm3.in>`_,               |
   |                    | `sm/con2dpm4.in                                                                              |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/con2dpm4.in>`_,               |
   |                    | `sm/con2dpm5.in                                                                              |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/con2dpm5.in>`_,               |
   |                    | `sm/con2dpm6.in                                                                              |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/con2dpm6.in>`_,               |
   |                    | `sm/con2dpm7.in                                                                              |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/con2dpm7.in>`_,               |
   |                    | `sm/con2dpm8.in                                                                              |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/con2dpm8.in>`_,               |
   |                    | `sm/con2dpm9.in                                                                              |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/con2dpm9.in>`_ (and 2 more)   |
   +--------------------+----------------------------------------------------------------------------------------------+

.. _sec:cdpm2f:

CDPM2F
~~~~~~

CDPM2F is a fibre extension of CDPM2 in :ref:`sec:cdpm2`. The theory of CDPM2F is described in [ZhoMarGra24]_.

The main idea of CDPM2F is to adjust the damage part of CDPM2 to predict the effect of fibres.

Below in Table :numref:`cdpm2f_table` the input parameters of CDPM2F in addition to those required for CDPM2 are described.

.. table:: CDPM2F -- summary.
   :name: cdpm2f_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | CDPM2F                                                                                       |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`cdpm2f` :elemparam:`lf{rn}` :elemparam:`vf0{rn}` :elemparam:`df{rn}`              |
   |                    | :elemparam:`tau0{rn}` :elemparam:`beta{rn}` [:elemparam:`f{rn}`] :elemparam:`ef{rn}`         |
   |                    | :elemparam:`sm{rn}` [:elemparam:`alpha{rn}`] [:elemparam:`xi{rn}`]                           |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`lf` length of fibre                                                                |
   |                    | - :param:`vf0` base volume fraction of fibres                                                |
   |                    | - :param:`df` diameter of fibre                                                              |
   |                    | - :param:`tau0` initial bond strength of interface between fibre and matrix                  |
   |                    | - :param:`beta` hardening parameter of interface between fibre and matrix                    |
   |                    | - :param:`f` snubbing factor                                                                 |
   |                    | - :param:`ef` Young's modulus of fibre material                                              |
   |                    | - :param:`sm` staturated crack spacing                                                       |
   |                    | - :param:`alpha` parameter related to dispersion                                             |
   |                    | - :param:`xi` parameter controling hardening of fibre stress strain response                 |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat, PlaneStrain, 1DMat                                                                    |
   +--------------------+----------------------------------------------------------------------------------------------+

.. _concreteFCM:

Fixed crack model for concrete - ConcreteFCM
--------------------------------------------

Implementation of a fixed crack model. Uncracked material is modeled as isotropic linear elastic characterized by Young's modulus and Poisson's ratio. Cracking is initiated when principal stress reaches tensile strength. Further loading is governed by a softening law. Proper amount energy dissipation is guaranteed by the crack-band approach.

Multiple cracking is allowed; the maximum number of allowed cracks is controlled by :ncracks: parameter. Only mutually perpendicular cracks are supported. If cracking occurs in more directions, the behavior on the crack planes is considered to be independent.

The secant stiffness is used for unloading and reloading. In a compression regime, this model correspond to an isotropic linear elastic material.

The model supports 3 different options for stiffness reduction in shear after cracking, 2 options for the limit of the maximum shear stress on a crack plane and 6 different laws for postpeak behavior. The model parameters are summarized in :numref:`concrete_fcm_table`.

ConcreteFCM 1 d 24.e-3 talpha 12.e-6 E 20000.~n 0.2 Gf 100e-6 ft 2.0 softType 2 shearType 1 beta 0.05 shearStrengthType 2 fc 30 ag 0.01\\ lengthscale 1.~multipleCrackShear

.. table:: Fixed crack model for concrete -- summary
   :name: concrete_fcm_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Fixed crack model for concrete                                                               |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`ConcreteFCM` :elemparam:`in` :elemparam:`d{rn}` :elemparam:`tAlpha{rn}`           |
   |                    | :elemparam:`E{rn}` :elemparam:`n{rn}` :optelemparam:`ncracks{in}`                            |
   |                    | :optelemstring:`multipleCrackShear` :optelemparam:`crackSpacing{rn}`                         |
   |                    | :optelemparam:`softType{in}` :optelemparam:`shearType{in}`                                   |
   |                    | :optelemparam:`shearStrengthType{in}` :optelemparam:`ecsm{rn}` :optelemparam:`Gf{rn}`        |
   |                    | :optelemparam:`ft{rn}` :optelemparam:`beta{rn}` :optelemparam:`sf{rn}`                       |
   |                    | :optelemparam:`fc{rn}` :optelemparam:`ag{rn}` :optelemparam:`lengthscale{rn}`                |
   |                    | :optelemparam:`soft_w{ra}` :optelemparam:`soft(w){ra}` :optelemparam:`soft_eps{ra}`          |
   |                    | :optelemparam:`soft(eps){ra}` :optelemparam:`beta_w{ra}` :optelemparam:`beta(w){ra}`         |
   |                    | :optelemparam:`H{rn}` :optelemparam:`eps_f{rn}`                                              |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - material model number                                                                      |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`tAlpha` thermal dilatation coefficient                                             |
   |                    | - :param:`E` Young's modulus                                                                 |
   |                    | - :param:`n` Poisson's ratio                                                                 |
   |                    | - :param:`ncracks` maximum allowed number of cracks                                          |
   |                    | - :param:`crackSpacing` specified distance between parallel cracks                           |
   |                    | - :param:`multipleCrackShear` if not given, shear stiffness computed from the dominant       |
   |                    |   crack, otherwise all cracks contribute                                                     |
   |                    | - :param:`softType` allows to select suitable softening law: ; 0 - no softening (default) ;  |
   |                    |   1 - exponential softening with parameters :param:`Gf` and :param:`ft` ; 2 - linear         |
   |                    |   softening with parameters :param:`Gf` and :param:`ft` ; 3 - Hordijk softening with         |
   |                    |   parameters :param:`Gf` and :param:`ft` ; 4 - user-defined wrt crack opening with           |
   |                    |   parameters :param:`ft`, :param:`soft_w`, and :param:`soft(w)` ; 5 - linear hardening wrt   |
   |                    |   strain with parameters :param:`ft`, :param:`H`, and optionally :param:`eps_f` ; 6 - user-  |
   |                    |   defined wrt strain with parameters :param:`ft`, :param:`soft_eps`, and :param:`soft(eps)`  |
   |                    | - :param:`shearType` offers to choose from different approaches for shear stiffness          |
   |                    |   reduction of a cracked element ; 0 - no shear reduction (default) ; 1 - constant shear     |
   |                    |   retention factor with parameter :param:`beta` ; 2 - constant shear factor coefficient with |
   |                    |   parameter :param:`sf` ; 3 - user-defined shear retention factor with parameters            |
   |                    |   :param:`beta_w` and :param:`beta(w)`                                                       |
   |                    | - :param:`shearStrengthType` allows to select a shear stress limit on a crack plane ; 0 - no |
   |                    |   stress limit (default) ; 1 - constant strength = :math:`f_t` ; 2 - Collins interlock with  |
   |                    |   parameters :param:`fc`, :param:`ag`, and :param:`lengthscale`                              |
   |                    | - :param:`ecsm` method used for evaluation of characteristic element size :math:`L`: 1 =     |
   |                    |   square root of area, 2 = projection centered, 3 = Oliver, 4 = Oliver modified, 0 (default) |
   |                    |   = projection                                                                               |
   |                    | - :param:`Gf` fracture energy                                                                |
   |                    | - :param:`ft` tensile strength                                                               |
   |                    | - :param:`beta` shear retention factor                                                       |
   |                    | - :param:`sf` shear factor coefficient                                                       |
   |                    | - :param:`fc` compressive strength in MPa                                                    |
   |                    | - :param:`ag` aggregate size                                                                 |
   |                    | - :param:`lengthscale` factor to convert crack opening and aggregate size in case of Collins |
   |                    |   aggregate interlock; 1 = analysis in meters, 1000 = in millimeters, etc.                   |
   |                    | - :param:`soft_w` specified values of crack opening and                                      |
   |                    | - :param:`soft(w)` corresponding values of traction normalized to :param:`ft`                |
   |                    | - :param:`soft_eps` specified values of cracking strain and                                  |
   |                    | - :param:`soft(eps)` corresponding values of traction normalized to :param:`ft`              |
   |                    | - :param:`beta_w` specified values of crack opening and                                      |
   |                    | - :param:`beta(w)` corresponding values of shear retention factor                            |
   |                    | - :param:`H` hardening modulus (expressed wrt cracking strain)                               |
   |                    | - :param:`eps_f` threshold for cracking strain after which traction is zero (applicable for  |
   |                    |   linear hardening only)                                                                     |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat, PlaneStress, PlaneStrain                                                              |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `sm/concrete_fcm_shear.in                                                                    |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/concrete_fcm_shear.in>`_      |
   +--------------------+----------------------------------------------------------------------------------------------+

.. _concreteFCMViscoelastic:

ConcreteFCMViscoelastic
~~~~~~~~~~~~~~~~~~~~~~~

The material model “ConcreteFcmViscoelastic” has been developed in order to
realistically capture the time-dependent behavior of cementitious materials
with the possibility of tensile cracking.
The model allows to link an arbitrary
viscoelastic material model present in OOFEM with the model "ConcreteFCM" for
tensile failure with fixed (not rotating) orientation of cracks
developed for concrete. The material model assumes an additive split of the total deformation
into the viscoelastic and cracking part. The viscoelastic material handles the
stress-independent strains, i.e. shrinkage and thermal dilation.

Incremental stiffness matrix
----------------------------

In the case of an elastic (initial) stiffness or uncracked material, the
stiffness matrix :math:`\bar{\mbf{D}}_{ve}` is evaluated as the product of incremental stiffness of the
viscoelastic material :math:`\bar{E}` and the unit elastic
stiffness matrix :math:`\mbf{D}_v`.
The evaluation of the incremental modulus :math:`\bar{E}` is
computationally demanding, therefore within one time step it is stored
in the material status. The value of stiffness depends on the history
of relative humidity and temperature at given Gauss point and
depends entirely on the viscoelastic material.


Secant stiffness matrix is obtained by adding the secant compliance of
the cohesive cracks to the viscoelastic compliance matrix and by subsequent
inversion of the result. 

.. math::

   \mbf{D} =   \left( \bar{\mbf{D}}^{-1}_{ve} + \mbf{D}^{-1}_{cr}  \right)^{-1}

Tangent stiffness matrix is computed as 

.. math::

   \mbf{D} = \bar{\mbf{D}}_{ve} - \bar{\mbf{D}}_{ve} \left(\bar{\mbf{D}}_{ve} + \mbf{D}_{cr} \right)^{-1} \bar{\mbf{D}}_{ve}

where :math:`\mbf{D}_{cr}` is the tangent stiffness of the cracked material.

The secant and tangent matrix are constructed in the local
coordinate system corresponding to the orientation of the cracks and
then are rotated to the global
coordinate system. The individual components of :math:`\mbf{D}_{cr}` are managed by
the ``ConcreteFCM`` material model.

Algorithm for stress evaluation
-------------------------------

The stress evaluation algorithm can be decomposed into 3 subsequent steps:


- Call the viscoelastic material to compute the stress-dependent strain increment

.. math::

   \Delta \varepsilon_\sigma  = \Delta \varepsilon - \Delta \varepsilon'' - \Delta \varepsilon_{sh} - \Delta \varepsilon_T  

where :math:`\Delta \varepsilon` is the total strain increment, :math:`\Delta \varepsilon''` is the strain increment due to creep, :math:`\Delta \varepsilon_{sh}` is the shrinkage increment, and :math:`\Delta \varepsilon_T` is the thermal dilation increment.

- Call the fixed-crack material to compute the new stress vector :math:`\sigma_{i+1}` fulfilling the nonlinear equation

.. math::

   \sigma_{i+1} = \sigma_{i} + \bar{\mbf{D}}_{ve} \left ( \Delta \varepsilon_\sigma - \Delta \varepsilon_{cr} \right)

where :math:`\sigma_{i}` is the old stress vector and :math:`\Delta \varepsilon_{cr}` is the increment of crack strain, and update the new crack strain vector

.. math::

   \varepsilon_{cr,i+1} = \varepsilon_{cr,i} + \Delta \varepsilon_{cr}  

- Pass the difference of the total strain and the cracking strain to the viscoelastic material and perform a regular viscoelastic step.

.. math::

   \sigma_{i+1} = \sigma_{i} + \bar{\mbf{D}}_{ve} \left ( \Delta \varepsilon - \Delta \varepsilon_{cr}  - \Delta \varepsilon'' - \Delta \varepsilon_{sh} - \Delta \varepsilon_T  \right)

Before the onset of cracking, the material is modeled as isotropic linear viscoelastic characterized by its compliance function, Poisson's ratio, and the evolution of shrinkage which can depend both on equivalent time and change in relative humidity provided that the problem is solved as coupled.

Cracking is initiated when principal stress reaches tensile strength. Further loading is governed by a softening law. Proper amount energy dissipation is guaranteed by the crack-band approach, the width of the crack band is given by the size of the finite element projected in the direction of the principal stress.

Multiple cracking is allowed, the maximum number of cracks is controlled by :param:`ncracks` parameter. Only mutually perpendicular cracks are supported. If cracking occurs in more directions, the behavior on the crack planes is considered to be independent.
The secant stiffness is used for unloading and reloading. In compression regime full stiffness is restored.

It is worth emphasizing that the status of the viscoelastic material
tracks and keeps only a part of the total deformation equal to the
difference of the total strain and the cracking strain.

.. _evolution-of-tensile-strength-and-fracture-energy:

Evolution of tensile strength and fracture energy
-------------------------------------------------

The current implementation supports both constant and time-evolving
values of tensile strength and fracture energy.
**Constant** values are used when the keyword :param:`timeDepFracturing`
is missing in the definition of the material. Then, the values are taken directly from the
:emph:`ConcreteFCM` material (keywords :param:`ft` and :param:`gf`).
On the other hand, when the flag :samp:`timeDepFracturing` is present
in the material definition, the tensile strength and fracture energy
is **time-evolving**. In that case, the values of :param:`ft` and :param:`gf`
are ignored, and the value of parameter :math:`s` (:param:`fib_s`), and the mean
value of the compressive strength :math:`f_{cm,28}` (:param:`fcm28`) need to be
provided and serve for the prediction of :math:`f_t(t)` and :math:`G_F(t)`.
Yet, it is possible to
override the predicted 28-day values by providing the fields
:param:`ft28` and :param:`gf28`.

It must be noted that once the tensile stress exceeds the current
value of tensile strength, the time development of the properties is
terminated on the particular crack plane.

Mean value of the uniaxial compressive strength at specific time which
is used for evaluation of the tensile strength and fracture energy is
computed using expressions (5.1-50) and (5.1-51) in :emph:`fib` MC2010 [fib:2010]_:

.. math::

   f_{cm}(t_e) = \exp \left[ s \left( 1. - \sqrt{28/t_e} \right)   \right] f_{cm,28}

in which :math:`s` is a coefficient reflecting cement type (table
5.1-9), :math:`t_e` is the equivalent time from the viscoelastic model, and
:math:`f_{cm,28}` is the mean value of the uniaxial compressive strength at
the age of 28 days.

The current value of tensile strength :math:`f_{ctm}(t)` at time :math:`t` is computed from the compressive strength at the corresponding age, similarly to section 5.1.5.1 in MC2010. The difference is that in MC2010, both tensile strength and fracture energy are treated as non-aging and are evaluated from :math:`f_{cm,28}`.

In the following expressions, :math:`f_{cm}` and :math:`f_{ctm}` are in MPa.

.. math::

   f_{ctm}(t_e) = 2.12 \: \ln \left(1 + 0.1 f_{cm}(t_e) \right) \: \dots \:
   \mathrm{for \:} f_{cm} \geq 58\:\mathrm{MPa}

.. math::

   f_{ctm}(t_e) = 0.3 \left(f_{cm}(t_e)-8 \right)^{2/3} \: \dots \:
   \mathrm{for \:} 20 \leq f_{cm} \leq 58\:\mathrm{MPa}

.. math::

   f_{ctm}(t_e) = 0.07862 \: f_{cm} \: \dots \:
   \mathrm{for \:} f_{cm} \leq 20 \: \mathrm{MPa}

The value of fracture energy at the age of 28 days is computed according to section 5.1.5.2 in MC2010

.. math::

   G_F = 73 \left( f_{cm,28} \right)^{0.18}

In this equation, :math:`G_F` is in N/m and the compressive strength is in MPa.

Here, the time development of fracture energy is chosen to be the same as the tensile strength

.. math::

   G_F(t) = G_F \frac{f_{ctm}(t)}{f_{ctm}(28)}

.. table:: Viscoelastic fixed crack model for concrete -- summary.
   :name: concretefcmviscoelastic_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Viscoelastic fixed crack model for concrete                                                  |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`ConcreteFCMViscoelastic` :elemparam:`viscoMat{in}` :elemparam:`gf{rn}`            |
   |                    | :elemparam:`ft{rn}` :optelemstring:`timeDepFracturing` :optelemparam:`fib_s{rn}`             |
   |                    | :optelemparam:`fcm28{rn}` :optelemparam:`timefactor{rn}` :optelemparam:`stiffnessfactor{rn}` |
   |                    | :optelemparam:`gf28{rn}` :optelemparam:`ft28{rn}`                                            |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - material model number                                                                      |
   |                    | - :param:`viscoMat` number of the slave viscoelastic material                                |
   |                    | - :param:`ft` tensile strength - constant value                                              |
   |                    | - :param:`gf` fracture energy - constant value                                               |
   |                    | - :param:`timeDepFracturing` flag to activate evolving :math:`f_t` and :math:`G_f` according |
   |                    |   to *fib MC2010*                                                                            |
   |                    | - :param:`fib_s` strength class of cement coefficient (Table 5.1-9 in fib MC2010)            |
   |                    | - :param:`fcm28` 28-day mean compressive strength of concrete used for prediction and        |
   |                    |   evolution of :math:`f_t` and :math:`G_f`                                                   |
   |                    | - :param:`timefactor` time transformation factor, 1 for days, 86400 for seconds etc.         |
   |                    | - :param:`stiffnessfactor` stiffness transformation factor, 1 for Pa, 1.e6 for MPa etc.      |
   |                    | - :param:`ft28` overrides predicted tensile strength                                         |
   |                    | - :param:`gf28` overrides predicted fracture energy                                          |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat, PlaneStress, PlaneStrain                                                              |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `sm/concrete_fcm_visco.in                                                                    |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/concrete_fcm_visco.in>`_      |
   +--------------------+----------------------------------------------------------------------------------------------+

.. subsubsection*:: Examples

   The following lines show a sample syntax specifying the viscoelastic fixed crack model with the compliance defined by the MPS model. Volume density is 24 kN/m\ :sup:`3`, thermal dilation coefficient :math:`12\times10^{-6}` K\ :sup:`-1`, Poisson's ratio~0.2, fracture energy :math:`G_f` and tensile strength :math:`f_t` predicted from the mean compressive strength :math:`f_{cm} = 30` MPa (defined :samp:`ft` and :samp:`gf` are discarded, to override use keywords :samp:`ft28` and :samp:`gf28`), exponential softening, maximum 2 perpendicular tensile cracks, constant shear factor coefficient, :samp:`sf = 20.`. All cracks contribute to the shear stiffness, shear strength equal to :math:`f_t`, numerical stiffness factors 1e-3, number of viscoelastic material is 2. The analysis uses [m], [MPa], [MN], and [day]:

   .. code-block:: text

      ConcreteFcmViscoelastic 1 d 24.e-3 talpha 0. E 30000. n 0.2
      Gf 100.e-6 ft 2. softtype 1 sheartype 2 sf 20.
      multiplecrackshear shearstrengthtype 1 ncracks 2
      shearCoeffNumer 0.001 normalCoeffNumer 0.001 viscomat 2
      timeDepFracturing fib_s 0.20 fcm28 30. timeFactor 1.
      stiffnessFactor 1.e6

The viscoelastic material can be defined, for example, by the following listing. Note that the coefficient of thermal expansion is zero in `ConcreteFcmViscoelastic` and nonzero in the viscoelastic material that handles all stress-independent strains. Also note that the same values of the Poisson's ratio are used.

mps 2 d 24.e-3 n 0.2 talpha 12.e-6 referencetemperature  296.
mode 0  fc 30. cc 340. w/c 0.523 a/c 5.28 stiffnessfactor 1.e6
timefactor 1. lambda0 1. begoftimeofinterest 1.e-6
endoftimeofinterest 1.e4 relMatAge 7. CoupledAnalysisType 2
ksh 0.0 k3 35. p 1000.

.. _frcfcm:


Fixed crack model for fiber reinforced composites - FRCFCM
----------------------------------------------------------


This material model is an extension of the `ConcreteFCM` described in :ref:`concreteFCM`.

It is possible to choose from three different “classes” of fibers. This choice is controlled by the keyword :param:`fiberType`: 0 = continuous aligned fibers (CAF), 1 = short aligned fibers (SAF), and 2 = short random fibers (SRF). Currently, it is not possible to combine more classes of fibers in one material model.

All of the above-mentioned fiber types are further defined by the material properties and geometry. Fiber quantity is captured by the dimensionless volume fraction :param:`Vf` (as a decimal). All fibers are assumed to have a circular cross-section (shape factor for shear :param:`kfib` with a default value of 0.9) and to possess the same geometry characterized by the diameter :param:`Df` and length :param:`Lf`.

The overall elastic stiffness of the fiber-reinforced composite is calculated as a weighted average of the moduli of the matrix :math:`E_m` and the fibers :math:`E_f`, and the Poisson's ratio is considered to be equal to the Poisson's ratio of the matrix. Similarly to `ConcreteFCM`, cracking is initiated once the tensile stress :math:`\sigma` in the matrix reaches the tensile strength :math:`f_t`.

The nominal bridging stress for continuous aligned fibers which are not perpendicular to the crack plane :math:`\sigma_{b,f,\theta}` can be very easily obtained by multiplying the nominal bridging stress for perpendicular fibers :math:`\sigma_{b,f}` by two terms: the first one reflecting the lower volume of inclined fibers passing through the crack plane and the second one capturing the snubbing effect:

.. math::

   \sigma_{b,f,\theta} = \sigma_{b,f} \cos(\theta) \exp(\theta f)

where :math:`f` is a snubbing coefficient.

For the **CAF** perpendicular to the crack, the nominal bridging stress can be derived as

.. math::

   \sigma_{b,f} = 2 V_f \sqrt{ \frac{ E_f (1+\eta) \tau_0 } {D_f} \: \bar{w} }

where :math:`\eta = (E_f V_f)/[E_m (1 - V_f)]`.

For **SRF** the nominal bridging stress is 

.. math::

   \sigma_{b,f}(w) &=& 2 V_f \sqrt{ \frac{E_f (1+\eta) \tau_0 \bar{w}}{D_f} } - \frac{V_f E_f (1+\eta) \bar{w} }{L_f}  \quad \mathrm{for} \: \bar{w} < w^* \\
   \sigma_{b,f}(w) &=& \frac{V_f L_f \tau_s(w)}{D_f} \left( 1- \frac{2 \bar{w}}{L_f} \right)^2 \quad \mathrm{for} \: w^* \leq \bar{w} < L_f/2 \\
   \sigma_{b,f}(w) &=& 0 \quad \mathrm{for} \: \bar{w} > L_f/2 

where :math:`w^* = \left(L_f^2 \tau_0 \right)/[(1+\eta) E_f D_f]`; :math:`\tau_0` is the bond shear strength between the fiber and matrix for small crack openings, :math:`w<w^\ast`.

Larger pull-out displacements can lead to significant physical changes in the fiber surface which can result into changes in the bond shear stress. This phenomenon is captured by function :math:`\tau_s(w)` relating the frictional bond to the crack opening and is implemented in three alternative formulations. (In order to keep :math:`\tau_s(w) = \tau_0` use :param:`fssType` = 0.) In conventional FRC with ordinary concrete matrix, the frictional bond usually decreases with increasing slip. To capture this type of behavior we adopt the function proposed by Sajdlová (activated with :param:`fssType` = 1) reads

.. math::

   \tau_s(w) = \tau_0 \left[ 1 + \mathrm{sign}(b_0) \left( 1 - \exp \left( -\frac{|b_0| w}{D_f}  \right) \right)  \right]

where :math:`b_0` is a micromechanical parameter.

In composites with high-strength matrix and coated high-strength steel fibers (HSFRC, UHPFRC) as well as in SHCC materials with polymeric fibers, the frictional bond-slip relation often exhibits hardening; this phenomenon can be well approximated by a cubic function (activated with :param:`fssType` = 2) proposed by Kabele

.. math::

   \tau_s(w) = \tau_0 \left[ 1 + b_1 \frac{w}{D_f} + b_2 \left( \frac{w}{D_f} \right)^2 + b_3 \left( \frac{w}{D_f} \right)^3 \right]

or an alternative formulation which results in smooth changes in the bridging stress (activated with :param:`fssType` = 3)

.. math::

   \tau_s(w) = \tilde{\tau}_0 + \tau_0 \left[ b_1 \frac{ \tilde{w} }{D_f} + b_2 \left( \frac{ \tilde{w} }{D_f} \right)^2 + b_3 \left( \frac{ \tilde{w} }{D_f} \right)^3 \right]

In the last two equations :math:`b_1`, :math:`b_2` and :math:`b_3` are micromechanical parameters and additionally in the last equation :math:`\tilde{w} = w - w^{\ast}`, :math:`\tilde{\tau_0} = \tau_0 (1.-w^{\ast}/L_f)^{-2}` for SRF and :math:`\tilde{\tau_0} = \tau_0 E_f (1+\eta) D_f / [ E_f (1+\eta) D_f - 2 L_f \tau_0 ]` for SAF.

The bridging stress for **SRF** is defined as

.. math::
   :label: eq:bridging-stress-srf

   \sigma_{b,f}(w) &= \frac{g V_f L_f \tau_0}{2 D_f} \left( 2 \sqrt{ \frac{\bar{w} }{w^{\ast}} } - \frac{\bar{w}}{ w^{\ast} } \right) \quad \text{for} \: \bar{w} < w^* \\
   \sigma_{b,f}(w) &= \frac{g V_f L_f \tau_s(w)}{2 D_f} \left( 1- \frac{2 w}{L_f} \right)^2 \quad \text{for} \: w^* \leq \bar{w} < L_f/2 \\
   \sigma_{b,f}(w) &= 0 \quad \text{for} \: \bar{w} > L_f/2 

where :math:`g` is the snubbing factor defined as

.. math::
   :label: eq:snubbing-factor

   g = 2 \frac{ 1 + \exp(\pi f / 2) }{ 4 + f^2}

It is possible to delay the activation of the stress in fibers using parameter :param:`fibreActivationOpening`, which can be imagined as a “lag” of the fiber-related crack opening behind the matrix-related crack opening.

During unloading, the stress in fibers does not decrease linearly to origin. The current implementation uses a power function

.. math::
   :label: eq:power-function

   \sigma_{b,f}(w) = \sigma_{b,f}(w_{max}) \left(\frac{\bar{w}}{\bar{w}_{max}}\right)^M

where :math:`w_{\text{max}}` is the maximum crack width reached in the entire previous history and :math:`M` is a positive constant, its default value is :math:`M = 4`.

The influence of crack opening and sliding on the bridging shear stress only due to fibers is expressed as

.. math::
   :label: eq:bridging-shear-stress

   \tau_{b,f} = \bar{V}_f k G_{f} \frac{u}{w_{max}} = \frac{\bar{V}_f k G_{f}}{\varepsilon_{cr,max}} \gamma_{cr} 

where :math:`\bar{V}_f` is the effective volume of fibers crossing a crack plane (:math:`V_f/2` for SRF and :math:`V_f \cos (\theta)` for CAF and SAF), and :math:`G_{f}` is the fiber shear modulus. This expression is motivated by the assumption that the fibers bridging the crack planes behave as Timoshenko beams subjected to shear. Note that the shear stiffness of fibers is not recovered upon unloading.

It has been found that in some high-performance fiber-reinforced cement composites, fibers rupture when cracks are exposed to shearing. This phenomenon is modeled by the damage parameter :math:`\omega`, which accounts for the ratio of ruptured fibers and varies between the values of 0 and 1. It is assumed that :math:`\omega` depends on the maximum shear strain sustained by the protruding portions of bridging fibers throughout the loading history. This crack shear strain can be expressed as:

.. math::
   :label: eq:crack-shear-strain

   \gamma_{f,max} = \max \left( \frac{\left| u_i(t) \right|}{ \max \left(w_i(t) \right)} \right) \quad \dots w(t) > \Delta w

where :math:`u_i` is the crack sliding displacement (CSD) and :math:`w_i` is the maximum value of the crack opening displacement of the i-th crack. This means that the damage does not grow if the crack closes (crack opening decreases). If more cracks exist, the maximum contribution is considered.

Two different one-parameter damage evolution laws are currently implemented. For :param:`fDamType` = 0 the damage is deactivated, with :param:`fDamType` = 1 damage is described by

.. math::

   \omega(\gamma_f) = \min \left( \frac{\gamma_f}{\gamma_{fc}},\: 1 \right)

and finally with :param:`fDamType` = 2

.. math::

   \omega(\gamma_f) = 1 - \exp \left( - \frac{\gamma_f}{\gamma_{fc}} \right)

where :math:`\gamma_{fc}` (:param:`gammaCrack` in the input record) is a parameter.

Since damage reduces the number of crack-bridging fibers, which is proportional to the fiber volume fraction, its effect can be suitably implemented by introducing the effective volume fraction

.. math::


V_f^{\ast} = V_f (1-\omega)

The material parameters are summarized in Tables :numref:`concrete_fcm_table` (matrix) and :numref:`frcfcm_table` (fiber extension).

Sample syntax for a fixed crack model reinforced with fibers with volume density 24 kN/m\ :sup:`3`, thermal dilation coefficient 12\ :sup:`-6` K\ :sup:`-1`, Young's modulus of the matrix 20 GPa, Poisson's ratio of matrix 0.2, fracture energy of matrix 100 N/m, tensile strength of matrix 2 MPa, linear tension softening, constant shear retention factor :math:`\beta = 0.05`, unlimited shear strength (:texttt{shearStrengthType = 0}), continuous aligned fibers, fiber volume 2\%, fiber diameter 0.04 mm, Young's modulus of fibers 20 GPa, shear modulus of fibers 1 GPa, fiber-matrix bond strength 1 MPa, snubbing coefficient 0.7, shear correction coefficient 0.9, deactivated fiber damage, fiber act if COD exceeds 10 :math:`\mu`\ m (with smoothing from :math:`w = 8` to :math:`11` :math:`\mu`\ m), fiber orientation at 45 degrees in x-y plane, automatic evaluation of crack spacing from composition; the analysis uses [m], [MPa] and [MN]:

.. code-block:: text

   FRCFCM 1 d 24.e-3 talpha 12.e-6 E 20000. n 0.2 Gf 100e-6 ft 2.0
   softType 2 shearType 1 beta 0.05 FiberType 0 Vf 0.02 Df 0.04e-3
   Ef 20000. Gfib 1000. tau_0 1. FSStype 0 f 0.7 kfib 0.9 fDamType 0
   fibreactivationopening 10.e-6 dw0 2.e-6 dw1 1.e-6 orientationVector 3 1. 1. 0. computeCrackSpacing

.. table:: Fixed crack model for fiber reinforced concrete -- summary.
   :name: frcfcm_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Fixed crack model for FRC                                                                    |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`FRCFCM` input record of ConcreteFCM :elemparam:`Vf{rn}` :elemparam:`Lf{rn}`       |
   |                    | :elemparam:`Df{rn}` :elemparam:`Ef{rn}` :optelemparam:`nuf{rn}` :optelemparam:`Gfib{rn}`     |
   |                    | :optelemparam:`kfib{rn}` :elemparam:`tau_0{rn}` :elemparam:`b0{rn}` :elemparam:`b1{rn}`      |
   |                    | :elemparam:`b2{rn}` :elemparam:`b3{rn}` :elemparam:`f{rn}` :optelemparam:`M{in}`             |
   |                    | :optelemparam:`fibreOrientationVector{ra}` :optelemparam:`fssType{in}`                       |
   |                    | :optelemparam:`fDamType{in}` :optelemparam:`fiberType{in}` :optelemparam:`gammaCrack{rn}`    |
   |                    | :optelemstring:`computeCrackSpacing` :optelemparam:`fibreActivationOpening{rn}`              |
   |                    | :optelemparam:`dw0{rn}` :optelemparam:`dw1{rn}`                                              |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`Vf` fiber content expressed as decimal                                             |
   |                    | - :param:`Lf` fiber length                                                                   |
   |                    | - :param:`Df` fiber diameter                                                                 |
   |                    | - :param:`Ef` fiber Young's modulus                                                          |
   |                    | - :param:`nuf` fiber Poisson's ratio                                                         |
   |                    | - :param:`Gfib` fiber shear modulus (read when :param:`nuf` is not provided)                 |
   |                    | - :param:`kfib` fiber cross-sectional shape correction factor                                |
   |                    | - :param:`tau_0` bond shear strength at zero slip                                            |
   |                    | - :param:`b0` micromechanical parameter for fiber shear according to Sajdlová                |
   |                    | - :param:`b1`, :param:`b2`, :param:`b3` micromechanical parameter for fiber shear according  |
   |                    |   to Kabele                                                                                  |
   |                    | - :param:`f` snubbing friction coefficient                                                   |
   |                    | - :param:`M` exponent related to fiber unloading                                             |
   |                    | - :param:`fibreOrientationVector` vector specifying orientation for CAF and SAF fibers       |
   |                    | - :param:`fssType` type of Fiber bond Shear Strength (bond shear strength vs. crack opening) |
   |                    |   ; 0 - constant shear strength ; 1 - bond shear strength with parameter :param:`b0` ; 2 -   |
   |                    |   bond shear strength with parameters :param:`b1`, :param:`b2`, :param:`b3` ; 3 - bond shear |
   |                    |   strength with parameters :param:`b1`, :param:`b2`, :param:`b3` which leads to smooth       |
   |                    |   traction-separation law                                                                    |
   |                    | - :param:`fDamType` type of damage law for fibers ; 0 - no damage ; 1 - damage controlled by |
   |                    |   shear slip deformation of the crack (with :param:`gammaCrack`), linear law ; 2 - damage    |
   |                    |   controlled by shear slip deformation of the crack (with :param:`gammaCrack`), exponential  |
   |                    |   law                                                                                        |
   |                    | - :param:`fiberType` type of reinforcing fibers ; 0 - CAF (continuous aligned fibers) ; 1 -  |
   |                    |   SAF (short aligned fibers) ; 2 - SRF (short randomly oriented fibers)                      |
   |                    | - :param:`gammaCrack` crack shear strain parameter applicable with fDamType = 1 or 2 (here   |
   |                    |   the crack shear strain is understood as the crack slip :math:`u` divided by the crack      |
   |                    |   opening :math:`w`)                                                                         |
   |                    | - :param:`computeCrackSpacing` crack spacing is evaluated automatically based on provided    |
   |                    |   composition                                                                                |
   |                    | - :param:`fibreActivationOpening` crack opening at which the fibers begin transferring       |
   |                    |   bridging stress                                                                            |
   |                    | - :param:`dw0`, :param:`dw1` applicable only if :param:`fibreActivationOpening` :math:`\neq  |
   |                    |   0`, then it allows to smooth the traction-separation law for fibers; lower bound is        |
   |                    |   :param:`fibreActivationOpening` - :param:`dw0` and the upper bound is                      |
   |                    |   :param:`fibreActivationOpening` - :param:`dw1`                                             |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat, PlaneStress, PlaneStrain                                                              |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `sm/frcfcm_shear.in                                                                          |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/frcfcm_shear.in>`_,           |
   |                    | `sm/frcfcm_tension.in                                                                        |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/frcfcm_tension.in>`_          |
   +--------------------+----------------------------------------------------------------------------------------------+

.. _frcfcmnl_table:

Nonlocal fixed crack model for fiber reinforced concrete -- summary.

.. list-table:: Nonlocal fixed crack model for FRC
   :header-rows: 1
   :widths: 20 80

   * - Description

     - Nonlocal fixed crack model for FRC
   * - Record Format

     - :descitem:`FRCFCMNL` input record of ConcreteFCM and FRCFCM
       :elemparam:`r` :rn
       :elemparam:`wft` :in

   * - Parameters

     - - :param:`r` nonlocal radius (reasonable value is several millimeters and its maximum is :math:`L_f/2` for short fibers)

       - :param:`wft` nonlocal averaging function, must be set to 4 (constant function)

   * - Supported modes

     - PlaneStress
