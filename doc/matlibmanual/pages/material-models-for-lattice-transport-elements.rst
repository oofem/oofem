Material models for lattice transport elements
==============================================


.. _latticetransmat:

Material for unsaturated flow in lattice models -- LatticeTransMat
------------------------------------------------------------------

Lattice transport elements have to be used with this material model. A positive sign is assumed for liquid tension, unlike the convention of soil mechanics which assumes compression positive.

These transport elements are idealised as one-dimensional conductive pipes. The gradient of hydraulic head, which governs flow rate along each transport element, is determined from the capillary pressures :math:`P_{c}` at the two nodes.

The mass balance equation describes the change in moisture inside a porous element as a consequence of liquid flow and solid-liquid retention. It leads to the following partial differential equation

.. math::
   :label: latticetransmat1

   c\frac{\partial P_{c}}{\partial t} + k\mbox{div}(\nabla(\frac{P_{c}}{g}-\rho z)) =0

where :math:`P_{c}` is the capillary pressure, :math:`c` is the mass capacity function(:math:`s^{2}m^{-2}`), :math:`k` is the Darcy hydraulic conductivity (:math:`ms^{-1}`), :math:`\rho` is the fluid mass density, :math:`g` is the acceleration of gravity, :math:`z` is the capillary height and :math:`t` is the time.

The hydraulic conductivity :math:`k` consists of

.. math::
   :label: latticetransmat11

   k=k_{0} + k_{c}

where :math:`k_{0}` is the hydraulic conductivity of the intact material and :math:`k_{c}` is the additional conductivity due to cracking.

Darcy hydraulic conductivity :math:`k_{0}` is defined as

.. math::
   :label: latticetransmat12

   k_{0}=\frac{\rho g}{\mu}\kappa

where :math:`\mu` is the dynamic viscosity (:math:`Pa.s`), :math:`\kappa` is the permeability also called intrinsic conductivity(:math:`m^{2}`), and :math:`\kappa_{r}` is the relative conductivity. :math:`\kappa_{r}` is a function of the effective degree of saturation.

The cracking part is

.. math::
   :label: latticetransmat13

   k_{c}=\xi\frac{\rho g}{\mu}\frac{w^{3}_{c}}{12h} \kappa_{r}

where :math:`\xi` is a tortuosity factor taking into account the roughness of the crack surface, :math:`w_{c}` is the equivalent crack opening of the dual mechanical lattice and :math:`h` is the length of the dual mechanical element.

.. _one-dimensional-transport-element:

One-dimensional transport element
---------------------------------

The discrete form of the differential equation for mass transport for a one-dimensional transport element is

.. math::
   :label: latticetransmat3

   \alpha_{e}P_{c} + C_{e}\frac{\partial P_{c}}{\partial t}  = f_{e}

where :math:`P_{c}` is a vector containing the nodal values of the capillary pressure, :math:`\alpha_{e}` is the conductivity matrix, :math:`C_{e}` is the capacity matrix and :math:`f_{e}` is the nodal flow rate vector (:math:`kg.m^{3}`). 

The capacity matrix is

.. math::
   :label: latticetransmat4

   C_{e} = \frac{Al}{12}c \left(  \begin{array}{ c c }  2 & 1 \\  1 & 2 \end{array}  \right)

where :math:`c` is the capacity of the material, :math:`l` is the length of the transport element and :math:`A` is the cross-sectional area of the transport element. 

The conductivity matrix is defined as

.. math::
   :label: latticetransmat5

   \alpha_{e} = \frac{A}{l}\frac{k}{g} \left(  \begin{array}{ c c }  1 & -1 \\  -1 & 1 \end{array}  \right)

The mass transport equation is based on the constitutive laws for the capacity :math:`c` and the hydraulic conductivity :math:`k`.

The capacity :math:`c` is defined as 

.. math::
   :label: latticetransmat6

   c=-\rho\frac{\partial \theta}{\partial P_{c}}

where :math:`\theta` is the volumetric water content (:math:`\theta=\frac{V_{w}}{V_{T}}` with :math:`V_{w}` the volume of water, and :math:`V_{T}` the total volume) which is calculated by a modified version of van Genuchten’s retention model. Note that the presence of a crack in an element does not influence the capacity in the present model.

The volumetric water content is 

.. math::
   :label: latticetransmat7

   \theta = S_{e} (\theta_{s} - \theta_{r})+\theta_{r}

where :math:`\theta_{r}` and :math:`\theta_{s}` are the residual and saturated water contents corresponding to effective saturation values of :math:`S_{e} = 0` and :math:`S_{e} = 1`, respectively.

The effective degree of saturation :math:`S_{e}` is defined as

.. math::
   :label: latticetransmat8

   S_{e}  = \left\{
       \begin{array}{ll}
            \frac{\theta_{m}-\theta_{r}}{\theta_{s}-\theta_{r}}\left(1+\left(\frac{ P_{c}}{a}\right)^{\frac{1}{1-m}}\right)^{-m}  & \mbox{if } P_{c} \ge P_{c(aev)} \\
            1 & \mbox{if } P_{c}<P_{c(aev)}
       \end{array}
   \right.

where :math:`\theta_{m}` is an additional model parameter and :math:`P_{c(aev)}` is the air-entry value of capillary pressure which separates saturated (:math:`P_{c}<P_{c(aev)}`) from unsaturated states (:math:`P_{c} \ge P_{c(aev)}`). It is intuitive that the smaller the pore size of the material, the larger the value of :math:`P_{c(aev)}` will be.

The relative conductivity :math:`\kappa_{r}` is a function of the effective degree of saturation and is defined as

.. math::
   :label: latticetransmat9

   \kappa_{r}  = \sqrt{S_{e}}\left(\frac{1-\left[1-\left(\frac{S_{e}}{  \frac{\theta_{m}-\theta_{r}}{\theta_{s}-\theta_{r}}}\right)^{\frac{1}{m}}\right]^{m}}{1-\left[1-\left(\frac{1}{  \frac{\theta_{m}-\theta_{r}}{\theta_{s}-\theta_{r}}}\right)^{\frac{1}{m}}\right]^{m}}\right)^{2}

If :math:`\theta_{m} = \theta_{s}`, we have :math:` \frac{\theta_{m}-\theta_{r}}{\theta_{s}-\theta_{r}}=1` : the equation reduces to the expression of the relative conductivity of the original van Genuchten model.

The model parameters are summarized

.. table:: Model Parameters
   :widths: auto

   +---------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | Description   | Material for fluid transport in lattice models                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
   +---------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | Record Format |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
   +---------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | Parameters    | - :param:`num` material model number - :param:`d` fluid mass density - :param:`k` permeability (:math:`m^{2}`) - :param:`vis` dynamic viscosity (:math:`Pa.s`) - :param:`contype` unsaturated flow allowed when :param:`contype`\ =1 - :param:`thetas` saturated water content - :param:`thetar` residual water content - :param:`paev` air-entry value of capillary pressure - :param:`m` van Genuchten parameter - :param:`a` van Genuchten parameter - :param:`thetam` additional model parameter for the modified version of van Genuchten’s retention model - :param:`ctor` coefficient of tortuosity (:math:`ctor=\frac{1}{\tau}\le1`) |
   +---------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

.. table:: Material for unsaturated flow in lattice models - summary.
   :name: Iatticetransmat_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Material for fluid transport in lattice models                                               |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`latticetransmat` :elemparam:`num{in}` :elemparam:`d{rn}` :elemparam:`k{rn}`       |
   |                    | :elemparam:`vis{rn}` :elemparam:`contype{in}` :elemparam:`thetas{rn}`                        |
   |                    | :elemparam:`thetar{rn}` :elemparam:`paev{rn}` :elemparam:`m{rn}` :elemparam:`a{rn}`          |
   |                    | :elemparam:`thetam{rn}` :optelemparam:`ctor{rn}`                                             |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`num` material model number                                                         |
   |                    | - :param:`d` fluid mass density                                                              |
   |                    | - :param:`k` permeability (:math:`m^{2}`)                                                    |
   |                    | - :param:`vis` dynamic viscosity (:math:`Pa.s`)                                              |
   |                    | - :param:`contype` unsaturated flow allowed when contype=1                                   |
   |                    | - :param:`thetas` saturated water content                                                    |
   |                    | - :param:`thetar` residual water content                                                     |
   |                    | - :param:`paev` air-entry value of capillary pressure                                        |
   |                    | - :param:`m` van Genuchten parameter                                                         |
   |                    | - :param:`a` van Genuchten parameter                                                         |
   |                    | - :param:`thetam` additional model parameter for the modified version of van Genuchten’s     |
   |                    |   retention model                                                                            |
   |                    | - :param:`ctor` coefficient of tortuosity (:math:`ctor=\frac{1}{\tau}\le1`)                  |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 2dMassLatticeTransport                                                                       |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `lm/lattice2dcrackinput.in                                                                   |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice2dcrackinput.in>`_,    |
   |                    | `lm/lattice3d_mt1.in                                                                         |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice3d_mt1.in>`_,          |
   |                    | `lm/lattice3d_mt2.in                                                                         |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice3d_mt2.in>`_,          |
   |                    | `lm/latticetransmat.in                                                                       |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/latticetransmat.in>`_         |
   +--------------------+----------------------------------------------------------------------------------------------+

References
~~~~~~~~~~

.. [Baehr:06] H.D. Baehr and K. Stephan: Heat and Mass Transfer, Springer, 2006.

.. [Bazant:72] Z.P. Bažant, L.J. Najjar. Nonlinear water diffusion in nonsaturated concrete. *Materials and Structures*, 5:3--20, 1972.

.. [Bazant:98] Z.P. Bažant, J. Planas: Fracture and Size Effect in Concrete and Other Quasibrittle Materials, CRC Press, 1998.

.. [BSB] S. Brunauer, J. Skalny, E.E. Bodor: J. Colloid Interface Sci, 30, 1969.

.. [NISTIR7232] D.P. Bentz: CEMHYD3D: A Three-Dimensional Cement Hydration and Microstructure Development Modeling Package. Version 3.0., NIST Building and Fire Research Laboratory, Gaithersburg, Maryland, Technical report, 2005.

.. [Cervera:99] M. Cervera, J. Oliver, and T. Prato: Thermo-chemo-mechanical model for concrete. I: Hydration and aging. *Journal of Engineering Mechanics ASCE*, 125(9):1018--1027, 1999.

.. [Gawin:06a] D. Gawin, F. Pesavento, and B. A. Schrefler: Hygro-thermo-chemo-mechanical modelling of concrete at early ages and beyond. Part I: Hydration and hygro-thermal phenomena. *International Journal for Numerical Methods in Engineering*, 67(3):299--331, 2006.

.. [GraJir] P. Grassl and M. Jirásek. Damage-plastic model for concrete failure. *International Journal of Solids and Structures*, 43:7166--7196, 2006.

.. [GraXenNys13] P. Grassl, D. Xenos, U. Nyström, R. Rempling and K. Gylltoft. CDPM2: A damage-plasticity approach to modelling the failure of concrete. *International Journal of Solids and Structures*, 50:3805--3816, 2013.

.. [GraJir10] P. Grassl and M. Jirásek. Meso-scale approach to modelling the fracture process zone of concrete subjected to uniaxial tension. *International Journal of Solids and Structures*, 47: 957--968, 2010.

.. [AlbGal17] M.G. Alberti, A. Enfedaque, J.C. Gálvez and E. Reyes. Numerical modelling of the fracture of polyolefin fibre reinforced concrete by using a cohesive fracture approach. *Composites Part B: Engineering*, 111: 200--210, 2017.

.. [Gra14] P. Grassl, D. Xenos, M. Jirásek and M. Horák. Evaluation of nonlocal approaches for modelling fracture near nonconvex boundaries. *International Journal of Solids and Structures*, 51: 3239–-3251, 2014.

.. [inel] M. Jirásek, Z.P. Bažant: Inelastic analysis of structures, John Wiley, 2001.

.. [Hansen] P.F. Hansen: Coupled Moisture/Heat Transport in Cross Sections of Structures, Beton og Konstruktionsinstituttet, 1985.

.. [Hoek] E. Hoek and Z.T. Bieniawski: Brittle Rock Fracture Propagation In Rock Under Compression, International Journal of Fracture Mechanics 1(3), 137-155, 1965.

.. [Kuenzel] H.M. K"unzel, H.M.: Simultaneous heat and moisture transport in building components, Ph.D. thesis, IRB-Verlag, 1995.

.. [mdm] M. Jirásek: Comments on microplane theory, Mechanics of Quasi-Brittle Materials and Structures, ed. G. Pijaudier-Cabot, Z. Bittnar, and B. Gérard, Hermès Science Publications, Paris, 1999, pp. 57-77.

.. [Rots] B. Lourenco, J.G. Rots: Multisurface Interface Model for Analysis of Masonry Structures, Journal of Engng Mech, vol. 123, No. 7, 1997.

.. [ortiz] M. Ortiz, E.P. Popov: Accuracy and stability of integration algorithms for elasto-plastic constitutive relations, Int. J. Numer. Methods Engrg, 21, 1561-1576, 1985.

.. [oofem] B. Patzák: OOFEM home page, http://www.oofem.org, 2003.

.. [SimoHughes] J.C. Simo, T.J.R. Hughes: Computational Inelasticity, Springer, 1998.

.. [Ruiz:01] J. Ruiz, A. Schindler, R. Rasmussen, P. Kim, G. Chang: Concrete temperature modeling and strength prediction using maturity concepts in the FHWA HIPERPAV software, 7th international conference on concrete pavements, Orlando (FL), USA, 2001.

.. [simo] J.C. Simo, J.G. Kennedy, S. Govindjee: Non-smooth multisurface plasticity and viscoplasticity. Loading/unloading conditions and numerical algorithms, Int. J. Numer. Methods Engrg, 26, 2161-2185, 1988.

.. [Schindler:2005] A. K. Schindler and K. J. Folliard: Heat of Hydration Models for Cementitious Materials, ACI Materials Journal, 102, 24 - 33, 2005.

.. [SimoPister] J.C. Simo, K.S. Pister: Remarks on rate constitutive equations for finite deformation problems: computational implications, Comp Methods in Applied Mech and Engng, 46, 201-215, 1984.

.. [Smilauer:09] V. Šmilauer and T. Krejčí, Multiscale Model for Temperature Distribution in Hydrating Concrete, International Journal for Multiscale Computational Engineering, 7 (2), 135-151, 2009.

.. [Vree:95] J.H.P. de Vree, W.A.M. Brekelmans, and M.A.J. van Gils: Comparison of nonlocal approaches in continuum damage mechanics. Computers and Structures 55(4), 581–588, 1995.

.. [Xi] Y. Xi, Z.P. Bažant, H.M. Jennings: Moisture Diffusion in Cementitious Materials, Advn Cem Bas Mat, 1994.

.. [Bazant-89-I] Z.P. Bažant, S. Prasannan: Solidification theory for concrete creep. I: Formulation. Journal of Engineering Mechanics 115(8), 1691–1703, 1989.

.. [Bazant-97-I] Z.P. Bažant, A.B. Hauggaard, F. Ulm: Microprestress-solidification theory for concrete creep. I: Aging and drying effects. Journal of Engineering Mechanics 123(11), 1188–1194, 1997.

.. [JirHav14a] M. Jirásek, P. Havlásek: Microprestress-Solidification Theory of Concrete Creep: Reformulation and Improvement. Cement and Concrete Research 60, 51–62, 2014.

.. [fib:2010] International Federation for Structural Concrete (*fib*): The *\ fib* Model Code for Concrete Structures 2010.

.. [GraDav11] P. Grassl and T. Davies: Lattice modelling of corrosion induced cracking and bond in reinforced concrete. Cement and Concrete Composites, 33, 918-924, 2011.

.. [AthWheGra18] I. Athanasiadis, S. Wheeler and P. Grassl: Hydro-mechanical network modelling of particulate composites. International Journal of Solids and Structures 130-131, 49-60, 2018.

.. [SciGraLarRun20] A. Sciegaj and P. Grassl and F. Larsson and K. Runesson and K. Lundgren. “Upscaling of three-dimensional reinforced concrete representative volume elements to effective beam and plate models”, International Journal of Solids and Structures, vol. 202, pp. 835-853, 2020.

.. [ZhoMarGra24] C. Zhou and A. Marlot and P. Grassl: CDPM2F: A Damage-Plasticity Approach to Modelling the Failure of Engineered Cementitious Composites, Available at SSRN 4819173, 2024. \endthebibliography

   :labelprefix: ref
   :style: plain
