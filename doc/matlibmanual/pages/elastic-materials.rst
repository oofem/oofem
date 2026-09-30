Elastic materials
=================


.. _IsoLE:

Isotropic linear elastic material - IsoLE
-----------------------------------------

Linear isotropic material model. The model parameters are summarized
in :numref:`IsoLE_table`.

.. table:: Linear Isotropic Material - summary.
   :name: IsoLE_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Linear isotropic elastic material                                                            |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`IsoLE` :elemparam:`num{in}` :elemparam:`d{rn}` :elemparam:`E{rn}`                 |
   |                    | :elemparam:`n{rn}` :elemparam:`tAlpha{rn}`                                                   |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`num` material model number                                                         |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`E` Young modulus                                                                   |
   |                    | - :param:`n` Poisson ratio                                                                   |
   |                    | - :param:`tAlpha` thermal dilatation coefficient                                             |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat, PlaneStress, PlaneStrain, 1dMat, 2dPlateLayer, 2dBeamLayer, 3dShellLayer, 2dPlate,    |
   |                    | 2dBeam, 3dShell, 3dBeam, PlaneStressRot                                                      |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Features           | Adaptivity support                                                                           |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `benchmark/sm/tr2shell7.in                                                                   |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/benchmark/sm/tr2shell7.in>`_,    |
   |                    | `mpm/mpms01.in <https://github.com/oofem/oofem/blob/devel/tests/regression/mpm/mpms01.in>`_, |
   |                    | `partests/brazil_2d_nl2/brazil_2d_nl.oofem.in <https://github.com/oofem/oofem/blob/devel/tes |
   |                    | ts/regression/partests/brazil_2d_nl2/brazil_2d_nl.oofem.in>`_, `sm/bondceb01.in              |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/bondceb01.in>`_,              |
   |                    | `sm/linkslip02.in                                                                            |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/linkslip02.in>`_,             |
   |                    | `sm/patch011.in                                                                              |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/patch011.in>`_,               |
   |                    | `sm/patch140.in                                                                              |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/patch140.in>`_,               |
   |                    | `sm/patch150.in                                                                              |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/patch150.in>`_,               |
   |                    | `sm/spring05.in                                                                              |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/spring05.in>`_,               |
   |                    | `sm/trshell02_test.in                                                                        |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/trshell02_test.in>`_ (and 199 |
   |                    | more)                                                                                        |
   +--------------------+----------------------------------------------------------------------------------------------+

.. _OrthoLE:

Orthotropic linear elastic material - OrthoLE
---------------------------------------------

Orthotropic, linear elastic material model. The model parameters are summarized
in :numref:`OrthoLE_table`.

Local coordinate system, which determines axes of material orthotropy
can be specified using :param:`lcs` array. This array contains six numbers,
where the first three numbers represent directional vector of a local x-axis,
and next three numbers represent directional vector of a local y-axis.
The local z-axis is determined using the vector product.
The right-hand coordinate system is assumed.

Local coordinate system
can also be specified using :param:`scs` parameter. Then local coordinate
system is specified in so called “shell”
coordinate system, which is defined locally on each particular element
and its definition is as follows: principal z-axis is perpendicular to
shell mid-section, x-axis is perpendicular to z-axis and normal to
user specified vector (so x-axis is parallel to plane, with  being
normal to this plane) and y-axis is perpendicular both to x and z
axes. This definition of coordinate system can be used only with plates
and shells elements.
When vector  is parallel to z-axis an error occurs. The :param:`scs` array contain three numbers
defining direction vector . If no local coordinate system is
specified, by default a global coordinate system is used.

For 3D case the material compliance matrix has the following form

.. math::

   \mbf{C}=\left[\begin{array}{cccccc}
   1/E_X & -\nu_{xy}/E_x& -\nu_{xz}/E_x& 0& 0&0\\
   -\nu_{yx}/E_y & 1/E_y& -\nu_{yz}/E_y& 0& 0&0\\
   -\nu_{zx}/E_z & -\nu_{zy}/E_z& 1/E_z& 0& 0&0\\
   0 & 0 & 0 & 1/G_{yz} & 0 & 0\\
   0 & 0 & 0 & 0 & 1/G_{xz} & 0\\
   0 & 0 & 0 & 0 & 0 & 1/G_{xy}
   \end{array}\right]

By inversion, the material stiffness matrix has the form

.. math::

   \mbf{D}=\left[\begin{array}{cccccc}
   d_{xx} & d_{xy} & d_{xz} & 0 & 0 & 0\\
   & d_{yy} & d_{yz} & 0 & 0 & 0\\
   \rm {sym} & & d_{zz} & 0 & 0 & 0\\
   0 & 0 & 0 & G_{yz} & 0 & 0\\
   0 & 0 & 0 & 0 & G_{xz} & 0\\
   0 & 0 & 0 & 0 & 0 & G_{xy}
   \end{array}\right]

where :math:`\xi=1-(\nu_{xy}*\nu_{yx}+\nu_{yz}*\nu_{zy}+\nu_{zx}*\nu_{xz})-(\nu_{xy}*\nu_{yz}*\nu_{zx}+\nu_{yx}*\nu_{zy}*\nu_{xz})` and

.. math::

   \begin{eqnarray}
   d_{xx}&=&E_X(1-\nu_{yz}*\nu_{zy})/\xi,\\
   d_{xy}&=&E_y*(\nu_{xy}+\nu_{xz}*\nu_{zy})/\xi,\\
   d_{xz}&=&E_z*(\nu_{xz}+\nu_{yz}*\nu_{xy})/\xi,\\
   d_{yy}&=&E_y*(1-\nu_{xz}*\nu_{zx})/\xi,\\
   d_{yz}&=&E_z*(\nu_{yz}+\nu_{yx}*\nu_{xz})/\xi,\\
   d_{zz}&=&E_z*(1-\nu_{yx}*\nu_{xy})/\xi.
   \end{eqnarray}

:math:`E_i` is Young's modulus in the :math:`i`-th direction, :math:`G_{ij}` is the shear modulus in :math:`ij` plane, :math:`\nu_{ij}` is the major Poisson ratio, and :math:`\nu_{ji}` is the minor Poisson ratio. Assuming that :math:`E_x>E_y>E_z`, :math:`\nu_{xy} > \nu_{yx}` etc., then :math:`\nu_{xy}` is referred to as the major Poisson ratio, while :math:`\nu_{yx}` is referred as the minor Poisson ratio.

Note that there are only nine independent material parameters,
because of symmetry conditions. The symmetry condition yields

.. math::

   \nu_{xy}E_y=\nu_{yx}E_x,\ \ \nu_{yz}E_z=\nu_{zy}E_y,\ \ \nu_{zx}E_x=\nu_{xz}E_z

The model description and parameters are summarized
in :numref:`OrthoLE_table`.

.. table:: Orthotropic, linear elastic material -- summary.
   :name: OrthoLE_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Orthotropic, linear elastic material                                                         |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`OrthoLE` :elemparam:`num{in}` :elemparam:`d{rn}` :elemparam:`Ex{rn}`              |
   |                    | :elemparam:`Ey{rn}` :elemparam:`Ez{rn}` :elemparam:`NYyz{rn}` :elemparam:`NYxz{rn}`          |
   |                    | :elemparam:`NYxy{rn}` :elemparam:`Gyz{rn}` :elemparam:`Gxz{rn}` :elemparam:`Gxy{rn}`         |
   |                    | :elemparam:`tAlphax{rn}` :elemparam:`tAlphay{rn}` :elemparam:`tAlphaz{rn}`                   |
   |                    | :optelemparam:`lcs{ra}` :optelemparam:`scs{ra}`                                              |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`num` material model number                                                         |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`Ex`, :param:`Ey`, :param:`Ez` Young moduli for x,y, and z directions               |
   |                    | - :param:`NYyz`, :param:`NYxz`, :param:`NYxy` major Poisson's ratio coefficients             |
   |                    | - :param:`Gyz`, :param:`Gxz`, :param:`Gxy` shear moduli                                      |
   |                    | - :param:`tAlphax`, :param:`tAlphay`, :param:`tAlphaz` thermal dilatation coefficients in    |
   |                    |   x,y,z directions                                                                           |
   |                    | - :param:`lcs` Array defining local material x and y axes of orthotrophy                     |
   |                    | - :param:`scs` Array defining a normal vector n. The local x axis is parallel to plane with  |
   |                    |   n being plane normal. The material local z-axis is perpendicular to shell mid-section.     |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat, PlaneStress, PlaneStrain, 1dMat, 2dPlateLayer, 2dBeamLayer, 3dShellLayer, 2dPlate,    |
   |                    | 2dBeam, 3dShell, 3dBeam, PlaneStressRot                                                      |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `benchmark/Tr2shell7XFEM_adaptiveDelam.in <https://github.com/oofem/oofem/blob/devel/tests/r |
   |                    | egression/benchmark/Tr2shell7XFEM_adaptiveDelam.in>`_,                                       |
   |                    | `benchmark/Tr2shell7XFEM_stressRecovery.in <https://github.com/oofem/oofem/blob/devel/tests/ |
   |                    | regression/benchmark/Tr2shell7XFEM_stressRecovery.in>`_, `sm/layered_cube.in                 |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/layered_cube.in>`_,           |
   |                    | `sm/layered_cube_lcs.in                                                                      |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/layered_cube_lcs.in>`_,       |
   |                    | `sm/materOrient01.in                                                                         |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/materOrient01.in>`_,          |
   |                    | `sm/plate_a1000.0_n10_h100.0.in <https://github.com/oofem/oofem/blob/devel/tests/regression/ |
   |                    | sm/plate_a1000.0_n10_h100.0.in>`_                                                            |
   +--------------------+----------------------------------------------------------------------------------------------+

.. _AnisoLE:

General anisotropic linear elastic material - AnisoLE
-----------------------------------------------------

Linear elastic material model with completely general material stiffness (21 independent elastic constants). The model parameters are summarized in :numref:`AnisoLE_table`. The material stiffness matrix is a completely arbitrary symmetric :math:`6 \times 6` matrix. The input line must contain an array with the upper triangle of this matrix. This array has length 21 and the stiffness coefficients are listed in each row from the diagonal to the last column.

.. table:: Anisotropic, linear elastic material -- summary.
   :name: AnisoLE_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Anisotropic, linear elastic material                                                         |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`OrthoLE` :elemparam:`num{in}` :elemparam:`d{rn}` :elemparam:`stiff{ra}`           |
   |                    | :elemparam:`tAlphax{ra}`                                                                     |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`num` material model number                                                         |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`stiff` real array of length 21 with stiffness coefficients :math:`D_{11}`,         |
   |                    |   :math:`D_{12}`, :math:`D_{13}\ldots D_{66}`                                                |
   |                    | - :param:`tAlpha` real array of length 0 or 3 with thermal dilatation coefficients in x,y,z  |
   |                    |   directions                                                                                 |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat, PlaneStress, PlaneStrain, 1dMat                                                       |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `benchmark/Tr2shell7XFEM_adaptiveDelam.in <https://github.com/oofem/oofem/blob/devel/tests/r |
   |                    | egression/benchmark/Tr2shell7XFEM_adaptiveDelam.in>`_,                                       |
   |                    | `benchmark/Tr2shell7XFEM_stressRecovery.in <https://github.com/oofem/oofem/blob/devel/tests/ |
   |                    | regression/benchmark/Tr2shell7XFEM_stressRecovery.in>`_, `sm/layered_cube.in                 |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/layered_cube.in>`_,           |
   |                    | `sm/layered_cube_lcs.in                                                                      |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/layered_cube_lcs.in>`_,       |
   |                    | `sm/materOrient01.in                                                                         |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/materOrient01.in>`_,          |
   |                    | `sm/plate_a1000.0_n10_h100.0.in <https://github.com/oofem/oofem/blob/devel/tests/regression/ |
   |                    | sm/plate_a1000.0_n10_h100.0.in>`_                                                            |
   +--------------------+----------------------------------------------------------------------------------------------+

.. _isoaxysymm1d:


1D linear elastic material with different tension and compression stiffness - isoAxysymm1D
------------------------------------------------------------------------------------------


Bi-linear material model for 1D elasticity, with different elastic moduli in tension and compression. The model parameters are summarized in :numref:`isoaxysymm1d_table`.

.. table:: Bi-Linear Elastic Material - summary.
   :name: isoaxysymm1d_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Bi-Linear elastic material                                                                   |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`IsoAxysymm1D` :elemparam:`num{in}` :elemparam:`d{rn}` :elemparam:`Et{rn}`         |
   |                    | :elemparam:`Ec{rn}` :elemparam:`tAlpha{rn}` [:elemparam:`m{rn}`]                             |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`num` material model number                                                         |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`Et` Modulus of elasticity in tension                                               |
   |                    | - :param:`Ec` Modulus of elasticity in compression                                           |
   |                    | - :param:`tAlpha` thermal dilatation coefficient                                             |
   |                    | - :param:`m` optional regularization coefficient, default value set to 15, higher value      |
   |                    |   makes response close to trully bilinear                                                    |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 1dMat                                                                                        |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Features           | Adaptivity support                                                                           |
   +--------------------+----------------------------------------------------------------------------------------------+

.. _hyperelastic_material_simo_pister_material:

Hyperelastic material - Simo-Pister Material
--------------------------------------------

This material model can describe elastic behavior at large strains. A hyperelastic model postulates the existence of a free energy potential. The existence of the potential implies reversibility of deformations and no energy dissipation during the loading process. Here, we use the free energy function introduced in [SimoHughes]_.

.. math::
   :label: freeEnergy

   \rho_0 \psi = \frac{1}{4} \left(K - \frac{2}{3}G \right) \left( J^2 - 2 \ln J - 1 \right) + G \left(\mathbf{E} : \mathbf{I} - \ln J \right)

where :math:`K` is the bulk modulus, :math:`G` is the shear modulus, :math:`J` is the Jacobian (determinant of the deformation gradient, corresponding to the ratio of the current and initial volume), and :math:`\mathbf{E}` is the Green-Lagrange strain. Then, the stress-strain law can be derived from (:eq:`freeEnergy`) as

.. math::

   \mathbf{S} = \rho_0 \frac{\partial \psi}{\partial \mathbf{E}} = \frac{1}{2}\left(K - \frac{2}{3}G \right) \left(J^2 - 1 \right)\mathbf{C}^{-1} + G \left(\mathbf{I} -\mathbf{C}^{-1} \right)

where :math:`\mathbf{S}` is the second Piola-Kirchhoff stress, :math:`\mathbf{E}` is the Green-Lagrange strain, and :math:`\mathbf{C}` is the right Cauchy-Green tensor.

The model description and parameters are summarized in :numref:`hyperElMat_table`.

.. table:: Hyperelastic material - summary.
   :name: hyperElMat_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Hyperelastic material                                                                        |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`SimoPisterMat` :elemparam:`in` :elemparam:`d{rn}` :elemparam:`K{rn}`              |
   |                    | :elemparam:`G{rn}`                                                                           |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - material number                                                                            |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`K` bulk modulus                                                                    |
   |                    | - :param:`G` shear modulus                                                                   |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat                                                                                        |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `benchmark/sm/contact/friction_wip/visual/contact3d_two_bricks_press_rotate_visual.in <https |
   |                    | ://github.com/oofem/oofem/blob/devel/tests/regression/benchmark/sm/contact/friction_wip/visu |
   |                    | al/contact3d_two_bricks_press_rotate_visual.in>`_                                            |
   +--------------------+----------------------------------------------------------------------------------------------+

.. _hyperelastic_material_compressible_mooney_rivlin:

Hyperelastic material - Compressible Mooney-Rivlin
--------------------------------------------------

The Mooney-Rivlin strain energy function is expressed by

.. math::
   :label: freeEnergyMR

   \rho_0 \psi = C_1(\bar{I}_1-3) + C_2(\bar{I}_2-3) + \frac{1}{2}K (\ln J)^2

where :math:`C_1` and :math:`C_2` are material constants, :math:`K` is the bulk modulus, :math:`J` is the Jacobian (determinant of the deformation gradient, corresponding to the ratio of the current and initial volume),  :math:`\bar{I}_1 = J^{-\frac{2}{3}} I_1`, :math:`\bar{I}_2 = J^{-\frac{2}{3}} I_2`, 
where :math:`I_1` and :math:`I_2` are the first and the second principal invariants of the right Cauchy-Green deformation tensor :math:`\mathbf{C}`.
Compressible neo-Hookean material model is obtained by setting :math:`C_2 = 0`.
Then stress-strain law can be derived from (:eq:`freeEnergyMR`) as

.. math::

   \mathbf{P} = \rho_0 \frac{\partial \psi}{\partial \mathbf{F}} = C_1 \frac{\partial\bar{I}_1}{\partial\mathbf{F}} + C_2\frac{\partial\bar{I}_2}{\partial\mathbf{F}} + K \ln J \mathbf{F}^{-T}

where
:math:`\mathbf{P}` is the first Piola-Kirchhoff stress,

.. math::

   \frac{\partial \bar{I}_1}{\partial\mathbf{F}} = \frac{2}{J^{\frac{2}{3}}}\mathbf{F} - \frac{2}{3}\bar{I}_1\mathbf{F}^{-T}

and

.. math::

   \frac{\partial \bar{I}_2}{\partial\mathbf{F}} = 2\bar{I}_1\mathbf{F}-\frac{4}{3}\bar{I}_2\mathbf{F}^{-T}-\frac{2}{J^\frac{4}{3}}\mathbf{F}\cdot \mathbf{C}

The first elasticity tensor is derived as

.. math::

   A_{ijkl} = \frac{\partial P_{ij}}{\partial F_{kl}} = C_1 A^1_{ijkl} + C_2  A^2_{ijkl} + K(F_{ji}^{-1}F_{lk}^{-1} - \ln J F_{jk}^{-1}F_{li}^{-1})

where 
:math:`A^1_{ijkl} = \frac{2}{3} J^{-\frac{2}{3}}\left[3\delta_{ik}\delta_{jl} + I_{1}F_{jk}^{-1}F_{li}^{-1} - 2F_{lk}^{-1}F_{ij}+\frac{2}{3}I_1 F_{ji}^{-1}F_{lk}^{-1}-2F_{ji}^{-1}F_{kl} \right]`

and

.. math::
   :label: A2ijkl

   A^2_{ijkl} = 2 J^{-\frac{4}{3}}\left[ I_1 \delta_{ik} \delta_{jl} + 2 F_{ij}F_{kl}-\frac{4}{3}I_1 F_{ij} F_{lk}^{-1} -\frac{8}{9}I_2 F_{ji}^{-1}F_{lk}^{-1}-\frac{4}{3}I_1 F_{ji}^{-1}F_{kl} \nonumber\right.\\ 
   \left.  + \frac{4}{3} F_{kn}C_{nl}F_{ji}^{-1}  +\frac{2}{3}I_2 F_{li}^{-1}F_{jk}^{-1} + \frac{4}{3}F_{kl}^{-1}F_{im}C_{mj} - \delta_{ik}C_{lj} - F_{il}F_{kj} + F_{km}F_{im}\delta_{jl}\right]

The model description and parameters are summarized in :numref:`MooneyRivlin_table`.

.. table:: Compressible Mooney-Rivlin - summary.
   :name: MooneyRivlin_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Mooney-Rivlin                                                                                |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`MooneyRivlinCompressibleMat` :elemparam:`in` :elemparam:`d{rn}`                   |
   |                    | :elemparam:`K{rn}` :elemparam:`C1{rn}` :elemparam:`C2{rn}`                                   |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - material number                                                                            |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`K` bulk modulus                                                                    |
   |                    | - :param:`C1` material constant                                                              |
   |                    | - :param:`C2` material constant                                                              |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat, PlaneStrain                                                                           |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `benchmark/sm/contact/friction_wip/visual/contact2d_inclined_plane_slip.in <https://github.c |
   |                    | om/oofem/oofem/blob/devel/tests/regression/benchmark/sm/contact/friction_wip/visual/contact2 |
   |                    | d_inclined_plane_slip.in>`_, `sm/contact2d.in                                                |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/contact2d.in>`_,              |
   |                    | `sm/contact2d_linesearch.in                                                                  |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/contact2d_linesearch.in>`_,   |
   |                    | `sm/contact2d_projection_boundary.in <https://github.com/oofem/oofem/blob/devel/tests/regres |
   |                    | sion/sm/contact2d_projection_boundary.in>`_, `sm/contact2d_projection_outside.in <https://gi |
   |                    | thub.com/oofem/oofem/blob/devel/tests/regression/sm/contact2d_projection_outside.in>`_,      |
   |                    | `sm/contact2d_sliding.in                                                                     |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/contact2d_sliding.in>`_,      |
   |                    | `sm/contact3d_generalized_edge_direct.in <https://github.com/oofem/oofem/blob/devel/tests/re |
   |                    | gression/sm/contact3d_generalized_edge_direct.in>`_, `sm/friction_wip/contact2d_friction.in  |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/friction_wip/contact2d_fricti |
   |                    | on.in>`_, `sm/mooneyrivlin1.in                                                               |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/mooneyrivlin1.in>`_,          |
   |                    | `sm/mooneyrivlin2.in                                                                         |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/mooneyrivlin2.in>`_ (and 16   |
   |                    | more)                                                                                        |
   +--------------------+----------------------------------------------------------------------------------------------+


Hyperelastic material - Compressible Ogden
------------------------------------------


The Ogden strain energy function is expressed by

.. math::
   :label: freeEnergyOgdenDeviatoricStretch

   \rho_0 \psi = \sum_{I=1}^N \frac{\mu_I}{\alpha_I} \left( \bar{\lambda}_1^{\alpha_I}+\bar{\lambda}_2^{\alpha_I}+\bar{\lambda}_3^{\alpha_I}-3\right)+ U(J)

where :math:`\mu_I` and :math:`\alpha_I` are material constants, :math:`U(J)` is the volumetric part of energy, :math:`J` is the Jacobian (determinant of the deformation gradient, corresponding to the ratio of the current and initial volume),  :math:`\bar{\lambda}_i = J^{-\frac{1}{3}} {\lambda}_i` are the deviatoric stretches.

The model description and parameters are summarized in :numref:`CompressibleOgden_table`.

.. table:: Compressible Ogden material - summary.
   :name: CompressibleOgden_table

   +-------------------------------------------+----------------------------------------------------------------------------------------------+
   | Description                               | Ogden material                                                                               |
   +-------------------------------------------+----------------------------------------------------------------------------------------------+
   | Record Format                             | :descitem:`OgdenCompressibleMat` :elemparam:`in` :elemparam:`d{rn}` :elemparam:`K{rn}`       |
   |                                           | :elemparam:`alpha{ra}` :elemparam:`mu{ra}`                                                   |
   +-------------------------------------------+----------------------------------------------------------------------------------------------+
   | Parameters                                | - material number                                                                            |
   |                                           | - :param:`d` material density                                                                |
   |                                           | - :param:`K` bulk modulus                                                                    |
   |                                           | - :param:`alpha` array of material constants                                                 |
   |                                           | - :param:`mu` array of material constants                                                    |
   +-------------------------------------------+----------------------------------------------------------------------------------------------+
   | alpha and mu have to have the same size   |                                                                                              |
   +-------------------------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples                            | `sm/ogden1.in <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/ogden1.in>`_,   |
   |                                           | `sm/ogden2.in <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/ogden2.in>`_,   |
   |                                           | `sm/timestepreduction.in                                                                     |
   |                                           | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/timestepreduction.in>`_       |
   +-------------------------------------------+----------------------------------------------------------------------------------------------+


Hyperelastic material - Blatz-Ko
--------------------------------


The Mooney-Rivlin strain energy function is expressed by

.. math::
   :label: freeEnergyBlatzKo

   \rho_0 \psi = \frac{\mu}{2}\left(\frac{I_2}{I_3} + 2\sqrt{I_3}-5 \right)

where :math:`\mu` is initial shear modulus, :math:`I_2` and :math:`I_3` are the second and third invariants of the Cauchy-Green tensor :math:`\mbf{C}`.

The model description and parameters are summarized in :numref:`BlatzKo_table`.

.. table:: Blatz-Ko material - summary.
   :name: BlatzKo_table

   +--------------------------------+----------------------------------------------------------------------------------------------+
   | Description                    | Blatz-Ko material                                                                            |
   +--------------------------------+----------------------------------------------------------------------------------------------+
   | Record Format                  | :descitem:`blatzkomat` :elemparam:`in` :elemparam:`d{rn}` :elemparam:`mu{rn}`                |
   +--------------------------------+----------------------------------------------------------------------------------------------+
   | Parameters                     | - material number                                                                            |
   |                                | - :param:`d` material density                                                                |
   |                                | - :param:`mu` shear modulus                                                                  |
   +--------------------------------+----------------------------------------------------------------------------------------------+
   | :math:`\nu` is fixed to 0.25   |                                                                                              |
   +--------------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples                 | `sm/blatzko1.in                                                                              |
   |                                | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/blatzko1.in>`_                |
   +--------------------------------+----------------------------------------------------------------------------------------------+
