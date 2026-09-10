Orthotropic damage model with fixed crack orientations for composites - CompDamMat
==================================================================================


.. _CompDamMat:

CompDamMat - Orthotropic damage model with fixed crack orientations for composites
----------------------------------------------------------------------------------

The model is designed for transversely isotropic elastic material defined by five elastic material constants. Typical example is a carbon fiber tow. Axis 1 represents the material principal direction. The orthotropic material constants are defined as

.. math::

   \nu_{12} &= \nu_{13},~\nu_{21}=\nu_{31},~\nu_{23}=\nu_{32},~E_{22}=E_{33}\\
   G_{12} &= G_{13}=G_{21}=G_{31},G_{23}=G_{32}\\
   \frac{\nu_{12}=\nu_{13}}{E_{11}} &= \frac{\nu_{21}=\nu_{31}}{E_{22}},~\frac{\nu_{31}=\nu_{21}}{E_{33}} = \frac{\nu_{13}=\nu_{12}}{E_{11}}

Material orientation on a finite element can be specified with the :param:`lcs` optional parameter. If unspecified, material orientation is the same as the global coordinate system. The :param:`lcs` array contains six numbers, where the first three numbers represent a directional vector of the local x-axis, and the next three numbers represent a directional vector of the local y-axis with the reference to the global coordinate system. The composite material is extended to 1D and is also suitable for beams and trusses. In such particular case, the :param:`lcs` has no effect and the 1D element orientation is aligned with the global xx component.

The index :math:`p,~p\in\{11,22,33,23,31,12\}` symbolizes six components of stress or strain vectors. The linear softening occurs after reaching a critical stress :math:`f_{p,0}`, see :numref:`comp_softening`. Orientation of cracks is assumed to be orthogonal and aligned with an orientation of material axes \cite[pp.236]{Bazant:98}. The transverse isotropy is generally lost upon fracture, material becomes orthotropic and six damage parameters :math:`d_p` are introduced.

.. figure:: /figures/Compodamagemat_diag.png
   :width: 70%
   :alt: Implemented stress-strain evolution with damage for 1D case. Tension and compression are separated, but sharing the same damage parameter.
   :align: center
   :name: comp_softening

   Implemented stress-strain evolution with damage for 1D case. Tension and compression are separated, but sharing the same damage parameter.

The compliance material matrix :math:`\mathbf{H}`, in the secant form and including damage parameters, reads

.. math::
   :label: comp_eq_H

   \mathbf{H}=
   \left[ \begin{array}{cccccc}
   \frac{1}{(1-d_{11})E_{11}} & -\frac{\nu_{21}}{E_{22}} & -\frac{\nu_{31}}{E_{33}} &0 &0 & 0\\
   -\frac{\nu_{12}}{E_{11}} & \frac{1}{(1-d_{22})E_{22}} & -\frac{\nu_{32}}{E_{33}} &0& 0&0\\
   -\frac{\nu_{13}}{E_{11}} & -\frac{\nu_{23}}{E_{22}} & \frac{1}{(1-d_{33})E_{33}} &0 &0 &0 \\
   0&0 &0 & \frac{1}{(1-d_{23})G_{23}} & 0&0 \\
   0& 0& 0& 0& \frac{1}{(1-d_{31})G_{31}} & 0\\
   0& 0&0 &0& 0& \frac{1}{(1-d_{12})G_{12}}
   \end{array} \right]

Damage occurs when any out of six stress tensor components exceeds a given strength :math:`f_{p,0}`

.. math::

   |\sigma_p| \geq |f_{p,0}|

Positive and negative ultimate strengths can be generally different but share the same damage variable.
At the point of damage initiation, see :numref:`comp_softening`, one evaluates :math:`\varepsilon_{p,E}` and characteristic element length :math:`l_p`, generally different for each damage mode. Given the fracture energy :math:`G_{F,p}`, the maximum strain at zero stress :math:`\varepsilon_{p,0}` is computed

.. math::
   :label: comp_eq_epsilon_p0

   \varepsilon_{p,0} =  \varepsilon_{p,E} + \frac{2G_{F,p}}{f_{p,0} l_p}

The point of damage initiation is never reached exactly, one needs to interpolate between the previous equilibrated step and current step to achieve objectivity.

The evolution of damage :math:`d_p` is based on the evolution of corresponding strain :math:`\varepsilon_p`. A maximum achieved strain is stored in the variable :math:`\kappa_p`. If :math:`\varepsilon_p > \kappa_p` the damage may grow so the corresponding damage variable :math:`d_p` may increase. Desired stress :math:`\sigma'_p` is evaluated from the actual strain :math:`\varepsilon_p`

.. math::
   :label: comp_eq_sigma_prime_p

   \sigma'_p =  f_{p,0} \frac{\varepsilon_{p,0} - \varepsilon_p}{\varepsilon_{p,0}-\varepsilon_{p,E}}

and the calculation of damage variables :math:`d_p` stems from :eq:`comp_eq_H`, for example

.. math::
   :label: comp_eq_dp11

   d_{11} &=& 1 - \frac{\sigma'_{11}}{E_{11}\left(\varepsilon_{11}+\frac{\nu_{21}}{E_{22}}\sigma_{22}+ \frac{\nu_{31}}{E_{33}}\sigma_{33} \right)}\\
   d_{12} &=& 1 - \frac{\sigma'_{12}}{G_{12}\varepsilon_{12}}

Damage is always controled not to decrease. :numref:`comp_performance` shows a typical performance for this damage model in one direction.

The damage initiation is based on a trial stress. It becomes necessary for higher precision to skip a few first iterations, typically 5, and then to introduce damage. A parameter :param:`afterIter` is designed for this purpose, and :param:`MinIter` forces a solver always to proceed a certain amount of iterations.

:param:`allowSnapBack` skips the checking of sufficient fracture energy for each direction. If not specified, all directions are checked to prevent snap-back, which dissipates an incorrect amount of energy.

.. figure:: /figures/Compodamagemat_test.png
   :width: 70%
   :alt: Typical loading/unloading material performance for homogenized stress and strain in the direction `:math:`{2}`'. Note that one damage parameter is common for both tension and compression.
   :name: comp_performance

   Typical loading/unloading material performance for homogenized stress and strain in the direction `:math:`{2}`'. Note that one damage parameter is common for both tension and compression.

.. table:: Orthotropic damage model with fixed crack orientations for composites -- summary
   :name: compdammat_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Orthotropic damage model with fixed crack orientations for composites                        |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`CompDamMat` :elemparam:`num{in}` :elemparam:`d{rn}` :elemparam:`Exx{rn}`          |
   |                    | :elemparam:`EyyEzz{rn}` :elemparam:`nuxynuxz{rn}` :elemparam:`nuyz{rn}`                      |
   |                    | :elemparam:`GxyGxz{rn}` :elemparam:`Tension_f0_Gf{ra}` :elemparam:`Compres_f0_Gf{ra}`        |
   |                    | [:elemparam:`afterIter{in}`] [:elemparam:`allowSnapBack`\ (ia)]                              |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`num` material model number                                                         |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`Exx` Young's modulus for principal direction :math:`xx`                            |
   |                    | - :param:`EyyEzz` Young's modulus in orthogonal directions to the principal direction        |
   |                    |   :math:`xx`                                                                                 |
   |                    | - :param:`nuxynuxz` Poisson's ratio in :math:`xy` and :math:`xz` directions                  |
   |                    | - :param:`nuyz` Poisson's ratio in :math:`yz` direction                                      |
   |                    | - :param:`GxyGxz` shear modulus in :math:`xy` and :math:`xz` directions                      |
   |                    | - :param:`Tension_f0_Gf` array with six pairs of positive numbers. Each pair describes       |
   |                    |   maximum stress in tension and fracture energy for each direction (:math:`xx`, :math:`yy`,  |
   |                    |   :math:`zz`, :math:`yz`, :math:`zx`, :math:`xy`)                                            |
   |                    | - :param:`Compres_f0_Gf` array with six pairs of numbers. Each pair describes maximum stress |
   |                    |   in compression (given as a negative number) and positive fracture energy for each          |
   |                    |   direction (:math:`xx`, :math:`yy`, :math:`zz`, :math:`yz`, :math:`zx`, :math:`xy`)         |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat, 1dMat                                                                                 |
   |                    |                                                                                              |
   |                    | - :param:`afterIter` how many iterations must pass until damage is computed from strains,    |
   |                    |   zero is default. User must ensure that the solver proceeds the minimum number of           |
   |                    |   iterations.                                                                                |
   |                    | - :param:`allowSnapBack` array to skip checking for snap-back. The array members are 1-6 for |
   |                    |   tension and 7-12 for compression components.                                               |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `sm/compoDamMat.in                                                                           |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/compoDamMat.in>`_             |
   +--------------------+----------------------------------------------------------------------------------------------+

.. _TrabBone3d:

TrabBone3d
~~~~~~~~~~

This model combines orthotropic elastoplasticity with isotropic damage.
Material orthotropy is described by the fabric tensor, i.e., a symmetric second-order
tensor with principal directions aligned with the axes of orthotropy and principal
values normalized such that their sum is 3. Elastic constants as well as coefficients
that appear in the yield condition are linked to the principal values of the
fabric tensor and to porosity. The yield condition is piecewise quadratic,
with different parameters in the regions of positive and negative volumetric strain.
