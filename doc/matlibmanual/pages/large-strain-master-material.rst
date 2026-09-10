
.. _large-strain-master-material:

Large-strain master material
============================

In this section, a large-strain master model based on generalized stress-strain measures is described. 
In the first step, strain measure is computed from equation :eq:`generalizedStrainMeasures`.

.. math::
   :label: generalizedStrainMeasures

   \boldsymbol{E}^{(m)} = \begin{cases}
    \displaystyle{\frac{1}{2m}}\left(\boldsymbol{C}^m-\boldsymbol{I}\right), & \text{if }m \neq 0 \\
    \\
    \displaystyle{\frac{1}{2}}\ln\boldsymbol{C}, & \text{if }m = 0
   \end{cases} 
   where $\boldsymbol{I}$ is the second-order unit tensor and $\boldsymbol{C} = \boldsymbol{F}^T\boldsymbol{F}$ is Cauchy-Green strain tensor. In the special cases when $m = 0$ and $m = 0.5$ we obtain the so-called Hencky (logarithmic) and Biot tensor, while for $m = 1$ we obtain the right Green-Lagrange strain tensor.   

In the second step, this strain measure enters a constitutive law of slave material and the stress measure conjugated to the strain measure defined in step one and appropriate stiffness matrix are computed. In the third step, the generalized stress tensor and stiffness matrix are transformed into the second Piola-Kirchhoff stress and the appropriate stiffness tensor. 
The model description and parameters are summarized in :numref:`LSmasterMat_table`.

.. table:: Large-strain master material material - summary.
   :name: LSmasterMat_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Large-strain master material material                                                        |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`LSmasterMat` :elemparam:`in` :elemparam:`m{rn}` :elemparam:`slavemat{in}`         |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - material number                                                                            |
   |                    | - :param:`m` parameter defining the strain measure                                           |
   |                    | - :param:`slavemat` number of slave material                                                 |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMatF                                                                                       |
   +--------------------+----------------------------------------------------------------------------------------------+
