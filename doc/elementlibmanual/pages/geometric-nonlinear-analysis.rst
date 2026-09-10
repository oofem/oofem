
Geometric nonlinear analysis
============================

To take int account geometric nonlinearity for a specific element, the keyword
:param:`nlgeo` must be specified. The :param:`nlgeo` parameter defines which formulation
of the momentum balance is solved and what deformation measure that is computed and sent
to the consitiutive models (see :numref:`strain_tensor_table`).  If
:param:`nlgeo`\ =:math:`1`, then the momentum balance is set up in the reference
configuration in terms of the First Piola-Kirchhoff stress tensor :math:`\mathbf{P}` and
the deformation tensor :math:`\mathbf{F}` as energy conjugates. This is also refered to
a Total Lagrangian formulation. The balance equation in weak form reads

.. math::

   \int_{\Omega} \delta \mathbf{F} : \mathbf{P} \mathrm{d} \Omega = \int_{\Gamma} \delta \mathbf{x} \cdot \mathbf{t}_P  \mathrm{d} \Gamma
   + \int_{\partial \Omega} \delta \mathbf{x} \cdot \mathbf{b}_P  \mathrm{d} \Omega

This equation can be rewritten in terms of the displacement :math:`\mathbf{u}` and the
displacement gradient :math:`\mathbf{H}`

.. math::

   \int_{\Omega} \delta \mathbf{H} : \mathbf{P} \mathrm{d} \Omega = \int_{\Gamma} \delta \mathbf{u} \cdot \mathbf{t}_P  \mathrm{d} \Gamma
   + \int_{\partial \Omega} \delta \mathbf{u} \cdot \mathbf{b}_P  \mathrm{d} \Omega

This equation is nearly identical to the one for small strains except that another
stress measure is used and we have the virtual displacement gradient  instead of the
virtual strains.

The corresponding FE-formulation is obtained as

.. math::

   \int_{\Omega} \mathbf{B}_H^{\mathrm{T}} \cdot \mathbf{P} \mathrm{d} \Omega = \int_{\Gamma} \mathbf{N}^{\mathrm{T}} \cdot \mathbf{t}_P  \mathrm{d} \Gamma
   + \int_{\partial \Omega} \delta \mathbf{N}^{\mathrm{T}} \cdot \mathbf{b}_P  \mathrm{d} \Omega

with the tangent stiffness

.. math::

   \mathbf{K}_{\mathrm{T}} = \int_{\Omega} \mathbf{B}_H^{\mathrm{T}} \cdot \frac{\partial \mathbf{P}}{\partial \mathbf{F} } \cdot \mathbf{B}_H \mathrm{d} \Omega

Thus, for an element to support large deformations (in addition to small deformation) it
needs only to implement the :math:`\mathbf{B}_H` matrix.  Similar to the regular **B**
matrix, which gives the strains in Voigt form when multiplied with the solution vector
:math:`\textbf{a}`,  :math:`\textbf{B}_H` should give the displacement gradient in Voigt
form with 9 components for a full 3D state.

.. table:: Nonlinear geometry modes
   :name: strain_tensor_table

   +-------------+------------------------------+
   | nlgeo       | strain tensor                |
   +=============+==============================+
   | 0 (default) | Small-strain tensor          |
   +-------------+------------------------------+
   | 1           | Green-Lagrange strain tensor |
   +-------------+------------------------------+
   | 2           | Deformation gradient         |
   +-------------+------------------------------+
