.. _errorestimators:

Error estimators and indicators
===============================

The type of error estimator is selected by :param:`eetype`, and the
:optelemparam:`errorestimatorparams{...}` group in the same record then carries
the parameters belonging to it.

.. table:: Supported error estimators and indicators.
   :name: eetypestable

   +---------------------------+--------------+
   | Error estimator/indicator | ``eetype``   |
   +===========================+==============+
   | EET_SEI                   | 0            |
   +---------------------------+--------------+
   | EET_ZZEE                  | 1            |
   +---------------------------+--------------+
   | EET_CZZSI                 | 2            |
   +---------------------------+--------------+

``EET_SEI``
    A scalar error indicator.  It indicates the element error from the value of
    some suitable scalar — damage level or plastic strain level, for example —
    obtained from the element integration points and the corresponding material
    model.

``EET_ZZEE``
    The Zienkiewicz-Zhu error estimator.  It requires special element
    algorithms, which may not be available for every element type.

``EET_CZZSI``
    A combined criterion: the Zienkiewicz-Zhu error estimator in the elastic
    regime and the scalar error indicator in the non-linear regime.

.. important::

   The error on the element level is evaluated using the element's default
   integration rule.  For the Zienkiewicz-Zhu estimator the error, in the L2 or
   the energy norm, is evaluated from the difference between the computed and
   the "recovered" stresses, which are approximated with the same interpolation
   functions as the displacements.  In many cases the default integration order
   is therefore not sufficient and a higher integration order must be used on
   the elements — see the Element Library Manual and the element's ``nip``
   parameter.

Parameters of EET_SEI
---------------------

.. record:: Syntax

   :elemparam:`vartype{in}` :elemparam:`minlim{rn}`
   :elemparam:`maxlim{rn}` :elemparam:`mindens{rn}`
   :elemparam:`maxdens{rn}` :elemparam:`defdens{rn}`
   :optelemparam:`regionskipmap{ia}`
   :optelemparam:`remeshingdensityratio{rn}`

**Parameters**

:optparam:`regionskipmap{ia}`
    Allows some regions to be skipped: the error is not evaluated there and the
    default mesh density is used.  Its size should equal the number of regions,
    and a nonzero entry marks a region to skip.

:param:`vartype{in}`
    Type of internal variable used as the error indicator.  The currently
    supported value is ``1``, a damage-based indicator.

:param:`minlim{rn}`, :param:`maxlim{rn}`, :param:`mindens{rn}`, :param:`maxdens{rn}`, :param:`defdens{rn}`
    If the indicator value lies in the range (:param:`minlim`,
    :param:`maxlim`), the proposed mesh density is interpolated linearly within
    the range (:param:`mindens`, :param:`maxdens`).  If the indicator value is
    less than :param:`minlim`, :param:`defdens` is used as the required
    density; if it is larger than :param:`maxlim`, :param:`maxdens` is used.

:optparam:`remeshingdensityratio{rn}`
    Allowed ratio between the proposed and the actual density.  Remeshing is
    forced whenever the actual ratio falls below this value.  Default ``0.80``.

Parameters of EET_ZZEE
----------------------

.. record:: Syntax

   :elemparam:`normtype{in}` :elemparam:`requirederror{rn}`
   :elemparam:`minelemsize{rn}` :optelemparam:`regionskipmap{ia}`

**Parameters**

:optparam:`regionskipmap{ia}`
    As for ``EET_SEI``.

:param:`normtype{in}`
    Type of norm used when the error is evaluated: ``0``, the default, uses the
    L2 norm, ``1`` uses the energy norm.

:param:`requirederror{rn}`
    Required error to reach, as a fraction — that is, percent divided by 100.

:param:`minelemsize{rn}`
    Minimum limit on the element size.

Parameters of EET_CZZSI
-----------------------

The union of the ``EET_SEI`` and ``EET_ZZEE`` parameters.  The inelastic
regions are driven by ``EET_SEI``, the elastic ones by ``EET_ZZEE``.
