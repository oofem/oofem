.. _TimeFunctionsRecords:

(Time) function records
=======================

These records describe the functions that generally describe the variation of
components during the solution.  They are counted by :param:`nltf` in the
:ref:`ComponentsSizeRecord`, and the numbering must start at one and end at
*n*, where *n* is the number of records.

.. record:: Record syntax common to all functions

   :descitem:`TimeFunctType` :elemparam:`num{in}`
   :optelemparam:`initialValue{rn}` :optelemparam:`paramtype{in}`

**Common parameters**

:optparam:`initialValue{rn}`
    Controls how the increment of the receiver is evaluated for the *first*
    solution step: that increment is the difference between the value of the
    receiver at the first step and this initial value, which defaults to zero.
    In subsequent steps the increment is the difference between the receiver
    evaluated at the given solution step and at the previous one.

    .. note::

       Deprecated in the sources, but still read.

:optparam:`paramtype{in}`
    What the function's argument is.  ``0``, the default, evaluates the
    function at the current solution time.  ``1`` evaluates it at a parameter
    supplied by whatever uses the function, so that the abscissa is not time at
    all.

    This is what lets a function describe a dependence on something other than
    time — a gap-dependent conductivity, for instance; see
    :ref:`thermal-s2s-contact`.

.. _ConstantFunctionRecord:

Constant function
-----------------

A constant function.

.. record::

   :descitem:`ConstantFunction` :elemparam:`num{in}` :elemparam:`f(t){rn}`

.. _PeakFunctionRecord:

Peak function
-------------

A peak function: if the time equals :param:`t`, the value of the function is
:param:`f(t)`; otherwise zero is returned.

.. record::

   :descitem:`PeakFunction` :elemparam:`num{in}` :elemparam:`t{rn}`
   :elemparam:`f(t){rn}`

.. _PiecewiseLinFunctionRecord:

Piecewise linear function
-------------------------

A piecewise linear function.  The time values in :param:`t` must be sorted
along the time scale, and the corresponding function values are in
:param:`f(t)`.  The value for a time not present in :param:`t` is computed by
linear interpolation.

.. record::

   :descitem:`PiecewiseLinFunction` :elemparam:`num{in}`
   ((:elemparam:`t{ra}` :elemparam:`f(t){ra}`) |
   :elemparam:`datafile{s}`)
   :optelemparam:`npoints{in}` :optelemparam:`paramtype{in}`

**Parameters**

:param:`t{ra}`, :param:`f(t){ra}`
    The abscissa and ordinate arrays.

:param:`datafile{s}`
    Read the data from an external ASCII file instead.  The data should be
    stored in two columns, or simply interleaved; lines starting with ``#`` are
    skipped.  The file name must be enclosed in double quotes.  An empty name
    is permitted, for a function that is intentionally never used.

:optparam:`npoints{in}`
    Number of time-value pairs.

    .. note::

       Accepted but ignored: like every array attribute, :param:`t` and
       :param:`f(t)` carry their own size (see :ref:`attribute-types`), and
       that is what is used.  It appears in many existing input files, so it is
       still tolerated.

.. _StepFunctionRecord:

Step function
-------------

Like the piecewise linear function, but the value is held constant between the
given points instead of being interpolated.

.. record::

   :descitem:`StepFunction` :elemparam:`num{in}`
   ((:elemparam:`t{ra}` :elemparam:`f(t){ra}`) |
   :elemparam:`datafile{s}`)

.. _PeriodicPiecewiseLinFunctionRecord:

Periodic piecewise linear function
----------------------------------

A piecewise linear function repeated with a given period.

.. record::

   :descitem:`PeriodicPiecewiseLinFunction` :elemparam:`num{in}`
   :elemparam:`t{ra}` :elemparam:`f(t){ra}`
   :elemparam:`period{rn}` :optelemparam:`addtf{in}`

**Parameters**

:param:`period{rn}`
    Length of one period.  After it the function repeats.

:optparam:`addtf{in}`
    Number of another function whose value is added to this one, which allows a
    periodic variation to be superposed on a trend.

.. _HeavisideLTFRecord:

Heaviside-like function
-----------------------

Zero up to the time given by :param:`origin`; equal to :param:`value` for
times greater than :param:`origin`.

.. record::

   :descitem:`HeavisideLTF` :elemparam:`num{in}` :elemparam:`origin{rn}`
   :elemparam:`value{rn}`

.. _PulseFunctionRecord:

Pulse function
--------------

A rectangular pulse: equal to :param:`value` between :param:`tmin` and
:param:`tmax`, and zero outside that interval.

.. record::

   :descitem:`PulseFunction` :elemparam:`num{in}` :elemparam:`tmin{rn}`
   :elemparam:`tmax{rn}` :elemparam:`value{rn}`
   :optelemparam:`tstime{rn}`

**Parameters**

:param:`tmin{rn}`, :param:`tmax{rn}`
    Start and end of the pulse.

:param:`value{rn}`
    Value of the function inside the interval.

:optparam:`tstime{rn}`
    Time-step time used to resolve the pulse.

.. _UsrDefLTFRecord:

User-defined function
---------------------

A user-defined function given as an expression, which may depend on the
parameter ``t``, for which the actual time is substituted before the expression
is evaluated.  See :ref:`attribute-types` for the expression syntax.

The first and second time derivatives may be required, depending on the type of
analysis.

.. record::

   :descitem:`UsrDefLTF` :elemparam:`num{in}` :elemparam:`f(t){expr}`
   :optelemparam:`dfdt(t){expr}` :optelemparam:`d2fdt2(t){expr}`

.. _PythonExpressionRecord:

User-defined Python expression
------------------------------

.. note::

   Requires a build configured with ``USE_PYTHON_EXTENSION=ON``.

A user-defined function given as a Python expression or script.  The code may
depend on the parameter ``t``, for which the actual time is substituted, and on
the array ``x``, which holds the position.

The expression or script must assign to a variable named ``ret``, whose value is
returned as the value of the function.

.. important::

   Text between double quotes is not lower-cased by the parser (see
   :ref:`anatomy-of-an-input-file`), so Python identifiers keep their case.

.. record::

   :descitem:`PythonExpression` :elemparam:`num{in}`
   (:elemparam:`f{s}` | :elemparam:`fpath{s}`)
   (:optelemparam:`dfdt{s}` | :optelemparam:`dfdtpath{s}`)
   (:optelemparam:`d2fdt2{s}` | :optelemparam:`d2fdt2path{s}`)

**Parameters**

:param:`f{s}`
    A string containing the Python expression.

:param:`fpath{s}`
    Alternatively, a string containing the path to a Python script file
    defining the code for the function.

:optparam:`dfdt{s}`, :optparam:`dfdtpath{s}`
    Expression, or path to a script, evaluating the first derivative of the
    function.

:optparam:`d2fdt2{s}`, :optparam:`d2fdt2path{s}`
    Expression, or path to a script, evaluating the second derivative.

The derivatives may be required depending on the context in which the function
is used.

.. _InterpolatingFunctionRecord:

Interpolating function
----------------------

Interpolates values from a data file over a structured grid of one, two or
three dimensions.

.. record::

   :descitem:`InterpolatingFunction` :elemparam:`num{in}`
   :elemparam:`name{s}` :elemparam:`dim{in}`

**Parameters**

:param:`name{s}`
    Path to the file holding the data.

:param:`dim{in}`
    Number of dimensions of the data.

.. _LocalGaussRandomFunctionRecord:

Local Gaussian random function
------------------------------

Returns a normally distributed random value, used to introduce randomness into
material properties over a domain.

.. record::

   :descitem:`LocalGaussRandomFunction` :elemparam:`num{in}`
   :elemparam:`mean{rn}` :elemparam:`variance{rn}`
   :optelemparam:`seed{in}`

**Parameters**

:param:`mean{rn}`, :param:`variance{rn}`
    Mean and variance of the distribution.

:optparam:`seed{in}`
    Seed of the random generator, so that a run can be reproduced.
