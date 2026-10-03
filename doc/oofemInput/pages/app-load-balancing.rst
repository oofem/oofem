.. _dynamicloadbalancing:

Dynamic load balancing parameters
=================================

.. note::

   Parallel build only.  Load balancing support requires the ParMETIS module
   to be configured and compiled.

Two factors generally cause a load imbalance between the individual subdomains:

#. factors coming from the nature of the application, such as a switch from
   linear to nonlinear response in certain regions, or local adaptive
   refinement;
#. external factors caused by resource reallocation, typical of non-dedicated
   cluster environments where the individual processors are shared between
   different applications and users, which makes the allocated processing power
   vary in time.

Load balance is recovered by repartitioning the problem domain and transferring
the work — typically represented by finite elements — from one subdomain to
another.

The corresponding part of the analysis record has the following syntax.

.. record:: Syntax

   :optelemparam:`lbflag{in}` :optelemparam:`forcelb1{in}`
   :optelemparam:`wtp{ia}` :optelemparam:`lbstep{in}`
   :optelemparam:`relwct{rn}` :optelemparam:`abswct{rn}`
   :optelemparam:`minwct{rn}`

**Parameters**

:optparam:`lbflag{in}`
    A nonzero value activates dynamic load balancing.  Default ``0``.

:optparam:`forcelb1{in}`
    A nonzero value forces load rebalancing after the first solution step.

:optparam:`wtp{ia}`
    Activates optional load balancing plugins.  At present the only supported
    value is ``1``, which activates the nonlocal plugin, necessary for nonlocal
    averaging to work properly while dynamic load balancing is active.

:optparam:`lbstep{in}`
    Rebalancing, if needed, is performed only every :param:`lbstep`-th
    solution step.  Default ``5``.

:optparam:`relwct{rn}`
    Relative wall-clock imbalance threshold.  When the relative imbalance
    between the wall-clock solution times of the individual processors exceeds
    it, the rebalancing procedure is activated.  Default ``0.1``.

:optparam:`abswct{rn}`
    Absolute wall-clock imbalance threshold.  When the absolute imbalance
    between the wall-clock solution times of the individual processors exceeds
    it, the rebalancing procedure is activated.  Default ``10.0``.

:optparam:`minwct{rn}`
    Minimum absolute imbalance for the relative check against :param:`relwct`
    to be performed at all; below it, only the absolute check is done.
    Default ``0.0``.
