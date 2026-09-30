Other problem types
===================

.. _DummyEngngModel:

DummyEngngModel
---------------

A dummy model that cannot perform any analysis.  Its purpose is to invoke the
configured export modules, so that the problem geometry can be exported without
actually solving anything.

.. record::

   :descitem:`Dummy` :optelemparam:`nmodules{in}`

**Parameters**

:optparam:`nmodules{in}`
    Number of export module records that follow; see :ref:`ExportModulesSec`.
    Without at least one export module the record does nothing.

.. _other-engineering-models:

Other registered engineering models
-----------------------------------

The models below are registered by the solver and can be used in an input file,
but are not yet documented in this manual.  Their input attributes can be read
off the ``_IFT_`` definitions in the header next to the listed source file, and
off the model's ``initializeFrom`` method.

.. list-table::
   :header-rows: 1
   :widths: 34 26 40

   * - Record keyword
     - Module
     - Source
   * - ``nonlineardynamic``
     - structural
     - ``src/sm/EngineeringModels/nlineardynamic.C``
   * - ``linearstability``
     - structural
     - ``src/sm/EngineeringModels/linearstability.C``
   * - ``qclinearstatic``
     - structural (quasicontinuum)
     - ``src/sm/EngineeringModels/qclinearstatic.C``
   * - ``structuralmaterialevaluator``
     - structural
     - ``src/sm/EngineeringModels/structuralmaterialevaluator.C``
   * - ``stokesflow``
     - fluid
     - ``src/fm/stokesflow.C``
   * - ``stokesflowvelocityhomogenization``
     - fluid
     - ``src/fm/stokesflowvelocityhomogenization.C``
   * - ``fluidmaterialevaluator``
     - fluid
     - ``src/fm/fluidmaterialevaluator.C``
   * - ``dgproblem``
     - mpm (discontinuous Galerkin)
     - ``src/mpm/dg.C``
   * - ``problemsequence``
     - core
     - ``src/core/problemsequence.C``
   * - ``test``
     - mpm (prototype)
     - ``src/mpm/prototype2.C``

The two ``*materialevaluator`` models are not boundary-value problems: they
exercise a single material model over a prescribed strain or strain-rate
history, which is useful for testing and calibrating constitutive models.
