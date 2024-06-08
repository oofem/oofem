XFEM module
================

An entity subject to XFEM modeling is described by an `EnrichmentItem`.
Such an entity can be e.g. a crack or a material interface. The
`EnrichmentItem` is responsible for the geometrical description as well as
the enrichment functions. The `EnrichmentItem` is also responsible for the
necessary updates if the interface is moving (e.g. if a crack is propagating).
The `EnrichmentItem` has different data members to fulfill these
requirements:

- The `EnrichmentItem`has an `EnrichmentDomain` that is
  responsible for the geometrical description. The `EnrichmentDomain` can be e.g. a
  polygon line or a list of node numbers.
- The `EnrichmentItem` has an `EnrichmentFunction` that
  describes the enrichment along the crack or interface.
- The `EnrichmentItem` has an `EnrichmentFront` that specifies
  a special treatment of the ends of the interface if desired. It may e.g. apply
  branch functions within a specified radius around the crack tips.
- The `EnrichmentItem` has a `PropagationLaw` that describes
  how the interface evolves in time. It may e.g. propagate a crack in a specified
  direction.



The enrichment functions are represented by classes derived from the pure
virtual `EnrichmentFunction` class. In particular, the following methods
have to be overloaded by all subclasses of `EnrichmentFunction`:

- `evaluateEnrFuncAt`: Computes the value of the enrichment
  function in a given point for a given level set value.
- `evaluateEnrFuncDerivAt`: Computes the gradient of the enrichment
  function in a given point for a given level set value.
- `giveJump`: Returns the jump of the enrichment function when the
  level set field changes sign, e.g. 1 for Heaviside, 2 for Sign and 0 for abs
  enrichment. This method is needed for a generic implementation of cohesive
  zones.


The enrichment items are stored in an array in the `XfemManager`
class, that also provides initialization and access to individual
enrichment items. The purpose of the `XfemManager` is to encapsulate XFEM
specific functionality. This avoids to a certain extent the need to modify and
penetrate the original code.

When a certain node in the problem domain is subjected to enrichment, the
corresponding enrichment item has to introduce additional degree(s) of freedom
in that node. Since several enrichment items can be active in a node, all DOFs
are assigned with a label, that allows to distinguish between them (in
traditional FE codes labels are also used to distinguish physical meaning of
particular DOF representing displacement along x,y,z axes, rotations, etc.).
The `XfemManager` class maintains a pool of available labels and is
responsible for their assignment to individual enrichment items.

Each enrichment item keeps its assigned DOF labels (provided by `XfemManager`).
These are used by individual enrichment items to introduce and later identify
related DOFs in individual nodes. The physical meaning of these DOFs depends on
the enrichment item. For example, a single DOF is needed to represent slip along
a slip line and two DOFs are required for a displacement discontinuity in 2D
modeled with Heaviside enrichment.


If the enrichment items evolve in time, the geometrical description needs to be
updated and new enriched DOFs may need to be created. This is done by the
`updateGeometry` method in the `XfemManager` class. The
sub-triangulation of cut elements may have to be updated as a result of the
geometrical update (e.g. previously uncut elements may be cut by a crack and
therefore need to be subdivided). The solver class `XFEMStatic` handles
this problem by creating a new domain with a new sub-triangulation and mapping
state variables as well as primary variables from the old domain to the new
domain.


The sub-triangulation of cut elements as well as computation of B-matrices for
2D XFEM elements are to a large extent encapsulated in the
`XfemElementInterface` class. Methods needed for cohesive zones are also
implemented here. The 2D elements with XFEM functionality are therefore derived
from the standard `Element` base class, as well as from the
`XfemElementInterface`. 

When a cut element has been subdivided by the `XfemElementInterface`,
a `PatchIntegrationRule` class is used to represent the integration scheme,
where Gauss integration is employed on individual triangles. When
initializing the integration points, their local coordinates (related to the
triangle geometry) are mapped into parent element coordinates and the
corresponding integration weights are scaled accordingly. There is no need to
store individual patch geometries, as the parent element geometry is used to
correctly evaluate transformation Jacobians (we may however wish to store the
them for debugging or visualization). 

If a cohesive crack model is employed, additional integration points along
the discontinuity are needed. These are stored separately, see
`XfemElementInterface` class. 

The following test cases involve the XFEM module:

- sm/xFemCrackVal.in
- sm/xFemCrackValBranch.in
- sm/xfemCohesiveZone1.in
- benchmark/xfem01.in
- sm/plasticRemap1.in
