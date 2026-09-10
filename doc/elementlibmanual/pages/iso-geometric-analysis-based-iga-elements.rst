
Iso Geometric Analysis based (IGA) elements
===========================================

The following record describes the common part of IGA element record:

\recentry**IGAElement**:elemparam:`num{in}`  \recentry\ :elemparam:`mat{in}`
:elemparam:`crossSect{in}` :elemparam:`nodes{ia}`    \recentry\
:elemparam:`knotvectoru{ra}` :elemparam:`knotvectorv{ra}` :elemparam:`knotvectorw{ra}`
\recentry\ :optelemparam:`knotmultiplicityu{ia}` :optelemparam:`knotmultiplicityv{ia}`
\recentry\ :optelemparam:`knotmultiplicityw{ia}`  \recentry\ :elemparam:`degree{ia}`
:elemparam:`nip{ia}`    \recentry⟨:optelemparam:`partitions{ia}`\ ⟩
⟨:optelemparam:`remote`\ ⟩

The :param:`knotvectoru`, :param:`knotvectorv`, and :param:`knotvectorw` parameters
specify knot vectors in individual parametric directions, considering only distinct
knots. Open knot vector is always assumed, so the multiplicity of the first and last
knot should be equal to :math:`p+1`, where :math:`p` is polynomial degree in
coresponding direction (determined by :param:`degree` parameter, see further).  The knot
multiplicity can be set using optional parameters :param:`knotmultiplicityu`,
:param:`knotmultiplicityv`, and :param:`knotmultiplicityw`. By default, the open knot
vector is assumed and multiplicity of internal knots is assumed to be equal to one.
Note, that total number of knots in particular direction (including multiplicity) must
be equal to number of control points in this direction increased by degree in this
direction plus 1.  The degree of approximation for each parametric direction is
determined from :param:`degree` array, dimension of which is equal to number of spatial
dimensions of the problem.  In case of elements with BSpline or Nurbs interpolation, the
nodes forming the rectangular array of control points of the element are ordered in a
such way, that u-index is changing most quickly, and w-index (or v-index in case of 2d
problems) most slowly. In case of elements with T-spline interpolation, the nodes
forming the T-mesh of the element are ordered arbitrarily.

The supported **IGAElement** values are following:  **Keyword**:
:param:`bsplineplanestresselement`  **Parameters**: None.  **Keyword**:
:param:`nurbsplanestresselement`  **Parameters**: None.  **Keyword**:
:param:`nurbs3delement`  **Parameters**: None.  **Keyword**:
:param:`tsplineplanestresselement`  **Parameters**:
:elemparam:`localindexknotvectoru{in}` :elemparam:`localindexknotvectorv{in}`
:elemparam:`localindexknotvectorw{in}`  The parameters :param:`localindexknotvectoru`,
:param:`localindexknotvectorv`,  and :param:`localindexknotvectorw` defined by the
indices to global knot vectors (given by :param:`knotvectoru`, :param:`knotvectorv`, and
:param:`knotvectorw` parameters) specify the local knot vectors for each control point
of T-mesh (node) in the same order as the nodes have been specified for the element. The
local knot vector in a particular direction has :math:`p+2` entries, where the :math:`p`
is the polynomial degree in that direction.

Tests/Examples: `tests/regression/sm/ex-bspline-01.in
<https://github.com/oofem/oofem/blob/devel/tests/regression/sm/ex-bspline-01.in>`_,
`tests/regression/sm/ex-bspline-02.in
<https://github.com/oofem/oofem/blob/devel/tests/regression/sm/ex-bspline-02.in>`_,
`tests/regression/sm/ex3d-nurbs-01.in
<https://github.com/oofem/oofem/blob/devel/tests/regression/sm/ex3d-nurbs-01.in>`_,
`tests/regression/sm/ex3d-nurbs-02.in
<https://github.com/oofem/oofem/blob/devel/tests/regression/sm/ex3d-nurbs-02.in>`_.
