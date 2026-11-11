# Adaptive analyses with Gmsh remeshing

These examples demonstrate the adaptive loops of `adaptlinearstatic` (linear analysis,
Zienkiewicz-Zhu error estimator) and `adaptnlinearstatic` (nonlinear analysis, damage
indicator with state mapping), both using the Gmsh mesher interface (`meshpackage 4`).

At every adaptive step:

1. the error is estimated and the required element size is computed;
2. the required size is written as a Gmsh background mesh (`<output>.bgm.<n>.pos`);
3. Gmsh remeshes the geometry, controlled by the background mesh, and the new mesh is
   converted to an oofem domain with `unv2oofem` (`<output>.domain.1.<n>.din`);
4. the analysis continues on the new mesh, as a new solution step.

   * **Linear analysis:** the problem is reanalyzed on the new mesh, and the loop stops
     when the error is acceptable or after `maxadaptsteps` remeshings.
   * **Nonlinear analysis:** the solution (primary unknowns and internal variables) is
     mapped to the new mesh and the loading continues.

Material, cross sections, boundary conditions, loads and sets are assigned by the names
of the Gmsh physical groups, as described in the unv2oofem ctrl file.

| case    | geometry       | ctrl file           | analysis               | elements          |
|---------|----------------|---------------------|------------------------|-------------------|
| `2d`    | `beam.geo`     | `beam.ctrl`         | linear, ZZ             | `TrPlaneStress2d` |
| `tet`   | `beam3d.geo`   | `beam3d_tet.ctrl`   | linear, ZZ             | `LTRSpace`        |
| `hex`   | `beam3d.geo`   | `beam3d_hex.ctrl`   | linear, ZZ             | `LSpace`          |
| `nl2d`  | `beam.geo`     | `beam_nl.ctrl`      | nonlocal damage        | `TrPlaneStress2d` |
| `nl3d`  | `beam3d.geo`   | `beam3d_nl.ctrl`    | nonlocal damage        | `LTRSpace`        |
| `ls2d`  | `lshape.geo`   | `lshape.ctrl`       | linear, ZZ             | `TrPlaneStress2d` |
| `lstet` | `lshape3d.geo` | `lshape3d_tet.ctrl` | linear, ZZ             | `LTRSpace`        |
| `lshex` | `lshape3d.geo` | `lshape3d_hex.ctrl` | linear, ZZ             | `LSpace`          |

## Geometries

**Notched three-point bending beam** (`beam.geo`, `beam3d.geo`)

- Dimensions and material constants are taken from the benchmark
  `tests/regression/benchmark/sm/concrete_3point_direct.in` (units mm, MPa):
  - span 2000, depth 500;
  - a notch 25 × 200 at midspan;
  - the beam rests on the bottom corners;
  - a vertical displacement is prescribed at the top midspan point.
- The 3D beam is 200 mm wide. Supports and load act along lines across the width.
- The zone x ∈ [900, 1100] around the notch (group `Concrete`) uses the damage material
  in the nonlinear analyses. The rest (group `Elastic`) is elastic.

**L-shaped body** (`lshape.geo`, `lshape3d.geo`)

- The left edge/face is clamped, and a traction acts on the bottom edge/face.
- The stress singularity at the re-entrant corner (2D) or edge (3D) drives the
  refinement.

Both geometries declare the initial mesh size with `DefineConstant`, so it can be
changed with `gmsh2oofem.py --setnumber lc <value>` (the beam also has `lcc` for the
concrete zone).

## Running

You need `gmsh` and `python3` on `PATH`, plus an oofem executable:

    OOFEM=/path/to/oofem ./run.sh 2d     # or tet, hex, nl2d, nl3d, ls2d, lstet, lshex

`run.sh` generates the initial mesh from the geometry with
`tools/unv2oofem/gmsh2oofem.py` and then runs oofem. A release build is recommended
for the 3D cases. Results go to `*.out`, and VTK files go to `*.vtu`/`*.pvd`, one step
per mesh (linear) or per load step (nonlinear).

## Adaptive linear analysis

The ZZ remeshing criterion asks for a new mesh when the global relative error exceeds
`requirederror`, or when the required element size is substantially (< 0.8×) smaller
than the current one. `minelemsize` limits the minimum element size.

Singularities produce a floor in the global error, because elements there cannot be
refined below the minimum size:
- in the beam: the point supports, the point load and the notch corners;
- in the L-shape: the re-entrant corner or edge.

Elements already at the minimum size do not trigger remeshing, so `requirederror` has to
be set above this floor.

Typical results (Gmsh 4.15, oofem release build). Reactions are in equilibrium with the
load on every mesh.

| case  | step 0 (initial mesh) | step 1             | last step                                 | time    |
|-------|-----------------------|--------------------|-------------------------------------------|---------|
| 2d    | 1348 elems, 19.7%     | 3414 elems, 15.9%  | step 2: 4880 elems, 14.9% (accepted)      | ~2 s    |
| tet   | 4209 elems, 29.9%     | 46804 elems, 19.8% | step 2: 100352 elems, 18.0% (accepted)    | ~1 min  |
| hex   | 3776 elems, 22.8%     | 14172 elems, 19.0% | step 3: 16760 elems, 18.6% (maxadaptsteps)| ~2.5 min|
| ls2d  | 48 elems, 43.5%       | 4531 elems, 4.4%   | step 2: 5011 elems, 3.7% (accepted)       | ~1 s    |
| lstet | 287 elems, 45.6%      | 24368 elems, 12.5% | step 2: 43936 elems, 10.8% (accepted)     | ~10 s   |
| lshex | 432 elems, 33.1%      | 20704 elems, 7.9%  | step 3: 14132 elems, 8.7% (maxadaptsteps) | ~3 min  |

Gmsh meshing is not strictly deterministic, so later steps may differ slightly between
runs.

## Adaptive nonlinear analysis

The cases `nl2d` and `nl3d` model the notched beam under displacement control. In the
concrete zone the material is the nonlocal isotropic damage model `idmnl1`, set up the
same way in 2D and 3D:
- E = 20 GPa, ν = 0.2;
- Rankine equivalent strain (`equivstraintype 4`) with e0 = ft / E = 1.25e-4
  (ft = 2.5 MPa);
- exponential softening with ef = 1e-3, which gives Gf ≈ 0.1 N/mm for the 40 mm nonlocal
  radius (bell function, `r 40 wft 1`);
- `regionmap 2 1 0` keeps the elastic zone out of the nonlocal averaging.

The nonlocal radius has to be resolved by the mesh in the damaged zone, with elements of
about r/3 or smaller. With r = 20 mm and the 3D mesh sizes used here (10–15 mm), the
averaging could not smooth the stress concentration under the load line. Spurious damage
then appeared at the top surface before the crack arrived. The larger radius avoids this
at an affordable 3D mesh size.

A crack starts at the notch tip and propagates upwards. At the end of each load
increment:

1. The damage-based error indicator gives the required element size at every node
   (`eetype 0 vartype 1`, ScalarErrorIndicator, with DirectErrorIndicatorRC).
   - Nodes with damage below `minlim` get `defdens`.
   - Damaged nodes get sizes interpolated between `mindens` (at `minlim`) and `maxdens`
     (at `maxlim`).
   - The sizes resolve the nonlocal radius: h = 5–12 mm in 2D and 10–15 mm in 3D (where
     `mindens` < `maxdens` refines most at the onset of damage, ahead of the crack).
2. When the required size drops below 0.8 of the current size anywhere, Gmsh creates a
   new mesh from the background size field.
3. The solution is mapped to the new mesh, both primary unknowns and damage state (`idm1`,
   `idmnl1` and `mdm` support state mapping).
4. With `equilmc 1`, the mapped configuration is brought back into equilibrium, and the
   analysis continues with the next increment.

Typical results (release build):

| case | initial mesh   | remeshings     | final mesh   | max. reaction   | time    |
|------|----------------|----------------|--------------|-----------------|---------|
| nl2d | 1348 triangles | 18 (25 steps)  | ~5900 elems  | 102.0           | ~18 s   |
| nl3d | 4209 tets      | 6 (25 steps)   | ~56000 tets  | 20100 (100/mm)  | ~6 min  |

- Both cases trace the full softening branch with no convergence problems, and the peak
  loads per unit width agree within 2% (102.0 in 2D, 100.5 per mm in 3D).
- The reaction changes by a few percent when the state is mapped to a new mesh.
- In 2D, the damage band stays above the notch tip (x ≈ 925–1075 mm) and does not reach
  the top surface.
- In 3D, every mapping of the state to a new mesh spreads damage by about one element
  size. Without grading control, the damage front entered the coarse mesh at almost every
  step: 20 remeshes, a band widening over the whole concrete zone, and damage reaching
  the top surface. `gmshgrading 0.1` limits how fast the required size may grow with
  distance (h_i <= h_j + 0.1 |x_i - x_j|). The fine zone then extends ahead of the crack:
  - 6 remeshes instead of 20;
  - the band (damage > 0.5) stays at x ≈ 966–1054;
  - damage near the top surface drops by about 70% and appears only late in the
    softening branch.

The model and parameters differ from the benchmark (nonlocal damage instead of
`Concrete3`), so the reactions are not expected to match it.

## Notes

* **Element types.** Gmsh UNV element types used in the ctrl files:

  | UNV type | element                 |
  |----------|-------------------------|
  | 91       | linear triangle         |
  | 94       | linear quadrilateral    |
  | 111      | linear tetrahedron      |
  | 115      | linear hexahedron       |
  | 21       | linear edge             |

  Boundary groups with an `etype[...]` line but no oofem element define
  `elementedges`/`elementboundaries` sets. Physical points and curves that are used only
  as node sets need no `etype` line.
* **Hexahedral meshes** are produced by subdividing tetrahedra (`gmshopts "--elemtype
  hex"`, which sets Gmsh `Mesh.SubdivisionAlgorithm=2`). Quad meshes in 2D are made the
  same way with `--elemtype quad`.
* **Mesh grading.** `gmshgrading g` (Gmsh interface) limits the gradation of the required
  size field before it is passed to Gmsh, so that the transition from fine to coarse mesh
  is gradual. Gmsh itself has no such control for background meshes. Values of 0.1–0.25
  are typical; smaller values give a smoother but larger mesh.
* **Linear solver.** All cases use the skyline solver with profile optimization
  (`profileopt 1`). On `nl3d`, the IML conjugate gradient solver was slower than this
  setup:
  - `lstype 1`, `stype 0`, with ILU or incomplete Cholesky;
  - tested at `lstol` 1e-4 to 1e-5.
* **Other mesher interfaces.** The same adaptive loops work with them. For example, T3d
  with an external remesh command:

      meshpackage 0 adapt 1 remeshcmd "t3d -i job.t3d -m %d -o job.t3d.out && t3d2oofem job.ctrl job.t3d.out job.in && tail -n +4 job.in > %o"

  The built-in subdivision (`meshpackage 3`) creates the new mesh in memory.
