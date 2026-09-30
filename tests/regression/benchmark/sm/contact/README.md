# Contact benchmark cases

Larger showcase examples for the structural penalty contact boundary
condition (`structuralpenaltycontactbc`), too big to run as part of the
regular test suite (see `tests/Guidelines.txt`). Each case is a single
runnable `.in` file with an embedded `#%BEGIN_CHECK%` block verifying a few
converged nodal displacements.

Open the matching `<name>.out.m0.pvd` file (generated on first run) in
ParaView to load the full VTU time series. Useful cell fields:

- `IST_ContactNormalGap`;
- `IST_ContactPressure`;
- `IST_ContactStatus` (`0` open, `1` active frictionless, `2` sticking,
  `3` sliding).

## Frictionless cases (active)

- `contact_test_01.in` — two touching cubes in frictionless contact,
  compressed in ten increments.
- `contact_test_02.in` — confined compression across a nonconforming
  two-facet/one-facet interface.
- `contact_test_03.in` — two 5x5-element blocks compressed against each
  other with two-pass contact (reversed master/slave penalty conditions).
- `contact_test_04.in` — two-pass "ironing" case: a stiff cylinder section
  presses down onto a softer block, then slides 3 length units along it.
  Uses automatic penalty selection and directional projection.
- `contact_test_05.in` — two-pass frictionless variant of the crossed
  cylindrical tubes case below (`contact_test_04_friction.in`), used to
  isolate the pure normal-contact search behavior from friction.

## Frictional cases (`friction_wip/`, not registered as tests yet)

Friction in `structuralpenaltycontactbc` is an experimental,
not-yet-fully-verified development feature (see the boundary condition's
own input-reference documentation and source comments). These cases are
kept for reference and future validation but deliberately live one
directory level below where CMake's test-discovery glob looks, so they are
not run automatically:

- `contact_test_01_friction.in` — smaller nonmatching cube compresses then
  slides with friction (mu=0.99).
- `contact_test_02_friction.in` — a frictional block presses onto a long
  plane then transitions from stick to sliding (mu=0.5).
- `contact_test_03_friction.in` — frictional variant of `contact_test_04.in`
  (two-pass ironing, mu=0.20).
- `contact_test_04_friction.in` — two crossed cylindrical tubes in frictional
  two-pass contact (mu=0.5). This case previously exhibited facet-ownership
  Newton chattering at a shared master-facet edge; it now completes all 50
  steps cleanly with `directionalprojection 1` (see the repository's
  contact-mechanics continuity notes for the full diagnosis and the not-yet-
  isolated fix). Kept here rather than in the active set purely because of
  the friction policy above, not because it is still failing.

`friction_wip/visual/` — simple 2D/3D Coulomb-friction cases meant for visual
inspection in ParaView rather than automated checking (stick/slip on an
inclined plane, a rotating two-brick stack) — no `#%BEGIN_CHECK%` block,
matching this project's existing convention for purely illustrative
benchmark cases.
