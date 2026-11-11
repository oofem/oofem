// L-shaped plate with a re-entrant corner at (1,1).
// The stress singularity at the corner drives the adaptive refinement.
DefineConstant[ lc = 0.25 ];   // initial (uniform) mesh size, can be changed by gmsh -setnumber lc value

Point(1) = {1.0, 0.0, 0, lc};
Point(2) = {1.5, 0.0, 0, lc};
Point(3) = {1.5, 1.5, 0, lc};
Point(4) = {0.0, 1.5, 0, lc};
Point(5) = {0.0, 1.0, 0, lc};
Point(6) = {1.0, 1.0, 0, lc};   // re-entrant corner

Line(1) = {1, 2};   // bottom edge  -> loaded
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};   // left edge    -> clamped
Line(5) = {5, 6};
Line(6) = {6, 1};

Curve Loop(1) = {1, 2, 3, 4, 5, 6};
Plane Surface(1) = {1};

// Physical groups are exported as UNV groups and referenced from the ctrl file
Physical Surface("Domain") = {1};
Physical Curve("Clamped") = {4};
Physical Curve("Load") = {1};
