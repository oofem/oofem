// Notched three-point bending beam (3D), dimensions taken from
// tests/regression/benchmark/sm/concrete_3point_direct.in [mm], width 200 mm:
// span 2000, depth 500, notch 25 x 200 (through the width) at midspan,
// line supports along the bottom edges at x = 0 and x = 2000,
// load (prescribed displacement) along the top line at midspan.
// The material zone x in [900, 1100] around the notch ("Concrete") is modelled
// by a (damage) material, the rest ("Elastic") is elastic.
SetFactory("OpenCASCADE");
DefineConstant[ lc  = 100 ];   // initial mesh size, elastic zone
DefineConstant[ lcc = 50 ];    // initial mesh size, concrete zone

DefineConstant[ W = 200 ];      // width
L = 2000; H = 500;             // span, depth
nw = 25; nd = 200;             // notch width and depth
x1 = 900; x2 = 1100;           // concrete zone

Box(1) = {0,  0, 0, x1, H, W};
Box(2) = {x1, 0, 0, L/2 - x1, H, W};
Box(3) = {L/2, 0, 0, x2 - L/2, H, W};      // split at x = L/2 -> load line
Box(4) = {x2, 0, 0, L - x2, H, W};
Box(5) = {L/2 - nw/2, 0, 0, nw, nd, W};    // notch
BooleanDifference{ Volume{2, 3}; Delete; }{ Volume{5}; Delete; }
BooleanFragments{ Volume{:}; Delete; }{}

e = 1e-3;
elastic()  = Volume In BoundingBox{-e, -e, -e, x1 + e, H + e, W + e};
elastic() += Volume In BoundingBox{x2 - e, -e, -e, L + e, H + e, W + e};
concrete() = Volume In BoundingBox{x1 - e, -e, -e, x2 + e, H + e, W + e};

MeshSize{ PointsOf{ Volume{:}; } } = lc;
MeshSize{ PointsOf{ Volume{concrete()}; } } = lcc;

// Physical groups are exported as UNV groups and referenced from the ctrl files
Physical Volume("Elastic") = {elastic()};
Physical Volume("Concrete") = {concrete()};
Physical Curve("SupportL") = Curve In BoundingBox{-e, -e, -e, e, e, W + e};
Physical Curve("SupportR") = Curve In BoundingBox{L - e, -e, -e, L + e, e, W + e};
Physical Curve("Load") = Curve In BoundingBox{L/2 - e, H - e, -e, L/2 + e, H + e, W + e};
