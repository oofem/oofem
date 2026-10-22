// Notched three-point bending beam (2D), dimensions taken from
// tests/regression/benchmark/sm/concrete_3point_direct.in [mm]:
// span 2000, depth 500, notch 25 x 200 at midspan, supports at bottom corners,
// load (prescribed displacement) at the top midspan point.
// The material zone x in [900, 1100] around the notch ("Concrete") is modelled
// by a (damage) material, the rest ("Elastic") is elastic.
SetFactory("OpenCASCADE");
DefineConstant[ lc  = 100 ];   // initial mesh size, elastic zone
DefineConstant[ lcc = 25 ];    // initial mesh size, concrete zone

L = 2000; H = 500;             // span, depth
nw = 25; nd = 200;             // notch width and depth
x1 = 900; x2 = 1100;           // concrete zone

Rectangle(1) = {0,  0, 0, x1, H};
Rectangle(2) = {x1, 0, 0, L/2 - x1, H};
Rectangle(3) = {L/2, 0, 0, x2 - L/2, H};   // split at x = L/2 -> node at load point
Rectangle(4) = {x2, 0, 0, L - x2, H};
Rectangle(5) = {L/2 - nw/2, 0, 0, nw, nd}; // notch
BooleanDifference{ Surface{2, 3}; Delete; }{ Surface{5}; Delete; }
BooleanFragments{ Surface{:}; Delete; }{}

e = 1e-3;
elastic()  = Surface In BoundingBox{-e, -e, -e, x1 + e, H + e, e};
elastic() += Surface In BoundingBox{x2 - e, -e, -e, L + e, H + e, e};
concrete() = Surface In BoundingBox{x1 - e, -e, -e, x2 + e, H + e, e};

MeshSize{ PointsOf{ Surface{:}; } } = lc;
MeshSize{ PointsOf{ Surface{concrete()}; } } = lcc;

// Physical groups are exported as UNV groups and referenced from the ctrl files
Physical Surface("Elastic") = {elastic()};
Physical Surface("Concrete") = {concrete()};
Physical Point("SupportL") = Point In BoundingBox{-e, -e, -e, e, e, e};
Physical Point("SupportR") = Point In BoundingBox{L - e, -e, -e, L + e, e, e};
Physical Point("Load") = Point In BoundingBox{L/2 - e, H - e, -e, L/2 + e, H + e, e};
