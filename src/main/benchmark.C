#include <benchmark/benchmark.h>
#include "fei3dhexalin.h"
#include "fei3dhexaquad.h"
#include "fei2dtrquad.h"
#include "floatarrayf.h"
#include "floatmatrixf.h"

using namespace oofem;


static void EigDyn(benchmark::State& state) {
#if 0
    FloatMatrix D = {
        {1., 2., 3., 0., 0., 0.},
        {2., -4., 5., 0., 1., 10.},
        {3., 5., 6., 4., 0., 0.},
        {0., 0., 4., 10., 0., 0.},
        {0., 1., 0., 0., 20., 0.},
        {0., 10., 0., 0., 0., 10.},
    };
#elif 0
    FloatMatrix D = {
        {1, 3},
        {3, 5},
    };
#else
    FloatMatrix D = {
        {1, 2, 3},
        {2, 4, 5},
        {3, 5, 6},
    };
#endif
    FloatArray e(3);
    FloatMatrix v(3,3);
    for (auto _ : state) {
        D.jaco_(e, v, 10);
        benchmark::DoNotOptimize(e);
        benchmark::DoNotOptimize(v);
    }
}
BENCHMARK(EigDyn);


static void EigFix(benchmark::State& state) {
#if 0
    FloatMatrixF<6,6> D = {
        1., 2., 3., 0., 0., 0.,
        2., -4., 5., 0., 1., 10.,
        3., 5., 6., 4., 0., 0.,
        0., 0., 4., 10., 0., 0.,
        0., 1., 0., 0., 20., 0.,
        0., 10., 0., 0., 0., 10.,
    };
#elif 0
    FloatMatrixF<2,2> D = {
        1., 3.,
        3., 5.,
    };
#else
    FloatMatrixF<3,3> D = {
        1., 2., 3.,
        2., 4., 5.,
        3., 5., 6.,
    };
#endif
    for (auto _ : state) {
        auto tmp = eig(D, 10);
        benchmark::DoNotOptimize(tmp);
    }
}
BENCHMARK(EigFix);


static void SolDyn(benchmark::State& state) {
#if 1
    FloatMatrix D = {
        {1, 2, 3, 7},
        {2, -4, 5, 0},
        {3, 5, 6, 9},
        {7, 0, 9, 10},
    };
    FloatArray x = {1., -2., 3., -4.};
#elif 0
    FloatMatrix D = {
        {1, 3},
        {3, 5},
    };
#else
    FloatMatrix D = {
        {1, 2, 3},
        {2, -4, 5},
        {3, 5, 6},
    };
    FloatArray x = {1., -2., 3.};
#endif
    FloatArray s;
    for (auto _ : state) {
        FloatMatrix D2 = D;
        D2.solveForRhs(x, s);
        benchmark::DoNotOptimize(s);
    }
}
BENCHMARK(SolDyn);


static void SolFix(benchmark::State& state) {
#if 0
    FloatMatrixF<4,4> D = {
        1., 2., 3., 7.,
        2., -4., 5., 0.,
        3., 5., 6., 9.,
        7., 0., 9., 10.,
    };
    FloatArrayF<4> x = {1., -2., 3., -4.};
#elif 0
    FloatMatrixF<2,2> D = {
        1., 3.,
        3., 5.,
    };
#else
    FloatMatrixF<3,3> D = {
        1., 2., 3.,
        2., -4., 5.,
        3., 5., 6.,
    };
    FloatArrayF<3> x = {1., -2., 3.};
#endif
    for (auto _ : state) {
        //auto tmp = solve(D, x);
        auto tmp = dot(inv(D), x);        
        benchmark::DoNotOptimize(tmp);
        benchmark::DoNotOptimize(D);
        benchmark::DoNotOptimize(x);
    }
}
BENCHMARK(SolFix);


static void InvDyn(benchmark::State& state) {
#if 1
    FloatMatrix D = {
        {1, 2, 3, 7},
        {2, -4, 5, 0},
        {3, 5, 6, 9},
        {7, 0, 9, 10},
    };
#elif 0
    FloatMatrix D = {
        {1, 3},
        {3, 5},
    };
#else
    FloatMatrix D = {
        {1, 2, 3},
        {2, -4, 5},
        {3, 5, 6},
    };
#endif
    FloatMatrix v(3,3);
    for (auto _ : state) {
        v.beInverseOf(D);
        benchmark::DoNotOptimize(v);
    }
}
BENCHMARK(InvDyn);


static void InvFix(benchmark::State& state) {
#if 1
    FloatMatrixF<4,4> D = {
        1., 2., 3., 7.,
        2., -4., 5., 0.,
        3., 5., 6., 9.,
        7., 0., 9., 10.,
    };
#elif 0
    FloatMatrixF<2,2> D = {
        1., 3.,
        3., 5.,
    };
#else
    FloatMatrixF<3,3> D = {
        1., 2., 3.,
        2., -4., 5.,
        3., 5., 6.,
    };
#endif
    for (auto _ : state) {
        auto tmp = inv(D);
        benchmark::DoNotOptimize(tmp);
    }
}
BENCHMARK(InvFix);



static void CopyD(benchmark::State& state) {
    FloatMatrix D(3,3);
    double E = 210;
    double nu = 0.3;
    double G = (E / ( 2.0 * ( 1. + nu ) ));
    //double K = E / ( 3.0 * ( 1. - 2. * nu ) );
    double ee = E / ( 1. - nu * nu );
    FloatMatrixF<3,3> tangent = {
        ee, nu*ee, 0.,
        nu*ee, ee, 0.,
        0., 0., G
    };

    for (auto _ : state) {
        auto xxx = tangent;
        benchmark::DoNotOptimize(xxx);
    }
}
BENCHMARK(CopyD);

static void ComputeD(benchmark::State& state) {
    //FloatMatrix D(3,3);
    double E = 210;
    double nu = 0.3;
    double G = (E / ( 2.0 * ( 1. + nu ) ));
    for (auto _ : state) {
        double e = E;
        double ee = e / ( 1. - nu * nu );
        double shear = G;
        FloatMatrixF<3,3> D;
        D.at(1, 1) = ee;
        D.at(1, 2) = nu * ee;
        D.at(2, 1) = nu * ee;
        D.at(2, 2) = ee;
        D.at(3, 3) = shear;
        benchmark::DoNotOptimize(D);
    }
}
BENCHMARK(ComputeD);

#if 0
const std::vector<FloatArrayF<3>> nodes_6 = {
    FloatArrayF<3>{0.,0.,0.},FloatArrayF<3>{2.,0.,0.},FloatArrayF<3>{0.,1.,0.},
    FloatArrayF<3>{1.,0.,0.},FloatArrayF<3>{1.,0.6,0.},FloatArrayF<3>{0.,0.5,0.},
};

const std::vector<FloatArrayF<3>> nodes_8 = {
    FloatArrayF<3>{-1.,-1., 1.},FloatArrayF<3>{-1.,1., 1.},FloatArrayF<3>{1.,1., 1.},FloatArrayF<3>{1.,-1., 1.},
    FloatArrayF<3>{-1.,-1.,-1.},FloatArrayF<3>{-1.,1.,-1.},FloatArrayF<3>{1.,1.,-1.},FloatArrayF<3>{1.,-1.,-1.},
};

const std::vector<FloatArrayF<3>> nodes_20 = {
    FloatArrayF<3>{-1.,-1., 1.},FloatArrayF<3>{-1.,1., 1.},FloatArrayF<3>{1.,1., 1.},FloatArrayF<3>{1.,-1., 1.},
    FloatArrayF<3>{-1.,-1.,-1.},FloatArrayF<3>{-1.,1.,-1.},FloatArrayF<3>{1.,1.,-1.},FloatArrayF<3>{1.,-1.,-1.},
    FloatArrayF<3>{-1., 0., 1.},FloatArrayF<3>{ 0.,1., 1.},FloatArrayF<3>{1.,0., 1.},FloatArrayF<3>{0.,-1., 1.},
    FloatArrayF<3>{-1., 0.,-1.},FloatArrayF<3>{ 0.,1.,-1.},FloatArrayF<3>{1.,0.,-1.},FloatArrayF<3>{0.,-1.,-1.},
    FloatArrayF<3>{-1.,-1., 0.},FloatArrayF<3>{-1.,1., 0.},FloatArrayF<3>{1.,1., 0.},FloatArrayF<3>{1.,-1., 0.}
};
#else
const std::vector<FloatArray> nodes_6 = {
    FloatArray{0.,0.,0.},FloatArray{2.,0.,0.},FloatArray{0.,1.,0.},
    FloatArray{1.,0.,0.},FloatArray{1.,0.6,0.},FloatArray{0.,0.5,0.},
};

const std::vector<FloatArray> nodes_8 = {
    FloatArray{-1.,-1., 1.},FloatArray{-1.,1., 1.},FloatArray{1.,1., 1.},FloatArray{1.,-1., 1.},
    FloatArray{-1.,-1.,-1.},FloatArray{-1.,1.,-1.},FloatArray{1.,1.,-1.},FloatArray{1.,-1.,-1.},
};

const std::vector<FloatArray> nodes_20 = {
    FloatArray{-1.,-1., 1.},FloatArray{-1.,1., 1.},FloatArray{1.,1., 1.},FloatArray{1.,-1., 1.},
    FloatArray{-1.,-1.,-1.},FloatArray{-1.,1.,-1.},FloatArray{1.,1.,-1.},FloatArray{1.,-1.,-1.},
    FloatArray{-1., 0., 1.},FloatArray{ 0.,1., 1.},FloatArray{1.,0., 1.},FloatArray{0.,-1., 1.},
    FloatArray{-1., 0.,-1.},FloatArray{ 0.,1.,-1.},FloatArray{1.,0.,-1.},FloatArray{0.,-1.,-1.},
    FloatArray{-1.,-1., 0.},FloatArray{-1.,1., 0.},FloatArray{1.,1., 0.},FloatArray{1.,-1., 0.}
};
#endif

const FEIVertexListGeometryWrapper tri_6 = nodes_6;
const FEIVertexListGeometryWrapper cube_8 = nodes_8;
const FEIVertexListGeometryWrapper cube_20 = nodes_20;
const FEIVoidCellGeometry void_cell;


static void BMatrixLin(benchmark::State& state) {
    FEI3dHexaLin interp;
    FloatArray lcoords = {0.2, 0.4, 0.3};
    for (auto _ : state) {
        FloatMatrix B;
        interp.evaldNdx(B, lcoords, cube_8);
        benchmark::DoNotOptimize(B);
    }
}
BENCHMARK(BMatrixLin);

static void BMatrixLinFixed(benchmark::State& state) {
    FEI3dHexaLin interp;
    FloatArrayF<3> lcoords = {0.2, 0.4, 0.3};
    for (auto _ : state) {
        auto x = interp.evaldNdx(lcoords, cube_8);
        benchmark::DoNotOptimize(x);
    }
}
BENCHMARK(BMatrixLinFixed);


static void HexLinN(benchmark::State& state) {
    FEI3dHexaLin interp;
    FloatArray lcoords = {0.2, 0.4, 0.3};
    FloatArray N(8);
    for (auto _ : state) {
        interp.evalN(N, lcoords, void_cell);
        benchmark::DoNotOptimize(N);
    }
}
BENCHMARK(HexLinN);

static void HexLinNFixed(benchmark::State& state) {
    FEI3dHexaLin interp;
    FloatArrayF<3> lcoords = {0.2, 0.4, 0.3};
    for (auto _ : state) {
        auto x = interp.evalN(lcoords);
        benchmark::DoNotOptimize(x);
    }
}
BENCHMARK(HexLinNFixed);


static void HexQuadB(benchmark::State& state) {
    FEI3dHexaQuad interp;
    FloatArray lcoords = {0.2, 0.4, 0.3};
    for (auto _ : state) {
        FloatMatrix B;
        interp.evaldNdx(B, lcoords, cube_20);
        benchmark::DoNotOptimize(B);
    }
}
BENCHMARK(HexQuadB);

static void HexQuadBFixed(benchmark::State& state) {
    FEI3dHexaQuad interp;
    FloatArrayF<3> lcoords = {0.2, 0.4, 0.3};
    for (auto _ : state) {
        auto x = interp.evaldNdx(lcoords, cube_20);
        benchmark::DoNotOptimize(x);
    }
}
BENCHMARK(HexQuadBFixed);


static void HexQuadN(benchmark::State& state) {
    FEI3dHexaQuad interp;
    FloatArray lcoords = {0.2, 0.4, 0.3};
    FloatArray N(20);
    for (auto _ : state) {
        interp.evalN(N, lcoords, void_cell);
        benchmark::DoNotOptimize(N);
    }
}
BENCHMARK(HexQuadN);

static void HexQuadNFixed(benchmark::State& state) {
    FEI3dHexaQuad interp;
    FloatArrayF<3> lcoords = {0.2, 0.4, 0.3};
    for (auto _ : state) {
        auto N = interp.evalN(lcoords);
        benchmark::DoNotOptimize(N);
    }
}
BENCHMARK(HexQuadNFixed);



static void TriQuadB(benchmark::State& state) {
    FEI2dTrQuad interp(1,2);
    FloatArray lcoords = {0.2, 0.4};
    for (auto _ : state) {
        FloatMatrix B;
        interp.evaldNdx(B, lcoords, tri_6);
        benchmark::DoNotOptimize(B);
    }
}
BENCHMARK(TriQuadB);

static void TriQuadBFixed(benchmark::State& state) {
    FEI2dTrQuad interp(1,2);
    FloatArrayF<2> lcoords = {0.2, 0.4};
    for (auto _ : state) {
        auto x = interp.evaldNdx(lcoords, tri_6);
        benchmark::DoNotOptimize(x);
    }
}
BENCHMARK(TriQuadBFixed);


static void TriQuadN(benchmark::State& state) {
    FEI2dTrQuad interp(1,2);
    FloatArray lcoords = {0.2, 0.4};
    for (auto _ : state) {
        FloatArray N;
        interp.evalN(N, lcoords, void_cell);
        benchmark::DoNotOptimize(N);
    }
}
BENCHMARK(TriQuadN);

static void TriQuadNFixed(benchmark::State& state) {
    FEI2dTrQuad interp(1,2);
    FloatArrayF<2> lcoords = {0.2, 0.4};
    for (auto _ : state) {
        auto x = interp.evalN(lcoords);
        benchmark::DoNotOptimize(x);
    }
}
BENCHMARK(TriQuadNFixed);


BENCHMARK_MAIN();
