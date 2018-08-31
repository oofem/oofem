#include <benchmark/benchmark.h>
#include "fei3dhexalin.h"
#include "fei3dhexaquad.h"
#include "fei2dtrquad.h"
#include "floatarrayf.h"
#include "floatmatrixf.h"

using namespace oofem;

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
    FloatArray N(6);
    for (auto _ : state) {
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
