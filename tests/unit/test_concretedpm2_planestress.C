#include "sm/Materials/ConcreteMaterials/concretedpm2.h"
#include "sm/Materials/isolinearelasticmaterial.h"
#include "material.h"
#include "dynamicinputrecord.h"
#include "floatarrayf.h"
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <memory>
#include <utility>
#include <vector>

using namespace oofem;

namespace {

// Constructs a ConcreteDPM2 with plausible concrete parameters. No Domain is
// needed: none of the plane-stress invariant/derivative helpers under test
// touch the domain, and FEMComponent's constructor only stores the pointer.
std::unique_ptr<ConcreteDPM2>
makeConcreteDPM2()
{
    auto ir = std::make_shared<DynamicInputRecord>();
    ir->setField(2400.e-9, _IFT_Material_density);
    ir->setField(30000., _IFT_IsotropicLinearElasticMaterial_e);
    ir->setField(0.2, _IFT_IsotropicLinearElasticMaterial_n);
    ir->setField(0., _IFT_IsotropicLinearElasticMaterial_talpha);
    ir->setField(30., _IFT_ConcreteDPM2_fc);
    ir->setField(3., _IFT_ConcreteDPM2_ft);
    ir->setField(5.e-5, _IFT_ConcreteDPM2_wf);

    auto mat = std::make_unique<ConcreteDPM2>(1, nullptr);
    mat->initializeFrom(ir);
    return mat;
}

bool
closeEnough(double a, double b, double relTol = 1.e-3, double absTol = 1.e-6)
{
    return std::abs(a - b) <= absTol + relTol * std::max(std::abs(a), std::abs(b));
}

} // namespace

TEST_CASE( "tests/unit/concretedpm2_planestress", "[ConcreteDPM2]" ) {
    auto mat = makeConcreteDPM2();

    const double h = 1.e-3; // MPa, principal-stress perturbation for stress-direction FDs
    const double hk = 1.e-4; // perturbation for kappa-direction FDs

    // Generic (non-symmetric) states, chosen to stay away from the Lode angle
    // boundaries theta = 0 (uniaxial tension / equibiaxial compression) and
    // theta = pi/3 (uniaxial compression / equibiaxial tension), where the
    // derivatives of k are documented as non-unique.
    std::vector< std::pair< double, double > > genericStates = {
        { -5.0, 2.0 }, { 1.0, -7.0 }, { 4.0, 1.5 }, { -3.0, -8.0 }
    };

    // Boundary / symmetric states. Only used to check quantities that are
    // proven smooth there (the yield function f and the invariants themselves);
    // see algoCDPM.tex, "Non-uniqueness of the Jacobi matrix".
    std::vector< std::pair< double, double > > boundaryStates = {
        { 3.0, 0.0 },      // uniaxial tension,        theta = 0
        { 0.0, -3.0 },     // uniaxial compression,    theta = pi/3
        { 2.0, 2.0 },      // equibiaxial tension,     theta = pi/3
        { -10.0, -10.0 }   // equibiaxial compression, theta = 0
    };

    std::vector< double > kappas = { 0.0, 0.3, 0.6, 1.0, 1.4 };
    // Kappa values used only when kappa itself is the finite-difference
    // direction: avoid landing exactly on the hardening-law kinks at kappa = 0
    // (computeHardeningOne is constant for kappa <= 0, cubic for 0 < kappa < 1)
    // and kappa = 1 (see algoCDPM.tex / computeHardeningOnePrime).
    std::vector< double > kappasForKappaFD = { 0.3, 0.6, 0.8, 1.4, 2.0 };

    SECTION( "computeInvariants_PlaneStress is internally consistent" ) {
        auto allStates = genericStates;
        allStates.insert(allStates.end(), boundaryStates.begin(), boundaryStates.end() );
        for ( auto &st : allStates ) {
            auto inv = mat->computeInvariants_PlaneStress(st.first, st.second);
            CHECK( inv.rho >= 0. );
            CHECK( closeEnough(inv.sig, ( st.first + st.second ) / 3., 1.e-9, 1.e-9) );
            CHECK( inv.theta >= -1.e-9 );
            CHECK( inv.theta <= M_PI / 3. + 1.e-9 );
            CHECK( std::isfinite(inv.dRho[0]) );
            CHECK( std::isfinite(inv.dRho[1]) );
            CHECK( std::isfinite(inv.dTheta[0]) );
            CHECK( std::isfinite(inv.dTheta[1]) );
        }
    }

    SECTION( "Lode angle derivative uses the sgn(0) = 0 convention at degenerate states" ) {
        // The hand-off note reports that sgn(0) = +-1 fails to converge for an
        // equibiaxial-tension trial state; sgn(0) = 0 is used instead. Check
        // that the degenerate cases (s1 = s2, or one of them exactly zero)
        // produce finite, vanishing one-sided-average derivatives rather than
        // propagating a +-1 sign jump.
        for ( double s : { 3.0, -5.0, 10.0 } ) {
            auto inv = mat->computeInvariants_PlaneStress(s, s);
            CHECK( std::isfinite(inv.dTheta[0]) );
            CHECK( std::isfinite(inv.dTheta[1]) );
            CHECK( std::abs(inv.dTheta[0]) < 1.e-9 );
            CHECK( std::abs(inv.dTheta[1]) < 1.e-9 );
        }
        auto inv = mat->computeInvariants_PlaneStress(0.0, -4.0);
        CHECK( std::isfinite(inv.dTheta[0]) );
        CHECK( std::isfinite(inv.dTheta[1]) );
    }

    SECTION( "f1, f2 match finite differences of the yield function, including at Lode angle boundaries" ) {
        auto F = [&](double kappa, double s1, double s2) {
                     auto i = mat->computeInvariants_PlaneStress(s1, s2);
                     return mat->computeYieldValue(i.sig, i.rho, i.theta, kappa);
                 };

        auto allStates = genericStates;
        allStates.insert(allStates.end(), boundaryStates.begin(), boundaryStates.end() );

        for ( double kappa : kappas ) {
            for ( auto &st : allStates ) {
                double s1 = st.first, s2 = st.second;
                auto inv = mat->computeInvariants_PlaneStress(s1, s2);
                auto f = mat->computeDFDPrincipal_PlaneStress(inv, kappa);

                double df1 = ( F(kappa, s1 + h, s2) - F(kappa, s1 - h, s2) ) / ( 2. * h );
                double df2 = ( F(kappa, s1, s2 + h) - F(kappa, s1, s2 - h) ) / ( 2. * h );

                CHECK( closeEnough(f[0], df1) );
                CHECK( closeEnough(f[1], df2) );
            }
        }
    }

    SECTION( "k1, k2 match finite differences of k = dKappa/dDeltaLambda, away from Lode angle boundaries" ) {
        auto K = [&](double kappa, double s1, double s2) {
                     auto i = mat->computeInvariants_PlaneStress(s1, s2);
                     return mat->computeDKappaDDeltaLambda(i.sig, i.rho, i.theta, kappa);
                 };

        for ( double kappa : kappas ) {
            for ( auto &st : genericStates ) {
                double s1 = st.first, s2 = st.second;
                auto inv = mat->computeInvariants_PlaneStress(s1, s2);
                auto k = mat->computeDKDPrincipal_PlaneStress(inv, kappa);

                double dk1 = ( K(kappa, s1 + h, s2) - K(kappa, s1 - h, s2) ) / ( 2. * h );
                double dk2 = ( K(kappa, s1, s2 + h) - K(kappa, s1, s2 - h) ) / ( 2. * h );

                CHECK( closeEnough(k[0], dk1) );
                CHECK( closeEnough(k[1], dk2) );
            }
        }
    }

    SECTION( "G_ij matches finite differences of g1, g2 (plastic potential has no Lode angle dependence)" ) {
        auto g1fun = [&](double kappa, double s1, double s2) {
                         auto i = mat->computeInvariants_PlaneStress(s1, s2);
                         return mat->computeDGDPrincipal_PlaneStress(i, kappa)[0];
                     };
        auto g2fun = [&](double kappa, double s1, double s2) {
                         auto i = mat->computeInvariants_PlaneStress(s1, s2);
                         return mat->computeDGDPrincipal_PlaneStress(i, kappa)[1];
                     };

        auto allStates = genericStates;
        allStates.insert(allStates.end(), boundaryStates.begin(), boundaryStates.end() );

        for ( double kappa : kappas ) {
            for ( auto &st : allStates ) {
                double s1 = st.first, s2 = st.second;
                auto inv = mat->computeInvariants_PlaneStress(s1, s2);
                auto G = mat->computeDDGDDPrincipal_PlaneStress(inv, kappa);

                double G11 = ( g1fun(kappa, s1 + h, s2) - g1fun(kappa, s1 - h, s2) ) / ( 2. * h );
                double G12 = ( g1fun(kappa, s1, s2 + h) - g1fun(kappa, s1, s2 - h) ) / ( 2. * h );
                double G21 = ( g2fun(kappa, s1 + h, s2) - g2fun(kappa, s1 - h, s2) ) / ( 2. * h );
                double G22 = ( g2fun(kappa, s1, s2 + h) - g2fun(kappa, s1, s2 - h) ) / ( 2. * h );

                CHECK( closeEnough(G.at(1, 1), G11) );
                CHECK( closeEnough(G.at(1, 2), G12) );
                CHECK( closeEnough(G.at(2, 1), G21) );
                CHECK( closeEnough(G.at(2, 2), G22) );
            }
        }
    }

    SECTION( "G_i-kappa matches finite differences of g1, g2 with respect to kappa" ) {
        for ( double kappa : kappasForKappaFD ) {
            for ( auto &st : genericStates ) {
                auto inv = mat->computeInvariants_PlaneStress(st.first, st.second);
                auto Gk = mat->computeDDGDPrincipalDKappa_PlaneStress(inv, kappa);

                auto gPlus = mat->computeDGDPrincipal_PlaneStress(inv, kappa + hk);
                auto gMinus = mat->computeDGDPrincipal_PlaneStress(inv, kappa - hk);

                double dG1k = ( gPlus[0] - gMinus[0] ) / ( 2. * hk );
                double dG2k = ( gPlus[1] - gMinus[1] ) / ( 2. * hk );

                CHECK( closeEnough(Gk[0], dG1k, 1.e-2) );
                CHECK( closeEnough(Gk[1], dG2k, 1.e-2) );
            }
        }
    }

    SECTION( "computeDFDKappa_NoCutoff matches finite difference of the yield function with respect to kappa" ) {
        for ( double kappa : kappasForKappaFD ) {
            for ( auto &st : genericStates ) {
                auto inv = mat->computeInvariants_PlaneStress(st.first, st.second);
                double fk = mat->computeDFDKappa_NoCutoff(inv.sig, inv.rho, inv.theta, kappa);

                double fPlus = mat->computeYieldValue(inv.sig, inv.rho, inv.theta, kappa + hk);
                double fMinus = mat->computeYieldValue(inv.sig, inv.rho, inv.theta, kappa - hk);
                double dfk = ( fPlus - fMinus ) / ( 2. * hk );

                CHECK( closeEnough(fk, dfk, 1.e-2) );
            }
        }
    }

    SECTION( "Principal stress <-> global transform is a round trip" ) {
        std::vector< FloatArrayF< 3 > > cartesian = {
            { 3., 0., 0. }, { 0., -3., 0. }, { 1., -7., 2.5 },
            { -5., 2., -1.2 }, { 2., 2., 0. }, { -10., -10., 0. }
        };
        for ( auto &s : cartesian ) {
            double s1, s2, c, sn;
            mat->computePrincipalStress_PlaneStress(s, s1, s2, c, sn);
            CHECK( s1 >= s2 );
            auto back = mat->transformPrincipalToGlobal_PlaneStress(s1, s2, c, sn);
            CHECK( closeEnough(back[0], s[0], 1.e-9, 1.e-9) );
            CHECK( closeEnough(back[1], s[1], 1.e-9, 1.e-9) );
            CHECK( closeEnough(back[2], s[2], 1.e-9, 1.e-9) );
        }
    }

    SECTION( "performRegularReturn_PlaneStress converges to an admissible point on the yield surface" ) {
        // Integration-level check of the actual 4x4 Newton solve (not just its
        // building blocks): for a trial state outside the yield surface, the
        // regular return must converge, land (approximately) on f = 0, and
        // produce a non-negative plastic multiplier (admissibility, see
        // "Admissibility of the solution" in algoCDPM.tex).
        //
        // The trial excursions below are deliberately modest (a realistic
        // single load increment, not an arbitrarily large jump): with the
        // default AHard/BHard/CHard/DHard, x_h(sigV) -> 0 in tension, so
        // k = dKappa/dDeltaLambda becomes very large there and the basin of
        // convergence of a single, undamped Newton solve shrinks accordingly.
        // Recovering from larger jumps is the job of the strain-subdivision
        // loop in performPlasticityReturn_PlaneStress (not tested here, since
        // it requires a GaussPoint/status), consistent with "Admissibility of
        // the solution" in algoCDPM.tex, which prescribes damping or strain
        // subdivision as the fallback for exactly this situation.
        std::vector< std::pair< double, double > > trialStates = {
            { 3.05, 0.0 },     // uniaxial tension, just past ft = 3
            { 0.0, -32.0 },    // uniaxial compression, just past fc = 30
            { -32.0, -32.0 },  // equibiaxial compression, just past fc
            { 2.7, 2.7 }       // equibiaxial tension, just past its (lower) threshold
        };

        for ( auto &st : trialStates ) {
            double s1 = st.first, s2 = st.second;
            double deltaLambda = -1.;
            ConcreteDPM2::ConcreteDPM2_ReturnResult result = ConcreteDPM2::RR_Unknown;

            double kappa = mat->performRegularReturn_PlaneStress(s1, s2, deltaLambda, result, 0.0);

            CHECK( result == ConcreteDPM2::RR_Converged );
            CHECK( deltaLambda >= 0. );
            CHECK( kappa >= 0. );

            auto inv = mat->computeInvariants_PlaneStress(s1, s2);
            double residual = mat->computeYieldValue(inv.sig, inv.rho, inv.theta, kappa);
            CHECK( std::abs(residual) < 1.e-4 );
        }
    }

    SECTION( "g3 out-of-plane chain rule factors match finite differences of the full 3D invariants" ) {
        // g3 = (1/3) dg/dsigV - (s1+s2)/(3 rho) dg/drho is derived from the full
        // 3D invariants differentiated w.r.t. sigma_3 and evaluated at sigma_3 = 0
        // (not from the sigma_3 = 0 reduced rho/theta formulas used elsewhere).
        // Check the two chain-rule factors directly against finite differences
        // of the general 3D sig(s1,s2,s3), rho(s1,s2,s3).
        double hz = 1.e-3;
        auto allStates = genericStates;
        allStates.insert(allStates.end(), boundaryStates.begin(), boundaryStates.end() );

        for ( auto &st : allStates ) {
            double s1 = st.first, s2 = st.second;

            auto sigRho = [&](double s3) {
                             double sig = ( s1 + s2 + s3 ) / 3.;
                             double rho = std::sqrt( ( 1. / 3. ) * ( ( s1 - s2 ) * ( s1 - s2 ) +
                                                                     ( s2 - s3 ) * ( s2 - s3 ) +
                                                                     ( s3 - s1 ) * ( s3 - s1 ) ) );
                             return std::make_pair(sig, rho);
                         };

            auto plus = sigRho(hz);
            auto minus = sigRho(-hz);
            auto zero = sigRho(0.);

            double dSigDSig3 = ( plus.first - minus.first ) / ( 2. * hz );
            CHECK( closeEnough(dSigDSig3, 1. / 3., 1.e-9, 1.e-9) );

            if ( zero.second > 1.e-6 ) {
                double dRhoDSig3 = ( plus.second - minus.second ) / ( 2. * hz );
                CHECK( closeEnough(dRhoDSig3, -( s1 + s2 ) / ( 3. * zero.second ) ) );
            }
        }
    }
}
