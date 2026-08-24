/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2026   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "latticedamage.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatmatrixf.h"
#include "floatarray.h"
#include "floatarrayf.h"
#include "CrossSections/latticecrosssection.h"
#include "engngm.h"
#include "mathfem.h"
#include "Elements/LatticeElements/latticestructuralelement.h"
#include "integrationrule.h"
#include "datastream.h"
#include "staggeredproblem.h"
#include "contextioerr.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(LatticeDamage);

LatticeDamage :: LatticeDamage(int n, Domain *d) : LatticeLinearElastic(n, d)
{}


bool
LatticeDamage :: hasMaterialModeCapability(MaterialMode mode) const
{
    return ( mode == _3dLattice );
}


void
LatticeDamage :: initializeFrom(const std::shared_ptr<InputRecord> &ir)
{
    LatticeLinearElastic :: initializeFrom(ir);

    softeningType = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, softeningType, _IFT_LatticeDamage_softeningType); // Macro

    // Peak: tensile strength as ft (=> e0 = ft/eNormalMean) or e0 directly.
    const bool ftGiven = ir->hasField(_IFT_LatticeDamage_ft);
    if ( ftGiven && ir->hasField(_IFT_LatticeDamage_e0Mean) ) {
        OOFEM_ERROR("give either ft or e0, not both (e0 is derived from ft)");
    }
    if ( ftGiven ) {
        double ftInput = 0.;
        IR_GIVE_FIELD(ir, ftInput, _IFT_LatticeDamage_ft);
        if ( ftInput <= 0. ) {
            OOFEM_ERROR("ft must be positive");
        }
        this->e0Mean = ftInput / this->eNormalMean; // eNormalMean set by the elastic base above
    } else {
        IR_GIVE_OPTIONAL_FIELD(ir, e0Mean, _IFT_LatticeDamage_e0Mean); // Macro
    }
    const double ft = this->eNormalMean * this->e0Mean; // effective tensile strength

    // Optional crack spacing. When given (> 0), the softening is regularized by sm (a fixed
    // crack spacing) instead of the element length, so cracking cannot localize below sm
    // (distributed cracking, e.g. reinforced concrete); crack widths then use sm as the band.
    IR_GIVE_OPTIONAL_FIELD(ir, sm, _IFT_LatticeDamage_sm);
    if ( sm < 0. ) {
        OOFEM_ERROR("sm must be positive");
    }

    // Bilinear (stype 2) kink ratios: kink stress = e0Ratio*ft, kink opening = wfRatio*wf.
    if ( softeningType == 2 ) {
        IR_GIVE_OPTIONAL_FIELD(ir, wfRatio, _IFT_LatticeDamage_wfRatio);
        IR_GIVE_OPTIONAL_FIELD(ir, e0Ratio, _IFT_LatticeDamage_e0Ratio);
        if ( wfRatio <= 0. || wfRatio >= 1. || e0Ratio <= 0. || e0Ratio >= 1. ) {
            OOFEM_ERROR("wfratio and e0ratio must lie in (0,1)");
        }
    } else if ( ir->hasField(_IFT_LatticeDamage_wfRatio) || ir->hasField(_IFT_LatticeDamage_e0Ratio) ) {
        OOFEM_ERROR("wfratio/e0ratio apply only to stype 2 (bilinear softening)");
    }

    // Softening extent: exactly one of wf (direct) or gf (derived from fracture energy).
    const bool wfGiven = ir->hasField(_IFT_LatticeDamage_wf);
    const bool gfGiven = ir->hasField(_IFT_LatticeDamage_gf);
    if ( wfGiven == gfGiven ) {
        OOFEM_ERROR("give exactly one softening extent: wf or gf");
    }
    if ( gfGiven ) { // derived from fracture energy
        if ( ft <= 0. ) {
            OOFEM_ERROR("gf requires a positive ft (give ft or e0)");
        }
        double gf = 0.;
        IR_GIVE_FIELD(ir, gf, _IFT_LatticeDamage_gf);
        if ( gf <= 0. ) {
            OOFEM_ERROR("gf must be positive");
        }
        if ( softeningType == 1 ) { // Gf = ft*wf/2
            this->wf = 2. * gf / ft;
        } else if ( softeningType == 2 ) { // Gf = ft*wf*(wfRatio + e0Ratio)/2
            this->wf = 2. * gf / ( ft * ( wfRatio + e0Ratio ) );
            this->wfOne = wfRatio * this->wf;
        } else { // Gf = ft*wf
            this->wf = gf / ft;
        }
    } else { // wf given directly
        IR_GIVE_FIELD(ir, wf, _IFT_LatticeDamage_wf); // Macro
        if ( softeningType == 2 ) {
            this->wfOne = wfRatio * this->wf;
        }
    }

    this->coh = 2.;
    IR_GIVE_OPTIONAL_FIELD(ir, coh, _IFT_LatticeDamage_coh); // Macro

    this->ec = 10.;
    IR_GIVE_OPTIONAL_FIELD(ir, ec, _IFT_LatticeDamage_ec); // Macro

    // biotCoefficient ("bio") is read by LatticeLinearElastic::initializeFrom (called above).
    this->biotType = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->biotType, _IFT_LatticeDamage_btype);

    this->calDissFlag = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->calDissFlag, _IFT_LatticeDamage_calDiss);
}


double
LatticeDamage :: computeEquivalentStrain(const FloatArrayF< 6 > &strain, GaussPoint *gp) const
{
    const double e0 = this->give(e0_ID, gp) * this->e0Mean;
    double paramA = 0.5 * ( e0 + ec * e0 );
    double paramB = ( coh * e0 ) / sqrt(1. - pow( ( ec * e0 - e0 ) / ( e0 + ec * e0 ), 2. ) );
    double paramC = 0.5 * ( this->ec * e0 - e0 );

    double shearNorm = norm(strain [ { 1, 2 } ]);
    return norm(FloatArrayF< 2 >{ this->alphaOne * shearNorm / paramB, ( strain.at(1) + paramC ) / paramA }) * paramA - paramC;
}


double
LatticeDamage :: computeDamageParamExplicit(double tempKappa, double e0, double wfParam, double eNormal, double le) const
{
    if ( softeningType == 1 ) {
        if ( tempKappa >= e0 && tempKappa < wfParam / le ) {
            if ( wfParam / le <= e0 ) {
                OOFEM_ERROR("e0>wf/Le \n Possible solutions: Increase fracture energy or reduce element size\n");
            }
            return ( 1. - e0 / tempKappa ) / ( 1. - e0 / ( wfParam / le ) );
        } else if ( tempKappa >= wfParam / le ) {
            return 1.;
        } else {
            return 0.;
        }
    } else if ( softeningType == 2 ) {
        if ( e0 > wfOne / le ) {
            OOFEM_ERROR("parameter wf1 is too small");
        } else if ( wfOne / le > wfParam / le ) {
            OOFEM_ERROR("parameter wf is too small");
        }
        if ( tempKappa > e0 ) {
            double helpStrain = this->e0Ratio * e0;
            double omega = ( 1 - e0 / tempKappa ) / ( ( helpStrain - e0 ) / this->wfOne * le + 1. );
            if ( omega * tempKappa * le > 0 && omega * tempKappa * le < this->wfOne ) {
                return omega;
            } else {
                omega = ( 1. - helpStrain / tempKappa - helpStrain * this->wfOne / ( tempKappa * ( wfParam - this->wfOne ) ) ) / ( 1. - helpStrain * le / ( wfParam - this->wfOne ) );
                if ( omega * tempKappa * le >= this->wfOne && omega * tempKappa * le < wfParam ) {
                    return omega;
                }
            }
            return clamp(omega, 0., 1.);
        } else {
            return 0.;
        }
    } else if ( softeningType == 3 ) {
        if ( tempKappa <= e0 ) {
            return 0.0;
        } else {
            int nite = 0;
            double R = 0., omega = 0.;
            double Ft = eNormal * e0;
            do {
                nite++;
                double help = le * omega * tempKappa / wfParam;
                double Lhs = eNormal * tempKappa - Ft * exp(-help) * le * tempKappa / wfParam;
                R = ( 1. - omega ) * eNormal * tempKappa - Ft * exp(-help);
                omega += R / Lhs;
                if ( nite > 40 ) {
                    OOFEM_ERROR("computeDamageParam: algorithm not converging");
                }
            } while ( fabs(R) >= 1.e-4 );
            if ( ( omega > 1.0 ) || ( omega < 0.0 ) ) {
                OOFEM_ERROR("computeDamageParam: internal error\n");
            }
            return omega;
        }
    } else {
        OOFEM_ERROR("computeDamageParam: unknown softening type");
    }
    return 0.;
}


double
LatticeDamage :: computeDamageParam(double tempKappa, GaussPoint *gp) const
{
    double e0      = this->give(e0_ID, gp) * this->e0Mean;
    double eNormal = this->give(eNormal_ID, gp) * this->eNormalMean;
    // Characteristic length: the crack spacing sm if given (so cracking cannot localize
    // below sm), otherwise the element length.
    double le = ( this->sm > 0. ) ? this->sm
              : static_cast< LatticeStructuralElement * >( gp->giveElement() )->giveLength();
    return computeDamageParamExplicit(tempKappa, e0, this->wf, eNormal, le);
}


std::unique_ptr< MaterialStatus >
LatticeDamage :: CreateStatus(GaussPoint *gp) const
{
    return std::make_unique< LatticeDamageStatus >(gp);
}


FloatMatrixF< 6, 6 >
LatticeDamage :: give3dLatticeStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    auto elastic = LatticeLinearElastic :: give3dLatticeStiffnessMatrix(mode, gp, tStep);

    if ( mode == ElasticStiffness ) {
        return elastic;
    } else if ( ( mode == SecantStiffness ) || ( mode == TangentStiffness ) ) {
        auto status = static_cast< LatticeDamageStatus * >( this->giveStatus(gp) );

        double omega = min(status->giveTempDamage(), 0.99999);
        return elastic * ( 1. - omega );
    } else {
        OOFEM_ERROR("Unsupported stiffness mode\n");
        return elastic;
    }
}


FloatMatrixF< 3, 3 >
LatticeDamage :: give2dLatticeStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{

    auto elastic = LatticeLinearElastic :: give2dLatticeStiffnessMatrix(mode, gp, tStep);

    if ( mode == ElasticStiffness ) {
        return elastic;
    } else if ( ( mode == SecantStiffness ) || ( mode == TangentStiffness ) ) {
        auto status = static_cast< LatticeDamageStatus * >( this->giveStatus(gp) );
        double omega = status->giveTempDamage();

        if ( omega > 0.99999 ) {
            omega = 0.99999;
        }

        return elastic * ( 1. - omega );
    } else {
        OOFEM_ERROR("Unsupported stiffness mode\n");
        return elastic;
    }
}


FloatArrayF< 6 >
LatticeDamage :: giveLatticeStress3d(const FloatArrayF< 6 > &strain, GaussPoint *gp, TimeStep *tStep)
{
    auto status = static_cast< LatticeDamageStatus * >( this->giveStatus(gp) );

    const double e0 = this->give(e0_ID, gp) * this->e0Mean;

    status->setE0(e0);
    this->initTempStatus(gp);

    // substract stress independent part
    auto reducedStrain = strain;
    auto thermalStrain = this->computeStressIndependentStrainVector(gp, tStep, VM_Total);
    if ( thermalStrain.giveSize() ) {
        reducedStrain -= FloatArrayF< 6 >(thermalStrain);
    }

    double omega = 0.;
    this->performDamageEvaluation(gp, reducedStrain);
    omega = status->giveTempDamage();

    auto stiffnessMatrix = LatticeLinearElastic :: give3dLatticeStiffnessMatrix(ElasticStiffness, gp, tStep);

    FloatArrayF< 6 >answer;
    for ( int i = 1; i <= 6; i++ ) { // only diagonal terms matter
        answer.at(i) = stiffnessMatrix.at(i, i) * reducedStrain.at(i) * ( 1. - omega );
    }

    //Read in fluid pressures from structural element if this is not a slave problem
    FloatArray pressures;
    if ( !domain->giveEngngModel()->giveMasterEngngModel() ) {
        static_cast< LatticeStructuralElement * >( gp->giveElement() )->givePressures(pressures);
    }


    double waterPressure = 0.;
    for ( int i = 0; i < pressures.giveSize(); i++ ) {
        waterPressure += 1. / pressures.giveSize() * pressures [ i ];
    }
    answer.at(1) += waterPressure;

    double tempDeltaDissipation = calDissFlag ? computeDeltaDissipation3d(omega, reducedStrain, gp, tStep) : 0.;
    double tempDissipation = status->giveDissipation() + tempDeltaDissipation;

    //Set all temp values
    status->setTempDissipation(tempDissipation);
    status->setTempDeltaDissipation(tempDeltaDissipation);

    status->letTempLatticeStrainBe(strain);
    status->letTempReducedLatticeStrainBe(reducedStrain);
    status->letTempLatticeStressBe(answer);
    status->setTempNormalLatticeStress(answer.at(1) );

    return answer;
}


void
LatticeDamage :: performDamageEvaluation(GaussPoint *gp, FloatArrayF< 6 > &reducedStrain) const
{
    auto status = static_cast< LatticeDamageStatus * >( this->giveStatus(gp) );
    // compute equivalent strain
    double equivStrain = this->computeEquivalentStrain(reducedStrain, gp);

    // compute value of loading function if strainLevel crit apply
    double f = equivStrain - status->giveKappa();

    double tempKappa, omega = 0.;
    if ( f <= 0.0 ) {
        // damage does not grow
        tempKappa = status->giveKappa();
        omega = status->giveDamage();
        if ( status->giveCrackFlag() != 0 ) {
            status->setTempCrackFlag(2);
        } else {
            status->setTempCrackFlag(0);
        }
    } else {
        // damage grows
        tempKappa = equivStrain;

        // evaluate damage parameter

        omega = this->computeDamageParam(tempKappa, gp);
        if ( omega > 0 ) {
            status->setTempCrackFlag(1);
        }
    }

    //Create history variables for the damage strain
    FloatArrayF< 6 >tempDamageLatticeStrain = omega * reducedStrain;
    status->letTempDamageLatticeStrainBe(tempDamageLatticeStrain);

    //Compute crack width; use the crack spacing sm as the band width when given.
    double length = ( this->sm > 0. ) ? this->sm
                  : ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveLength();
    double crackWidth = omega * norm(reducedStrain [ { 0, 1, 2 } ]) * length;

    status->setTempEquivalentStrain(equivStrain);
    status->setTempKappa(tempKappa);
    status->setTempDamage(omega);
    status->setTempCrackWidth(crackWidth);
}

double
LatticeDamage :: computeBiot(double omega, double kappa, double le) const
{
    if ( this->softeningType == 1 || this->softeningType == 3 ) {
        if ( omega == 0 ) {
            return this->biotCoefficient;
        } else if ( omega * kappa * le > 0 && omega * kappa * le < this->wf ) {
            return this->biotCoefficient + ( 1. - biotCoefficient ) * omega * kappa * le / this->wf;
        } else {
            return 1.;
        }
    } else {
        OOFEM_ERROR("Wrong stype for btype=1. Only linear and exponential softening considered so far\n");
        return 0.;
    }
}


double
LatticeDamage :: computeReferenceGf(GaussPoint *gp) const
{
    const double e0 = this->give(e0_ID, gp) * this->e0Mean;
    const double eNormal = this->give(eNormal_ID, gp) * this->eNormalMean;
    const double ft = e0 * eNormal;
    // Per-element reference dissipation. With a crack spacing sm the softening spans more
    // strain than one element, so the element sees a fraction he/sm of the crack energy.
    double factor = 1.;
    if ( this->sm > 0. ) {
        double he = static_cast< LatticeStructuralElement * >( gp->giveElement() )->giveLength();
        factor = he / this->sm;
    }
    if ( softeningType == 1 ) { // linear
        return 0.5 * ft * this->wf * factor;
    } else if ( softeningType == 2 ) { // bilinear
        return 0.5 * ft * ( this->wfOne + this->e0Ratio * this->wf ) * factor;
    } else { // exponential
        return ft * this->wf * factor;
    }
}


double
LatticeDamage :: computeIntervals(double testDissipation, double referenceGf) const
{
    if ( testDissipation / ( referenceGf ) > 0.01 ) {
        return 1000. * min(testDissipation / referenceGf, 1.);
    } else {
        return 1.;
    }
}


double
LatticeDamage :: computeDeltaDissipation2d(double omega,
                                           const FloatArrayF< 3 > &reducedStrain,
                                           GaussPoint *gp,
                                           TimeStep *tStep) const
{
    auto status = static_cast< LatticeDamageStatus * >( this->giveStatus(gp) );
    const double length = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveLength();
    const double eNormal = this->give(eNormal_ID, gp) * this->eNormalMean;
    const double eShear = this->alphaOne * eNormal;
    const double eTorsion = this->alphaTwo * eNormal;
    const FloatArrayF< 3 >tangent = { eNormal, eShear, eTorsion };

    const auto reducedStrainOld = status->giveReducedLatticeStrain() [ { 0, 1, 5 } ];
    const double omegaOld = status->giveDamage();

    auto oldIntermediateStrain = reducedStrainOld;
    double oldIntermediateOmega = omegaOld;
    double deltaOmega = ( omega - omegaOld );
    auto testStrain = ( reducedStrain + reducedStrainOld ) * 0.5;
    auto testStress = mult(tangent, testStrain);
    double testDissipation = 0.5 * length * norm(testStress) * deltaOmega;

    double referenceGf = computeReferenceGf(gp);
    double intervals = computeIntervals(testDissipation, referenceGf);

    double oldKappa = status->giveKappa();
    if ( deltaOmega > 0 ) {
        double tempDeltaDissipation = 0.;
        for ( int k = 0; k < intervals; k++ ) {
            auto intermediateStrain = reducedStrainOld + ( k + 1 ) / intervals * ( reducedStrain - reducedStrainOld );
            double equivStrain = this->computeEquivalentStrain(assemble< 6 >(intermediateStrain, { 0, 1, 5 }), gp);
            double f = equivStrain - oldKappa;
            if ( f > 0 ) {
                auto intermediateOmega = this->computeDamageParam(equivStrain, gp);
                auto deltaOmega = ( intermediateOmega - oldIntermediateOmega );
                auto midStrain = ( intermediateStrain + oldIntermediateStrain ) / 2.;
                auto midStress = mult(tangent, midStrain);
                auto deltaTempDeltaDissipation = 0.5 * length * norm(midStress) * deltaOmega;
                oldKappa = equivStrain;
                oldIntermediateOmega = intermediateOmega;
                tempDeltaDissipation += deltaTempDeltaDissipation;
            }

            oldIntermediateStrain = intermediateStrain;
        }
        return max(tempDeltaDissipation, 2. * referenceGf);
    } else {
        return 0.;
    }
}


double
LatticeDamage :: computeDeltaDissipation3d(double omega,
                                           const FloatArrayF< 6 > &reducedStrain,
                                           GaussPoint *gp,
                                           TimeStep *atTime) const
{
    auto status = static_cast< LatticeDamageStatus * >( this->giveStatus(gp) );
    const double length = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveLength();
    const double eNormal = this->give(eNormal_ID, gp) * this->eNormalMean;
    const double eShear = this->alphaOne * eNormal;
    const double eTorsion = this->alphaTwo * eNormal;
    const FloatArrayF< 6 >tangent = { eNormal, eShear, eShear, eTorsion, eTorsion, eTorsion };

    const auto &reducedStrainOld = status->giveReducedLatticeStrain();
    const double omegaOld = status->giveDamage();

    auto oldIntermediateStrain = reducedStrainOld;
    double oldIntermediateOmega = omegaOld;
    double deltaOmega = ( omega - omegaOld );
    auto testStrain = ( reducedStrain + reducedStrainOld ) * 0.5;
    auto testStress = mult(tangent, testStrain);
    double testDissipation = 0.5 * length * norm(testStress) * deltaOmega;

    double referenceGf = computeReferenceGf(gp);
    double intervals = computeIntervals(testDissipation, referenceGf);

    double oldKappa = status->giveKappa();
    if ( deltaOmega > 0 ) {
        double tempDeltaDissipation = 0.;
        for ( int k = 0; k < intervals; k++ ) {
            auto intermediateStrain = reducedStrainOld + ( k + 1 ) / intervals * ( reducedStrain - reducedStrainOld );
            double equivStrain = this->computeEquivalentStrain(intermediateStrain, gp);
            double f = equivStrain - oldKappa;
            if ( f > 0 ) {
                auto intermediateOmega = this->computeDamageParam(equivStrain, gp);
                auto deltaOmega = ( intermediateOmega - oldIntermediateOmega );
                auto midStrain = ( intermediateStrain + oldIntermediateStrain ) / 2.;
                auto midStress = mult(tangent, midStrain);
                auto deltaTempDeltaDissipation = 0.5 * length * norm(midStress) * deltaOmega;
                oldKappa = equivStrain;
                oldIntermediateOmega = intermediateOmega;
                tempDeltaDissipation += deltaTempDeltaDissipation;
            }

            oldIntermediateStrain = intermediateStrain;
        }
        return max(tempDeltaDissipation, 2. * referenceGf);
    } else {
        return 0.;
    }
}


double
LatticeDamage :: give(int aProperty, GaussPoint *gp) const
{
  this->giveStatus(gp);
  
    double answer;
    if ( RandomMaterialExtensionInterface :: give(aProperty, gp, answer) ) {
        if ( answer < 0.1 ) { //Introduce cut off to avoid numerical problems
            answer = 0.1;
        } else if ( answer > 10 ) {
            answer = 10;
        }
        return answer;
    } else if ( aProperty == e0_ID ) {
        return 1.;
    } else if ( aProperty == ef_ID ) {
        return 1.;
    } else {
        return LatticeLinearElastic :: give(aProperty, gp);
    }
}

int
LatticeDamage :: giveIPValue(FloatArray &answer,
                             GaussPoint *gp,
                             InternalStateType type,
                             TimeStep *atTime)
{
    LatticeDamageStatus *status = static_cast< LatticeDamageStatus * >( this->giveStatus(gp) );
    if ( type == IST_CrackedFlag ) {
        answer.resize(1);
        answer.zero();
        answer.at(1) = status->giveCrackFlag();
        return 1;
    } else if ( type == IST_DamageScalar ) {
        answer.resize(1);
        answer.zero();
        answer.at(1) = status->giveDamage();
        return 1;
    } else if ( type == IST_DamageTensor ) {
        answer.resize(6);
        answer.zero();
        answer.at(1) = answer.at(2) = answer.at(3) = status->giveDamage();
        return 1;
    } else if ( type == IST_DissWork ) {
        answer.resize(1);
        answer.zero();
        answer.at(1) = status->giveDissipation();
        return 1;
    } else if ( type == IST_DeltaDissWork ) {
        answer.resize(1);
        answer.zero();
        answer.at(1) = status->giveDeltaDissipation();
        return 1;
    } else if ( type == IST_CrackWidth ) {
        answer.resize(1);
        answer.zero();
        answer.at(1) = status->giveCrackWidth();
        return 1;
    } else if ( type == IST_NormalStress ) {
        answer.resize(1);
        answer.zero();
        answer.at(1) = status->giveNormalLatticeStress();
        return 1;
    } else if ( type == IST_CharacteristicLength ) {
        answer.resize(1);
        answer.zero();
        answer.at(1) = static_cast< LatticeStructuralElement * >( gp->giveElement() )->giveLength();
        return 1;
    } else if ( type == IST_TensileStrength ) {
        answer.resize(1);
        answer.zero();
        answer.at(1) = this->give(e0_ID, gp) * this->e0Mean;
        return 1;
    } else {
        return LatticeLinearElastic :: giveIPValue(answer, gp, type, atTime);
    }
}

LatticeDamageStatus :: LatticeDamageStatus(GaussPoint *g) :
    LatticeMaterialStatus(g)
{ }

void
LatticeDamageStatus :: initTempStatus()
{
    LatticeMaterialStatus :: initTempStatus();

    this->tempKappa = this->kappa;
    this->tempEquivStrain = this->equivStrain;
    this->tempDamage = this->damage;
}

void
LatticeDamageStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
{
    LatticeMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "kappa %f, equivStrain %f, damage %f, dissipation %f, deltaDissipation %f, e0 %f, crackFlag %d\n", this->kappa, this->equivStrain, this->damage, this->dissipation, this->deltaDissipation, this->e0, this->crackFlag);
}


void
LatticeDamageStatus :: updateYourself(TimeStep *atTime)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reached equilibrium.
//
{
    LatticeMaterialStatus :: updateYourself(atTime);

    this->kappa = this->tempKappa;
    this->equivStrain = this->tempEquivStrain;
    this->damage = this->tempDamage;
}

void
LatticeDamageStatus :: saveContext(DataStream &stream, ContextMode mode)
//
// saves full information stored in this Status
// no temp variables stored
//
{
    LatticeMaterialStatus :: saveContext(stream, mode);

    if ( !stream.write(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(equivStrain) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(damage) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(e0) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(biot) ) {
        THROW_CIOERR(CIO_IOERR);
    }
}

void
LatticeDamageStatus :: restoreContext(DataStream &stream, ContextMode mode)
//
// restores full information stored in stream to this Status
//
{
    LatticeMaterialStatus :: restoreContext(stream, mode);

    if ( !stream.read(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(equivStrain) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(damage) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(e0) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(biot) ) {
        THROW_CIOERR(CIO_IOERR);
    }
}
} // end namespace oofem
