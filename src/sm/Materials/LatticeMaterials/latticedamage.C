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
 *               Copyright (C) 1993 - 2019   Borek Patzak
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
LatticeDamage :: initializeFrom(InputRecord &ir)
{
    LatticeLinearElastic :: initializeFrom(ir);

    softeningType = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, softeningType, _IFT_LatticeDamage_softeningType); // Macro

    IR_GIVE_FIELD(ir, wf, _IFT_LatticeDamage_wf); // Macro

    if ( softeningType == 2 ) { //bilinear softening
        wfOne = 0.15 * wf;
        IR_GIVE_OPTIONAL_FIELD(ir, wfOne, _IFT_LatticeDamage_wfOne); // Macro
        e0OneMean = 0.3 * e0Mean;
        IR_GIVE_OPTIONAL_FIELD(ir, e0OneMean, _IFT_LatticeDamage_e0OneMean);
    }

    IR_GIVE_OPTIONAL_FIELD(ir, e0Mean, _IFT_LatticeDamage_e0Mean); // Macro

    this->coh = 2.;
    IR_GIVE_OPTIONAL_FIELD(ir, coh, _IFT_LatticeDamage_coh); // Macro

    this->ec = 10.;
    IR_GIVE_OPTIONAL_FIELD(ir, ec, _IFT_LatticeDamage_ec); // Macro

    this->biotCoefficient = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, this->biotCoefficient, _IFT_LatticeDamage_bio);

    this->biotType = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->biotType, _IFT_LatticeDamage_btype);
}


double
LatticeDamage :: computeEquivalentStrain(const FloatArrayF< 6 > &strain, GaussPoint *gp) const
{
    const double e0 = this->give(e0_ID, gp) * this->e0Mean;
    double paramA = 0.5 * ( e0 + ec * e0 );
    double paramB = ( coh * e0 ) / sqrt(1. - pow( ( ec * e0 - e0 ) / ( e0 + ec * e0 ), 2. ) );
    double paramC = 0.5 * ( this->ec * e0 - e0 );

    double shearNorm = norm(strain [ { 1, 2 } ]);
    return norm({ this->alphaOne * shearNorm / paramB, ( strain.at(1) + paramC ) / paramA }) * paramA - paramC;
}


double
LatticeDamage :: computeDamageParam(double tempKappa, GaussPoint *gp) const
{
    double le = static_cast< LatticeStructuralElement * >( gp->giveElement() )->giveLength();
    const double e0 = this->give(e0_ID, gp) * this->e0Mean;
    double eNormal = this->give(eNormal_ID, gp) * this->eNormalMean;

    if ( softeningType == 1 ) { //linear
        if ( tempKappa >= e0 && tempKappa < this->wf / le ) {
            //linear stress-crack opening relation
            //check if input parameter make sense
            if ( this->wf / le <= e0 ) {
                OOFEM_ERROR("e0>wf/Le \n Possible solutions: Increase fracture energy or reduce element size\n");
            }

            return ( 1. - e0 / tempKappa ) / ( 1. - e0 / ( this->wf / le ) );
        } else if ( tempKappa >= this->wf / le ) {
            return 1.;
        } else {
            return 0.;
        }
    } else if ( softeningType == 2 ) {      //bilinear softening
        //Check if input parameter make sense
        if ( e0 > wfOne / le ) {
            OOFEM_ERROR("parameter wf1 is too small");
        } else if ( wfOne / le >  this->wf / le ) {
            OOFEM_ERROR("parameter wf is too small");
        }

        if ( tempKappa > e0 ) {
            double helpStrain = 0.3 * e0;
            double omega = ( 1 - e0 / tempKappa ) / ( ( helpStrain - e0 ) / this->wfOne * le + 1. );

            if ( omega * tempKappa * le > 0 && omega * tempKappa * le < this->wfOne ) {
                return omega;
            } else {
                omega = ( 1. - helpStrain / tempKappa - helpStrain * this->wfOne / ( tempKappa * ( this->wf - this->wfOne ) ) ) / ( 1. - helpStrain * le / ( this->wf - this->wfOne ) );

                if ( omega * tempKappa * le >= this->wfOne  && omega * tempKappa * le < this->wf  ) {
                    return omega;
                }
            }

            return clamp(omega, 0., 1.);
        } else {
            return 0.;
        }
    } else if ( softeningType == 3 ) {
        //exponential softening
        //  iteration to achieve objectivity
        //   we are finding state, where elastic stress is equal to
        //   stress from crack-opening relation (wf = wf characterizes the carc opening diagram)
        if ( tempKappa <= e0 ) {
            return 0.0;
        } else {
            int nite = 0;
            double R = 0., omega = 0.;
            double Ft = eNormal * e0;
            do {
                nite++;
                double help = le * omega * tempKappa / this->wf;
                double Lhs = eNormal * tempKappa - Ft * exp(-help) * le * tempKappa / this->wf;
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
}


MaterialStatus *
LatticeDamage :: CreateStatus(GaussPoint *gp) const
{
    return new LatticeDamageStatus(gp);
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

    double tempDeltaDissipation = computeDeltaDissipation3d(omega, reducedStrain, gp, tStep);
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

    //Compute crack width
    double length = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveLength();
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
    if ( softeningType == 1 ) {
        return e0 * eNormal * this->wf / 2.;
    } else {   //This is for the exponential law. Should also implement it for the bilinear one.
        return e0 * eNormal * this->wf;
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
