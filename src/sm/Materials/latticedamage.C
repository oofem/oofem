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
 *               Copyright (C) 1993 - 2013   Borek Patzak
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
#include "floatarray.h"
#include "CrossSections/structuralcrosssection.h"
#include "engngm.h"
#include "mathfem.h"
#include "Elements/latticestructuralelement.h"
#include "datastream.h"
#include "staggeredproblem.h"
#include "contextioerr.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(LatticeDamage);

LatticeDamage :: LatticeDamage(int n, Domain *d) : LatticeLinearElastic(n, d)
{}


LatticeDamage :: ~LatticeDamage()
//
// destructor
//
{}

int
LatticeDamage :: hasMaterialModeCapability(MaterialMode mode)
//
// returns whether receiver supports given mode
//
{
    if ( ( mode == _2dLattice ) || ( mode == _3dLattice ) ) {
        return 1;
    }

    return 0;
}


IRResultType
LatticeDamage :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                             // Required by IR_GIVE_FIELD macro

    softeningType = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, softeningType, _IFT_LatticeDamage_softeningType); // Macro

    IR_GIVE_FIELD(ir, wf, _IFT_LatticeDamage_wf); // Macro

    if ( softeningType == 2 ) { //bilinear softening
        wfOne = 0.15 * wf;
        IR_GIVE_OPTIONAL_FIELD(ir, wfOne, _IFT_LatticeDamage_wfOne); // Macro
        e0OneMean = 0.3 * e0Mean;
        IR_GIVE_OPTIONAL_FIELD(ir, e0OneMean, _IFT_LatticeDamage_e0OneMean);
    }


    IR_GIVE_FIELD(ir, e0Mean, _IFT_LatticeDamage_e0Mean); // Macro

    this->coh = 2.;
    IR_GIVE_OPTIONAL_FIELD(ir, coh, _IFT_LatticeDamage_coh); // Macro

    this->ec = 10.;
    IR_GIVE_OPTIONAL_FIELD(ir, ec, _IFT_LatticeDamage_ec); // Macro

    this->biotCoefficient = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, this->biotCoefficient, _IFT_LatticeDamage_bio);

    this->biotType = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->biotType, _IFT_LatticeDamage_btype);

    return LatticeLinearElastic :: initializeFrom(ir);
}


void
LatticeDamage :: computeEquivalentStrain(double &tempEquivStrain, const FloatArray &strain, GaussPoint *gp, TimeStep *atTime)
{
    if ( strain.isEmpty() ) {
        tempEquivStrain = 0.;
        return;
    }

    const double e0 = this->give(e0_ID, gp) * this->e0Mean;
    double paramA, paramB, paramC;
    double shearNorm;

    paramA = 0.5 * ( e0 + ec * e0 );
    paramB = ( coh * e0 ) / sqrt( 1. - pow( ( ec * e0 - e0 ) / ( e0 + ec * e0 ), 2. ) );
    paramC = 0.5 * ( this->ec * e0 - e0 );


    if ( gp->giveMaterialMode() == _2dLattice ) {
        shearNorm = strain.at(2);
    } else {
        shearNorm = sqrt( pow(strain.at(2), 2.) + pow(strain.at(3), 2.) );
    }

    tempEquivStrain =  sqrt( pow(this->alphaOne * shearNorm / paramB, 2.) + pow( ( strain.at(1) + paramC ) / paramA, 2. ) ) * paramA - paramC;

    return;
}


void
LatticeDamage :: computeDamageParam(double &omega, double tempKappa, GaussPoint *gp)
{
    double le = static_cast< LatticeStructuralElement * >( gp->giveElement() )->giveLength();
    const double e0 = this->give(e0_ID, gp) * this->e0Mean;
    double eNormal = this->give(eNormal_ID, gp) * this->eNormalMean;

    omega = 0.0;

    int nite = 0;
    double R, Lhs, Ft, help;

    if ( softeningType == 1 ) { //linear
        if ( tempKappa >= e0 && tempKappa < this->wf / le ) {
            //linear stress-crack opening relation
            //check if input parameter make sense
            if ( this->wf / le <= e0 ) {
                OOFEM_ERROR("e0>wf/Le \n Possible solutions: Increase fracture energy or reduce element size\n");
            }

            omega = ( 1. - e0 / tempKappa ) / ( 1. - e0 / ( this->wf / le ) );
        } else if ( tempKappa >= this->wf / le ) {
            omega = 1.;
        } else {
            omega = 0.;
        }
    } else if ( softeningType == 2 ) {      //bilinear softening
        double helpStrain = 0.3 * e0;

        //Check if input parameter make sense
        if ( e0 > wfOne / le ) {
            OOFEM_ERROR("parameter wf1 is too small");
        } else if ( wfOne / le >  this->wf / le ) {
            OOFEM_ERROR("parameter wf is too small");
        }

        //
        if ( tempKappa > e0 ) {
            omega = ( 1 - e0 / tempKappa ) / ( ( helpStrain - e0 ) / this->wfOne * le + 1. );

            if ( omega * tempKappa * le > 0 && omega * tempKappa * le < this->wfOne ) {
                return;
            } else {
                omega = ( 1. - helpStrain / tempKappa - helpStrain * this->wfOne / ( tempKappa * ( this->wf - this->wfOne ) ) ) / ( 1. - helpStrain * le / ( this->wf - this->wfOne ) );

                if ( omega * tempKappa * le >= this->wfOne  && omega * tempKappa * le < this->wf  ) {
                    return;
                }
            }

            if ( omega > 1. ) {
                omega = 1.;
            } else if ( omega < 0 ) {
                omega = 0;
            }
        } else {
            omega = 0.;
        }
    } else if ( softeningType == 3 ) {      //exponential softening
        //  iteration to achieve objectivity
        //   we are finding state, where elastic stress is equal to
        //   stress from crack-opening relation (wf = wf characterizes the carc opening diagram)

        if ( tempKappa <= e0 ) {
            omega = 0.0;
        } else {
            omega = 0.0;
            Ft = eNormal * e0;
            do {
                nite++;
                help = le * omega * tempKappa / this->wf;
                R = ( 1. - omega ) * eNormal * tempKappa - Ft *exp(-help);
                Lhs = eNormal * tempKappa - Ft *exp(-help) * le * tempKappa / this->wf;
                omega += R / Lhs;
                if ( nite > 40 ) {
                    OOFEM_ERROR("computeDamageParam: algorithm not converging");
                }
            } while ( fabs(R) >= 1.e-4 );

            if ( ( omega > 1.0 ) || ( omega < 0.0 ) ) {
                OOFEM_ERROR("computeDamageParam: internal error\n");
            }
        }
    } else {
        OOFEM_ERROR("Unknown softening type");
    }
}



MaterialStatus *
LatticeDamage :: CreateStatus(GaussPoint *gp) const
{
    return new LatticeDamageStatus(gp);
}



void
LatticeDamage :: giveRealStressVector(FloatArray &answer,
                                      GaussPoint *gp,
                                      const FloatArray &totalStrain,
                                      TimeStep *atTime)
//
// returns real stress vector in 3d stress space of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
{
    LatticeDamageStatus *status = static_cast< LatticeDamageStatus * >( this->giveStatus(gp) );

    const double e0 = this->give(e0_ID, gp) * this->e0Mean;

    status->setE0(e0);
    FloatArray reducedStrain, reducedStrainOld;

    double f, equivStrain, tempKappa, omega = 0.;

    this->initTempStatus(gp);

    FloatArray testStrainOld( status->giveStrainVector() );

    // substract stress independent part
    this->giveStressDependentPartOfStrainVector(reducedStrain, gp, totalStrain, atTime, VM_Total);

    // compute equivalent strain
    this->computeEquivalentStrain(equivStrain, reducedStrain, gp, atTime);

    // compute value of loading function if strainLevel crit apply
    f = equivStrain - status->giveKappa();

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

        this->computeDamageParam(omega, tempKappa, gp);
        if ( omega > 0 ) {
            status->setTempCrackFlag(1);
        }
    }


    FloatMatrix stiffnessMatrix;
    this->giveStiffnessMatrix(stiffnessMatrix, ElasticStiffness, gp, atTime);

    int rsize = StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() );


    answer.resize(rsize);
    answer.zero();
    for ( int i = 1; i <= rsize; i++ ) { // only diagonal terms matter
        answer.at(i) = stiffnessMatrix.at(i, i) * reducedStrain.at(i) * ( 1. - omega );
    }



    //Compute crack width
    double length = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveLength();


    double crackWidth;
    if ( gp->giveMaterialMode() == _2dLattice ) {
        crackWidth = omega * sqrt( pow(reducedStrain.at(1), 2.) + pow(reducedStrain.at(2), 2.) ) * length;
    } else {
        crackWidth = omega * sqrt( pow(reducedStrain.at(1), 2.) + pow(reducedStrain.at(2), 2.) + pow(reducedStrain.at(3), 2.) ) * length;
    }


    // todo - resolve and unify !!!

    double waterPressure = 0.;

    if ( gp->giveMaterialMode() == _2dLattice ) {
        IntArray coupledModels;

        //Calculate the the bio coefficient;
        double biot = 0.;
        if ( this->biotType == 0 ) {
            biot = this->biotCoefficient;
        } else if ( this->biotType == 1 ) {
            biot = computeBiot(omega, tempKappa, length);
        } else {
            OOFEM_ERROR("Unknown biot type\n");
        }

        answer.at(1) = answer.at(1) + biot * waterPressure;

        status->setBiotCoefficientInStatus(biot);
    } else { //3d_lattice
        //Read in fluid pressures from structural element if this is not a slave problem
        FloatArray pressures;
        if ( !domain->giveEngngModel()->giveMasterEngngModel() ) {
            static_cast< LatticeStructuralElement * >( gp->giveElement() )->givePressures(pressures);
        }

        waterPressure = 0.;
        for ( int i = 0; i < pressures.giveSize(); i++ ) {
            waterPressure += 1. / pressures.giveSize() * pressures.at(i + 1);
        }

        answer.at(1) += waterPressure;
    }


    double tempDissipation = status->giveDissipation();
    double tempDeltaDissipation;


    if ( gp->giveMaterialMode() == _2dLattice ) {
        tempDeltaDissipation = computeDeltaDissipation2d(omega, reducedStrain, gp, atTime);
    } else {
        tempDeltaDissipation = computeDeltaDissipation3d(omega, reducedStrain, gp, atTime);
    }

    tempDissipation += tempDeltaDissipation;




    //Set all temp values
    status->setTempDissipation(tempDissipation);
    status->setTempDeltaDissipation(tempDeltaDissipation);

    status->setTempEquivalentStrain(equivStrain);
    status->letTempStrainVectorBe(totalStrain);
    status->letTempReducedStrainBe(reducedStrain);
    status->letTempStressVectorBe(answer);
    status->setTempKappa(tempKappa);
    status->setTempDamage(omega);

    status->setTempNormalStress( answer.at(1) );
    status->setTempCrackWidth(crackWidth);

    return;
}



double
LatticeDamage :: computeBiot(double omega,
                             double kappa,
                             double le)
{
    double biot = 0.;

    if ( this->softeningType == 1 || this->softeningType == 3 ) {
        if ( omega == 0 ) {
            biot = this->biotCoefficient;
        } else if ( omega * kappa * le > 0 && omega * kappa * le < this->wf ) {
            biot = this->biotCoefficient + ( 1. - biotCoefficient ) * omega * kappa * le / this->wf;
        } else {
            biot = 1.;
        }
    } else {
        OOFEM_ERROR("Wrong stype for btype=1. Only linear and exponential softening considered so far\n");
    }

    return biot;
}



double
LatticeDamage :: computeDeltaDissipation2d(double omega,
                                           FloatArray &reducedStrain,
                                           GaussPoint *gp,
                                           TimeStep *tStep)
{
    LatticeDamageStatus *status = static_cast< LatticeDamageStatus * >( this->giveStatus(gp) );
    double length = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveLength();
    const double e0 = this->give(e0_ID, gp) * this->e0Mean;
    double eNormal = this->give(eNormal_ID, gp) * this->eNormalMean;

    const double eShear =  this->alphaOne * eNormal;
    const double eTorsion =  this->alphaTwo * eNormal;


    FloatArray reducedStrainOld;

    reducedStrainOld = status->giveReducedStrain();
    double omegaOld = status->giveDamage();
    double deltaOmega;

    FloatArray crackOpeningOld(3);
    crackOpeningOld.times(omegaOld);
    crackOpeningOld.times(length);
    FloatArray stressOld( status->giveStressVector() );
    FloatArray intermediateStrain(3);

    double tempDeltaDissipation = 0.;
    double deltaTempDeltaDissipation = 0.;

    double intermediateOmega = 0;
    FloatArray oldIntermediateStrain(3);
    oldIntermediateStrain = reducedStrainOld;
    double oldIntermediateOmega = omegaOld;
    deltaOmega = ( omega - omegaOld );
    double testDissipation =
        0.5 * length * ( pow( ( reducedStrain(0) + reducedStrainOld(0) ) / 2., 2. ) * eNormal +
                         pow( ( reducedStrain(1) + reducedStrainOld(1) ) / 2., 2. ) * eShear +
                         pow( ( reducedStrain(2) + reducedStrainOld(2) ) / 2., 2. ) * eTorsion ) * deltaOmega;


    double intervals = 0.;

    double referenceGf = 0;

    if ( softeningType == 1 ) {
        referenceGf = e0 * eNormal * this->wf / 2.;
    } else {   //This is for the exponential law. Should also implement it for the bilinear one.
        referenceGf = e0 * eNormal * this->wf;
    }

    if ( testDissipation / ( referenceGf ) > 0.01 ) {
        intervals = 1000. * testDissipation / referenceGf;
    } else {
        intervals = 1.;
    }

    if ( intervals > 1000. ) {
        intervals = 1000.;
    }

    double oldKappa = status->giveKappa();
    double f, equivStrain;
    if ( deltaOmega > 0 ) {
        for ( int k = 0; k < intervals; k++ ) {
            intermediateStrain(0) = reducedStrainOld(0) + ( k + 1 ) / intervals * ( reducedStrain(0) - reducedStrainOld(0) );
            intermediateStrain(1) = reducedStrainOld(1) + ( k + 1 ) / intervals * ( reducedStrain(1) - reducedStrainOld(1) );
            intermediateStrain(2) = reducedStrainOld(2) + ( k + 1 ) / intervals * ( reducedStrain(2) - reducedStrainOld(2) );
            this->computeEquivalentStrain(equivStrain, intermediateStrain, gp, tStep);
            f = equivStrain - oldKappa;
            if ( f > 0 ) {
                this->computeDamageParam(intermediateOmega, equivStrain, gp);
                deltaOmega = ( intermediateOmega - oldIntermediateOmega );
                deltaTempDeltaDissipation =
                    0.5 * length * ( pow( ( intermediateStrain(0) + oldIntermediateStrain(0) ) / 2., 2. ) * eNormal +
                                     pow( ( intermediateStrain(1) + oldIntermediateStrain(1) ) / 2., 2. ) * eShear +
                                     pow( ( intermediateStrain(2) + oldIntermediateStrain(2) ) / 2., 2. ) * eTorsion ) * deltaOmega;

                oldKappa = equivStrain;
                oldIntermediateOmega = intermediateOmega;
            } else {
                deltaTempDeltaDissipation = 0.;
            }

            tempDeltaDissipation += deltaTempDeltaDissipation;
            oldIntermediateStrain = intermediateStrain;
        }
    } else {
        tempDeltaDissipation = 0.;
    }

    if ( tempDeltaDissipation >= 2. * referenceGf ) {
        tempDeltaDissipation = 2. * referenceGf;
    }

    return tempDeltaDissipation;
}



double
LatticeDamage :: computeDeltaDissipation3d(double omega,
                                           FloatArray &reducedStrain,
                                           GaussPoint *gp,
                                           TimeStep *atTime)
{
    LatticeDamageStatus *status = static_cast< LatticeDamageStatus * >( this->giveStatus(gp) );
    double length = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveLength();
    const double e0 = this->give(e0_ID, gp) * this->e0Mean;
    double eNormal = this->give(eNormal_ID, gp) * this->eNormalMean;

    const double eShear =  this->alphaOne * eNormal;

    const double eTorsion =  this->alphaTwo * eNormal;

    FloatArray reducedStrainOld;

    reducedStrainOld = status->giveReducedStrain();
    double omegaOld = status->giveDamage();
    double deltaOmega;

    FloatArray crackOpeningOld(6);
    crackOpeningOld.times(omegaOld);
    crackOpeningOld.times(length);
    FloatArray stressOld( status->giveStressVector() );
    FloatArray intermediateStrain(6);

    double tempDeltaDissipation = 0.;
    double deltaTempDeltaDissipation = 0.;

    double intermediateOmega = 0;
    FloatArray oldIntermediateStrain(6);
    oldIntermediateStrain = reducedStrainOld;
    double oldIntermediateOmega = omegaOld;
    deltaOmega = ( omega - omegaOld );
    double testDissipation  = 0.5 * length * ( pow( ( reducedStrain.at(1) + reducedStrainOld.at(1) ) / 2., 2. ) * eNormal +
                                               pow( ( reducedStrain.at(2) + reducedStrainOld.at(2) ) / 2., 2. ) * eShear +
                                               pow( ( reducedStrain.at(3) + reducedStrainOld.at(3) ) / 2., 2. ) * eShear +
                                               pow( ( reducedStrain.at(4) + reducedStrainOld.at(4) ) / 2., 2. ) * eTorsion +
                                               pow( ( reducedStrain.at(5) + reducedStrainOld.at(5) ) / 2., 2. ) * eTorsion +
                                               pow( ( reducedStrain.at(6) + reducedStrainOld.at(6) ) / 2., 2. ) * eTorsion ) * deltaOmega;
    double intervals = 0.;

    double referenceGf = 0;

    if ( softeningType == 1 ) {
        referenceGf = e0 * eNormal * this->wf / 2.;
    } else {   //This is for the exponential law. Should also implement it for the bilinear one.
        referenceGf = e0 * eNormal * this->wf;
    }

    if ( testDissipation / ( referenceGf ) > 0.001 ) {
        intervals = 1000. * testDissipation / referenceGf;
    } else {
        intervals = 1.;
    }

    if ( intervals > 1000. ) {
        intervals = 1000.;
    }

    double oldKappa = status->giveKappa();
    double f, equivStrain;
    if ( deltaOmega > 0 ) {
        for ( int k = 0; k < intervals; k++ ) {
            intermediateStrain(0) = reducedStrainOld(0) + ( k + 1 ) / intervals * ( reducedStrain(0) - reducedStrainOld(0) );
            intermediateStrain(1) = reducedStrainOld(1) + ( k + 1 ) / intervals * ( reducedStrain(1) - reducedStrainOld(1) );
            intermediateStrain(2) = reducedStrainOld(2) + ( k + 1 ) / intervals * ( reducedStrain(2) - reducedStrainOld(2) );
            intermediateStrain(3) = reducedStrainOld(3) + ( k + 1 ) / intervals * ( reducedStrain(3) - reducedStrainOld(3) );
            intermediateStrain(4) = reducedStrainOld(4) + ( k + 1 ) / intervals * ( reducedStrain(4) - reducedStrainOld(4) );
            intermediateStrain(5) = reducedStrainOld(5) + ( k + 1 ) / intervals * ( reducedStrain(5) - reducedStrainOld(5) );

            this->computeEquivalentStrain(equivStrain, intermediateStrain, gp, atTime);
            f = equivStrain - oldKappa;
            if ( f > 0 ) {
                this->computeDamageParam(intermediateOmega, equivStrain, gp);
                deltaOmega = ( intermediateOmega - oldIntermediateOmega );
                deltaTempDeltaDissipation =
                    0.5 * length * ( pow( ( intermediateStrain(0) + oldIntermediateStrain(0) ) / 2., 2. ) * eNormal +
                                     pow( ( intermediateStrain(1) + oldIntermediateStrain(1) ) / 2., 2. ) * eShear +
                                     pow( ( intermediateStrain(2) + oldIntermediateStrain(2) ) / 2., 2. ) * eShear +
                                     pow( ( intermediateStrain(3) + oldIntermediateStrain(3) ) / 2., 2. ) * eTorsion +
                                     pow( ( intermediateStrain(4) + oldIntermediateStrain(4) ) / 2., 2. ) * eTorsion +
                                     pow( ( intermediateStrain(5) + oldIntermediateStrain(5) ) / 2., 2. ) * eTorsion ) * deltaOmega;

                oldKappa = equivStrain;
                oldIntermediateOmega = intermediateOmega;
            } else {
                deltaTempDeltaDissipation = 0.;
            }

            tempDeltaDissipation += deltaTempDeltaDissipation;
            oldIntermediateStrain = intermediateStrain;
        }
    } else {
        tempDeltaDissipation = 0.;
    }

    if ( tempDeltaDissipation >= 2. * referenceGf ) {
        tempDeltaDissipation = 2. * referenceGf;
    }

    return tempDeltaDissipation;
}


void
LatticeDamage :: give2dLatticeStiffMtrx(FloatMatrix &answer, MatResponseMode rmode, GaussPoint *gp, TimeStep *atTime)
{
    LatticeLinearElastic :: give2dLatticeStiffMtrx(answer, rmode, gp, atTime);

    if ( rmode == ElasticStiffness ) {
        return;
    } else if ( ( rmode == SecantStiffness ) || ( rmode == TangentStiffness ) ) {
        LatticeDamageStatus *status = static_cast< LatticeDamageStatus * >( this->giveStatus(gp) );
        double omega = status->giveTempDamage();

        if ( omega > 0.99999 ) {
            omega = 0.99999;
        }

        answer.times(1. - omega);

        return;
    } else {
        OOFEM_ERROR("Unsupported stiffness mode\n");
    }
}

void
LatticeDamage :: give3dLatticeStiffMtrx(FloatMatrix &answer, MatResponseMode rmode, GaussPoint *gp, TimeStep *atTime)
{
    LatticeLinearElastic :: give3dLatticeStiffMtrx(answer, rmode, gp, atTime);

    if ( rmode == ElasticStiffness ) {
        return;
    } else if ( ( rmode == SecantStiffness ) || ( rmode == TangentStiffness ) ) {
        LatticeDamageStatus *status = static_cast< LatticeDamageStatus * >( this->giveStatus(gp) );
        double omega = status->giveTempDamage();

        if ( omega > 0.99999 ) {
            omega = 0.99999;
        }

        answer.times(1. - omega);

        return;
    } else {
        OOFEM_ERROR("Unsupported stiffness mode\n");
    }
}


double
LatticeDamage :: give(int aProperty, GaussPoint *gp)
{
    double answer;
    if ( static_cast< LatticeDamageStatus* >( this->giveStatus(gp) )->_giveProperty(aProperty, answer) ) {
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
    if ( type == IST_CrackStatuses ) {
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
        answer.at(1) = status->giveNormalStress();
        return 1;
    }
    else if ( type == IST_CharacteristicLength ) {
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
{
    e0 = 0.;
    damage = tempDamage = 0.;
    equivStrain = tempEquivStrain = 0.;
    kappa = tempKappa = 0.;
}

void
LatticeDamageStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
// builds new crackMap
//
{
    LatticeMaterialStatus :: initTempStatus();
    this->tempKappa = this->kappa;
    this->tempEquivStrain = this->equivStrain;
    this->tempDamage = this->damage;
}

void
LatticeDamageStatus :: printOutputAt(FILE *file, TimeStep *tStep)
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

    // write a raw data
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

    // read raw data
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
}     // end namespace oofem
