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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#include "latticedamage2d.h"
#include "isolinearelasticmaterial.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "structuralcrosssection.h"
#include "engngm.h"
#include <math.h>
#include "datastream.h"
#include "contextioerr.h"

namespace oofem {
LatticeDamage2d :: LatticeDamage2d(int n, Domain *d) : StructuralMaterial(n, d), RandomMaterialExtensionInterface()
{}


LatticeDamage2d :: ~LatticeDamage2d()
{}

int
LatticeDamage2d :: hasMaterialModeCapability(MaterialMode mode)
{
    if ( mode == _2dLattice ) {
        return 1;
    }

    return 0;
}


IRResultType
LatticeDamage2d :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                             // Required by IR_GIVE_FIELD macro

    StructuralMaterial :: initializeFrom(ir);
    RandomMaterialExtensionInterface :: initializeFrom(ir);

    double value = 0.;
    IR_GIVE_FIELD(ir, value, IFT_IsotropicLinearElasticMaterial_talpha, "talpha"); // Macro
    propertyDictionary->add(tAlpha, value);

    IR_GIVE_FIELD(ir, eNormal, IFT_LatticeDamage2d_eNormal, "e"); // Macro

    //factor which relates the shear stiffness to the normal stiffness. Default is 1
    alphaOne = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, alphaOne, IFT_LatticeDamage2d_alphaOne, "a1"); // Macro
    eShear = alphaOne * eNormal;

    //Parameter which is used for the definition of the moment.
    alphaTwo = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, alphaTwo, IFT_LatticeDamage2d_alphaTwo, "a2"); // Macro

    eTorsion = alphaTwo * eNormal / 12.;

    softeningType = 0;
    IR_GIVE_FIELD(ir, softeningType, IFT_LatticeDamage2d_softeningType, "stype"); // Macro

    if ( softeningType == 1 || softeningType == 3 ) { //linear or exponential softening
        IR_GIVE_FIELD(ir, wf, IFT_LatticeDamage2d_wf, "wf"); // Macro
    } else if ( softeningType == 2 ) {      //bilinear softening
        IR_GIVE_FIELD(ir, wf, IFT_LatticeDamage2d_wf, "wf"); // Macro
        wfOne = 0.15 * wf;
        IR_GIVE_OPTIONAL_FIELD(ir, wfOne, IFT_LatticeDamage2d_wfOne, "wf1"); // Macro
    } else {
        _error("Unknown softening type");
    }

    localRandomType = 0; //Default: No local random field
    IR_GIVE_OPTIONAL_FIELD(ir, localRandomType, IFT_LatticeDamage2d_localrandomtype, "randomtype"); // Macro
    if ( localRandomType == 1 ) { //Gaussian random generator
        coefficientOfVariation = 0.;
        IR_GIVE_FIELD(ir, coefficientOfVariation, IFT_LatticeDamage2d_coefficientOfVariation, "cov"); // Macro
    }

    int equivType = 0;
    paramDuct = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, equivType, IFT_LatticeDamage2d_equivType, "equivtype"); // Macro
    e0Mean = 0;
    IR_GIVE_FIELD(ir, e0Mean, IFT_LatticeDamage2d_e0, "e0"); // Macro

    IR_GIVE_FIELD(ir, coh, IFT_LatticeDamage2d_coh, "coh"); // Macro
    IR_GIVE_FIELD(ir, ec, IFT_LatticeDamage2d_ec, "ec"); // Macro
    IR_GIVE_FIELD(ir, paramDuct, IFT_LatticeDamage2d_paramDuct, "duct"); // Macro

    return IRRT_OK;
}

void
LatticeDamage2d :: computeEquivalentStrain(double &tempEquivStrain, const FloatArray &strain, GaussPoint *gp, TimeStep *atTime)
{
    if ( strain.isEmpty() ) {
        tempEquivStrain = 0.;
        return;
    }

    LatticeDamage2dStatus *status = ( LatticeDamage2dStatus * ) this->giveStatus(gp);
    const double e0 = this->give(e0_ID, gp) * this->e0Mean;

    double paramA = 0.;
    double paramB = 0.;
    double paramC = 0.;
    double equivStrain;

    paramA = 0.5 * ( e0 + ec * e0 );
    paramB = ( coh * e0 ) / sqrt( 1. - pow( ( ec * e0 - e0 ) / ( e0 + ec * e0 ), 2. ) );
    paramC = 0.5 * ( this->ec * e0 - e0 );

    equivStrain = status->giveEquivalentStrain();
    tempEquivStrain =  sqrt( pow(this->alphaOne * strain.at(2) / paramB, 2.) + pow( ( strain.at(1) + paramC ) / paramA, 2. ) ) * paramA - paramC;

    return;
}

void
LatticeDamage2d :: computeDamageParam(double &omega, double tempKappa, const FloatArray &strain, GaussPoint *gp)
{
    LatticeDamage2dStatus *status = ( LatticeDamage2dStatus * ) this->giveStatus(gp);
    const double e0 = this->give(e0_ID, gp) * this->e0Mean;
    omega = 0.0;

    int nite = 0;
    double R, Lhs, Ft, help;

    if ( softeningType == 1 ) { //linear
        if ( tempKappa >= e0 && tempKappa < this->wf / status->giveLe() ) {
            //linear stress-crack opening relation
            //check if input parameter make sense
            if ( this->wf / status->giveLe() <= e0 ) {
                _error("e0>wf/Le \n Possible solutions: Increase fracture energy or reduce element size\n");
            }

            omega = ( 1. - e0 / tempKappa ) / ( 1. - e0 / ( this->wf / status->giveLe() ) );
        } else if ( tempKappa >= this->wf / status->giveLe() ) {
            omega = 1.;
        } else {
            omega = 0.;
        }
    } else if ( softeningType == 2 ) {      //bilinear softening
        double helpStrain = 0.3 * e0;

        //Check if input parameter make sense
        if ( e0 > wfOne / status->giveLe() ) {
            _error("parameter wf1 is too small");
        } else if ( wfOne / status->giveLe() >  this->wf / status->giveLe() ) {
            _error("parameter wf is too small");
        }

        //
        if ( tempKappa > e0 ) {
            omega = ( this->wfOne / status->giveLe() - helpStrain ) / ( this->wfOne / status->giveLe() - e0 ) * ( 1 - e0 / tempKappa );

            if ( omega * tempKappa > 0 && omega * tempKappa < this->wfOne / status->giveLe() ) {
                return;
            } else {
                omega = 1. - helpStrain / ( this->wf / status->giveLe() - this->wfOne / status->giveLe() ) * ( this->wf / status->giveLe() / tempKappa - 1 );

                if ( omega * tempKappa > this->wfOne / status->giveLe() && omega * tempKappa < this->wf / status->giveLe() ) {
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
            Ft = this->eNormal * e0;
            do {
                nite++;
                help = status->giveLe() * omega * tempKappa / this->wf;
                R = ( 1. - omega ) * eNormal * tempKappa - Ft *exp(-help);
                Lhs = eNormal * tempKappa - Ft *exp(-help) * status->giveLe() * tempKappa / this->wf;
                omega += R / Lhs;
                if ( nite > 40 ) {
                    _error("computeDamageParam: algorithm not converging");
                }
            } while ( fabs(R) >= 1.e-4 );

            if ( ( omega > 1.0 ) || ( omega < 0.0 ) ) {
                _error("computeDamageParam: internal error\n");
            }
        }
    } else {
        _error("Unknown softening type");
    }
}


void
LatticeDamage2d :: initDamaged(double kappa, FloatArray &strainVector, GaussPoint *gp)
{
    int i, indx = 1;
    double le;
    FloatArray principalStrains, crackPlaneNormal(3), fullstrain;
    FloatMatrix principalDir(3, 3);
    LatticeDamage2dStatus *status = ( LatticeDamage2dStatus * ) this->giveStatus(gp);

    //get the random variable from the status
    const double e0 = this->give(e0_ID, gp) * this->e0Mean;

    StructuralCrossSection *crossSection = ( StructuralCrossSection * ) gp->giveElement()->giveCrossSection();

    crossSection->giveFullCharacteristicVector(fullstrain, gp, strainVector);

    if ( ( kappa > e0 ) && ( status->giveDamage() == 0. ) ) {
        this->computePrincipalValDir(principalStrains, principalDir, fullstrain, principal_strain);
        // finfd index of max positive principal strain
        for ( i = 2; i <= 3; i++ ) {
            if ( principalStrains.at(i) > principalStrains.at(indx) ) {
                indx = i;
            }
        }

        for ( i = 1; i <= 3; i++ ) {
            crackPlaneNormal.at(i) = principalDir.at(i, indx);
        }

        le = gp->giveElement()->giveCharacteristicLenght(gp, crackPlaneNormal);
        // remember le in cooresponding status
        status->setLe(le);
    }
}


double
LatticeDamage2d :: computeDuctilityMeasure(const FloatArray &strain, GaussPoint *gp)
{
    double ductilityMeasure = 1.;
    if ( strain.at(1) < 0 ) {
        ductilityMeasure = sqrt( pow(strain.at(1), 2.) * pow(this->paramDuct, 2.) / pow(this->ec, 2.) + pow(strain.at(2), 2.) ) /
                           sqrt( pow(strain.at(1), 2.) + pow(strain.at(2), 2.) );
    }

    return ductilityMeasure;
}

MaterialStatus *
LatticeDamage2d :: CreateStatus(GaussPoint *gp) const
{
    LatticeDamage2dStatus *answer = new LatticeDamage2dStatus(1, LatticeDamage2d :: domain, gp);
    return answer;
}

MaterialStatus *
LatticeDamage2d :: giveStatus(GaussPoint *gp) const
{
    MaterialStatus *status;

    status = gp->giveMaterialStatus();
    if ( status == NULL ) {
        // create a new one
        status = this->CreateStatus(gp);

        if ( status != NULL ) {
            gp->setMaterialStatus(status);
            this->_generateStatusVariables(gp);
        }
    }

    return status;
}



void
LatticeDamage2d :: giveRealStressVector(FloatArray &answer,
                                        MatResponseForm form,
                                        GaussPoint *gp,
                                        const FloatArray &totalStrain,
                                        TimeStep *atTime)
//
// returns real stress vector in 3d stress space of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
{
    LatticeDamage2dStatus *status = ( LatticeDamage2dStatus * ) this->giveStatus(gp);

    if ( gp->giveElement()->giveNumber() == 3299 ) {
        printf("Debug\n");
    }

    const double e0 = this->give(e0_ID, gp) * this->e0Mean;

    FloatArray strainVector, reducedStrain, reducedStrainOld;

    double f, equivStrain, tempKappa, omega = 0.;

    this->initGpForNewStep(gp);
    reducedStrain = totalStrain;

    FloatArray testStrainOld( status->giveStrainVector() );

    // subtract stress independent part
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
        this->initDamaged(tempKappa, reducedStrain, gp);

        // evaluate damage parameter
        this->computeDamageParam(omega, tempKappa, reducedStrain, gp);
        if ( omega > 0 ) {
            status->setTempCrackFlag(1);
        }
    }

    //Compute the stress
    answer.resize(3);
    answer.zero();
    answer.at(1) = ( 1. - omega ) * eNormal * reducedStrain.at(1);
    answer.at(2) = ( 1. - omega ) * eShear * reducedStrain.at(2);
    answer.at(3) = ( 1. - omega ) * eTorsion * reducedStrain.at(3);

    //Compute dissipation for post-processing reasons

    FloatArray crackPlaneNormal(3); //Dummy floatArray so that I can use giveCharacteristiclenght. This needs to be improved
    double length = gp->giveElement()->giveCharacteristicLenght(gp, crackPlaneNormal);
    status->giveReducedStrain(reducedStrainOld);
    double omegaOld = status->giveDamage();
    double deltaOmega;

    FloatArray crackOpeningOld(3);
    crackOpeningOld.times(omegaOld);
    crackOpeningOld.times(length);
    FloatArray stressOld( status->giveStressVector() );
    FloatArray intermediateStrain(3);
    double tempDissipation = status->giveDissipation();
    double tempDeltaDissipation = 0.;
    double deltaTempDeltaDissipation = 0.;

    //Need to do this only if f is greater than zero. Otherwise dissipation is zero.
    //Loop over 100 steps to get the excact dissipation

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

    if ( intervals > 1000 ) {
        intervals = 1000;
    }

    double oldKappa = status->giveKappa();

    if ( deltaOmega > 0 ) {
        for ( int k = 0; k < intervals; k++ ) {
            intermediateStrain(0) = reducedStrainOld(0) + ( k + 1 ) / intervals * ( reducedStrain(0) - reducedStrainOld(0) );
            intermediateStrain(1) = reducedStrainOld(1) + ( k + 1 ) / intervals * ( reducedStrain(1) - reducedStrainOld(1) );
            intermediateStrain(2) = reducedStrainOld(2) + ( k + 1 ) / intervals * ( reducedStrain(2) - reducedStrainOld(2) );
            this->computeEquivalentStrain(equivStrain, intermediateStrain, gp, atTime);
            f = equivStrain - oldKappa;
            if ( f > 0 ) {
                this->computeDamageParam(intermediateOmega, equivStrain, intermediateStrain, gp);
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

    tempDissipation += tempDeltaDissipation;

    if ( tempDissipation >= 2. * referenceGf ) {
        tempDissipation = 2. * referenceGf;
    }

    //Set all the temp values
    status->setTempDissipation(tempDissipation);
    status->setTempDeltaDissipation(tempDeltaDissipation);

    status->setTempEquivalentStrain(equivStrain);
    status->letTempStrainVectorBe(totalStrain);
    status->letTempReducedStrainBe(reducedStrain);
    status->letTempStressVectorBe(answer);
    status->setTempKappa(tempKappa);
    status->setTempDamage(omega);

    double le = status->giveLe();
    double crackWidth = omega * sqrt( pow(reducedStrain.at(1), 2.) + pow(reducedStrain.at(2), 2.) + pow(reducedStrain.at(3), 2.) ) * le;

    status->setTempCrackWidth(crackWidth);

    return;
}

void
LatticeDamage2d :: giveStressStrainMask(IntArray &answer, MatResponseForm form,
                                        MaterialMode mmode) const
{
    int i;

    if ( mmode == _2dLattice ) {
        answer.resize(3);
        for ( i = 1; i <= 3; i++ ) {
            answer.at(i) = i;
        }
    } else {
        _error("Error: Unknown materialmode");
    }
}

void LatticeDamage2d :: giveRandomParameters(FloatArray &param)
{
    param.resize(3);
    param.at(1) = localRandomType;

    if ( localRandomType == 1 ) { //Gaussian
        param.at(2) = coefficientOfVariation;
    } else {
        _error("Error: Unknown local random type:\n randomtype 1 = Gaussian\n");
    }

    return;
}


Interface *
LatticeDamage2d :: giveInterface(InterfaceType type)
{
    return NULL;
}


int
LatticeDamage2d :: giveStressStrainComponentIndOf(MatResponseForm form, MaterialMode mmode, int ind)
{
    if ( mmode == _2dLattice ) {
        return ind;
    } else {
        _error("Unknown material mode\n");
    }

    return 0;
}

int
LatticeDamage2d :: giveSizeOfReducedStressStrainVector(MaterialMode mode)
{
    switch ( mode ) {
    case _2dLattice:
        return 3;

    default:
        _error("Unknown material mode \n");
        return 0;
    }
}

void
LatticeDamage2d :: giveCharacteristicMatrix(FloatMatrix &answer,
                                            MatResponseForm form, MatResponseMode rMode,
                                            GaussPoint *gp, TimeStep *atTime)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _2dLattice:

        if ( rMode == ElasticStiffness ) {
            this->giveElasticStiffnessMatrix(answer, gp, atTime);
        } else if ( rMode == SecantStiffness ) {
            this->giveSecantStiffnessMatrix(answer, gp, atTime);
        } else if ( rMode == TangentStiffness ) {
            this->giveSecantStiffnessMatrix(answer, gp, atTime);
        } else {
            _error("Unsupported stiffness mode\n");
        }

        break;
    default:
        _error("unknown material mode");
    }
}


void
LatticeDamage2d :: giveSecantStiffnessMatrix(FloatMatrix &answer,
                                             GaussPoint *gp,
                                             TimeStep *atTime)
{
    LatticeDamage2dStatus *status = ( LatticeDamage2dStatus * ) this->giveStatus(gp);

    double omega = status->giveTempDamage();

    if ( omega > 1. - 1.e-9 ) { //To avoid convergence problems
        omega = 1. - 1.e-9;
    }

    /* Returns elastic moduli in reduced stress-strain space*/
    answer.resize(3, 3);
    answer.zero();
    answer.at(1, 1) = ( 1 - omega ) * eNormal;
    answer.at(2, 2) = ( 1 - omega ) * eShear;
    answer.at(3, 3) = ( 1 - omega ) * eTorsion;
    return;
}


void
LatticeDamage2d :: giveTangentStiffnessMatrix(FloatMatrix &answer,
                                              GaussPoint *gp,
                                              TimeStep *atTime)
{
    _error("tangent stiffness not implemented\n");
    return;
}



void
LatticeDamage2d :: giveElasticStiffnessMatrix(FloatMatrix &answer,
                                              GaussPoint *gp,
                                              TimeStep *atTime)
{
    /* Returns elastic moduli in reduced stress-strain space*/
    answer.resize(3, 3);
    answer.zero();
    answer.at(1, 1) = eNormal;
    answer.at(2, 2) = eShear;
    answer.at(3, 3) = eTorsion;
    return;
}



void
LatticeDamage2d :: giveReducedCharacteristicVector(FloatArray &answer, GaussPoint *gp,
                                                   const FloatArray &charVector3d)
//
// returns reduced stressVector or strainVector from full 3d vector reduced
// to vector required by gp->giveStressStrainMode()
//
{
    MaterialMode mode = gp->giveMaterialMode();

    if ( mode == _2dLattice ) {
        answer = charVector3d;
        return;
    } else {
        _error("Unknown material mode\n");
    }
}



void
LatticeDamage2d :: giveThermalDilatationVector(FloatArray &answer,
                                               GaussPoint *gp,  TimeStep *tStep)
{
    answer.resize(3);
    answer.zero();
    answer.at(1) = this->give(tAlpha, gp);

    return;
}


double
LatticeDamage2d :: give(int aProperty, GaussPoint *gp)
{
    double answer;
    if ( RandomMaterialExtensionInterface :: give(aProperty, gp, answer) ) {
        return answer;
    } else if ( aProperty == e0_ID ) {
        return 1.;
    } else if ( aProperty == ef_ID ) {
        return 1.;
    } else {
        return LatticeDamage2d :: give(aProperty, gp);
    }
}


void
LatticeDamage2d :: giveFullCharacteristicVector(FloatArray &answer,
                                                GaussPoint *gp,
                                                const FloatArray &strainVector)

{
    MaterialMode mode = gp->giveMaterialMode();
    if ( mode == _2dLattice ) {
        answer = strainVector;
        return;
    } else {
        _error("Unknown material model\n");
    }
}

int
LatticeDamage2d::giveIPValue(FloatArray& answer,
                         GaussPoint* gp,
                         InternalStateType type,
                         TimeStep* atTime)
{
    LatticeDamage2dStatus *status = (LatticeDamage2dStatus*) this -> giveStatus (gp);
    if (type == IST_CrackStatuses){
        answer.resize(1);
        answer(0) = status->giveCrackFlag();
        return 1;
    }
    else return StructuralMaterial::giveIPValue(answer, gp, type, atTime);
}

int
LatticeDamage2d::giveIPValueSize(InternalStateType type,
                              GaussPoint* gp)
{
    if (type == IST_CrackStatuses) return 1;
    else return StructuralMaterial::giveIPValueSize (type, gp);
}

int
LatticeDamage2d::giveIntVarCompFullIndx(IntArray& answer,
                                    InternalStateType type, MaterialMode mmode)
{
    if(type == IST_CrackStatuses) {
        answer.resize (1);
        answer.at(1) = 1;
        return 1;
    } else
        return StructuralMaterial::giveIntVarCompFullIndx (answer, type, mmode);
}

InternalStateValueType
LatticeDamage2d::giveIPValueType (InternalStateType type)
{
    if(type==IST_CrackStatuses) {
        return ISVT_SCALAR;
    }
    else{
        return StructuralMaterial::giveIPValueType (type);
    }
}

LatticeDamage2dStatus :: LatticeDamage2dStatus(int n, Domain *d, GaussPoint *g) :
    LatticeMaterialStatus(n, d, g), RandomMaterialStatusExtensionInterface(), reducedStrain(3), tempReducedStrain(3)
{
    le = 0.0;
    crack_flag = temp_crack_flag = 0;
    crackWidth = tempCrackWidth = 0;
    e0 = 0.;
    damage = tempDamage = 0.;
    equivStrain = tempEquivStrain = 0.;
    kappa = tempKappa = 0.;
    dissipation = tempDissipation = 0.;
    deltaDissipation = tempDeltaDissipation = 0.;
}

void
LatticeDamage2dStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
// builds new crackMap
//
{
    StructuralMaterialStatus :: initTempStatus();
    this->tempReducedStrain = this->reducedStrain;
    this->tempKappa = this->kappa;
    this->tempEquivStrain = this->equivStrain;
    this->tempDamage = this->damage;
    this->tempDissipation = this->dissipation;
    this->tempDeltaDissipation = this->deltaDissipation;
    this->temp_crack_flag = this->crack_flag;
    this->tempCrackWidth = this->crackWidth;
}

void
LatticeDamage2dStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { ");
    fprintf(file, "reduced strains ");
    int rSize = reducedStrain.giveSize();
    for ( int k = 1; k <= rSize; k++ ) {
        fprintf( file, "% .4e ", reducedStrain.at(k) );
    }

    fprintf(file, "kappa %f, equivStrain %f, damage %f, dissipation %f, deltaDissipation %f, e0 %f, crack_flag %d ", this->kappa, this->equivStrain, this->damage, this->dissipation, this->deltaDissipation, this->e0, this->crack_flag);
    fprintf(file, "}\n");
}


void
LatticeDamage2dStatus :: updateYourself(TimeStep *atTime)
{
    StructuralMaterialStatus :: updateYourself(atTime);

    this->reducedStrain = this->tempReducedStrain;
    this->kappa = this->tempKappa;
    this->equivStrain = this->tempEquivStrain;
    this->damage = this->tempDamage;
    this->dissipation = this->tempDissipation;
    this->deltaDissipation = this->tempDeltaDissipation;
    this->crack_flag = this->temp_crack_flag;
    this->crackWidth = this->tempCrackWidth;
}

Interface *
LatticeDamage2dStatus :: giveInterface(InterfaceType type)
{
    if ( type == RandomMaterialStatusExtensionInterfaceType ) {
        return ( RandomMaterialStatusExtensionInterface * ) this;
    } else {
        return NULL;
    }
}


contextIOResultType
LatticeDamage2dStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves full information stored in this Status
// no temp variables stored
//
{
    contextIOResultType iores;
    // save parent class status
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }



    // write a raw data
    if ( ( iores = reducedStrain.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write a raw data
    if ( !stream->write(& kappa, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& equivStrain, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& damage, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& le, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& dissipation, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& deltaDissipation, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& crack_flag, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& e0, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}

void
LatticeDamage2dStatus :: setTempCrackFlag(int val)
{
    temp_crack_flag = val;
}



void
LatticeDamage2dStatus :: setTempCrackWidth(double val)
{
    tempCrackWidth = val;
}


void
LatticeDamage2dStatus :: setVariableInStatus(double variable) {
    e0 = variable;
}

int
LatticeDamage2dStatus :: giveCrackFlag()
{
    if ( crack_flag != 0 ) {}

    return crack_flag;
}


contextIOResultType
LatticeDamage2dStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restores full information stored in stream to this Status
//
{
    contextIOResultType iores;
    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = reducedStrain.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
    if ( !stream->read(& kappa, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& equivStrain, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& damage, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& le, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& dissipation, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& deltaDissipation, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& crack_flag, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& e0, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}
}     // end namespace oofem
