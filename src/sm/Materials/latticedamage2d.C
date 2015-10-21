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
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "latticedamage2d.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "datastream.h"
#include "contextioerr.h"
#include "mathfem.h"
#include "Elements/latticestructuralelement.h"
#include "Materials/isolinearelasticmaterial.h"
#include "staggeredproblem.h"
#include "classfactory.h"
#ifdef __TM_MODULE
 #include "../tm/latticetransportelement.h"
#endif

#include <cstdlib>

namespace oofem {
REGISTER_Material(LatticeDamage2d);

LatticeDamage2d :: LatticeDamage2d(int n, Domain *d) : StructuralMaterial(n, d), RandomMaterialExtensionInterface()
{ }


LatticeDamage2d :: ~LatticeDamage2d()
{ }

int
LatticeDamage2d :: hasMaterialModeCapability(MaterialMode mode)
{
    return mode == _2dLattice;
}


IRResultType
LatticeDamage2d :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                             // Required by IR_GIVE_FIELD macro

    result = StructuralMaterial :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;
    
    result = RandomMaterialExtensionInterface :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    double value = 0.;
    IR_GIVE_FIELD(ir, value, _IFT_IsotropicLinearElasticMaterial_talpha);
    propertyDictionary.add(tAlpha, value);

    IR_GIVE_FIELD(ir, eNormal, _IFT_LatticeDamage2d_eNormal);

    cAlpha = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, cAlpha, _IFT_LatticeDamage2d_calpha);


    //factor which relates the shear stiffness to the normal stiffness. Default is 1
    alphaOne = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, alphaOne, _IFT_LatticeDamage2d_alphaOne);
    eShear = alphaOne * eNormal;

    //Parameter which is used for the definition of the moment.
    alphaTwo = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, alphaTwo, _IFT_LatticeDamage2d_alphaTwo);

    eTorsion = alphaTwo * eNormal;

    localRandomType = 0; //Default: No local random field
    IR_GIVE_OPTIONAL_FIELD(ir, localRandomType, _IFT_LatticeDamage2d_localrandomtype);
    if ( localRandomType == 1 ) { //Gaussian random generator
        coefficientOfVariation = 0.;
        IR_GIVE_FIELD(ir, coefficientOfVariation, _IFT_LatticeDamage2d_coefficientOfVariation);
    }

    e0Mean = 0;
    IR_GIVE_FIELD(ir, e0Mean, _IFT_LatticeDamage2d_e0Mean);

    IR_GIVE_FIELD(ir, coh, _IFT_LatticeDamage2d_coh);
    IR_GIVE_FIELD(ir, ec, _IFT_LatticeDamage2d_ec);

    softeningType = 0;
    IR_GIVE_FIELD(ir, softeningType, _IFT_LatticeDamage2d_softeningType);

    if ( softeningType == 1 || softeningType == 3 ) { //linear or exponential softening
        IR_GIVE_FIELD(ir, wf, _IFT_LatticeDamage2d_wf);
    } else if ( softeningType == 2 ) {      //bilinear softening
        IR_GIVE_FIELD(ir, wf, _IFT_LatticeDamage2d_wf);
        wfOne = 0.15 * wf;
        IR_GIVE_OPTIONAL_FIELD(ir, wfOne, _IFT_LatticeDamage2d_wfOne);
        e0OneMean = 0.3 * e0Mean;
        IR_GIVE_OPTIONAL_FIELD(ir, e0OneMean, _IFT_LatticeDamage2d_e0OneMean);
    } else {
        OOFEM_WARNING("Unknown softening type");
        return IRRT_BAD_FORMAT;
    }

    this->biotCoefficient = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, this->biotCoefficient, _IFT_LatticeDamage2d_bio);

    this->biotType = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->biotType, _IFT_LatticeDamage2d_btype);

    return IRRT_OK;
}

void
LatticeDamage2d :: computeEquivalentStrain(double &tempEquivStrain, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    if ( strain.isEmpty() ) {
        tempEquivStrain = 0.;
        return;
    }

    //LatticeDamage2dStatus *status = static_cast< LatticeDamage2dStatus * >( this->giveStatus(gp) );
    const double e0 = this->give(e0_ID, gp) * this->e0Mean;

    double paramA, paramB, paramC;

    paramA = 0.5 * ( e0 + ec * e0 );
    paramB = ( coh * e0 ) / sqrt( 1. - pow( ( ec * e0 - e0 ) / ( e0 + ec * e0 ), 2. ) );
    paramC = 0.5 * ( this->ec * e0 - e0 );

    //double equivStrain = status->giveEquivalentStrain();
    tempEquivStrain =  sqrt( pow(this->alphaOne * strain.at(2) / paramB, 2.) + pow( ( strain.at(1) + paramC ) / paramA, 2. ) ) * paramA - paramC;
}

void
LatticeDamage2d :: computeDamageParam(double &omega, double tempKappa, const FloatArray &strain, GaussPoint *gp)
{
    double le = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveLength();
    const double e0 = this->give(e0_ID, gp) * this->e0Mean;
    omega = 0.0;

    int nite = 0;
    double R, Lhs, Ft, help;

    if ( softeningType == 1 ) { //linear
        if ( tempKappa >= e0 && tempKappa < this->wf / le ) {
            //linear stress-crack opening relation
            //check if input parameter make sense
            if ( this->wf / le <= e0 ) {
                OOFEM_ERROR("e0>wf/Le \n Possible solutions: Increase fracture energy or reduce element size");
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
            Ft = this->eNormal * e0;
            do {
                nite++;
                help = le * omega * tempKappa / this->wf;
                R = ( 1. - omega ) * eNormal * tempKappa - Ft *exp(-help);
                Lhs = eNormal * tempKappa - Ft *exp(-help) * le * tempKappa / this->wf;
                omega += R / Lhs;
                if ( nite > 40 ) {
                    OOFEM_ERROR("algorithm not converging");
                }
            } while ( fabs(R) >= 1.e-4 );

            if ( ( omega > 1.0 ) || ( omega < 0.0 ) ) {
                OOFEM_ERROR("internal error");
            }
        }
    } else {
        OOFEM_ERROR("Unknown softening type");
    }
}


void
LatticeDamage2d :: computeStressIndependentStrainVector(FloatArray &answer,
                                                        GaussPoint *gp,
                                                        TimeStep *tStep,
                                                        ValueModeType mode)
{
    FloatArray et;
    LatticeStructuralElement *lselem = static_cast< LatticeStructuralElement * >( gp->giveElement() );

    answer.resize(3);
    answer.zero();

    if ( tStep->giveIntrinsicTime() < this->castingTime ) {
        return;
    }

    if ( lselem ) {
        lselem->computeResultingIPTemperatureAt(et, tStep, gp, mode);
    }

    if ( et.giveSize() == 0 ) {
        answer.clear();
        return;
    }

    if ( et.giveSize() < 1 ) {
        OOFEM_ERROR("Bad format of TemperatureLoad");
        exit(1);
    }

    if ( mode == VM_Total ) {
        // compute temperature difference
        double deltaTemperature = et.at(1) - this->referenceTemperature;
        answer.at(1) = this->give(tAlpha, gp) * deltaTemperature;
    } else {
        answer.at(1) = this->give(tAlpha, gp) * et.at(1);
    }

    double length = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveLength();

    answer.at(1) += this->cAlpha * et.at(1) / length;
    return;
}

MaterialStatus *
LatticeDamage2d :: CreateStatus(GaussPoint *gp) const
{
    return new LatticeDamage2dStatus(1, LatticeDamage2d :: domain, gp);
}

MaterialStatus *
LatticeDamage2d :: giveStatus(GaussPoint *gp) const
{
    MaterialStatus *status = static_cast< MaterialStatus * >( gp->giveMaterialStatus() );
    if ( status == NULL ) {
        // create a new one
        status = this->CreateStatus(gp);

        if ( status != NULL ) {
            gp->setMaterialStatus( status, this->giveNumber() );
            this->_generateStatusVariables(gp);
        }
    }

    return status;
}


void
LatticeDamage2d :: giveRealStressVector(FloatArray &answer,
                                        GaussPoint *gp,
                                        const FloatArray &totalStrain,
                                        TimeStep *tStep)
//
// returns real stress vector in 3d stress space of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
{
    LatticeDamage2dStatus *status = static_cast< LatticeDamage2dStatus * >( this->giveStatus(gp) );

    const double e0 = this->give(e0_ID, gp) * this->e0Mean;
    status->setE0(e0);

    FloatArray reducedStrain;

    double f, equivStrain, tempKappa, omega = 0.;

    this->initTempStatus(gp);
    reducedStrain = totalStrain;

    // subtract stress independent part
    this->giveStressDependentPartOfStrainVector(reducedStrain, gp, totalStrain, tStep, VM_Total);

    // compute equivalent strain
    this->computeEquivalentStrain(equivStrain, reducedStrain, gp, tStep);


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

    //Add the water pressure to the normal component of the stress
    IntArray coupledModels;
    double waterPressure = 0.;

#ifdef __TM_MODULE
    if ( domain->giveEngngModel()->giveMasterEngngModel() ) {
        ( static_cast< StaggeredProblem * >( domain->giveEngngModel()->giveMasterEngngModel() ) )->giveCoupledModels(coupledModels);
        int couplingFlag = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveCouplingFlag();

        if ( couplingFlag == 1 && coupledModels.at(2) != 0 && !tStep->isTheFirstStep() ) {
            IntArray couplingNumbers;

            static_cast< LatticeStructuralElement * >( gp->giveElement() ) ->giveCouplingNumbers(couplingNumbers);
            LatticeTransportElement *coupledElement;
            coupledElement  = static_cast< LatticeTransportElement * >(domain->giveEngngModel()->giveMasterEngngModel()->giveSlaveProblem( coupledModels.at(2) )->giveDomain(1)->giveElement(couplingNumbers.at(1)) );
            waterPressure = coupledElement->givePressure();
        }
    }
#endif
    
    //Compute dissipation
    double tempDissipation = status->giveDissipation();
    double tempDeltaDissipation = computeDeltaDissipation(omega, reducedStrain, gp, tStep);

    tempDissipation += tempDeltaDissipation;

    //Compute crack width
    double le = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveLength();
    double crackWidth = omega * sqrt( pow(reducedStrain.at(1), 2.) + pow(reducedStrain.at(2), 2.) ) * le;

    //Calculate the the bio coefficient;
    double biot = 0.;
    if ( this->biotType == 0 ) {
        biot = this->biotCoefficient;
    } else if ( this->biotType == 1 ) {
        biot = computeBiot(omega, tempKappa, le);
    } else {
      OOFEM_ERROR("Unknown biot type\n");
    }
    
    answer.at(1) = answer.at(1) + biot * waterPressure;


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
    
    status->setBiotCoefficientInStatus(biot);
}

double
LatticeDamage2d :: computeBiot(double omega,
                               double kappa,
                               double le)
{
  
    double biot = 0.;
    
    if ( this->softeningType == 1 || this->softeningType == 3 ) {
        if ( omega == 0 ) {
            biot = this->biotCoefficient;
        } else if ( omega*kappa*le > 0 && omega*kappa*le < this->wf ) {
            biot = this->biotCoefficient + (1.-biotCoefficient)*omega*kappa*le/this->wf;
        } else {
            biot = 1.;
        }
    } else {
        OOFEM_ERROR("Wrong stype for btype=1. Only linear and exponential softening considered so far\n");
    } 
    
    return biot;
}


double
LatticeDamage2d :: computeDeltaDissipation(double omega,
                                           FloatArray &reducedStrain,
                                           GaussPoint *gp,
                                           TimeStep *tStep)
{
    LatticeDamage2dStatus *status = static_cast< LatticeDamage2dStatus * >( this->giveStatus(gp) );
    double length = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveLength();
    const double e0 = this->give(e0_ID, gp) * this->e0Mean;

    FloatArray reducedStrainOld;

    reducedStrainOld = status->giveReducedStrain();
    double omegaOld = status->giveDamage();
    double deltaOmega;

    FloatArray crackOpeningOld(3);
    crackOpeningOld.times(omegaOld);
    crackOpeningOld.times(length);
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

    if ( tempDeltaDissipation >= 2. * referenceGf ) {
        tempDeltaDissipation = 2. * referenceGf;
    }

    return tempDeltaDissipation;
}


void LatticeDamage2d :: giveRandomParameters(FloatArray &param)
{
    param.resize(3);
    param.at(1) = localRandomType;

    if ( localRandomType == 1 ) { //Gaussian
        param.at(2) = coefficientOfVariation;
    } else {
        OOFEM_ERROR("Unknown local random type:\n randomtype 1 = Gaussian");
    }
}


Interface *
LatticeDamage2d :: giveInterface(InterfaceType type)
{
    return NULL;
}


void
LatticeDamage2d :: giveStiffnessMatrix(FloatMatrix &answer,
                                       MatResponseMode rMode,
                                       GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _2dLattice:

        if ( rMode == ElasticStiffness ) {
            this->giveElasticStiffnessMatrix(answer, gp, tStep);
        } else if ( rMode == SecantStiffness ) {
            this->giveSecantStiffnessMatrix(answer, gp, tStep);
        } else if ( rMode == TangentStiffness ) {
            this->giveSecantStiffnessMatrix(answer, gp, tStep);
        } else {
            OOFEM_ERROR("Unsupported stiffness mode");
        }

        break;
    default:
        OOFEM_ERROR("unknown material mode");
    }
}


void
LatticeDamage2d :: giveSecantStiffnessMatrix(FloatMatrix &answer,
                                             GaussPoint *gp,
                                             TimeStep *tStep)
{
    LatticeDamage2dStatus *status = static_cast< LatticeDamage2dStatus * >( this->giveStatus(gp) );

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
}

void
LatticeDamage2d :: giveElasticStiffnessMatrix(FloatMatrix &answer,
                                              GaussPoint *gp,
                                              TimeStep *tStep)
{
    /* Returns elastic moduli in reduced stress-strain space*/
    answer.resize(3, 3);
    answer.zero();
    answer.at(1, 1) = eNormal;
    answer.at(2, 2) = eShear;
    answer.at(3, 3) = eTorsion;
}


void
LatticeDamage2d :: giveThermalDilatationVector(FloatArray &answer,
                                               GaussPoint *gp,  TimeStep *tStep)
{
    answer.resize(3);
    answer.zero();
    answer.at(1) = this->give(tAlpha, gp);
}


double
LatticeDamage2d :: give(int aProperty, GaussPoint *gp)
{
    double answer;
    if ( static_cast< LatticeDamage2dStatus * >( this->giveStatus(gp) )->_giveProperty(aProperty, answer) ) {
        return answer;
    } else if ( aProperty == e0_ID ) {
        return 1.;
    } else if ( aProperty == ef_ID ) {
        return 1.;
    } else {
        return StructuralMaterial :: give(aProperty, gp);
    }
}


int
LatticeDamage2d :: giveIPValue(FloatArray &answer,
                               GaussPoint *gp,
                               InternalStateType type,
                               TimeStep *tStep)
{
    LatticeDamage2dStatus *status = static_cast< LatticeDamage2dStatus * >( this->giveStatus(gp) );
    if ( type == IST_CrackStatuses ) {
        answer.resize(1);
        answer.at(1) = status->giveCrackFlag();
        return 1;
    } else if ( type == IST_DamageScalar ) {
        answer.resize(1);
        answer.at(1) = status->giveDamage();
        return 1;
    } else if ( type == IST_DamageTensor ) {
        answer.resize(6);
        answer.at(1) = answer.at(2) = answer.at(3) = status->giveDamage();
        return 1;
    } else if ( type == IST_DissWork ) {
        answer.resize(1);
        answer.at(1) = status->giveDissipation();
        return 1;
    } else if ( type == IST_DeltaDissWork ) {
        answer.resize(1);
        answer.at(1) = status->giveDeltaDissipation();
        return 1;
    } else if ( type == IST_CrackWidth ) {
        answer.resize(1);
        answer.at(1) = status->giveCrackWidth();
        return 1;
    } else {
        return StructuralMaterial :: giveIPValue(answer, gp, type, tStep);
    }
}

LatticeDamage2dStatus :: LatticeDamage2dStatus(int n, Domain *d, GaussPoint *g) :
    LatticeMaterialStatus(n, d, g), RandomMaterialStatusExtensionInterface(), reducedStrain(3), tempReducedStrain(3)
{
    le = 0.0;
    crack_flag = temp_crack_flag = 0;
    crackWidth = tempCrackWidth = 0;
    normalStress = tempNormalStress = 0;
    oldNormalStress = 0;
    e0 = 0.;
    damage = tempDamage = 0.;
    equivStrain = tempEquivStrain = 0.;
    kappa = tempKappa = 0.;
    dissipation = tempDissipation = 0.;
    deltaDissipation = tempDeltaDissipation = 0.;
    biot = 0.;
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
    this->tempNormalStress = this->normalStress;    
    this->updateFlag = 0;
}

void
LatticeDamage2dStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { ");
    fprintf(file, "reduced strains ");
    for ( auto &val : reducedStrain ) {
        fprintf( file, "%.4e ", val );
    }

    fprintf(file, "kappa %f, equivStrain %f, damage %f, dissipation %f, deltaDissipation %f, e0 %f, crack_flag %d, crackWidth %.8e, biotCoeff %f", this->kappa, this->equivStrain, this->damage, this->dissipation, this->deltaDissipation, this->e0, this->crack_flag, this->crackWidth, this->biot);
    fprintf(file, "}\n");
}


void
LatticeDamage2dStatus :: updateYourself(TimeStep *tStep)
{
    StructuralMaterialStatus :: updateYourself(tStep);

    this->reducedStrain = this->tempReducedStrain;
    this->kappa = this->tempKappa;
    this->equivStrain = this->tempEquivStrain;
    this->damage = this->tempDamage;
    this->dissipation = this->tempDissipation;
    this->deltaDissipation = this->tempDeltaDissipation;
    this->crack_flag = this->temp_crack_flag;
    this->crackWidth = this->tempCrackWidth;
    //Need to preserve oldValues for coupling;
    this->oldNormalStress = this->normalStress;
    this->normalStress = this->tempNormalStress;
    //set updateflag
    this->updateFlag = 1;
}

Interface *
LatticeDamage2dStatus :: giveInterface(InterfaceType type)
{
    if ( type == RandomMaterialStatusExtensionInterfaceType ) {
        return static_cast< RandomMaterialStatusExtensionInterface * >(this);
    } else {
        return NULL;
    }
}

void
LatticeDamage2dStatus :: setTempNormalStress(double val)
{
    tempNormalStress = val;
}

void
LatticeDamage2dStatus :: setOldNormalStress(double val)
{
    oldNormalStress = val;
}


contextIOResultType
LatticeDamage2dStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
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
    if ( ( iores = reducedStrain.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

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

    if ( !stream.write(le) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(dissipation) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(deltaDissipation) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(crack_flag) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(crackWidth) ) {
        THROW_CIOERR(CIO_IOERR);
    }


    if ( !stream.write(e0) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(biot) ) {
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
LatticeDamage2dStatus :: setVariableInStatus(double variable)
{
    e0 = variable;
}


void
LatticeDamage2dStatus :: setBiotCoefficientInStatus(double variable)
{
    biot = variable;
}

int
LatticeDamage2dStatus :: giveCrackFlag()
{
    if ( crack_flag != 0 ) { }

    return crack_flag;
}


contextIOResultType
LatticeDamage2dStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
//
// restores full information stored in stream to this Status
//
{
    contextIOResultType iores;
    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = reducedStrain.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

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

    if ( !stream.read(le) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(dissipation) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(deltaDissipation) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(crack_flag) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(crackWidth) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(e0) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(biot) ) {
        THROW_CIOERR(CIO_IOERR);
    }


    return CIO_OK;
}
}     // end namespace oofem
