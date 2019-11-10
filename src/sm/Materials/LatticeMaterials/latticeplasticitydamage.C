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

#include "latticeplasticitydamage.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatmatrixf.h"
#include "floatarray.h"
#include "floatarrayf.h"
#include "engngm.h"
#include "mathfem.h"
#include <math.h>
#include "latticematstatus.h"
#include "intarray.h"
#include "Elements/LatticeElements/latticestructuralelement.h"
#include "datastream.h"
#include "contextioerr.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(LatticePlasticityDamage);

LatticePlasticityDamage :: LatticePlasticityDamage(int n, Domain *d) : LatticeLinearElastic(n, d)
    //
    // constructor
    //
{
    yieldTol = 0.;
    myPi = 3.141592653;
}


bool
LatticePlasticityDamage :: hasMaterialModeCapability(MaterialMode mode) const
{
    return mode == _3dLattice;
}

void
LatticePlasticityDamage :: initializeFrom(InputRecord &ir)
{

    LatticeLinearElastic :: initializeFrom(ir);

    yieldTol = 1.e-6;
    IR_GIVE_OPTIONAL_FIELD(ir, yieldTol, _IFT_LatticePlasticityDamage_tol);

    newtonIter = 100;
    IR_GIVE_OPTIONAL_FIELD(ir, newtonIter, _IFT_LatticePlasticityDamage_iter);

    numberOfSubIncrements = 10;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfSubIncrements, _IFT_LatticePlasticityDamage_sub);

    IR_GIVE_FIELD(ir, ft, _IFT_LatticePlasticityDamage_ft);

    IR_GIVE_FIELD(ir, fc, _IFT_LatticePlasticityDamage_fc);

    this->frictionAngleOne = 0.5;
    IR_GIVE_OPTIONAL_FIELD(ir, frictionAngleOne, _IFT_LatticePlasticityDamage_angle1);
    this->frictionAngleTwo = 0.5;
    IR_GIVE_OPTIONAL_FIELD(ir, frictionAngleTwo, _IFT_LatticePlasticityDamage_angle2);
    this->flowAngleOne = 0.25;
    IR_GIVE_OPTIONAL_FIELD(ir, flowAngleOne, _IFT_LatticePlasticityDamage_flow);
    this->flowAngleTwo = this->frictionAngleTwo;

    softeningType = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, softeningType, _IFT_LatticePlasticityDamage_stype);

    if ( softeningType > 1 ) {
        OOFEM_ERROR("Unknown softening type");
    }

    IR_GIVE_FIELD(ir, wf, _IFT_LatticePlasticityDamage_wf);


    if ( softeningType == 1 ) { //bilinear softening
        this->ftOne = 0.15 * this->ft;
        IR_GIVE_OPTIONAL_FIELD(ir, ftOne, _IFT_LatticePlasticityDamage_ft1);
        this->wfOne = 0.1 * this->wf;
        IR_GIVE_OPTIONAL_FIELD(ir, wfOne, _IFT_LatticePlasticityDamage_wf1);
    }

    this->aHard = 0.1;
    IR_GIVE_OPTIONAL_FIELD(ir, aHard, _IFT_LatticePlasticityDamage_ahard);

    this->damageFlag = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, damageFlag, _IFT_LatticePlasticityDamage_damage);

}

void
LatticePlasticityDamage :: computeDamageParam(double &omega, double kappaDOne, double kappaDTwo, GaussPoint *gp)
{
    omega = 0.0;

    int nite = 0;
    double R = 0.;
    double Lhs = 0.;
    double help = 0.;
    double helpOne = 0.;
    double helpTwo = 0;
    double omegaOne = 0;
    double omegaTwo = 0;

    double le = static_cast< LatticeStructuralElement * >( gp->giveElement() )->giveLength();

    double strength;

    if ( this->ft == 0 ) {
        strength = ( this->fc - this->frictionAngleOne * this->frictionAngleTwo * this->ft ) / ( 1 + this->frictionAngleOne * this->frictionAngleTwo );
    } else {
        strength = this->ft;
    }

    if ( softeningType == 0 ) { //Default: Exponential
        omega = 0;
        do {
            nite++;
            help = le * ( kappaDOne + omega * kappaDTwo ) / this->wf;
            R = ( 1 - omega ) * kappaDTwo * this->eNormalMean - strength * exp(-help);
            Lhs = le * kappaDTwo / this->wf * strength * exp(-help) - kappaDTwo * this->eNormalMean;
            omega -= R / Lhs;
            if ( nite > 40 ) {
                OOFEM_ERROR("computeDamageParam: algorithm not converging");
            }
        } while ( fabs(R) >= 1.e-6 || omega < 0.0 );
    } else if ( softeningType == 1 ) {      //bilinear: Calculate all damage parameters and check which makes sense
        omegaOne = ( this->eNormalMean * kappaDTwo * this->wfOne - this->ft * this->wfOne - ( this->ftOne - this->ft ) * kappaDOne * le ) /
                   ( this->eNormalMean * kappaDTwo * this->wfOne + ( this->ftOne - this->ft ) * le * kappaDTwo );
        helpOne = le * kappaDOne + le * omegaOne * kappaDTwo;


        omegaTwo = ( this->eNormalMean * kappaDTwo * ( this->wf - this->wfOne ) - this->ftOne * ( this->wf - this->wfOne ) +
                     this->ftOne * kappaDOne * le  - this->ftOne * this->wfOne ) /
                   ( this->eNormalMean * kappaDTwo * ( this->wf - this->wfOne )  - this->ftOne * le * kappaDTwo );
        helpTwo = le * kappaDOne + le * omegaTwo * kappaDTwo;


        if ( helpOne >= 0. && helpOne < wfOne ) {
            omega = omegaOne;
        } else if ( helpTwo > wfOne && helpTwo < wf ) {
            omega = omegaTwo;
        } else if ( helpTwo > wf ) {
            omega = 0.99999;
        }
    } else {
        OOFEM_ERROR("Unknown softening type\n");
    }

    if ( omega > 1.0 ) {
        omega = 0.99999;
    } else if ( omega < 0.0 ) {
        omega = 0.;
    }

    return;
}


MaterialStatus *
LatticePlasticityDamage :: CreateStatus(GaussPoint *gp) const
{
    LatticePlasticityDamageStatus *answer = new LatticePlasticityDamageStatus(1, LatticePlasticityDamage :: domain, gp);
    return answer;
}


double
LatticePlasticityDamage :: computeYieldValue(const FloatArray &stress,
                                             const double kappa,
                                             GaussPoint *gp)
{
    double yieldValue = 0;
    double hardening = computeHardening(kappa, gp);

    double shearNorm = sqrt(pow(stress.at(2), 2.) + pow(stress.at(3), 2.) );


    double transition = -( this->fc - this->frictionAngleOne * this->frictionAngleTwo * this->ft ) / ( 1. + this->frictionAngleOne * this->frictionAngleTwo ) * hardening;

    if ( stress.at(1) >= transition ) { //main ellipse
        yieldValue = pow(shearNorm, 2.) + pow(this->frictionAngleOne, 2.) * pow(stress.at(1), 2.) +
                     2. * ( this->fc - this->frictionAngleOne * this->frictionAngleTwo * this->ft ) * pow(this->frictionAngleOne, 2.) / ( 1. + this->frictionAngleOne * this->frictionAngleTwo ) * hardening * stress.at(1) -
                     ( 2. * this->fc * this->ft * pow(this->frictionAngleOne, 2.) + ( 1. - this->frictionAngleOne * this->frictionAngleTwo ) * pow(this->ft, 2.) * pow(this->frictionAngleOne, 2.) ) / ( 1. + this->frictionAngleOne * this->frictionAngleTwo ) * pow(hardening, 2.);
    } else {   //cap ellipse
        yieldValue =  pow(shearNorm, 2.) + pow(stress.at(1) / this->frictionAngleTwo, 2.) +
                     2. * ( this->fc - this->frictionAngleOne * this->frictionAngleTwo * this->ft ) / ( pow(this->frictionAngleTwo, 2.) * ( 1. + this->frictionAngleOne * this->frictionAngleTwo ) ) * hardening * stress.at(1) +
                     ( pow(this->fc, 2.) * ( 1. - pow(this->frictionAngleOne * this->frictionAngleTwo, 2.) ) -
                       2. * this->fc * this->ft * this->frictionAngleOne * this->frictionAngleTwo * ( 1. + this->frictionAngleOne * this->frictionAngleTwo ) ) / ( pow(this->frictionAngleTwo, 2.) * pow(1. + this->frictionAngleOne * this->frictionAngleTwo, 2.) ) * pow(hardening, 2.);
    }

    return yieldValue;
}

double
LatticePlasticityDamage :: computeHardening(const double kappa, GaussPoint *gp)
{
    double hardening = exp(kappa / this->aHard);


    return hardening;
}

double
LatticePlasticityDamage :: computeDHardeningDKappa(const double kappa, GaussPoint *gp)
{
    double dHardeningDKappa = 1. / this->aHard * exp(kappa / this->aHard);

    return dHardeningDKappa;
}

double
LatticePlasticityDamage :: computeDDHardeningDDKappa(const double kappa, GaussPoint *gp)
{
    double dDHardeningDDKappa = 1. / pow(this->aHard, 2.) * exp(kappa / this->aHard);

    return dDHardeningDDKappa;
}

void
LatticePlasticityDamage :: computeFVector(FloatArray &answer,
                                          const FloatArray &stress,
                                          const double kappa,
                                          GaussPoint *gp)
{
    double hardening = computeHardening(kappa, gp);
    double dHardeningDKappa = computeDHardeningDKappa(kappa, gp);


    double shearNorm = sqrt(pow(stress.at(2), 2.) + pow(stress.at(3), 2.) );

    double transition = -( this->fc - this->frictionAngleOne * this->frictionAngleTwo * this->ft ) / ( 1. + this->frictionAngleOne * this->frictionAngleTwo ) * hardening;

    if ( stress.at(1) >= transition ) { //main ellipse
        answer.at(1) = 2. * pow(this->frictionAngleOne, 2.) * stress.at(1) + 2. * ( this->fc - this->frictionAngleOne * this->frictionAngleTwo * this->ft ) * pow(this->frictionAngleOne, 2.) / ( 1. + this->frictionAngleOne * this->frictionAngleTwo ) * hardening;
        answer.at(2) = 2. * shearNorm;
        answer.at(3) = 2. * ( this->fc - this->frictionAngleOne * this->frictionAngleTwo * this->ft ) * pow(this->frictionAngleOne, 2.) / ( 1. + this->frictionAngleOne * this->frictionAngleTwo ) * stress.at(1) * dHardeningDKappa -
                       2. * ( 2. * this->fc * this->ft * pow(this->frictionAngleOne, 2.) + ( 1. - this->frictionAngleOne * this->frictionAngleTwo ) * pow(this->ft, 2.) * pow(this->frictionAngleOne, 2.) ) / ( 1. + this->frictionAngleOne * this->frictionAngleTwo ) * hardening * dHardeningDKappa;
    } else {   //cap ellipse
        answer.at(1) = 2. * stress.at(1) / pow(this->frictionAngleTwo, 2.) + 2. * ( this->fc - this->frictionAngleOne * this->frictionAngleTwo * this->ft ) / ( pow(this->frictionAngleTwo, 2.) * ( 1. + this->frictionAngleOne * this->frictionAngleTwo ) ) * hardening;
        answer.at(2) = 2. * shearNorm;
        answer.at(3) = 2. * ( this->fc - this->frictionAngleOne * this->frictionAngleTwo * this->ft ) / ( pow(this->frictionAngleTwo, 2.) * ( 1. + this->frictionAngleOne * this->frictionAngleTwo ) ) * stress.at(1) * dHardeningDKappa +
                       2. * ( pow(this->fc, 2.) * ( 1. - pow(this->frictionAngleOne * this->frictionAngleTwo, 2.) ) -
                              2. * this->fc * this->ft * this->frictionAngleOne * this->frictionAngleTwo * ( 1. + this->frictionAngleOne * this->frictionAngleTwo ) ) /
                       ( pow(frictionAngleTwo, 2.) * pow(1. + this->frictionAngleOne * this->frictionAngleTwo, 2.) ) * hardening * dHardeningDKappa;
    }

    return;
}

void
LatticePlasticityDamage :: computeMVector(FloatArray &answer,
                                          const FloatArray &stress,
                                          const double kappa,
                                          GaussPoint *gp)
{
    double hardening = computeHardening(kappa, gp);

    double shearNorm = sqrt(pow(stress.at(2), 2.) + pow(stress.at(3), 2.) );

    double transition = -( this->fc - this->flowAngleOne * this->frictionAngleTwo * this->ft ) / ( 1. + this->flowAngleOne * this->frictionAngleTwo ) * hardening;

    if ( stress.at(1) >= transition ) { //ellipse
        answer.at(1) = 2. * pow(this->flowAngleOne, 2.) * stress.at(1) + 2. * ( this->fc - this->flowAngleOne * this->flowAngleTwo * this->ft ) * pow(this->flowAngleOne, 2.) / ( 1. + this->flowAngleOne * this->flowAngleTwo ) * hardening;
        answer.at(2) = 2. * shearNorm;
        answer.at(3) = fabs(answer.at(1) );
    } else {   //circle
        answer.at(1) = 2. * stress.at(1) / pow(this->flowAngleTwo, 2.) + 2. * ( this->fc - this->flowAngleTwo * this->flowAngleOne * this->ft ) / ( pow(this->flowAngleTwo, 2.) * ( 1. + this->flowAngleOne * this->flowAngleTwo ) ) * hardening;
        answer.at(2) = 2. * shearNorm;
        answer.at(3) = fabs(answer.at(1) );
    }

    return;
}


void
LatticePlasticityDamage :: computeDMMatrix(FloatMatrix &answer,
                                           const FloatArray &stress,
                                           const double kappa,
                                           GaussPoint *gp)
{
    double hardening = computeHardening(kappa, gp);
    double dHardeningDKappa = computeDHardeningDKappa(kappa, gp);

    FloatArray mVector(3);
    computeMVector(mVector, stress, kappa, gp);

    double transition = -( this->fc - this->flowAngleOne * this->flowAngleTwo * this->ft ) / ( 1. + this->flowAngleOne * this->flowAngleTwo ) * hardening;

    if ( stress.at(1) >= transition ) { //main ellipse
        //Derivatives of dGDSig
        answer.at(1, 1) = 2. * pow(this->flowAngleOne, 2.);
        answer.at(1, 2) = 0.;
        answer.at(1, 3) = 2. * ( this->fc - this->flowAngleOne * this->flowAngleTwo * this->ft ) * pow(this->flowAngleOne, 2.) / ( 1. + this->flowAngleOne * this->flowAngleTwo ) * dHardeningDKappa;

        //Derivatives of dGDTau
        answer.at(2, 1) = 0.;
        answer.at(2, 2) = 2.;
        answer.at(2, 3) = 0;

        //Derivates of evolution law
        answer.at(3, 1) = answer.at(1, 1);
        answer.at(3, 2) = answer.at(1, 2);
        answer.at(3, 3) = answer.at(1, 3);
    } else {   //cap ellipse
               //Derivatives of dGDSig
        answer.at(1, 1) = 2. / pow(this->flowAngleTwo, 2.);
        answer.at(1, 2) = 0.;
        answer.at(1, 3) = 2. * ( this->fc - this->flowAngleOne * this->flowAngleTwo * this->ft ) / ( pow(this->flowAngleTwo, 2.) * ( 1. + this->flowAngleOne * this->flowAngleTwo ) ) * dHardeningDKappa;

        //Derivatives of dGDTau
        answer.at(2, 1) = 0.;
        answer.at(2, 2) = 2.;
        answer.at(2, 3) = 0;

        //Derivates of evolution law
        answer.at(3, 1) = -answer.at(1, 1);
        answer.at(3, 2) = -answer.at(1, 2);
        answer.at(3, 3) = -answer.at(1, 3);
    }
    return;
}


void
LatticePlasticityDamage :: computeAMatrix(FloatMatrix &answer,
                                          const FloatArray &stress,
                                          const double kappa,
                                          const double deltaLambda,
                                          GaussPoint *gp)
{
    FloatMatrix dMMatrix(3, 3);
    computeDMMatrix(dMMatrix, stress, kappa, gp);
    /* Compute matrix*/
    answer.at(1, 1) = 1 / this->eNormalMean + deltaLambda * dMMatrix.at(1, 1);
    answer.at(1, 2) = deltaLambda * dMMatrix.at(1, 2);
    answer.at(1, 3) = deltaLambda * dMMatrix.at(1, 3);
    /**/
    answer.at(2, 1) = deltaLambda * dMMatrix.at(2, 1);
    answer.at(2, 2) = 1 / ( this->alphaOne * this->eNormalMean ) + deltaLambda * dMMatrix.at(2, 2);
    answer.at(2, 3) = deltaLambda * dMMatrix.at(2, 3);
    /**/
    answer.at(3, 1) = deltaLambda * dMMatrix.at(3, 1);
    answer.at(3, 2) = deltaLambda * dMMatrix.at(3, 2);
    answer.at(3, 3) = deltaLambda * dMMatrix.at(3, 3) - 1.;
}


void
LatticePlasticityDamage :: giveReducedStrain(FloatArray &answer,
                                             GaussPoint *gp,
                                             TimeStep *tStep)
{
    LatticePlasticityDamageStatus *status = ( LatticePlasticityDamageStatus * ) ( this->giveStatus(gp) );
    answer = status->giveReducedStrain();
}



void
LatticePlasticityDamage :: performPlasticityReturn(FloatArray &stress,
                                                   GaussPoint *gp,
                                                   const FloatArray &strain,
                                                   TimeStep *tStep)
{
    LatticePlasticityDamageStatus *status =
        ( LatticePlasticityDamageStatus * ) ( giveStatus(gp) );

    double tempKappa = 0.;
    //Get tempKappa from the status
    tempKappa = status->giveTempKappaP();
    /* Get plastic strain vector from status*/
    FloatArray tempPlasticStrain;

    tempPlasticStrain = status->giveTempPlasticStrain();

    /* Compute trial stress*/
    stress.resize(3);
    stress.zero();
    stress.at(1) = ( strain.at(1) - tempPlasticStrain.at(1) ) * this->eNormalMean;
    stress.at(2) = ( strain.at(2) - tempPlasticStrain.at(2) ) * this->alphaOne * this->eNormalMean;
    stress.at(3) = ( strain.at(3) - tempPlasticStrain.at(3) ) * this->alphaOne * this->eNormalMean;


    //Introduce variables for subincrementation
    //Only _3dLattice is possible

    FloatArray convergedStrain(3);

    FloatArray oldReducedStrain(6);
    FloatArray oldStrain(3);

    this->giveReducedStrain(oldReducedStrain, gp, tStep);

    for ( int i = 1; i <= 3; i++ ) {
        oldStrain.at(i) = oldReducedStrain.at(i);
    }

    FloatArray tempStrain(3);
    FloatArray deltaStrain(3);



    // introduce a subincrementation flag
    int subIncrementFlag = 0;

    /* Compute yield value*/
    double yieldValue = computeYieldValue(stress, tempKappa, gp);
    /* Check yield condition, i.e. if the yield value is less than the yield tolerance.
     * If yield condition is valid. Do perform regular return (closest point return)*/
    if ( yieldValue / pow(fc, 2.) > yieldTol ) {
        subIncrementFlag = 0;
        convergedStrain = oldStrain;
        tempStrain = strain;
        deltaStrain = oldStrain;
        deltaStrain.times(-1.);
        deltaStrain.add(strain);
        //To get into the loop
        returnResult = RR_NotConverged;
        while ( returnResult == RR_NotConverged || subIncrementFlag == 1 ) {
            stress.at(1) = ( tempStrain.at(1) - tempPlasticStrain.at(1) ) * eNormalMean;
            stress.at(2) = ( tempStrain.at(2) - tempPlasticStrain.at(2) ) * ( this->alphaOne * eNormalMean );
            stress.at(3) = ( tempStrain.at(3) - tempPlasticStrain.at(3) ) * ( this->alphaOne * eNormalMean );

            tempKappa = performRegularReturn(stress, yieldValue, gp);

            if ( returnResult == RR_NotConverged ) {
                subIncrementFlag = 1;
                deltaStrain.times(0.5);
                tempStrain = convergedStrain;
                tempStrain.add(deltaStrain);
            } else if ( returnResult == RR_Converged && subIncrementFlag == 1 ) {
                //	printf("Subincrementation sucessful\n");

                tempPlasticStrain.at(1) = tempStrain.at(1) - stress.at(1) / eNormalMean;
                tempPlasticStrain.at(2) = tempStrain.at(2) - stress.at(2) / ( this->alphaOne * eNormalMean );
                tempPlasticStrain.at(3) = tempStrain.at(3) - stress.at(3) / ( this->alphaOne * eNormalMean );

                status->letTempPlasticStrainBe(tempPlasticStrain);
                status->setTempKappaP(tempKappa);

                subIncrementFlag = 0;
                returnResult = RR_NotConverged;
                convergedStrain = tempStrain;
                deltaStrain = convergedStrain;
                deltaStrain.times(-1.);
                deltaStrain.add(strain);
                tempStrain = strain;
            } else {
                status->setTempKappaP(tempKappa);
            }
        }
    }


    tempPlasticStrain.at(1) = strain.at(1) - stress.at(1) / eNormalMean;
    tempPlasticStrain.at(2) = strain.at(2) - stress.at(2) / ( this->alphaOne * eNormalMean );
    tempPlasticStrain.at(3) = strain.at(3) - stress.at(3) / ( this->alphaOne * eNormalMean );

    status->letTempPlasticStrainBe(tempPlasticStrain);
    status->letTempStressVectorBe(stress);

    printf("plastic strain =\n");
    tempPlasticStrain.printYourself();
}


double
LatticePlasticityDamage :: performRegularReturn(FloatArray &stress,
                                                double yieldValue,
                                                GaussPoint *gp)
{
    LatticePlasticityDamageStatus *status = ( LatticePlasticityDamageStatus * ) ( giveStatus(gp) );
    FloatArray trialStress(3), tempStress(3), residuals(4), residualsNorm(4), mVector(3), fVector(3);
    FloatMatrix jacobian(4, 4), inverseOfJacobian(4, 4);
    FloatArray rVector(3), helpVector(3), helpVector2(3);

    MaterialMode matMode = gp->giveMaterialMode();

    FloatMatrix aMatrix(3, 3), aMatrixInv(3, 3);
    FloatArray deltaIncrement(3);
    deltaIncrement.zero();
    FloatArray deltaIncrementTemp(3);
    deltaIncrementTemp.zero();
    FloatArray tempPlasticStrain(matMode);
    FloatArray plasticStrain(matMode);
    double deltaLambda = 0.;
    double normOfResiduals = 0.;
    int iterationCount = 0;

    FloatArray jacobianTimesAnswerIncrement;
    FloatArray unknownsTrial;
    FloatArray residualsTrial;
    FloatArray deltaUnknowns;
    FloatArray jacobianTimesDeltaUnknowns;
    FloatArray unknowns(4);

    trialStress = stress;
    tempStress = trialStress;

    double trialShearStressNorm;
    trialShearStressNorm = sqrt(pow(trialStress.at(2), 2.) + pow(trialStress.at(3), 2.) );

    double tempShearStressNorm = trialShearStressNorm;

    plasticStrain = status->givePlasticStrain();
    tempPlasticStrain = plasticStrain;

    double thetaTrial;
    thetaTrial = atan2(stress.at(3), stress.at(2) );

    // Do the same for kappa
    double kappa = status->giveKappaP();
    double tempKappa = kappa;


    //initialise unknowns
    unknowns.at(1) = trialStress.at(1);
    unknowns.at(2) = trialShearStressNorm;
    unknowns.at(3) = tempKappa;
    unknowns.at(4) = 0.;


    // Look at the magnitudes of the residuals. You have to scale the yieldValue down.
    yieldValue = computeYieldValue(tempStress, tempKappa, gp);

    //initiate residuals
    residuals.zero();
    residuals.at(4) = yieldValue;

    normOfResiduals  = 1.; //just to get into the loop

    while ( normOfResiduals > yieldTol ) {
        iterationCount++;
        if ( iterationCount == newtonIter ) {
            returnResult = RR_NotConverged;
            return 0.;
        }

        //Normalize residuals. Think about it more.
        residualsNorm.at(1) = residuals.at(1) / this->fc;
        residualsNorm.at(2) = residuals.at(2) / this->fc;
        residualsNorm.at(3) = residuals.at(3);
        residualsNorm.at(4) = residuals.at(4) / pow(this->fc, 2.);

        normOfResiduals = residualsNorm.computeNorm();

        //First check if return has failed
        if ( isnan(normOfResiduals) ) {
            returnResult = RR_NotConverged;
            return 0.;
        }


        if ( normOfResiduals > yieldTol ) {
            // Test to run newton iteration using inverse of Jacobian
            computeJacobian(jacobian, tempStress, tempKappa, deltaLambda, gp);

            if ( computeInverseOfJacobian(inverseOfJacobian, jacobian) ) {
                returnResult = RR_NotConverged;
                return 0.;
            }


            deltaIncrement.beProductOf(inverseOfJacobian, residuals);
            deltaIncrement.times(-1.);


            //compute trial values
            unknownsTrial = unknowns;
            residualsTrial = residuals;

            //compute Unknowns
            for ( int i = 0; i < 4; i++ ) {
                unknowns(i) = unknownsTrial(i) + deltaIncrement(i);
            }
            if ( unknowns.at(4) <= 0. ) { //Keep deltaLambda greater than zero!
                unknowns.at(4) = 0.;
            }
            if ( unknowns.at(2) <= 0. ) { //Keep rho greater than zero!
                unknowns.at(2) = 0.;
            }
            if ( unknowns.at(3) - kappa <= 0. ) { //Keep deltaKappa greater than zero!
                unknowns.at(3) = kappa;
            }


            /* Update increments final values and DeltaLambda*/
            tempStress.at(1) = unknowns.at(1);
            tempShearStressNorm  = unknowns.at(2);

            tempStress.at(2) = tempShearStressNorm * cos(thetaTrial);
            tempStress.at(3) = tempShearStressNorm * sin(thetaTrial);

            tempKappa = unknowns.at(3);
            deltaLambda = unknowns.at(4);

            /* Compute the mVector holding the derivatives of the g function and the hardening function*/
            computeMVector(mVector, tempStress, tempKappa, gp);

            residuals.at(1) = tempStress.at(1) - trialStress.at(1) + this->eNormalMean * deltaLambda * mVector.at(1);
            residuals.at(2) = tempShearStressNorm - trialShearStressNorm + this->alphaOne * this->eNormalMean * deltaLambda * mVector.at(2);
            residuals.at(3) = -tempKappa + kappa + deltaLambda * mVector.at(3);
            residuals.at(4) = computeYieldValue(tempStress, tempKappa, gp);
        }
    }

    //  status-> setTempKappaP (tempKappa);

    returnResult = RR_Converged;

    stress = tempStress;

    return tempKappa;
}



void
LatticePlasticityDamage :: computeJacobian(FloatMatrix &answer,
                                           const FloatArray &stress,
                                           const double kappa,
                                           const double deltaLambda,
                                           GaussPoint *gp)
{
    //Variables
    FloatArray fVector(3);
    FloatArray mVector(3);
    FloatMatrix dMMatrix(3, 3);

    computeDMMatrix(dMMatrix, stress, kappa, gp);

    computeMVector(mVector, stress, kappa, gp);
    computeFVector(fVector, stress, kappa, gp);

    /* Compute matrix*/
    answer.at(1, 1) = 1. + this->eNormalMean * deltaLambda * dMMatrix.at(1, 1);
    answer.at(1, 2) = this->eNormalMean * deltaLambda * dMMatrix.at(1, 2);
    answer.at(1, 3) = this->eNormalMean * deltaLambda * dMMatrix.at(1, 3);
    answer.at(1, 4) = this->eNormalMean * mVector.at(1);
    /**/
    answer.at(2, 1) = this->alphaOne * this->eNormalMean * deltaLambda * dMMatrix.at(2, 1);
    answer.at(2, 2) = 1. + this->alphaOne * this->eNormalMean * deltaLambda * dMMatrix.at(2, 2);
    answer.at(2, 3) = this->alphaOne * this->eNormalMean * deltaLambda * dMMatrix.at(2, 3);
    answer.at(2, 4) = this->alphaOne * this->eNormalMean * mVector.at(2);
    /**/
    answer.at(3, 1) = deltaLambda * dMMatrix.at(3, 1);
    answer.at(3, 2) = deltaLambda * dMMatrix.at(3, 2);
    answer.at(3, 3) = deltaLambda * dMMatrix.at(3, 3) - 1.;
    answer.at(3, 4) = mVector.at(3);
    /**/
    answer.at(4, 1) = fVector.at(1);
    answer.at(4, 2) = fVector.at(2);
    answer.at(4, 3) = fVector.at(3);
    answer.at(4, 4) = 0.;

    return;
}



int
LatticePlasticityDamage :: computeInverseOfJacobian(FloatMatrix &answer, const FloatMatrix &src)
// Receiver becomes inverse of given parameter src. If necessary, size is adjusted.
{
#  ifdef DEBUG
    if ( !src.isSquare() ) {
        OOFEM_ERROR("FloatMatrix::beInverseOf : cannot inverse matrix since it is not square");
    }
#  endif

    int nRows = src.giveNumberOfRows();

    //gaussian elimination - slow but safe
    //
    int i, j, k;
    double piv, linkomb;
    FloatMatrix tmp = src;
    answer.zero();
    // initialize answer to be unity matrix;
    for ( i = 1; i <= nRows; i++ ) {
        answer.at(i, i) = 1.0;
    }

    // lower triangle elimination by columns
    for ( i = 1; i < nRows; i++ ) {
        piv = tmp.at(i, i);
        if ( fabs(piv) < 1.e-20 ) {
            //indication that return does not converge. Output flag for subincrementing
            return 1;
        }

        for ( j = i + 1; j <= nRows; j++ ) {
            linkomb = tmp.at(j, i) / tmp.at(i, i);
            for ( k = i; k <= nRows; k++ ) {
                tmp.at(j, k) -= tmp.at(i, k) * linkomb;
            }

            for ( k = 1; k <= nRows; k++ ) {
                answer.at(j, k) -= answer.at(i, k) * linkomb;
            }
        }
    }

    // upper triangle elimination by columns
    for ( i = nRows; i > 1; i-- ) {
        piv = tmp.at(i, i);
        for ( j = i - 1; j > 0; j-- ) {
            linkomb = tmp.at(j, i) / tmp.at(i, i);
            for ( k = i; k > 0; k-- ) {
                tmp.at(j, k) -= tmp.at(i, k) * linkomb;
            }

            for ( k = nRows; k > 0; k-- ) {
                answer.at(j, k) -= answer.at(i, k) * linkomb;
            }
        }
    }

    // diagonal scaling
    for ( i = 1; i <= nRows; i++ ) {
        for ( j = 1; j <= nRows; j++ ) {
            answer.at(i, j) /= tmp.at(i, i);
        }
    }

    return 0;
}

FloatArrayF< 6 >
LatticePlasticityDamage :: giveLatticeStress3d(const FloatArrayF< 6 > &originalStrain, GaussPoint *gp, TimeStep *tStep)
{
    FloatArray answer;
    answer.resize(6);
    answer.zero();


    LatticePlasticityDamageStatus *status = ( LatticePlasticityDamageStatus * ) this->giveStatus(gp);
    FloatArray reducedStrain;

    this->giveStressDependentPartOfStrainVector(reducedStrain, gp, originalStrain, tStep, VM_Total);

    gp->giveMaterialMode();

    //Subset of reduced strain.
    //Rotational components are not used
    FloatArray strain;
    strain.resize(3);
    strain.zero();

    strain.at(1) = reducedStrain.at(1);
    strain.at(2) = reducedStrain.at(2);
    strain.at(3) = reducedStrain.at(3);

    double omega = 0.;

    this->performPlasticityReturn(answer, gp, strain, tStep);

    if ( damageFlag == 1 ) {
        this->performDamageEvaluation(gp, strain);
        omega = status->giveTempDamage();
    }

    answer.resizeWithValues(6);

    answer.at(1) = ( 1. - omega ) * answer.at(1);
    answer.at(2) = ( 1. - omega ) * answer.at(2);
    answer.at(3) = ( 1. - omega ) * answer.at(3);
    answer.at(4) = ( 1. - omega ) * reducedStrain.at(4) * this->alphaTwo * this->eNormalMean;
    answer.at(5) = ( 1. - omega ) * reducedStrain.at(5) * this->alphaTwo * this->eNormalMean;
    answer.at(6) = ( 1. - omega ) * reducedStrain.at(6) * this->alphaTwo * this->eNormalMean;

    status->letTempStrainVectorBe(originalStrain);
    status->letTempReducedStrainBe(reducedStrain);
    status->letTempStressVectorBe(answer);

    return answer;
}




void
LatticePlasticityDamage :: performDamageEvaluation(GaussPoint *gp, FloatArray &reducedStrain)
{
    LatticePlasticityDamageStatus *status = ( LatticePlasticityDamageStatus * ) this->giveStatus(gp);
    // FloatArray tempStress, tempEffectiveStress;
    double tempKappaDOne, tempKappaDTwo, omega = 0.;

    double le = static_cast< LatticeStructuralElement * >( gp->giveElement() )->giveLength();

    FloatArray tempPlasticStrain = status->giveTempPlasticStrain();

    printf("plastic strain in damage =\n");
    tempPlasticStrain.printYourself();


    FloatArray plasticStrain = status->givePlasticStrain();
    FloatArray deltaPlasticStrain;
    deltaPlasticStrain = tempPlasticStrain;
    deltaPlasticStrain.subtract(plasticStrain);

    //Compute the damage history variable
    tempKappaDOne = status->giveKappaDOne();
    tempKappaDTwo = status->giveKappaDTwo();

    // Compute deltaKappaP from the plastic Strain
    double tempKappaP = status->giveTempKappaP();
    double kappaP = status->giveKappaP();
    double deltaKappaP = tempKappaP - kappaP;
    double hardening = computeHardening(tempKappaP, gp);
    //    FloatArray elasticStrain;
    //    FloatArray deltaElasticStrain;

    double crackWidth = 0.;

    double strength = ( this->fc - this->frictionAngleOne * this->frictionAngleTwo * this->ft ) / ( 1 + this->frictionAngleOne * this->frictionAngleTwo );

    if ( deltaKappaP <= 0. ) { //unloading or reloading
        tempKappaDOne = status->giveKappaDOne();
        tempKappaDTwo = status->giveKappaDTwo();
        omega = status->giveDamage();
        if ( status->giveCrackFlag() != 0 ) {
            status->setTempCrackFlag(2);
        } else {
            status->setTempCrackFlag(0);
        }
    } else {   //plastic loading
        if ( deltaPlasticStrain.at(1) <= 0. ) { //compressive side! No damage
            tempKappaDOne = status->giveKappaDOne();
            tempKappaDTwo = status->giveKappaDTwo();
            omega = status->giveDamage();
            status->setTempCrackFlag(3);
        } else {
            tempKappaDOne += deltaPlasticStrain.at(1);

            if ( this->ft == 0 ) {
                tempKappaDTwo = hardening * strength / this->eNormalMean;
            } else {
                tempKappaDTwo = hardening * this->ft / this->eNormalMean;
            }

            // evaluate damage parameter
            this->computeDamageParam(omega, tempKappaDOne, tempKappaDTwo, gp);

            //threshold for crack patterns
            if ( ( tempKappaDOne + omega * tempKappaDTwo ) * le > 0. ) {
                status->setTempCrackFlag(1);
            } else {
                status->setTempCrackFlag(0);
            }
        }
    }

    crackWidth = sqrt(pow(tempPlasticStrain.at(1) + omega * ( reducedStrain.at(1) - tempPlasticStrain.at(1) ), 2.) +
                      pow(tempPlasticStrain.at(2) + omega * ( reducedStrain.at(2) - tempPlasticStrain.at(2) ), 2.) +
                      pow(tempPlasticStrain.at(3) + omega * ( reducedStrain.at(3) - tempPlasticStrain.at(3) ), 2.) ) * le;

    //TODO: Compute dissipation
    // double tempDissipation = status->giveDissipation();
    // double tempDeltaDissipation = 0.;
    // tempDeltaDissipation = computeDeltaDissipation(omega, reducedStrain, gp, atTime);
    // tempDissipation += tempDeltaDissipation;

    status->setTempKappaDOne(tempKappaDOne);
    status->setTempKappaDTwo(tempKappaDTwo);
    status->setTempDamage(omega);
    status->setTempCrackWidth(crackWidth);

    printf("damage = %e\n", omega);

    return;
}

/*double LatticePlasticityDamage :: giveMaterialParameter()
 * {
 *  return this->ft;
 *  }*/


// double
// LatticePlasticityDamage :: computeDeltaDissipation(double omega,
//                                            FloatArray &reducedStrain,
//                                            GaussPoint *gp,
//                                            TimeStep *atTime)
// {
//     LatticeDamageStatus *status = static_cast< LatticeDamageStatus * >( this->giveStatus(gp) );
//     double length = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveLength();
//     const double e0 = this->give(e0_ID, gp) * this->e0Mean;
//     double eNormal = this->give(eNormal_ID, gp) * this->eNormalMean;

//     const double eShear =  this->alphaOne * eNormal;
//     const double eTorsion =  this->alphaTwo * eNormal / 12.;

//     FloatArray reducedStrainOld;

//     reducedStrainOld = status->giveReducedStrain();
//     double omegaOld = status->giveDamage();
//     double deltaOmega;

//     FloatArray crackOpeningOld(6);
//     crackOpeningOld.times(omegaOld);
//     crackOpeningOld.times(length);
//     FloatArray stressOld( status->giveStressVector() );
//     FloatArray intermediateStrain(6);

//     double tempDeltaDissipation = 0.;
//     double deltaTempDeltaDissipation = 0.;

//     double intermediateOmega = 0;
//     FloatArray oldIntermediateStrain(6);
//     oldIntermediateStrain = reducedStrainOld;
//     double oldIntermediateOmega = omegaOld;
//     deltaOmega = ( omega - omegaOld );
//     double testDissipation  = 0.5 * length * ( pow( ( reducedStrain.at(1) + reducedStrainOld.at(1) ) / 2., 2. ) * eNormal +
//                                                pow( ( reducedStrain.at(2) + reducedStrainOld.at(2) ) / 2., 2. ) * eShear +
//                                                pow( ( reducedStrain.at(3) + reducedStrainOld.at(3) ) / 2., 2. ) * eShear +
//                                                pow( ( reducedStrain.at(4) + reducedStrainOld.at(4) ) / 2., 2. ) * eTorsion +
//                                                pow( ( reducedStrain.at(5) + reducedStrainOld.at(5) ) / 2., 2. ) * eTorsion +
//                                                pow( ( reducedStrain.at(6) + reducedStrainOld.at(6) ) / 2., 2. ) * eTorsion ) * deltaOmega;
//     double intervals = 0.;

//     double referenceGf = 0;

//     if ( softeningType == 1 ) {
//         referenceGf = e0 * eNormal * this->wf / 2.;
//     } else {   //This is for the exponential law. Should also implement it for the bilinear one.
//         referenceGf = e0 * eNormal * this->wf;
//     }

//     if ( testDissipation / ( referenceGf ) > 0.001 ) {
//         intervals = 1000. * testDissipation / referenceGf;
//     } else {
//         intervals = 1.;
//     }

//     if ( intervals > 1000. ) {
//         intervals = 1000.;
//     }

//     double oldKappa = status->giveKappa();
//     double f, equivStrain;
//     if ( deltaOmega > 0 ) {
//         for ( int k = 0; k < intervals; k++ ) {
//             intermediateStrain(0) = reducedStrainOld(0) + ( k + 1 ) / intervals * ( reducedStrain(0) - reducedStrainOld(0) );
//             intermediateStrain(1) = reducedStrainOld(1) + ( k + 1 ) / intervals * ( reducedStrain(1) - reducedStrainOld(1) );
//             intermediateStrain(2) = reducedStrainOld(2) + ( k + 1 ) / intervals * ( reducedStrain(2) - reducedStrainOld(2) );
//             intermediateStrain(3) = reducedStrainOld(3) + ( k + 1 ) / intervals * ( reducedStrain(3) - reducedStrainOld(3) );
//             intermediateStrain(4) = reducedStrainOld(4) + ( k + 1 ) / intervals * ( reducedStrain(4) - reducedStrainOld(4) );
//             intermediateStrain(5) = reducedStrainOld(5) + ( k + 1 ) / intervals * ( reducedStrain(5) - reducedStrainOld(5) );

//             this->computeEquivalentStrain(equivStrain, intermediateStrain, gp, atTime);
//             f = equivStrain - oldKappa;
//             if ( f > 0 ) {
//                 this->computeDamageParam(intermediateOmega, equivStrain, intermediateStrain, gp);
//                 deltaOmega = ( intermediateOmega - oldIntermediateOmega );
//                 deltaTempDeltaDissipation =
//                     0.5 * length * ( pow( ( intermediateStrain(0) + oldIntermediateStrain(0) ) / 2., 2. ) * eNormal +
//                                      pow( ( intermediateStrain(1) + oldIntermediateStrain(1) ) / 2., 2. ) * eShear +
//                                      pow( ( intermediateStrain(2) + oldIntermediateStrain(2) ) / 2., 2. ) * eShear +
//                                      pow( ( intermediateStrain(3) + oldIntermediateStrain(3) ) / 2., 2. ) * eTorsion +
//                                      pow( ( intermediateStrain(4) + oldIntermediateStrain(4) ) / 2., 2. ) * eTorsion +
//                                      pow( ( intermediateStrain(5) + oldIntermediateStrain(5) ) / 2., 2. ) * eTorsion ) * deltaOmega;

//                 oldKappa = equivStrain;
//                 oldIntermediateOmega = intermediateOmega;
//             } else {
//                 deltaTempDeltaDissipation = 0.;
//             }

//             tempDeltaDissipation += deltaTempDeltaDissipation;
//             oldIntermediateStrain = intermediateStrain;
//         }
//     } else {
//         tempDeltaDissipation = 0.;
//     }

//     if ( tempDeltaDissipation >= 2. * referenceGf ) {
//         tempDeltaDissipation = 2. * referenceGf;
//     }

//     return tempDeltaDissipation;
// }


FloatMatrixF< 6, 6 >
LatticePlasticityDamage :: give3dLatticeStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    FloatMatrix answer;
    answer.resize(6, 6);
    answer.zero();

    answer = LatticeLinearElastic :: give3dLatticeStiffnessMatrix(mode, gp, tStep);

    if ( mode == ElasticStiffness ) {
        return answer;
    } else if ( ( mode == SecantStiffness ) || ( mode == TangentStiffness ) ) {
        LatticePlasticityDamageStatus *status = static_cast< LatticePlasticityDamageStatus * >( this->giveStatus(gp) );
        double omega = status->giveTempDamage();

        if ( omega > 0.99999 ) {
            omega = 0.99999;
        }

        answer.times(1. - omega);

        return answer;
    } else {
        OOFEM_ERROR("Unsupported stiffness mode\n");
    }

    return answer;
}

int
LatticePlasticityDamage :: giveIPValue(FloatArray &answer,
                                       GaussPoint *gp,
                                       InternalStateType type,
                                       TimeStep *atTime)
{
    LatticePlasticityDamageStatus *status = ( LatticePlasticityDamageStatus * ) this->giveStatus(gp);

    if ( type == IST_DamageTensor ) {
        answer.resize(6);
        answer.zero();
        answer.at(1) = answer.at(2) = answer.at(3) = status->giveDamage();
        return 1;
    } else if ( type == IST_DamageScalar ) {
        answer.resize(1);
        answer.zero();
        answer.at(1) = status->giveDamage();
        return 1;
    } else if ( type == IST_CrackStatuses ) {
        answer.resize(1);
        answer.zero();
        answer.at(1) = status->giveCrackFlag();
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
    } else if ( type == IST_CharacteristicLength ) {
        answer.resize(1);
        answer.zero();
        answer.at(1) = static_cast< LatticeStructuralElement * >( gp->giveElement() )->giveLength();
        return 1;
    } else {
        return StructuralMaterial :: giveIPValue(answer, gp, type, atTime);
    }
}


LatticePlasticityDamageStatus :: LatticePlasticityDamageStatus(int n, Domain *d, GaussPoint *g) :  LatticeMaterialStatus(g)
{
    damage = tempDamage = 0.;
    kappaDOne = tempKappaDOne = 0.;
    kappaDTwo = tempKappaDTwo = 0.;
    kappaP = tempKappaP = 0.;
    compressionFlag = 0;
}

void
LatticePlasticityDamageStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
// builds new crackMap
//
{
    LatticeMaterialStatus :: initTempStatus();
    this->tempPlasticStrain = this->plasticStrain;
    this->tempKappaP = this->kappaP;
    this->tempKappaDOne = this->kappaDOne;
    this->tempKappaDTwo = this->kappaDTwo;
    this->tempDamage = this->damage;
}

void
LatticePlasticityDamageStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    LatticeMaterialStatus :: printOutputAt(file, tStep);

    fprintf(file, "plasticStrains ");
    for ( int k = 1; k <= plasticStrain.giveSize(); k++ ) {
        fprintf(file, "% .8e ", plasticStrain.at(k) );
    }

    fprintf(file, "\n \t");
    fprintf(file, " kappaP %.8e, kappaDOne %.8e, kappaDTwo %.8e, damage %.8e, deltaDissipation %.8e, dissipation %.8e, crackFlag %d, crackWidth %.8e }\n", this->kappaP, this->kappaDOne,  this->kappaDTwo, this->damage, this->deltaDissipation, this->dissipation, this->crackFlag, this->crackWidth);
    return;
}


void
LatticePlasticityDamageStatus :: updateYourself(TimeStep *atTime)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reached equilibrium.
//
{
    LatticeMaterialStatus :: updateYourself(atTime);
    this->plasticStrain = this->tempPlasticStrain;
    this->kappaP = this->tempKappaP;
    this->kappaDOne = this->tempKappaDOne;
    this->kappaDTwo = this->tempKappaDTwo;
    this->damage = this->tempDamage;
}


void
LatticePlasticityDamageStatus :: saveContext(DataStream &stream, ContextMode mode)
//
// saves full information stored in this Status
// no temp variables stored
//
{
    contextIOResultType iores;
    // save parent class status
    LatticeMaterialStatus :: saveContext(stream, mode);


    if ( ( iores = plasticStrain.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream.write(& kappaP, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }
    if ( !stream.write(& kappaDOne, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }
    if ( !stream.write(& kappaDTwo, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }
    if ( !stream.write(& damage, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }
}


void
LatticePlasticityDamageStatus :: restoreContext(DataStream &stream, ContextMode mode)
//
// restores full information stored in stream to this Status
//
{
    contextIOResultType iores;
    // read parent class status
    LatticeMaterialStatus :: restoreContext(stream, mode);

    //FloatArrays
    if ( ( iores = plasticStrain.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
    if ( !stream.read(& kappaP, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }
    if ( !stream.read(& kappaDOne, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }
    if ( !stream.read(& kappaDTwo, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }
    if ( !stream.read(& damage, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }
}
}     // end namespace oofem
