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

#include "latticebondplasticity.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
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
REGISTER_Material(LatticeBondPlasticity);


LatticeBondPlasticity :: LatticeBondPlasticity(int n, Domain *d) : LatticeLinearElastic(n, d)
{
}


bool
LatticeBondPlasticity :: hasMaterialModeCapability(MaterialMode mode) const
{
    return mode == _3dLattice;
}

double
LatticeBondPlasticity :: computeHardening(double kappa)
{
    return exp(pow(kappa / this->ef, 2.) );
}

double
LatticeBondPlasticity :: computeDHardeningDKappa(double kappa)
{
    return 2. * kappa / ( pow(this->ef, 2.) ) * exp(pow(kappa / this->ef, 2.) );
}


void
LatticeBondPlasticity :: initializeFrom(InputRecord &ir)
{

    LatticeLinearElastic :: initializeFrom(ir);

    yieldTol = 1.e-6;
    IR_GIVE_OPTIONAL_FIELD(ir, yieldTol, _IFT_LatticeBondPlasticity_tol);

    newtonIter = 100;
    IR_GIVE_OPTIONAL_FIELD(ir, newtonIter, _IFT_LatticeBondPlasticity_iter);

    numberOfSubIncrements = 64;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfSubIncrements, _IFT_LatticeBondPlasticity_sub);

    IR_GIVE_FIELD(ir, this->fc, _IFT_LatticeBondPlasticity_fc);

    this->frictionAngleOne = 0.2;
    IR_GIVE_OPTIONAL_FIELD(ir, this->frictionAngleOne, _IFT_LatticeBondPlasticity_angle1);

    this->frictionAngleTwo = this->frictionAngleOne;
    //    IR_GIVE_OPTIONAL_FIELD(ir, this->frictionAngleTwo, _IFT_LatticeBondPlasticity_angle2);

    this->flowAngle = this->frictionAngleOne;
    //    IR_GIVE_OPTIONAL_FIELD(ir, this->flowAngle, _IFT_LatticeBondPlasticity_flow);


    this->ef = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, this->ef, _IFT_LatticeBondPlasticity_ef);

}

MaterialStatus *
LatticeBondPlasticity :: CreateStatus(GaussPoint *gp) const
{
    LatticeBondPlasticityStatus *answer = new LatticeBondPlasticityStatus(1, LatticeBondPlasticity :: domain, gp);
    return answer;
}

double
LatticeBondPlasticity :: computeShift(const double kappa)
{
    double hardening = computeHardening(kappa);

    double paramA = computeParamA(kappa);

    double shift = this->fc * hardening - paramA;

    return shift;
}


double
LatticeBondPlasticity :: computeParamA(const double kappa)
{
    double hardening = computeHardening(kappa);

    double paramA = this->frictionAngleTwo * this->frictionAngleOne * this->fc * hardening /
                    ( this->frictionAngleTwo * this->frictionAngleOne +
                      sqrt(1. + pow(this->frictionAngleTwo, 2.) * pow(this->frictionAngleOne, 2.) ) );

    return paramA;
}



double
LatticeBondPlasticity :: computeDShiftDKappa(const double kappa)
{
    double dHardeningDKappa = computeDHardeningDKappa(kappa);

    double dParamADKappa = computeDParamADKappa(kappa);

    double dShiftDKappa = this->fc * dHardeningDKappa - dParamADKappa;

    return dShiftDKappa;
}


double
LatticeBondPlasticity :: computeDParamADKappa(const double kappa)
{
    double dHardeningDKappa = computeDHardeningDKappa(kappa);

    double A2 = this->frictionAngleTwo * this->frictionAngleOne + sqrt(1. + pow(this->frictionAngleTwo, 2.) * pow(this->frictionAngleOne, 2.) );

    double dA1DKappa = this->frictionAngleTwo * this->frictionAngleOne * this->fc * dHardeningDKappa;

    double dParamADKappa = dA1DKappa / A2;

    return dParamADKappa;
}



void
LatticeBondPlasticity :: giveReducedStrain(FloatArray &answer,
                                           GaussPoint *gp,
                                           TimeStep *tStep)
{
    LatticeBondPlasticityStatus *status = ( LatticeBondPlasticityStatus * ) ( this->giveStatus(gp) );
    answer = status->giveReducedStrain();
}

void
LatticeBondPlasticity :: performPlasticityReturn(FloatArray &stress,
                                                 GaussPoint *gp,
                                                 const FloatArray &strain,
                                                 TimeStep *tStep)
{
    LatticeBondPlasticityStatus *status = static_cast< LatticeBondPlasticityStatus * >( this->giveStatus(gp) );

    double tempKappa = 0.;
    //Get tempKappa from the status
    tempKappa = status->giveTempKappaP();

    /* Get plastic strain vector from status*/
    FloatArray tempPlasticStrain;

    int transitionFlag;

    tempPlasticStrain = status->giveTempPlasticStrain();

    /* Compute trial stress*/
    stress.resize(3);
    stress.zero();
    stress.at(1) = ( strain.at(1) - tempPlasticStrain.at(1) ) * this->eNormalMean;
    stress.at(2) = ( strain.at(2) - tempPlasticStrain.at(2) ) * this->alphaOne * this->eNormalMean;
    stress.at(3) = ( strain.at(3) - tempPlasticStrain.at(3) ) * this->alphaOne * this->eNormalMean;

    //Introduce variables for subincrementation
    FloatArray convergedStrain(3);

    FloatArray oldReducedStrain(6);
    FloatArray oldStrain(3);

    this->giveReducedStrain(oldReducedStrain, gp, tStep);

    for ( int i = 1; i <= 3; i++ ) {
        oldStrain.at(i) = oldReducedStrain.at(i);
    }

    FloatArray tempStrain(3);
    FloatArray deltaStrain(3);
    FloatArray deltaStrainIncrement(3);

    /* Compute yield value*/
    transitionFlag = 0;
    double yieldValue = computeYieldValue(stress, tempKappa, transitionFlag, gp);

    double transition = computeTransition(tempKappa, gp);

    //If shear surface is violated then return to this surface.
    //No need to check situation with ellipse
    //If only ellipse is violated then try this.
    if ( yieldValue < yieldTol && stress.at(1) < transition - yieldTol ) {
        transitionFlag = 1;
        yieldValue = computeYieldValue(stress, tempKappa, 1, gp);
    }

    // introduce subincrementation flag
    int subIncrementFlag = 0;
    double subIncrementCounter = 1;

    /* Compute yield value*/
    double error = 0.;
    if ( transitionFlag == 0 ) {
        error = yieldValue / fc;
    } else {
        error = yieldValue / pow(fc, 2.);
    }

    if ( error > yieldTol ) {
        if ( checkForVertexCase(stress, gp) ) {
            performVertexReturn(stress, gp);
            surfaceType = ST_Vertex;
        } else {
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

                tempKappa = performRegularReturn(stress, yieldValue, transitionFlag, gp);

                if ( returnResult == RR_NotConverged ) {
                    subIncrementFlag = 1;


                    subIncrementCounter *= 2;
                    if ( subIncrementCounter >= numberOfSubIncrements ) {
                        OOFEM_ERROR("Subincrementation = %e did not help. Stop here with numberOfSubIncrements = %d.\n", subIncrementCounter, numberOfSubIncrements);
                    }

                    deltaStrain.times(0.5);
                    tempStrain = convergedStrain;
                    tempStrain.add(deltaStrain);
                } else if ( returnResult == RR_Converged && subIncrementFlag == 1 ) {
                    tempPlasticStrain.at(1) = tempStrain.at(1) - stress.at(1) / eNormalMean;
                    tempPlasticStrain.at(2) = tempStrain.at(2) - stress.at(2) / ( this->alphaOne * eNormalMean );
                    tempPlasticStrain.at(3) = tempStrain.at(3) - stress.at(3) / ( this->alphaOne * eNormalMean );

                    status->letTempPlasticStrainBe(tempPlasticStrain);

                    status->setTempKappaP(tempKappa);

                    subIncrementCounter -= 1;

                    convergedStrain = tempStrain;
                    deltaStrain = convergedStrain;
                    deltaStrain.times(-1.);
                    deltaStrain.add(strain);
                    deltaStrainIncrement = deltaStrain;
                    deltaStrainIncrement.times(1. / subIncrementCounter);
                    tempStrain = convergedStrain;
                    tempStrain.add(deltaStrainIncrement);
                    if ( subIncrementCounter > 1 ) {
                        subIncrementFlag = 1;
                    } else   {
                        subIncrementFlag = 0;
                    }
                    returnResult = RR_NotConverged;
                } else if ( returnResult == RR_Elastic && subIncrementFlag == 1 ) {
                    convergedStrain = tempStrain;
                    deltaStrain = convergedStrain;
                    deltaStrain.times(-1.);
                    deltaStrain.add(strain);
                    deltaStrainIncrement = deltaStrain;
                    deltaStrainIncrement.times(1. / subIncrementCounter);
                    tempStrain = convergedStrain;
                    tempStrain.add(deltaStrainIncrement);
                    if ( subIncrementCounter > 1 ) {
                        subIncrementFlag = 1;
                    } else   {
                        subIncrementFlag = 0;
                    }
                    returnResult = RR_NotConverged;
                } else if ( returnResult == RR_Converged && subIncrementFlag == 0 )      {
                    status->setTempKappaP(tempKappa);
                }
            }
        }
    }

    tempPlasticStrain.at(1) = strain.at(1) - stress.at(1) / eNormalMean;
    tempPlasticStrain.at(2) = strain.at(2) - stress.at(2) / ( this->alphaOne * eNormalMean );
    tempPlasticStrain.at(3) = strain.at(3) - stress.at(3) / ( this->alphaOne * eNormalMean );


    status->letTempPlasticStrainBe(tempPlasticStrain);
    status->letTempStressVectorBe(stress);
    status->setSurfaceValue(surfaceType);
}


double
LatticeBondPlasticity :: performRegularReturn(FloatArray &stress,
                                              double yieldValue,
                                              int transitionFlag,
                                              GaussPoint *gp)
{
    LatticeBondPlasticityStatus *status = ( LatticeBondPlasticityStatus * ) ( giveStatus(gp) );
    FloatArray trialStress(3), tempStress(3), residuals(4), residualsNorm(4), mVector(3), fVector(3);
    FloatMatrix jacobian(4, 4), inverseOfJacobian(4, 4);
    FloatArray rVector(3), helpVector(3), helpVector2(3);

    FloatMatrix aMatrix(3, 3), aMatrixInv(3, 3);
    FloatArray deltaIncrement(3);
    deltaIncrement.zero();
    FloatArray deltaIncrementTemp(3);
    deltaIncrementTemp.zero();
    FloatArray tempPlasticStrain;
    FloatArray plasticStrain;
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
    double kappa = status->giveTempKappaP();
    double tempKappa = kappa;


    //initialise unknowns
    unknowns.at(1) = trialStress.at(1);
    unknowns.at(2) = trialShearStressNorm;
    unknowns.at(3) = tempKappa;
    unknowns.at(4) = 0.;

    yieldValue = computeYieldValue(tempStress, tempKappa, transitionFlag, gp);

    //initiate residuals
    residuals.zero();
    residuals.at(4) = yieldValue;

    if ( transitionFlag == 0 ) {
        residualsNorm.at(4) = residuals.at(4) / this->fc;
    } else {
        residualsNorm.at(4) = residuals.at(4) / pow(this->fc, 2.);
    }

    if ( residualsNorm.at(4) < yieldTol ) {//regular return not needed
        returnResult = RR_Elastic;
        return kappa;
    }

    normOfResiduals = residualsNorm.computeNorm();
    //    normOfResiduals  = 1.; //just to get into the loop

    while ( normOfResiduals > yieldTol ) {
        iterationCount++;
        if ( iterationCount == newtonIter ) {
            returnResult = RR_NotConverged;
            return kappa;
        }

        residualsNorm.at(1) = residuals.at(1) / this->fc;
        residualsNorm.at(2) = residuals.at(2) / this->fc;
        residualsNorm.at(3) = residuals.at(3);
        if ( transitionFlag == 0 ) {
            residualsNorm.at(4) = residuals.at(4) / this->fc;
        } else {
            residualsNorm.at(4) = residuals.at(4) / pow(this->fc, 2.);
        }

        normOfResiduals = residualsNorm.computeNorm();

        if ( isnan(normOfResiduals) ) {
            returnResult = RR_NotConverged;
            return kappa;
        }

        if ( normOfResiduals > yieldTol ) {
            computeJacobian(jacobian, tempStress, tempKappa, deltaLambda, transitionFlag, gp);

            if ( computeInverseOfJacobian(inverseOfJacobian, jacobian) ) {
                returnResult = RR_NotConverged;
                return kappa;
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
            computeMVector(mVector, tempStress, tempKappa, transitionFlag, gp);

            residuals.at(1) = tempStress.at(1) - trialStress.at(1) + this->eNormalMean * deltaLambda * mVector.at(1);
            residuals.at(2) = tempShearStressNorm - trialShearStressNorm + this->alphaOne * this->eNormalMean * deltaLambda * mVector.at(2);
            residuals.at(3) = -tempKappa + kappa + deltaLambda * mVector.at(3);
            residuals.at(4) = computeYieldValue(tempStress, tempKappa, transitionFlag, gp);
        } else {
            //Check transition
            double transition = computeTransition(tempKappa, gp);
            if ( transitionFlag == 0 ) {
                //Check if the stress is on the correct side
                if ( tempStress.at(1) < transition - yieldTol ) {
                    transitionFlag = 1;
                    tempKappa = kappa;
                    deltaLambda = 0.;
                    tempStress = trialStress;
                    tempShearStressNorm = sqrt(pow(trialStress.at(2), 2.) + pow(trialStress.at(3), 2.) );

                    unknowns.at(1) = trialStress.at(1);
                    unknowns.at(2) = trialShearStressNorm;
                    unknowns.at(3) = tempKappa;
                    unknowns.at(4) = 0.;

                    residuals.zero();
                    residuals.at(4) = computeYieldValue(tempStress, tempKappa, transitionFlag, gp);
                    normOfResiduals = 1;

                    tempPlasticStrain = plasticStrain;
                    iterationCount = 0;
                } else {  //Accept the stress return
                    this->surfaceType = ST_Shear;
                    returnResult = RR_Converged;
                }
            } else if ( transitionFlag == 1 ) {
                if ( tempStress.at(1) > transition + yieldTol ) {
                    returnResult = RR_NotConverged;
                    return kappa;
                } else {  //accept stress return
                    this->surfaceType = ST_Compression;
                    returnResult = RR_Converged;
                }
            }
        }
    }

    stress = tempStress;

    return tempKappa;
}




void
LatticeBondPlasticity :: computeJacobian(FloatMatrix &answer,
                                         const FloatArray &stress,
                                         const double kappa,
                                         const double deltaLambda,
                                         int transitionFlag,
                                         GaussPoint *gp)
{
    //Variables
    FloatArray fVector(3);
    FloatArray mVector(3);
    FloatMatrix dMMatrix(3, 3);

    computeDMMatrix(dMMatrix, stress, kappa, transitionFlag, gp);

    computeMVector(mVector, stress, kappa, transitionFlag, gp);
    computeFVector(fVector, stress, kappa, transitionFlag, gp);

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


double
LatticeBondPlasticity :: computeYieldValue(const FloatArray &stress,
                                           const double kappa,
                                           const int transitionFlag,
                                           GaussPoint *gp)
{
    double yieldValue = 0;

    double shift = computeShift(kappa);
    double paramA = computeParamA(kappa);

    double shearNorm = sqrt(pow(stress.at(2), 2.) + pow(stress.at(3), 2.) );

    if ( transitionFlag == 0 ) {//friction
        yieldValue = shearNorm + this->frictionAngleOne * stress.at(1);
    } else {
        yieldValue = pow(shearNorm, 2.) + pow(stress.at(1) + shift, 2.) / pow(this->frictionAngleTwo, 2.) - pow(paramA, 2.) / pow(this->frictionAngleTwo, 2.);
    }

    return yieldValue;
}




double
LatticeBondPlasticity :: computeTransition(const double kappa,
                                           GaussPoint *gp)
{
    double paramA = computeParamA(kappa);

    double transition = -paramA / ( this->frictionAngleTwo * this->frictionAngleOne *
                                    sqrt(1. + pow(this->frictionAngleTwo, 2.) *
                                         pow(this->frictionAngleOne, 2.) ) );
    return transition;
}



void
LatticeBondPlasticity :: computeFVector(FloatArray &answer,
                                        const FloatArray &stress,
                                        const double kappa,
                                        const int transitionFlag,
                                        GaussPoint *gp)
{
    double shearNorm = sqrt(pow(stress.at(2), 2.) + pow(stress.at(3), 2.) );

    double paramA = computeParamA(kappa);
    double shift = computeShift(kappa);

    double dParamADKappa = computeDParamADKappa(kappa);
    double dShiftDKappa = computeDShiftDKappa(kappa);

    if ( transitionFlag == 0 ) {//line
        answer.at(1) = this->frictionAngleOne;
        answer.at(2) = 1.;
        answer.at(3) = 0.;
    } else {  //cap ellipse
        answer.at(1) = 2. * ( stress.at(1) + shift ) / pow(this->frictionAngleTwo, 2.);
        answer.at(2) = 2. * shearNorm;
        answer.at(3) = 2. * ( stress.at(1) + shift ) / pow(this->frictionAngleTwo, 2.) * dShiftDKappa -
                       2. * paramA / pow(this->frictionAngleTwo, 2.) * dParamADKappa;
    }

    return;
}

void
LatticeBondPlasticity :: computeMVector(FloatArray &answer,
                                        const FloatArray &stress,
                                        const double kappa,
                                        const int transitionFlag,
                                        GaussPoint *gp)
{
    double shearNorm = sqrt(pow(stress.at(2), 2.) + pow(stress.at(3), 2.) );

    double shift = computeShift(kappa);

    if ( transitionFlag == 0 ) {
        answer.at(1) = this->flowAngle;
        answer.at(2) = 1.;
        //No hardening for the Mohr-Coulomb
        answer.at(3) = 0.;
    } else {
        answer.at(1) = 2. * ( stress.at(1) + shift ) / pow(this->frictionAngleTwo, 2.);
        answer.at(2) = 2. * shearNorm;
        if ( stress.at(1) < -shift - yieldTol ) {
            answer.at(3) = sqrt(pow(answer.at(1), 2.) );
        } else   {
            answer.at(3) = 0;
        }
    }

    return;
}



void
LatticeBondPlasticity :: computeDMMatrix(FloatMatrix &answer,
                                         const FloatArray &stress,
                                         const double kappa,
                                         const int transitionFlag,
                                         GaussPoint *gp)
{
    FloatArray mVector(3);
    computeMVector(mVector, stress, kappa, transitionFlag, gp);
    double dShiftDKappa = computeDShiftDKappa(kappa);
    double shift = computeShift(kappa);

    if ( transitionFlag == 0 ) {
        answer.zero();
    } else {  //cap ellipse
              //Derivatives of dGDSig
        answer.at(1, 1) = 2. / pow(this->frictionAngleTwo, 2.);
        answer.at(1, 2) = 0.;
        answer.at(1, 3) = 2. * dShiftDKappa / pow(this->frictionAngleTwo, 2.);

        //Derivatives of dGDTau
        answer.at(2, 1) = 0.;
        answer.at(2, 2) = 2.;
        answer.at(2, 3) = 0;

        //Derivates of evolution law
        if ( stress.at(1) < -shift - yieldTol ) {
            answer.at(3, 1) = 1. / mVector.at(3) * mVector.at(1) * answer.at(1, 1);
            answer.at(3, 2) = 0.;
            answer.at(3, 3) = 1. / mVector.at(3) * ( mVector.at(1) * answer.at(1, 3) );
        } else   {
            answer.at(3, 1) = answer.at(3, 2) = answer.at(3, 3) = 0.;
        }
    }

    return;
}


void
LatticeBondPlasticity :: computeAMatrix(FloatMatrix &answer,
                                        const FloatArray &stress,
                                        const double kappa,
                                        const double deltaLambda,
                                        const int transitionFlag,
                                        GaussPoint *gp)
{
    FloatMatrix dMMatrix(3, 3);
    computeDMMatrix(dMMatrix, stress, kappa, transitionFlag, gp);
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

bool
LatticeBondPlasticity :: checkForVertexCase(FloatArray &stress,
                                            GaussPoint *gp)
{
    double shearNorm = sqrt(pow(stress.at(2), 2.) + pow(stress.at(3), 2.) );

    double ratio = 0.;

    if ( shearNorm / this->eNormalMean > 1.e-30 ) {
        ratio = this->alphaOne * stress.at(1) / shearNorm;
    } else {
        ratio = 5. * this->flowAngle;//set to arbitrary value which is greater than flowAngle
    }
    if ( stress.at(1) > 0 && this->flowAngle < ratio ) {
        return true;
    } else {
        return false;
    }
}


void
LatticeBondPlasticity :: performVertexReturn(FloatArray &stress,
                                             GaussPoint *gp)
{
    LatticeBondPlasticityStatus *status = ( LatticeBondPlasticityStatus * ) ( giveStatus(gp) );

    double shearStressNorm = sqrt(pow(stress.at(2), 2.) + pow(stress.at(3), 2.) );
    //Update kappa
    double kappa = status->giveKappaP();
    double deltaKappa = pow(stress.at(1) / this->eNormalMean, 2.) + pow(shearStressNorm / ( this->alphaOne * this->eNormalMean ), 2.);

    double tempKappa = kappa + deltaKappa;

    stress.zero();

    status->setTempKappaP(tempKappa);

    return;
}

int
LatticeBondPlasticity :: computeInverseOfJacobian(FloatMatrix &answer, const FloatMatrix &src)
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



void
LatticeBondPlasticity :: giveRealStressVector(FloatArray &answer,
                                              GaussPoint *gp,
                                              const FloatArray &originalStrain,
                                              TimeStep *atTime)
{
    LatticeBondPlasticityStatus *status = ( LatticeBondPlasticityStatus * ) this->giveStatus(gp);
    FloatArray reducedStrain;

    this->giveStressDependentPartOfStrainVector(reducedStrain, gp, originalStrain, atTime, VM_Total);



    FloatArray strain(3);
    strain.at(1) = reducedStrain.at(1);
    strain.at(2) = reducedStrain.at(2);
    strain.at(3) = reducedStrain.at(3);

    this->performPlasticityReturn(answer, gp, strain, atTime);

    answer.resizeWithValues(6);

    answer.at(4) = reducedStrain.at(4) * this->alphaTwo * this->eNormalMean;
    answer.at(5) = reducedStrain.at(5) * this->alphaTwo * this->eNormalMean;
    answer.at(6) = reducedStrain.at(6) * this->alphaTwo * this->eNormalMean;

    status->letTempStrainVectorBe(originalStrain);
    status->letTempReducedStrainBe(reducedStrain);
    status->letTempStressVectorBe(answer);

    return;
}


LatticeBondPlasticityStatus :: LatticeBondPlasticityStatus(int n, Domain *d, GaussPoint *gp)
    : LatticeMaterialStatus(gp)
{
    kappaP = tempKappaP = 0.;
}

void
LatticeBondPlasticityStatus :: initTempStatus()
{
    LatticeMaterialStatus :: initTempStatus();

    this->tempKappaP = this->kappaP;
}

void
LatticeBondPlasticityStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    LatticeMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { ");
    fprintf(file, "plasticStrains ");
    for ( int k = 1; k <= plasticStrain.giveSize(); k++ ) {
        fprintf(file, "% .8e ", plasticStrain.at(k) );
    }

    fprintf(file, "\n \t");
    fprintf(file, " kappaP %.8e, surfaceType %d\n", this->kappaP, this->surfaceValue);
    fprintf(file, "}\n");
    return;
}


void
LatticeBondPlasticityStatus :: updateYourself(TimeStep *atTime)
{
    LatticeMaterialStatus :: updateYourself(atTime);
    this->kappaP = this->tempKappaP;
}

void
LatticeBondPlasticityStatus :: saveContext(DataStream &stream, ContextMode mode)
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
}



void
LatticeBondPlasticityStatus :: restoreContext(DataStream &stream, ContextMode mode)
{
    contextIOResultType iores;

    LatticeMaterialStatus :: restoreContext(stream, mode);

    if ( ( iores = plasticStrain.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream.read(& kappaP, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }
}
}     // end namespace oofem
