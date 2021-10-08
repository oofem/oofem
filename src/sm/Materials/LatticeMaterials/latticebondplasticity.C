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


LatticeBondPlasticity::LatticeBondPlasticity(int n, Domain *d) : LatticeLinearElastic(n, d)
{}


bool
LatticeBondPlasticity::hasMaterialModeCapability(MaterialMode mode) const
{
    return mode == _3dLattice;
}

double
LatticeBondPlasticity::computeHardening(double kappa) const
{
    return exp( pow(kappa / this->ef, 2.) );
}

double
LatticeBondPlasticity::computeDHardeningDKappa(double kappa) const
{
    return 2. * kappa / ( pow(this->ef, 2.) ) * exp( pow(kappa / this->ef, 2.) );
}


void
LatticeBondPlasticity::initializeFrom(InputRecord &ir)
{
    LatticeLinearElastic::initializeFrom(ir);

    this->yieldTol = 1.e-6;
    IR_GIVE_OPTIONAL_FIELD(ir, this->yieldTol, _IFT_LatticeBondPlasticity_tol);

    newtonIter = 100;
    IR_GIVE_OPTIONAL_FIELD(ir, newtonIter, _IFT_LatticeBondPlasticity_iter);

    numberOfSubIncrements = 10;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfSubIncrements, _IFT_LatticeBondPlasticity_sub);

    IR_GIVE_FIELD(ir, this->fc, _IFT_LatticeBondPlasticity_fc);

    this->frictionAngleOne = 0.2;
    IR_GIVE_OPTIONAL_FIELD(ir, this->frictionAngleOne, _IFT_LatticeBondPlasticity_angle1);

    this->frictionAngleTwo = this->frictionAngleOne;

    this->flowAngle = this->frictionAngleOne;

    this->ef = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, this->ef, _IFT_LatticeBondPlasticity_ef);
}

MaterialStatus *
LatticeBondPlasticity::CreateStatus(GaussPoint *gp) const
{
    return new LatticeBondPlasticityStatus(1, LatticeBondPlasticity::domain, gp);
}

double
LatticeBondPlasticity::computeShift(const double kappa) const
{
    double hardening = computeHardening(kappa);

    double paramA = computeParamA(kappa);

    double shift = this->fc * hardening - paramA;

    return shift;
}


double
LatticeBondPlasticity::computeParamA(const double kappa) const
{
    double hardening = computeHardening(kappa);

    double paramA = this->frictionAngleTwo * this->frictionAngleOne * this->fc * hardening /
                    ( this->frictionAngleTwo * this->frictionAngleOne +
                      sqrt( 1. + pow(this->frictionAngleTwo, 2.) * pow(this->frictionAngleOne, 2.) ) );

    return paramA;
}


double
LatticeBondPlasticity::computeDShiftDKappa(const double kappa) const
{
    double dHardeningDKappa = computeDHardeningDKappa(kappa);

    double dParamADKappa = computeDParamADKappa(kappa);

    double dShiftDKappa = this->fc * dHardeningDKappa - dParamADKappa;

    return dShiftDKappa;
}


double
LatticeBondPlasticity::computeDParamADKappa(const double kappa) const
{
    double dHardeningDKappa = computeDHardeningDKappa(kappa);

    double A2 = this->frictionAngleTwo * this->frictionAngleOne + sqrt( 1. + pow(this->frictionAngleTwo, 2.) * pow(this->frictionAngleOne, 2.) );

    double dA1DKappa = this->frictionAngleTwo * this->frictionAngleOne * this->fc * dHardeningDKappa;

    double dParamADKappa = dA1DKappa / A2;

    return dParamADKappa;
}


FloatArrayF< 6 >
LatticeBondPlasticity::giveReducedStrain(GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< LatticeBondPlasticityStatus * >( this->giveStatus(gp) );
    return status->giveReducedLatticeStrain();
}


FloatArrayF< 3 >
LatticeBondPlasticity::performPlasticityReturn(GaussPoint *gp,
                                               const FloatArrayF< 3 > &strain,
                                               TimeStep *tStep) const
{
    auto status = static_cast< LatticeBondPlasticityStatus * >( this->giveStatus(gp) );

    LatticeBondPlasticity_SurfaceType surfaceType = ST_Unknown;

    //Get tempKappa from the status
    double tempKappa = status->giveTempKappaP();

    /* Get plastic strain vector from status*/
    auto tempPlasticStrain = status->giveTempPlasticLatticeStrain() [ { 0, 1, 2 } ];

    FloatArrayF< 3 >tangent = { this->eNormalMean, this->alphaOne * this->eNormalMean, this->alphaOne * this->eNormalMean };

    /* Compute trial stress*/
    auto stress = mult(tangent, strain - tempPlasticStrain);

    //Introduce variables for subincrementation
    auto oldStrain = this->giveReducedStrain(gp, tStep) [ { 0, 1, 2 } ];

    /* Compute yield value*/
    int transitionFlag = 0;
    double yieldValue = computeYieldValue(stress, tempKappa, transitionFlag, gp);

    double transition = computeTransition(tempKappa, gp);

    //If shear surface is violated then return to this surface.
    //No need to check situation with ellipse
    //If only ellipse is violated then try this.
    if ( yieldValue < this->yieldTol && stress.at(1) < transition - this->yieldTol ) {
        transitionFlag = 1;
        yieldValue = computeYieldValue(stress, tempKappa, 1, gp);
    }

    /* Compute yield value*/
    double error = 0.;
    if ( transitionFlag == 0 ) {
        error = yieldValue / fc;
    } else {
        error = yieldValue / pow(fc, 2.);
    }

    if ( error > this->yieldTol ) {
        if ( checkForVertexCase(stress, gp) ) {
            performVertexReturn(stress, gp);
            surfaceType = ST_Vertex;
        } else {
            int subIncrementFlag = 0;
            double subIncrementCounter = 1;
            auto convergedStrain = oldStrain;
            auto tempStrain = strain;
            auto deltaStrain = strain - oldStrain;
            //To get into the loop
            LatticeBondPlasticity_ReturnResult returnResult = RR_NotConverged;
            while ( returnResult == RR_NotConverged || subIncrementFlag == 1 ) {
                stress = mult(tangent, tempStrain - tempPlasticStrain);
                tempKappa = performRegularReturn(stress, returnResult, surfaceType, yieldValue, transitionFlag, gp);

                if ( returnResult == RR_NotConverged ) {
                    subIncrementFlag = 1;


                    subIncrementCounter *= 2;
                    if ( subIncrementCounter >= numberOfSubIncrements ) {
                        OOFEM_ERROR("Subincrementation = %e did not help. Stop here with numberOfSubIncrements = %d.\n", subIncrementCounter, numberOfSubIncrements);
                    }

                    deltaStrain *= 0.5;
                    tempStrain = convergedStrain + deltaStrain;
                } else if ( returnResult == RR_Converged && subIncrementFlag == 1 ) {
                    tempPlasticStrain.at(1) = tempStrain.at(1) - stress.at(1) / eNormalMean;
                    tempPlasticStrain.at(2) = tempStrain.at(2) - stress.at(2) / ( this->alphaOne * eNormalMean );
                    tempPlasticStrain.at(3) = tempStrain.at(3) - stress.at(3) / ( this->alphaOne * eNormalMean );

                    status->letTempPlasticLatticeStrainBe( assemble< 6 >(tempPlasticStrain, { 0, 1, 2 }) );

                    status->setTempKappaP(tempKappa);

                    subIncrementCounter -= 1;

                    convergedStrain = tempStrain;
                    deltaStrain = strain - convergedStrain;
                    auto deltaStrainIncrement = deltaStrain * ( 1. / subIncrementCounter );
                    tempStrain = convergedStrain + deltaStrainIncrement;
                    if ( subIncrementCounter > 1 ) {
                        subIncrementFlag = 1;
                    } else {
                        subIncrementFlag = 0;
                    }
                    returnResult = RR_NotConverged;
                } else if ( returnResult == RR_Elastic && subIncrementFlag == 1 ) {
                    convergedStrain = tempStrain;
                    deltaStrain = strain - convergedStrain;
                    auto deltaStrainIncrement = deltaStrain * ( 1. / subIncrementCounter );
                    tempStrain = convergedStrain + deltaStrainIncrement;
                    if ( subIncrementCounter > 1 ) {
                        subIncrementFlag = 1;
                    } else {
                        subIncrementFlag = 0;
                    }
                    returnResult = RR_NotConverged;
                } else if ( returnResult == RR_Converged && subIncrementFlag == 0 ) {
                    status->setTempKappaP(tempKappa);
                }
            }
        }
    }

    tempPlasticStrain.at(1) = strain.at(1) - stress.at(1) / eNormalMean;
    tempPlasticStrain.at(2) = strain.at(2) - stress.at(2) / ( this->alphaOne * eNormalMean );
    tempPlasticStrain.at(3) = strain.at(3) - stress.at(3) / ( this->alphaOne * eNormalMean );

    status->letTempPlasticLatticeStrainBe( assemble< 6 >(tempPlasticStrain, { 0, 1, 2 }) );
    status->letTempLatticeStressBe( assemble< 6 >(stress, { 0, 1, 2 }) );
    status->setSurfaceValue(surfaceType);

    return stress;
}

int
LatticeBondPlasticity::giveIPValue(FloatArray &answer,
                                   GaussPoint *gp,
                                   InternalStateType type,
                                   TimeStep *atTime)
{
    auto status = static_cast< LatticeBondPlasticityStatus * >( this->giveStatus(gp) );

    if ( type == IST_PlasticLatticeStrain ) {
        answer = status->givePlasticLatticeStrain();
        return 1;
    } else {
        return LatticeLinearElastic::giveIPValue(answer, gp, type, atTime);
    }
}

double
LatticeBondPlasticity::performRegularReturn(FloatArrayF< 3 > &stress,
                                            LatticeBondPlasticity_ReturnResult &returnResult,
                                            LatticeBondPlasticity_SurfaceType &surfaceType,
                                            double yieldValue,
                                            int transitionFlag,
                                            GaussPoint *gp) const
{
    auto status = static_cast< LatticeBondPlasticityStatus * >( this->giveStatus(gp) );

    double deltaLambda = 0.;

    auto trialStress = stress;
    auto tempStress = trialStress;

    double trialShearStressNorm = norm(trialStress [ { 1, 2 } ]);
    double tempShearStressNorm = trialShearStressNorm;

    auto plasticStrain = status->givePlasticLatticeStrain() [ { 0, 1, 2 } ];
    auto tempPlasticStrain = plasticStrain;

    double thetaTrial = atan2( stress.at(3), stress.at(2) );

    // Do the same for kappa
    double kappa = status->giveTempKappaP();
    double tempKappa = kappa;


    //initialise unknowns
    FloatArrayF< 4 >unknowns;
    unknowns.at(1) = trialStress.at(1);
    unknowns.at(2) = trialShearStressNorm;
    unknowns.at(3) = tempKappa;
    unknowns.at(4) = 0.;

    yieldValue = computeYieldValue(tempStress, tempKappa, transitionFlag, gp);

    //initiate residuals
    FloatArrayF< 4 >residuals;
    residuals.at(4) = yieldValue;

    double normOfResiduals = this->yieldTol * 2 + 1.; //just to get into the loop

    int iterationCount = 0;
    while ( normOfResiduals > this->yieldTol ) {
        iterationCount++;
        if ( iterationCount == newtonIter ) {
            returnResult = RR_NotConverged;
            return kappa;
        }

        FloatArrayF< 4 >residualsNorm;
        residualsNorm.at(1) = residuals.at(1) / this->fc;
        residualsNorm.at(2) = residuals.at(2) / this->fc;
        residualsNorm.at(3) = residuals.at(3);
        if ( transitionFlag == 0 ) {
            residualsNorm.at(4) = residuals.at(4) / this->fc;
        } else {
            residualsNorm.at(4) = residuals.at(4) / pow(this->fc, 2.);
        }

        normOfResiduals = norm(residualsNorm);

        if ( std::isnan(normOfResiduals) ) {
            returnResult = RR_NotConverged;
            return kappa;
        }

        if ( normOfResiduals > this->yieldTol ) {
            auto jacobian = computeJacobian(tempStress, tempKappa, deltaLambda, transitionFlag, gp);

            auto solution = solve_check(jacobian, residuals);
            if ( solution.first ) {
                unknowns -= solution.second;
            } else {
                returnResult = RR_NotConverged;
                return kappa;
            }

            unknowns.at(2) = max(unknowns.at(2), 0.); //Keep rho greater than zero!
            unknowns.at(3) = max(unknowns.at(3), kappa); //Keep deltaKappa greater than zero!
            unknowns.at(4) = max(unknowns.at(4), 0.); //Keep deltaLambda greater than zero!

            /* Update increments final values and DeltaLambda*/
            tempStress.at(1) = unknowns.at(1);
            tempShearStressNorm = unknowns.at(2);
            tempStress.at(2) = tempShearStressNorm * cos(thetaTrial);
            tempStress.at(3) = tempShearStressNorm * sin(thetaTrial);
            tempKappa = unknowns.at(3);
            deltaLambda = unknowns.at(4);

            /* Compute the mVector holding the derivatives of the g function and the hardening function*/
            auto mVector = computeMVector(tempStress, tempKappa, transitionFlag, gp);

            residuals.at(1) = tempStress.at(1) - trialStress.at(1) + this->eNormalMean * deltaLambda * mVector.at(1);
            residuals.at(2) = tempShearStressNorm - trialShearStressNorm + this->alphaOne * this->eNormalMean * deltaLambda * mVector.at(2);
            residuals.at(3) = -tempKappa + kappa + deltaLambda * mVector.at(3);
            residuals.at(4) = computeYieldValue(tempStress, tempKappa, transitionFlag, gp);
        } else {
            //Check transition
            double transition = computeTransition(tempKappa, gp);
            if ( transitionFlag == 0 ) {
                //Check if the stress is on the correct side
                if ( tempStress.at(1) < transition - this->yieldTol ) {
                    transitionFlag = 1;
                    tempKappa = kappa;
                    deltaLambda = 0.;
                    tempStress = trialStress;
                    tempShearStressNorm = norm(trialStress [ { 1, 2 } ]);

                    unknowns.at(1) = trialStress.at(1);
                    unknowns.at(2) = trialShearStressNorm;
                    unknowns.at(3) = tempKappa;
                    unknowns.at(4) = 0.;

                    residuals = zeros< 4 >();
                    residuals.at(4) = computeYieldValue(tempStress, tempKappa, transitionFlag, gp);
                    normOfResiduals = 1;

                    tempPlasticStrain = plasticStrain;
                    iterationCount = 0;
                } else {  //Accept the stress return
                    surfaceType = ST_Shear;
                    returnResult = RR_Converged;
                }
            } else if ( transitionFlag == 1 ) {
                if ( tempStress.at(1) > transition + this->yieldTol ) {
                    returnResult = RR_NotConverged;
                    return kappa;
                } else {  //accept stress return
                    surfaceType = ST_Compression;
                    returnResult = RR_Converged;
                }
            }
        }
    }

    stress = tempStress;

    return tempKappa;
}


FloatMatrixF< 4, 4 >
LatticeBondPlasticity::computeJacobian(const FloatArrayF< 3 > &stress,
                                       const double kappa,
                                       const double deltaLambda,
                                       int transitionFlag,
                                       GaussPoint *gp) const
{
    auto dMMatrix = computeDMMatrix(stress, kappa, transitionFlag, gp);
    auto mVector = computeMVector(stress, kappa, transitionFlag, gp);
    auto fVector = computeFVector(stress, kappa, transitionFlag, gp);

    /* Compute matrix*/
    FloatMatrixF< 4, 4 >jacobian;
    jacobian.at(1, 1) = 1. + this->eNormalMean * deltaLambda * dMMatrix.at(1, 1);
    jacobian.at(1, 2) = this->eNormalMean * deltaLambda * dMMatrix.at(1, 2);
    jacobian.at(1, 3) = this->eNormalMean * deltaLambda * dMMatrix.at(1, 3);
    jacobian.at(1, 4) = this->eNormalMean * mVector.at(1);
    /**/
    jacobian.at(2, 1) = this->alphaOne * this->eNormalMean * deltaLambda * dMMatrix.at(2, 1);
    jacobian.at(2, 2) = 1. + this->alphaOne * this->eNormalMean * deltaLambda * dMMatrix.at(2, 2);
    jacobian.at(2, 3) = this->alphaOne * this->eNormalMean * deltaLambda * dMMatrix.at(2, 3);
    jacobian.at(2, 4) = this->alphaOne * this->eNormalMean * mVector.at(2);
    /**/
    jacobian.at(3, 1) = deltaLambda * dMMatrix.at(3, 1);
    jacobian.at(3, 2) = deltaLambda * dMMatrix.at(3, 2);
    jacobian.at(3, 3) = deltaLambda * dMMatrix.at(3, 3) - 1.;
    jacobian.at(3, 4) = mVector.at(3);
    /**/
    jacobian.at(4, 1) = fVector.at(1);
    jacobian.at(4, 2) = fVector.at(2);
    jacobian.at(4, 3) = fVector.at(3);
    jacobian.at(4, 4) = 0.;

    return jacobian;
}


double
LatticeBondPlasticity::computeYieldValue(const FloatArrayF< 3 > &stress,
                                         const double kappa,
                                         const int transitionFlag,
                                         GaussPoint *gp) const
{
    double shift = computeShift(kappa);
    double paramA = computeParamA(kappa);

    double shearNorm = norm(stress [ { 1, 2 } ]);

    if ( transitionFlag == 0 ) { //friction
        return shearNorm + this->frictionAngleOne * stress.at(1);
    } else {
        return pow(shearNorm, 2.) + pow(stress.at(1) + shift, 2.) / pow(this->frictionAngleTwo, 2.) - pow(paramA, 2.) / pow(this->frictionAngleTwo, 2.);
    }
}


double
LatticeBondPlasticity::computeTransition(const double kappa,
                                         GaussPoint *gp) const
{
    double paramA = computeParamA(kappa);

    double transition = -paramA / ( this->frictionAngleTwo * this->frictionAngleOne *
                                    sqrt( 1. + pow(this->frictionAngleTwo, 2.) *
                                          pow(this->frictionAngleOne, 2.) ) );
    return transition;
}


FloatArrayF< 3 >
LatticeBondPlasticity::computeFVector(const FloatArrayF< 3 > &stress,
                                      const double kappa,
                                      const int transitionFlag,
                                      GaussPoint *gp) const
{
    double shearNorm = norm(stress [ { 1, 2 } ]);

    double paramA = computeParamA(kappa);
    double shift = computeShift(kappa);

    double dParamADKappa = computeDParamADKappa(kappa);
    double dShiftDKappa = computeDShiftDKappa(kappa);

    FloatArrayF< 3 >f;
    if ( transitionFlag == 0 ) {//line
        f.at(1) = this->frictionAngleOne;
        f.at(2) = 1.;
        f.at(3) = 0.;
    } else {  //cap ellipse
        f.at(1) = 2. * ( stress.at(1) + shift ) / pow(this->frictionAngleTwo, 2.);
        f.at(2) = 2. * shearNorm;
        f.at(3) = 2. * ( stress.at(1) + shift ) / pow(this->frictionAngleTwo, 2.) * dShiftDKappa -
                  2. * paramA / pow(this->frictionAngleTwo, 2.) * dParamADKappa;
    }

    return f;
}


FloatArrayF< 3 >
LatticeBondPlasticity::computeMVector(const FloatArrayF< 3 > &stress,
                                      const double kappa,
                                      const int transitionFlag,
                                      GaussPoint *gp) const
{
    double shearNorm = norm(stress [ { 1, 2 } ]);

    double shift = computeShift(kappa);

    FloatArrayF< 3 >m;
    if ( transitionFlag == 0 ) {
        m.at(1) = this->flowAngle;
        m.at(2) = 1.;
        //No hardening for the Mohr-Coulomb
        m.at(3) = 0.;
    } else {
        m.at(1) = 2. * ( stress.at(1) + shift ) / pow(this->frictionAngleTwo, 2.);
        m.at(2) = 2. * shearNorm;
        if ( stress.at(1) < -shift - this->yieldTol ) {
            m.at(3) = fabs( m.at(1) );
        } else {
            m.at(3) = 0;
        }
    }

    return m;
}


FloatMatrixF< 3, 3 >
LatticeBondPlasticity::computeDMMatrix(const FloatArrayF< 3 > &stress,
                                       const double kappa,
                                       const int transitionFlag,
                                       GaussPoint *gp) const
{
    auto mVector = computeMVector(stress, kappa, transitionFlag, gp);
    double dShiftDKappa = computeDShiftDKappa(kappa);
    double shift = computeShift(kappa);

    if ( transitionFlag == 0 ) {
        return FloatMatrixF< 3, 3 >();
    } else {  //cap ellipse
        FloatMatrixF< 3, 3 >dm;
        //Derivatives of dGDSig
        dm.at(1, 1) = 2. / pow(this->frictionAngleTwo, 2.);
        dm.at(1, 2) = 0.;
        dm.at(1, 3) = 2. * dShiftDKappa / pow(this->frictionAngleTwo, 2.);

        //Derivatives of dGDTau
        dm.at(2, 1) = 0.;
        dm.at(2, 2) = 2.;
        dm.at(2, 3) = 0;

        //Derivates of evolution law
        if ( stress.at(1) < -shift - this->yieldTol ) {
            dm.at(3, 1) = 1. / mVector.at(3) * mVector.at(1) * dm.at(1, 1);
            dm.at(3, 2) = 0.;
            dm.at(3, 3) = 1. / mVector.at(3) * ( mVector.at(1) * dm.at(1, 3) );
        } else {
            dm.at(3, 1) = dm.at(3, 2) = dm.at(3, 3) = 0.;
        }
        return dm;
    }
}


FloatMatrixF< 3, 3 >
LatticeBondPlasticity::computeAMatrix(const FloatArrayF< 3 > &stress,
                                      const double kappa,
                                      const double deltaLambda,
                                      const int transitionFlag,
                                      GaussPoint *gp) const
{
    auto dMMatrix = computeDMMatrix(stress, kappa, transitionFlag, gp);
    /* Compute matrix*/
    FloatMatrixF< 3, 3 >a;
    a.at(1, 1) = 1 / this->eNormalMean + deltaLambda * dMMatrix.at(1, 1);
    a.at(1, 2) = deltaLambda * dMMatrix.at(1, 2);
    a.at(1, 3) = deltaLambda * dMMatrix.at(1, 3);
    /**/
    a.at(2, 1) = deltaLambda * dMMatrix.at(2, 1);
    a.at(2, 2) = 1 / ( this->alphaOne * this->eNormalMean ) + deltaLambda * dMMatrix.at(2, 2);
    a.at(2, 3) = deltaLambda * dMMatrix.at(2, 3);
    /**/
    a.at(3, 1) = deltaLambda * dMMatrix.at(3, 1);
    a.at(3, 2) = deltaLambda * dMMatrix.at(3, 2);
    a.at(3, 3) = deltaLambda * dMMatrix.at(3, 3) - 1.;
    return a;
}


bool
LatticeBondPlasticity::checkForVertexCase(FloatArrayF< 3 > &stress, GaussPoint *gp) const
{
    double shearNorm = norm(stress [ { 1, 2 } ]);

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
LatticeBondPlasticity::performVertexReturn(FloatArrayF< 3 > &stress, GaussPoint *gp) const
{
    auto status = static_cast< LatticeBondPlasticityStatus * >( giveStatus(gp) );

    double shearStressNorm = norm(stress [ { 1, 2 } ]);
    //Update kappa
    double kappa = status->giveKappaP();
    double deltaKappa = pow(stress.at(1) / this->eNormalMean, 2.) + pow(shearStressNorm / ( this->alphaOne * this->eNormalMean ), 2.);

    double tempKappa = kappa + deltaKappa;

    stress = zeros< 3 >();

    status->setTempKappaP(tempKappa);
}


FloatArrayF< 6 >
LatticeBondPlasticity::giveLatticeStress3d(const FloatArrayF< 6 > &originalStrain, GaussPoint *gp, TimeStep *tStep)
{
    auto status = static_cast< LatticeBondPlasticityStatus * >( this->giveStatus(gp) );
    status->initTempStatus();

    auto reducedStrain = originalStrain;
    auto thermalStrain = this->computeStressIndependentStrainVector(gp, tStep, VM_Total);
    if ( thermalStrain.giveSize() ) {
        reducedStrain -= FloatArrayF< 6 >(thermalStrain);
    }

    auto stress3 = this->performPlasticityReturn(gp, reducedStrain [ { 0, 1, 2 } ], tStep);

    auto stress = assemble< 6 >(stress3, { 0, 1, 2 });
    stress.at(4) = reducedStrain.at(4) * this->alphaTwo * this->eNormalMean;
    stress.at(5) = reducedStrain.at(5) * this->alphaTwo * this->eNormalMean;
    stress.at(6) = reducedStrain.at(6) * this->alphaTwo * this->eNormalMean;

    status->letTempLatticeStrainBe(originalStrain);
    status->letTempReducedLatticeStrainBe(reducedStrain);
    status->letTempLatticeStressBe(stress);

    return stress;
}


LatticeBondPlasticityStatus::LatticeBondPlasticityStatus(int n, Domain *d, GaussPoint *gp)
    : LatticeMaterialStatus(gp)
{ }


void
LatticeBondPlasticityStatus::initTempStatus()
{
    LatticeMaterialStatus::initTempStatus();
    this->tempKappaP = this->kappaP;
}

void
LatticeBondPlasticityStatus::printOutputAt(FILE *file, TimeStep *tStep) const
{
    LatticeMaterialStatus::printOutputAt(file, tStep);
    fprintf(file, "status { ");
    fprintf(file, "plasticStrains ");
    for ( double s : this->plasticLatticeStrain ) {
        fprintf(file, "% .8e ", s);
    }

    fprintf(file, "\n \t");
    fprintf(file, " kappaP %.8e, surfaceType %d\n", this->kappaP, this->surfaceValue);
    fprintf(file, "}\n");
}


void
LatticeBondPlasticityStatus::updateYourself(TimeStep *atTime)
{
    LatticeMaterialStatus::updateYourself(atTime);
    this->kappaP = this->tempKappaP;
}


void
LatticeBondPlasticityStatus::saveContext(DataStream &stream, ContextMode mode)
{
    LatticeMaterialStatus::saveContext(stream, mode);

    if ( !stream.write(& kappaP, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }
}


void
LatticeBondPlasticityStatus::restoreContext(DataStream &stream, ContextMode mode)
{
    LatticeMaterialStatus::restoreContext(stream, mode);

    if ( !stream.read(& kappaP, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }
}
}     // end namespace oofem
