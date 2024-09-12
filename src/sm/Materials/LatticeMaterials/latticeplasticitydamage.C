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

LatticePlasticityDamage::LatticePlasticityDamage(int n, Domain *d) : LatticeLinearElastic(n, d)
{}


bool
LatticePlasticityDamage::hasMaterialModeCapability(MaterialMode mode) const
{
    return mode == _3dLattice;
}

void
LatticePlasticityDamage::initializeFrom(InputRecord &ir)
{
    LatticeLinearElastic::initializeFrom(ir);

    yieldTol = 1.e-6;
    IR_GIVE_OPTIONAL_FIELD(ir, yieldTol, _IFT_LatticePlasticityDamage_tol);

    newtonIter = 100;
    IR_GIVE_OPTIONAL_FIELD(ir, newtonIter, _IFT_LatticePlasticityDamage_iter);

    numberOfSubIncrements = 10;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfSubIncrements, _IFT_LatticePlasticityDamage_sub);

    IR_GIVE_FIELD(ir, this->ft, _IFT_LatticePlasticityDamage_ft);

    IR_GIVE_FIELD(ir, this->fc, _IFT_LatticePlasticityDamage_fc);

    this->frictionAngleOne = 0.23; //Based on fc = 10*ft and fq= ft
    IR_GIVE_OPTIONAL_FIELD(ir, frictionAngleOne, _IFT_LatticePlasticityDamage_angle1);
    if ( this->frictionAngleOne > 0.5 ) {
        OOFEM_WARNING("Friction angle angle1 very large. Really intended? Default value is 0.23");
    }

    this->frictionAngleTwo = 0.5;
    IR_GIVE_OPTIONAL_FIELD(ir, frictionAngleTwo, _IFT_LatticePlasticityDamage_angle2);

    this->flowAngleOne = this->frictionAngleOne;
    IR_GIVE_OPTIONAL_FIELD(ir, flowAngleOne, _IFT_LatticePlasticityDamage_flow);
    if ( this->flowAngleOne > this->frictionAngleOne ) {
        OOFEM_WARNING("Flow angle flow exceeds friction angle angle1. Set flow equal to angle1.\n");
        this->flowAngleOne = this->frictionAngleOne;
    }

    this->flowAngleTwo = this->frictionAngleTwo;

    softeningType = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, softeningType, _IFT_LatticePlasticityDamage_stype);

    if ( softeningType > 1 ) {
        OOFEM_ERROR("Unknown softening type");
    }

    IR_GIVE_FIELD(ir, wf, _IFT_LatticePlasticityDamage_wf);

    double ftOne = 0.;
    if ( softeningType == 1 ) { //bilinear softening
        ftOne = 0.15 * this->ft;
        IR_GIVE_OPTIONAL_FIELD(ir, ftOne, _IFT_LatticePlasticityDamage_ft1);
        this->ftOneRatio = ftOne / this->ft;
        this->wfOne = 0.1 * this->wf;
        IR_GIVE_OPTIONAL_FIELD(ir, wfOne, _IFT_LatticePlasticityDamage_wf1);
    }

    this->aHard = 0.1;
    IR_GIVE_OPTIONAL_FIELD(ir, aHard, _IFT_LatticePlasticityDamage_ahard);

    this->damageFlag = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, damageFlag, _IFT_LatticePlasticityDamage_damage);
}

double
LatticePlasticityDamage::computeDamageParam(double kappaDOne, double kappaDTwo, GaussPoint *gp, TimeStep *tStep) const
{
    double ftLocal =  giveTensileStrength(gp, tStep);
    double fcLocal =  giveCompressiveStrength(gp, tStep);

    double le = static_cast< LatticeStructuralElement * >( gp->giveElement() )->giveLength();
    double strength;
    if ( ftLocal == 0 ) {
        strength = ( fcLocal - this->frictionAngleOne * this->frictionAngleTwo * ftLocal ) / ( 1 + this->frictionAngleOne * this->frictionAngleTwo );
    } else {
        strength = ftLocal;
    }

    double omega = 0.;
    if ( softeningType == 0 ) { //Default: Exponential
        omega = 0;
        double R = 0.;
        int nite = 0;
        do {
            nite++;
            double help = le * ( kappaDOne + omega * kappaDTwo ) / this->wf;
            double Lhs = le * kappaDTwo / this->wf * strength * exp(-help) - kappaDTwo * this->eNormalMean;
            R = ( 1 - omega ) * kappaDTwo * this->eNormalMean - strength * exp(-help);
            omega -= R / Lhs;
            if ( nite > 40 ) {
                OOFEM_ERROR("computeDamageParam: algorithm not converging");
            }
        } while ( fabs(R) >= 1.e-6 || omega < 0.0 );
    } else if ( softeningType == 1 ) {      //bilinear: Calculate all damage parameters and check which makes sense
        double omegaOne = ( this->eNormalMean * kappaDTwo * this->wfOne - ftLocal * this->wfOne - ( this->ftOneRatio * ftLocal - ftLocal ) * kappaDOne * le ) /
                          ( this->eNormalMean * kappaDTwo * this->wfOne + ( this->ftOneRatio * ftLocal - this->ftOneRatio * ftLocal ) * le * kappaDTwo );
        double helpOne = le * kappaDOne + le * omegaOne * kappaDTwo;


        double omegaTwo = ( this->eNormalMean * kappaDTwo * ( this->wf - this->wfOne ) - this->ftOneRatio * ftLocal * ( this->wf - this->wfOne ) +
                            this->ftOneRatio * ftLocal * kappaDOne * le  - this->ftOneRatio * ftLocal * this->wfOne ) /
                          ( this->eNormalMean * kappaDTwo * ( this->wf - this->wfOne )  - this->ftOneRatio * ftLocal * le * kappaDTwo );
        double helpTwo = le * kappaDOne + le * omegaTwo * kappaDTwo;


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

    return omega;
}


MaterialStatus *
LatticePlasticityDamage::CreateStatus(GaussPoint *gp) const
{
    return new LatticePlasticityDamageStatus(1, LatticePlasticityDamage::domain, gp);
}


double
LatticePlasticityDamage::computeYieldValue(const FloatArrayF< 3 > &stress,
                                           const double kappa,
                                           GaussPoint *gp,
                                           TimeStep *tStep) const
{
    double ftLocal =  giveTensileStrength(gp, tStep);
    double fcLocal =  giveCompressiveStrength(gp, tStep);

    double yieldValue = 0;
    double hardening = computeHardening(kappa, gp);
    double shearNorm = norm(stress [ { 1, 2 } ]);
    double transition = -( fcLocal - this->frictionAngleOne * this->frictionAngleTwo * ftLocal ) / ( 1. + this->frictionAngleOne * this->frictionAngleTwo ) * hardening;

    if ( stress.at(1) >= transition ) { //main ellipse
        yieldValue = pow(shearNorm, 2.) + pow(this->frictionAngleOne, 2.) * pow(stress.at(1), 2.) +
                     2. * ( fcLocal - this->frictionAngleOne * this->frictionAngleTwo * ftLocal ) * pow(this->frictionAngleOne, 2.) / ( 1. + this->frictionAngleOne * this->frictionAngleTwo ) * hardening * stress.at(1) -
                     ( 2. * fcLocal * ftLocal * pow(this->frictionAngleOne, 2.) + ( 1. - this->frictionAngleOne * this->frictionAngleTwo ) * pow(ftLocal, 2.) * pow(this->frictionAngleOne, 2.) ) / ( 1. + this->frictionAngleOne * this->frictionAngleTwo ) * pow(hardening, 2.);
    } else {   //cap ellipse
        yieldValue =  pow(shearNorm, 2.) + pow(stress.at(1) / this->frictionAngleTwo, 2.) +
                     2. * ( fcLocal - this->frictionAngleOne * this->frictionAngleTwo * ftLocal ) / ( pow(this->frictionAngleTwo, 2.) * ( 1. + this->frictionAngleOne * this->frictionAngleTwo ) ) * hardening * stress.at(1) +
                     ( pow(fcLocal, 2.) * ( 1. - pow(this->frictionAngleOne * this->frictionAngleTwo, 2.) ) -
                       2. * fcLocal * ftLocal * this->frictionAngleOne * this->frictionAngleTwo * ( 1. + this->frictionAngleOne * this->frictionAngleTwo ) ) / ( pow(this->frictionAngleTwo, 2.) * pow(1. + this->frictionAngleOne * this->frictionAngleTwo, 2.) ) * pow(hardening, 2.);
    }

    return yieldValue;
}

double
LatticePlasticityDamage::computeHardening(const double kappa, GaussPoint *gp) const
{
    return exp(kappa / this->aHard);
}

double
LatticePlasticityDamage::computeDHardeningDKappa(const double kappa, GaussPoint *gp) const
{
    return 1. / this->aHard * exp(kappa / this->aHard);
}

double
LatticePlasticityDamage::computeDDHardeningDDKappa(const double kappa, GaussPoint *gp) const
{
    return 1. / pow(this->aHard, 2.) * exp(kappa / this->aHard);
}

FloatArrayF< 3 >
LatticePlasticityDamage::computeFVector(const FloatArrayF< 3 > &stress,
                                        const double kappa,
                                        GaussPoint *gp,
                                        TimeStep *tStep) const
{
    double ftLocal =  giveTensileStrength(gp, tStep);
    double fcLocal =  giveCompressiveStrength(gp, tStep);

    double hardening = computeHardening(kappa, gp);
    double dHardeningDKappa = computeDHardeningDKappa(kappa, gp);
    double shearNorm = norm(stress [ { 1, 2 } ]);
    double transition = -( fcLocal - this->frictionAngleOne * this->frictionAngleTwo * ftLocal ) / ( 1. + this->frictionAngleOne * this->frictionAngleTwo ) * hardening;

    FloatArrayF< 3 >f;
    if ( stress.at(1) >= transition ) { //main ellipse
        f.at(1) = 2. * pow(this->frictionAngleOne, 2.) * stress.at(1) + 2. * ( fcLocal - this->frictionAngleOne * this->frictionAngleTwo * ftLocal ) * pow(this->frictionAngleOne, 2.) / ( 1. + this->frictionAngleOne * this->frictionAngleTwo ) * hardening;
        f.at(2) = 2. * shearNorm;
        f.at(3) = 2. * ( fcLocal - this->frictionAngleOne * this->frictionAngleTwo * ftLocal ) * pow(this->frictionAngleOne, 2.) / ( 1. + this->frictionAngleOne * this->frictionAngleTwo ) * stress.at(1) * dHardeningDKappa -
                  2. * ( 2. * fcLocal * ftLocal * pow(this->frictionAngleOne, 2.) + ( 1. - this->frictionAngleOne * this->frictionAngleTwo ) * pow(ftLocal, 2.) * pow(this->frictionAngleOne, 2.) ) / ( 1. + this->frictionAngleOne * this->frictionAngleTwo ) * hardening * dHardeningDKappa;
    } else {   //cap ellipse
        f.at(1) = 2. * stress.at(1) / pow(this->frictionAngleTwo, 2.) + 2. * ( fcLocal - this->frictionAngleOne * this->frictionAngleTwo * ftLocal ) / ( pow(this->frictionAngleTwo, 2.) * ( 1. + this->frictionAngleOne * this->frictionAngleTwo ) ) * hardening;
        f.at(2) = 2. * shearNorm;
        f.at(3) = 2. * ( fcLocal - this->frictionAngleOne * this->frictionAngleTwo * ftLocal ) / ( pow(this->frictionAngleTwo, 2.) * ( 1. + this->frictionAngleOne * this->frictionAngleTwo ) ) * stress.at(1) * dHardeningDKappa +
                  2. * ( pow(fcLocal, 2.) * ( 1. - pow(this->frictionAngleOne * this->frictionAngleTwo, 2.) ) -
                         2. * fcLocal * ftLocal * this->frictionAngleOne * this->frictionAngleTwo * ( 1. + this->frictionAngleOne * this->frictionAngleTwo ) ) /
                  ( pow(frictionAngleTwo, 2.) * pow(1. + this->frictionAngleOne * this->frictionAngleTwo, 2.) ) * hardening * dHardeningDKappa;
    }

    return f;
}

FloatArrayF< 3 >
LatticePlasticityDamage::computeMVector(const FloatArrayF< 3 > &stress,
                                        const double kappa,
                                        GaussPoint *gp,
                                        TimeStep *tStep) const
{
    double ftLocal =  giveTensileStrength(gp, tStep);
    double fcLocal =  giveCompressiveStrength(gp, tStep);

    double hardening = computeHardening(kappa, gp);
    double shearNorm = norm(stress [ { 1, 2 } ]);
    double transition = -( fcLocal - this->flowAngleOne * this->frictionAngleTwo * ftLocal ) / ( 1. + this->flowAngleOne * this->frictionAngleTwo ) * hardening;

    FloatArrayF< 3 >m;
    if ( stress.at(1) >= transition ) { //ellipse
        m.at(1) = 2. * pow(this->flowAngleOne, 2.) * stress.at(1) + 2. * ( fcLocal - this->flowAngleOne * this->flowAngleTwo * ftLocal ) * pow(this->flowAngleOne, 2.) / ( 1. + this->flowAngleOne * this->flowAngleTwo ) * hardening;
        m.at(2) = 2. * shearNorm;
        m.at(3) = fabs( m.at(1) );
    } else {   //circle
        m.at(1) = 2. * stress.at(1) / pow(this->flowAngleTwo, 2.) + 2. * ( fcLocal - this->flowAngleTwo * this->flowAngleOne * ftLocal ) / ( pow(this->flowAngleTwo, 2.) * ( 1. + this->flowAngleOne * this->flowAngleTwo ) ) * hardening;
        m.at(2) = 2. * shearNorm;
        m.at(3) = fabs( m.at(1) );
    }

    return m;
}


FloatMatrixF< 3, 3 >
LatticePlasticityDamage::computeDMMatrix(const FloatArrayF< 3 > &stress, const double kappa, GaussPoint *gp, TimeStep *tStep) const
{
    double ftLocal =  giveTensileStrength(gp, tStep);
    double fcLocal =  giveCompressiveStrength(gp, tStep);

    double hardening = computeHardening(kappa, gp);
    double dHardeningDKappa = computeDHardeningDKappa(kappa, gp);
    double transition = -( fcLocal - this->flowAngleOne * this->flowAngleTwo * this->ft ) / ( 1. + this->flowAngleOne * this->flowAngleTwo ) * hardening;

    FloatMatrixF< 3, 3 >dm;
    if ( stress.at(1) >= transition ) { //main ellipse
        //Derivatives of dGDSig
        dm.at(1, 1) = 2. * pow(this->flowAngleOne, 2.);
        dm.at(1, 2) = 0.;
        dm.at(1, 3) = 2. * ( fcLocal - this->flowAngleOne * this->flowAngleTwo * ftLocal ) * pow(this->flowAngleOne, 2.) / ( 1. + this->flowAngleOne * this->flowAngleTwo ) * dHardeningDKappa;

        //Derivatives of dGDTau
        dm.at(2, 1) = 0.;
        dm.at(2, 2) = 2.;
        dm.at(2, 3) = 0;

        //Derivates of evolution law
        dm.at(3, 1) = dm.at(1, 1);
        dm.at(3, 2) = dm.at(1, 2);
        dm.at(3, 3) = dm.at(1, 3);
    } else {   //cap ellipse
               //Derivatives of dGDSig
        dm.at(1, 1) = 2. / pow(this->flowAngleTwo, 2.);
        dm.at(1, 2) = 0.;
        dm.at(1, 3) = 2. * ( fcLocal - this->flowAngleOne * this->flowAngleTwo * ftLocal ) / ( pow(this->flowAngleTwo, 2.) * ( 1. + this->flowAngleOne * this->flowAngleTwo ) ) * dHardeningDKappa;

        //Derivatives of dGDTau
        dm.at(2, 1) = 0.;
        dm.at(2, 2) = 2.;
        dm.at(2, 3) = 0;

        //Derivates of evolution law
        dm.at(3, 1) = -dm.at(1, 1);
        dm.at(3, 2) = -dm.at(1, 2);
        dm.at(3, 3) = -dm.at(1, 3);
    }
    return dm;
}


FloatArrayF< 6 >
LatticePlasticityDamage::giveReducedStrain(GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< LatticePlasticityDamageStatus * >( this->giveStatus(gp) );
    return status->giveReducedLatticeStrain();
}


FloatArrayF< 6 >
LatticePlasticityDamage::performPlasticityReturn(GaussPoint *gp,
                                                 const FloatArrayF< 6 > &reducedStrain,
                                                 TimeStep *tStep) const
{
    LatticePlasticityDamage_ReturnResult returnResult = RR_Unknown;

    double fcLocal =  giveCompressiveStrength(gp, tStep);
    auto status = static_cast< LatticePlasticityDamageStatus * >( this->giveStatus(gp) );

    //Get kappa from the status
    double tempKappa = status->giveKappaP();

    //Subset of reduced strain.
    //Rotational components are not used for plasticity return
    auto strain = reducedStrain [ { 0, 1, 2 } ];


    /* Get plastic strain vector from status*/
    auto tempPlasticStrain = status->givePlasticLatticeStrain() [ { 0, 1, 2 } ];

    FloatArrayF< 3 >tangent = { this->eNormalMean, this->alphaOne * this->eNormalMean, this->alphaOne * this->eNormalMean };
    /* Compute trial stress*/
    auto stress = mult(tangent, strain - tempPlasticStrain);

    //Introduce variables for subincrementation
    //Only _3dLattice is possible

    auto oldStrain = this->giveReducedStrain(gp, tStep) [ { 0, 1, 2 } ];

    /* Compute yield value*/
    double yieldValue = computeYieldValue(stress, tempKappa, gp, tStep);
    int subIncrementCounter = 0;

    /* Check yield condition, i.e. if the yield value is less than the yield tolerance.
     * If yield condition is valid. Do perform regular return (closest point return)*/

    if ( yieldValue / pow(fcLocal, 2.) > yieldTol ) {
        // introduce a subincrementation flag
        int subIncrementFlag = 0;
        auto convergedStrain = oldStrain;
        auto tempStrain = strain;
        auto deltaStrain = strain - oldStrain;
        //To get into the loop
        returnResult = RR_NotConverged;
        while ( returnResult == RR_NotConverged || subIncrementFlag == 1 ) {
            stress = mult(tangent, tempStrain - tempPlasticStrain);

            tempKappa = performRegularReturn(stress, returnResult, yieldValue, gp, tStep);

            if ( returnResult == RR_NotConverged ) {
                subIncrementCounter++;
                if ( subIncrementCounter > numberOfSubIncrements ) {
                    OOFEM_LOG_INFO( "Unstable element %d \n", gp->giveElement()->giveGlobalNumber() );
                    OOFEM_LOG_INFO("Yield value %e \n", yieldValue);
                    OOFEM_LOG_INFO( "ConvergedStrain value %e %e %e\n", convergedStrain.at(1), convergedStrain.at(2), convergedStrain.at(3) );
                    OOFEM_LOG_INFO( "tempStrain value %e %e %e\n", tempStrain.at(1), tempStrain.at(2), tempStrain.at(3) );
                    OOFEM_LOG_INFO( "deltaStrain value %e %e %e\n", deltaStrain.at(1), deltaStrain.at(2), deltaStrain.at(3) );
                    OOFEM_LOG_INFO( "targetstrain value %e %e %e\n", strain.at(1), strain.at(2), strain.at(3) );

                    OOFEM_ERROR("LatticePlasticityDamage :: performPlasticityReturn - Could not reach convergence with small deltaStrain, giving up.");
                }

                subIncrementFlag = 1;
                deltaStrain *= 0.5;
                tempStrain = convergedStrain + deltaStrain;
            } else if ( returnResult == RR_Converged && subIncrementFlag == 1 ) {
                tempPlasticStrain.at(1) = tempStrain.at(1) - stress.at(1) / eNormalMean;
                tempPlasticStrain.at(2) = tempStrain.at(2) - stress.at(2) / ( this->alphaOne * eNormalMean );
                tempPlasticStrain.at(3) = tempStrain.at(3) - stress.at(3) / ( this->alphaOne * eNormalMean );

                status->letTempPlasticLatticeStrainBe( assemble< 6 >(tempPlasticStrain, { 0, 1, 2 }) );
                status->setTempKappaP(tempKappa);

                subIncrementFlag = 0;
                returnResult = RR_NotConverged;
                convergedStrain = tempStrain;
                deltaStrain = strain - convergedStrain;
                tempStrain = strain;
                subIncrementCounter = 0;
            } else {
                status->setTempKappaP(tempKappa);
            }
        }
    }


    tempPlasticStrain.at(1) = strain.at(1) - stress.at(1) / eNormalMean;
    tempPlasticStrain.at(2) = strain.at(2) - stress.at(2) / ( this->alphaOne * eNormalMean );
    tempPlasticStrain.at(3) = strain.at(3) - stress.at(3) / ( this->alphaOne * eNormalMean );

    status->letTempPlasticLatticeStrainBe( assemble< 6 >(tempPlasticStrain, { 0, 1, 2 }) );


    //    status->letTempLatticeStressBe(assemble< 6 >(stress, { 0, 1, 2 }) );

    auto answer = assemble< 6 >(stress, { 0, 1, 2 });
    answer.at(4) = this->alphaTwo * this->eNormalMean * reducedStrain.at(4);
    answer.at(5) = this->alphaTwo * this->eNormalMean * reducedStrain.at(5);
    answer.at(6) = this->alphaTwo * this->eNormalMean * reducedStrain.at(6);

    return answer;
}


double
LatticePlasticityDamage::performRegularReturn(FloatArrayF< 3 > &stress,
                                              LatticePlasticityDamage_ReturnResult &returnResult,
                                              double yieldValue,
                                              GaussPoint *gp,
                                              TimeStep *tStep) const
{
    double fcLocal =  giveCompressiveStrength(gp, tStep);

    auto status = static_cast< LatticePlasticityDamageStatus * >( this->giveStatus(gp) );

    double deltaLambda = 0.;

    auto trialStress = stress;
    auto tempStress = trialStress;

    double trialShearStressNorm = norm(trialStress [ { 1, 2 } ]);

    double tempShearStressNorm = trialShearStressNorm;

    double thetaTrial = atan2( stress.at(3), stress.at(2) );

    // Do the same for kappa
    double kappa = status->giveKappaP();
    double tempKappa = kappa;

    //initialise unknowns
    FloatArrayF< 4 >unknowns;
    unknowns.at(1) = trialStress.at(1);
    unknowns.at(2) = trialShearStressNorm;
    unknowns.at(3) = tempKappa;
    unknowns.at(4) = 0.;

    // Look at the magnitudes of the residuals. You have to scale the yieldValue down.
    yieldValue = computeYieldValue(tempStress, tempKappa, gp, tStep);

    //initiate residuals
    FloatArrayF< 4 >residuals;
    residuals.at(4) = yieldValue;

    double normOfResiduals  = 1.; //just to get into the loop

    int iterationCount = 0;
    while ( normOfResiduals > yieldTol ) {
        iterationCount++;
        if ( iterationCount == newtonIter ) {
            returnResult = RR_NotConverged;
            return 0.;
        }

        //Normalize residuals. Think about it more.
        FloatArrayF< 4 >residualsNorm;
        residualsNorm.at(1) = residuals.at(1) / fcLocal;
        residualsNorm.at(2) = residuals.at(2) / fcLocal;
        residualsNorm.at(3) = residuals.at(3);
        residualsNorm.at(4) = residuals.at(4) / pow(fcLocal, 2.);

        normOfResiduals = norm(residualsNorm);

        //First check if return has failed
        if ( std::isnan(normOfResiduals) ) {
            returnResult = RR_NotConverged;
            return 0.;
        }

        if ( normOfResiduals > yieldTol ) {
            // Test to run newton iteration using inverse of Jacobian
            auto jacobian = computeJacobian(tempStress, tempKappa, deltaLambda, gp, tStep);

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
            auto mVector = computeMVector(tempStress, tempKappa, gp, tStep);

            residuals.at(1) = tempStress.at(1) - trialStress.at(1) + this->eNormalMean * deltaLambda * mVector.at(1);
            residuals.at(2) = tempShearStressNorm - trialShearStressNorm + this->alphaOne * this->eNormalMean * deltaLambda * mVector.at(2);
            residuals.at(3) = -tempKappa + kappa + deltaLambda * mVector.at(3);
            residuals.at(4) = computeYieldValue(tempStress, tempKappa, gp, tStep);
        }
    }

    returnResult = RR_Converged;

    stress = tempStress;

    return tempKappa;
}



FloatMatrixF< 4, 4 >
LatticePlasticityDamage::computeJacobian(const FloatArrayF< 3 > &stress,
                                         const double kappa,
                                         const double deltaLambda,
                                         GaussPoint *gp,
                                         TimeStep *tStep) const
{
    auto dMMatrix = computeDMMatrix(stress, kappa, gp, tStep);
    auto mVector = computeMVector(stress, kappa, gp, tStep);
    auto fVector = computeFVector(stress, kappa, gp, tStep);

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
LatticePlasticityDamage::give(int aProperty, GaussPoint *gp) const
{
    this->giveStatus(gp);

    double answer;
    if ( RandomMaterialExtensionInterface::give(aProperty, gp, answer) ) {
        return answer;
    } else if ( aProperty == fc_strength ) {
        return 1.;
    } else if ( aProperty == ft_strength ) {
        return 1.;
    } else {
        return LatticeLinearElastic::give(aProperty, gp);
    }
}


FloatArrayF< 6 >
LatticePlasticityDamage::giveLatticeStress3d(const FloatArrayF< 6 > &originalStrain, GaussPoint *gp, TimeStep *tStep)
{
    auto status = static_cast< LatticePlasticityDamageStatus * >( this->giveStatus(gp) );
    status->initTempStatus();

    auto reducedStrain = originalStrain;
    auto thermalStrain = this->computeStressIndependentStrainVector(gp, tStep, VM_Total);
    if ( thermalStrain.giveSize() ) {
        reducedStrain -= FloatArrayF< 6 >(thermalStrain);
    }

    auto stress = this->performPlasticityReturn(gp, reducedStrain, tStep);

    double omega = 0.;
    if ( damageFlag == 1 ) {
        this->performDamageEvaluation(gp, reducedStrain, tStep);
        omega = status->giveTempDamage();
    }

    stress *= ( 1. - omega );

    status->letTempLatticeStrainBe(originalStrain);
    status->letTempReducedLatticeStrainBe(reducedStrain);
    status->letTempLatticeStressBe(stress);

    return stress;
}


void
LatticePlasticityDamage::performDamageEvaluation(GaussPoint *gp, FloatArrayF< 6 > &reducedStrain, TimeStep *tStep) const
{
    double ftLocal =  giveTensileStrength(gp, tStep);
    double fcLocal =  giveCompressiveStrength(gp, tStep);

    auto status = static_cast< LatticePlasticityDamageStatus * >( this->giveStatus(gp) );

    double le = static_cast< LatticeStructuralElement * >( gp->giveElement() )->giveLength();

    auto tempPlasticStrain = status->giveTempPlasticLatticeStrain();
    auto plasticStrain = status->givePlasticLatticeStrain();
    auto deltaPlasticStrain = tempPlasticStrain - plasticStrain;

    //Compute the damage history variable
    double omega = status->giveDamage();
    double tempKappaDOne = status->giveKappaDOne();
    double tempKappaDTwo = status->giveKappaDTwo();

    // Compute deltaKappaP from the plastic Strain
    double tempKappaP = status->giveTempKappaP();
    double kappaP = status->giveKappaP();
    double deltaKappaP = tempKappaP - kappaP;
    double hardening = computeHardening(tempKappaP, gp);

    double strength = ( fcLocal - this->frictionAngleOne * this->frictionAngleTwo * ftLocal ) / ( 1 + this->frictionAngleOne * this->frictionAngleTwo );

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

            if ( ftLocal == 0 ) {
                tempKappaDTwo = hardening * strength / this->eNormalMean;
            } else {
                tempKappaDTwo = hardening * ftLocal / this->eNormalMean;
            }

            // evaluate damage parameter
            omega = this->computeDamageParam(tempKappaDOne, tempKappaDTwo, gp, tStep);

            //threshold for crack patterns
            if ( ( tempKappaDOne + omega * tempKappaDTwo ) * le > 0. ) {
                status->setTempCrackFlag(1);
            } else {
                status->setTempCrackFlag(0);
            }
        }
    }

    //Create history variables for the damage strain
    FloatArrayF< 6 >elasticReducedStrain = reducedStrain - tempPlasticStrain;
    FloatArrayF< 6 >tempDamageLatticeStrain = omega * elasticReducedStrain;
    status->letTempDamageLatticeStrainBe(tempDamageLatticeStrain);

    double crackWidth = norm( tempPlasticStrain + omega * ( reducedStrain - tempPlasticStrain ) ) * le;

    //TODO: Compute dissipation
    // double tempDissipation = status->giveDissipation();
    // double tempDeltaDissipation = 0.;
    // tempDeltaDissipation = computeDeltaDissipation(omega, reducedStrain, gp, atTime);
    // tempDissipation += tempDeltaDissipation;

    status->setTempKappaDOne(tempKappaDOne);
    status->setTempKappaDTwo(tempKappaDTwo);
    status->setTempDamage(omega);
    status->setTempCrackWidth(crackWidth);
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
LatticePlasticityDamage::give3dLatticeStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    auto elastic = LatticeLinearElastic::give3dLatticeStiffnessMatrix(mode, gp, tStep);

    if ( mode == ElasticStiffness ) {
        return elastic;
    } else if ( ( mode == SecantStiffness ) || ( mode == TangentStiffness ) ) {
        auto status = static_cast< LatticePlasticityDamageStatus * >( this->giveStatus(gp) );
        double omega = min(status->giveTempDamage(), 0.99999);
        return elastic * ( 1. - omega );
    } else {
        OOFEM_ERROR("Unsupported stiffness mode\n");
    }
}

int
LatticePlasticityDamage::giveIPValue(FloatArray &answer,
                                     GaussPoint *gp,
                                     InternalStateType type,
                                     TimeStep *atTime)
{
    auto status = static_cast< LatticePlasticityDamageStatus * >( this->giveStatus(gp) );

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
    } else if ( type == IST_CrackedFlag ) {
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
    } else if ( type == IST_PlasticLatticeStrain ) {
        answer = status->givePlasticLatticeStrain();
        return 1;
    } else if ( type == IST_TensileStrength ) {
        answer.resize(1);
        answer.at(1);
        answer.at(1) = giveTensileStrength(gp, atTime);
        return 1;
    } else {
        return LatticeLinearElastic::giveIPValue(answer, gp, type, atTime);
    }
}


LatticePlasticityDamageStatus::LatticePlasticityDamageStatus(int n, Domain *d, GaussPoint *g) :  LatticeMaterialStatus(g)
{ }

void
LatticePlasticityDamageStatus::initTempStatus()
{
    LatticeMaterialStatus::initTempStatus();
    this->tempKappaP = this->kappaP;
    this->tempKappaDOne = this->kappaDOne;
    this->tempKappaDTwo = this->kappaDTwo;
    this->tempDamage = this->damage;
}

void
LatticePlasticityDamageStatus::printOutputAt(FILE *file, TimeStep *tStep) const
{
    LatticeMaterialStatus::printOutputAt(file, tStep);

    fprintf(file, "plasticStrains ");
    for ( double s : this->plasticLatticeStrain ) {
        fprintf(file, "% .8e ", s);
    }

    fprintf(file, ", kappaP %.8e, kappaDOne %.8e, kappaDTwo %.8e, damage %.8e, deltaDissipation %.8e, dissipation %.8e, crackFlag %d, crackWidth %.8e \n", this->kappaP, this->kappaDOne,  this->kappaDTwo, this->damage, this->deltaDissipation, this->dissipation, this->crackFlag, this->crackWidth);
}


void
LatticePlasticityDamageStatus::updateYourself(TimeStep *atTime)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reached equilibrium.
//
{
    LatticeMaterialStatus::updateYourself(atTime);
    this->kappaP = this->tempKappaP;
    this->kappaDOne = this->tempKappaDOne;
    this->kappaDTwo = this->tempKappaDTwo;
    this->damage = this->tempDamage;
}


void
LatticePlasticityDamageStatus::saveContext(DataStream &stream, ContextMode mode)
//
// saves full information stored in this Status
// no temp variables stored
//
{
    LatticeMaterialStatus::saveContext(stream, mode);

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
LatticePlasticityDamageStatus::restoreContext(DataStream &stream, ContextMode mode)
//
// restores full information stored in stream to this Status
//
{
    LatticeMaterialStatus::restoreContext(stream, mode);

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
