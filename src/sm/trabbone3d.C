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

#include "trabbone3d.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "structuralcrosssection.h"
#include "mathfem.h"
#include "internalstatetype.h"
#include "contextioerr.h"

namespace oofem {

TrabBone3D :: TrabBone3D(int n, Domain *d) : StructuralMaterial(n, d)
{}


int
TrabBone3D :: hasMaterialModeCapability(MaterialMode mode)
{
    if ( mode == _3dMat ) {
        return 1;
    }

    return 0;
}


void TrabBone3D :: computePlasStrainEnerDensity(GaussPoint *gp, const FloatArray &totalStrain, const FloatArray &totalStress)
{
    TrabBone3DStatus *status = ( TrabBone3DStatus * ) this->giveStatus(gp);


    double tsed, tempPSED, tempTSED, tempESED;
    FloatArray oldTotalDef, tempPlasDef, oldStress;

    oldTotalDef = status->giveStrainVector();

    tsed = status->giveTSED();
    tempPlasDef = * status->giveTempPlasDef();
    oldStress = status->giveTempStressVector();

    FloatArray elStrain, deltaStrain, tmpStress;

    elStrain.beDifferenceOf(totalStrain, tempPlasDef);
    deltaStrain.beDifferenceOf(totalStrain, oldTotalDef);

    ///@todo Is it really supposed to be totalStress + oldStress? That is what the old code did, but it doesn't seem logical.
    // deltaStress.beDifferenceOf(totalStress, oldStress); is what i would have expected.
    // I can't tell if it is correct, since tempTSED isn't documented at all. // Mikael
    tmpStress = totalStress;
    tmpStress.add(oldStress);

    tempESED = 0.5 * elStrain.dotProduct(totalStress);
    tempTSED = tsed + 0.5 * deltaStrain.dotProduct(tmpStress);
    tempPSED = tempTSED - tempESED;

    status->setTempTSED(tempTSED);
    status->setTempPSED(tempPSED);
}


void
TrabBone3D :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                            MatResponseForm form, MatResponseMode mode, GaussPoint *gp,
                                            TimeStep *atTime)
{
    double tempDam, beta, tempKappa, kappa;
    FloatArray tempEffectiveStress, tempTensor2, prodTensor, plasFlowDirec;
    FloatMatrix elasticity, compliance, SSaTensor, secondTerm, thirdTerm, tangentMatrix;
    TrabBone3DStatus *status = ( TrabBone3DStatus * ) this->giveStatus(gp);

    if ( mode == ElasticStiffness ) {
        this->constructAnisoComplTensor(compliance);
        elasticity.beInverseOf(compliance);
        answer = elasticity;
    } else if ( mode == SecantStiffness || ( (mode == TangentStiffness) && (status->giveNsubsteps()>1))) {
        if (printflag)
            printf("secant\n");
        this->constructAnisoComplTensor(compliance);
        elasticity.beInverseOf(compliance);
        tempDam = status->giveTempDam();
        answer = elasticity;
        answer.times(1.0 - tempDam);
    } else if ( mode == TangentStiffness ) {
        kappa = status->giveKappa();
        tempKappa = status->giveTempKappa();
        tempDam = status->giveTempDam();

        if ( tempKappa > kappa ) { // plastic loading
            // Imports
            tempEffectiveStress = * status->giveTempEffectiveStress();
            plasFlowDirec = * status->givePlasFlowDirec();
            SSaTensor = * status->giveSSaTensor();
            beta = status->giveBeta();
            // Construction of the dyadic product tensor
            prodTensor.beTProductOf(SSaTensor, plasFlowDirec);
            // Construction of the tangent stiffness third term
            thirdTerm.beDyadicProductOf(tempEffectiveStress, prodTensor);
            thirdTerm.times( -expDam * critDam * exp(-expDam * tempKappa) );
            thirdTerm.times(1. / beta);
            // Construction of the tangent stiffness second term
            tempTensor2.beProductOf(SSaTensor, plasFlowDirec);
            secondTerm.beDyadicProductOf(tempTensor2, prodTensor);
            secondTerm.times(-( 1.0 - tempDam ) / beta);
            // Construction of the tangent stiffness
            tangentMatrix = SSaTensor;
            tangentMatrix.times(1.0 - tempDam);
            tangentMatrix.add(secondTerm);
            tangentMatrix.add(thirdTerm);
            answer = tangentMatrix;
        } else { // elastic behavior with damage
            // Construction of the secant stiffness
            this->constructAnisoComplTensor(compliance);
            elasticity.beInverseOf(compliance);
            answer = elasticity;
            answer.times(1.0 - tempDam);
        }
    }

    status->setSmtrx(answer);
}


double
TrabBone3D :: evaluateCurrentYieldStress(const double kappa)
{
    return ( 1.0 + plasHardFactor * ( 1.0 - exp(-kappa * expPlasHard) ) );
}


double
TrabBone3D :: evaluateCurrentPlasticModulus(const double kappa)
{
    return ( plasHardFactor * expPlasHard * exp(-kappa * expPlasHard) );
}

// for given trialEffectiveStress and elasticity tensor computes:
// tempKappa
// tempEffectiveStress
// tempPlasDef
bool
TrabBone3D :: projectOnYieldSurface(double &tempKappa, FloatArray &tempEffectiveStress, FloatArray &tempPlasDef, const FloatArray &trialEffectiveStress, const FloatMatrix &elasticity, const FloatMatrix &compliance, TrabBone3DStatus *status)
{
    bool convergence;
    int flagLoop;
    double deltaKappa, incKappa, halfSpaceCriterion, beta, tempScalar, norm;
    double plasCriterion, plasModulus, toSolveScalar, errorF, errorR;
    FloatArray incTempEffectiveStress, tempStressDiff, errLoop;
    FloatArray toSolveTensor, plasFlowDirec, yieldDerivative, tempTensor2, tensorFF_S;
    FloatMatrix fabric, SSaTensor, tempTensor4, normAdjust, normedFFFF, derivPlasFlowDirec;

    // Selection of the plastic criterion half-space

    halfSpaceCriterion = ( pow(m1, -2.0 * expq) * trialEffectiveStress.at(1) + pow(m2, -2.0 * expq) * trialEffectiveStress.at(2) + pow( ( 3.0 - ( m1 + m2 ) ), -2.0 * expq ) * trialEffectiveStress.at(3) );

    if ( halfSpaceCriterion < 0.0 ) {
        this->constructAnisoFabricTensor(fabric, 0);
        if (printflag)
            printf(" - ");
    } else {
        this->constructAnisoFabricTensor(fabric, 1);
        if (printflag)
            printf(" + ");
    }

    // Evaluation of the yield function

    tempTensor2.beProductOf(fabric, trialEffectiveStress);
    plasCriterion = sqrt( trialEffectiveStress.dotProduct(tempTensor2) ) - evaluateCurrentYieldStress(tempKappa);

    if ( plasCriterion < rel_yield_tol ) {         // trial stress in elastic domain
        convergence = true;
    } else {         // return to the yield surface needed
        // Initial values

        toSolveTensor.resize(6);
        toSolveTensor.zero();
        toSolveScalar = plasCriterion;
        errorF = plasCriterion;
        errorR = 0;
        SSaTensor = elasticity;

        // evaluation of the norm

        // scaling matrix
        this->constructNormAdjustTensor(normAdjust);
        tempTensor4.beProductOf(normAdjust, fabric);
        normedFFFF.beProductOf(fabric, tempTensor4);
        tempTensor2.beProductOf(normedFFFF, tempEffectiveStress);
        norm = sqrt( tempEffectiveStress.dotProduct(tempTensor2) );

        // product FF.S
        tensorFF_S.beProductOf(fabric, tempEffectiveStress);

        // direction of plastic flow
        plasFlowDirec = tensorFF_S;
        plasFlowDirec.times(1.0 / norm);
        deltaKappa = 0.;

        // iteration loop - solution of a set of nonlinear equations
        flagLoop = 1;
        do {
            plasModulus = evaluateCurrentPlasticModulus(tempKappa + deltaKappa);

            //*************************************
            //Evaluation of the Recursive Equations
            //*************************************

            tempTensor2.beProductOf(SSaTensor, plasFlowDirec);
            beta = plasFlowDirec.dotProduct(tempTensor2);
            beta += sqrt( tempEffectiveStress.dotProduct(tensorFF_S) ) / norm * plasModulus;

            // Construction of the equation of Delta Kappa

            tempTensor2.beProductOf(SSaTensor, toSolveTensor);
            tempScalar = plasFlowDirec.dotProduct(tempTensor2);
            incKappa = ( sqrt( tempEffectiveStress.dotProduct(tensorFF_S) ) / norm * toSolveScalar - tempScalar ) / beta;
            deltaKappa += incKappa;

            // Construction of the equation of stress

            // this is the abaqus-style approach
            if (abaqus){
                tempTensor2.beProductOf(normedFFFF, tempEffectiveStress);
                tempTensor4.beDyadicProductOf(tensorFF_S, tempTensor2);
                tempTensor4.times( -1.0 / ( norm * norm * norm ) );
                derivPlasFlowDirec = fabric;
                derivPlasFlowDirec.times(1.0 / norm);
                derivPlasFlowDirec.add(tempTensor4);

                tempTensor4 = derivPlasFlowDirec;
                tempTensor4.times(deltaKappa);
                tempTensor4.add(compliance);
                SSaTensor.beInverseOf(tempTensor4);
            }

            tempTensor2 = plasFlowDirec;
            tempTensor2.times(incKappa);
            tempTensor2.add(toSolveTensor);
            incTempEffectiveStress.beProductOf(SSaTensor, tempTensor2);
            tempEffectiveStress.subtract(incTempEffectiveStress);

            if ( printflag ) {
                printf("   %d %g %g %g %g\n", flagLoop, tempEffectiveStress.at(1), tempEffectiveStress.at(3), incKappa, deltaKappa);
            }

            //*************************************
            // Evaluation of the f and R
            //*************************************

            // Evaluation of the norm

            tempTensor2.beProductOf(normedFFFF, tempEffectiveStress);
            norm = sqrt( tempEffectiveStress.dotProduct(tempTensor2) );

            // Construction of the direction of plastic flow

            tensorFF_S.beProductOf(fabric, tempEffectiveStress);
            plasFlowDirec = tensorFF_S;
            plasFlowDirec.times(1.0 / norm);

            // Construction of the derivative of the plastic flow

            tempTensor2.beProductOf(normedFFFF, tempEffectiveStress);
            tempTensor4.beDyadicProductOf(tensorFF_S, tempTensor2);
            tempTensor4.times( -1.0 / ( norm * norm * norm ) );
            derivPlasFlowDirec = fabric;
            derivPlasFlowDirec.times(1.0 / norm);
            derivPlasFlowDirec.add(tempTensor4);

            // Construction of the gradient Nabla_S of R and SSa tensor

            tempTensor4 = derivPlasFlowDirec;
            tempTensor4.times(deltaKappa);
            tempTensor4.add(compliance);
            SSaTensor.beInverseOf(tempTensor4);

            // Evaluation of R
            tempStressDiff.beDifferenceOf(tempEffectiveStress, trialEffectiveStress);
            tempTensor2.beProductOf( compliance, tempStressDiff );
            toSolveTensor = tempTensor2;
            toSolveTensor.add(deltaKappa, plasFlowDirec);

            // Evaluation of f

            tempTensor2.beProductOf(fabric, tempEffectiveStress);
            toSolveScalar = sqrt( tempEffectiveStress.dotProduct(tempTensor2) ) - evaluateCurrentYieldStress(tempKappa + deltaKappa);

            //*************************************
            // Evaluation of the error
            //*************************************

            errLoop = toSolveTensor;
            errorR = errLoop.computeNorm();
            errorF = toSolveScalar;
            flagLoop++;
            convergence = ( fabs(errorF) < rel_yield_tol && errorR < strain_tol );
        } while ( flagLoop <= max_num_iter && !convergence );         // end of iterative loop

        if ( convergence ) {
            plasModulus = evaluateCurrentPlasticModulus(tempKappa + deltaKappa);
            tempTensor2.beProductOf(SSaTensor, plasFlowDirec);
            beta = plasFlowDirec.dotProduct(tempTensor2);
            beta += sqrt( tempEffectiveStress.dotProduct(tensorFF_S) ) / norm * plasModulus;
            tempPlasDef.add(deltaKappa, plasFlowDirec);
            tempKappa += deltaKappa;

            status->setBeta(beta);
            status->setPlasFlowDirec(plasFlowDirec);
            status->setSSaTensor(SSaTensor);
        }
    }         // end of else statement

    return convergence;
}


void
TrabBone3D :: performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain)
{
    bool convergence;
    int isubstep, nsubstep;
    double tempKappa;
    FloatArray tempPlasDef, tempEffectiveStress, trialEffectiveStress, tempStrainDiff;
    FloatArray strainAfterSubstep, strainIncrement, strainSubIncrement;
    FloatMatrix elasticity, compliance;

    TrabBone3DStatus *status = ( TrabBone3DStatus * ) this->giveStatus(gp);

    // elastic compliance
    this->constructAnisoComplTensor(compliance);
    // elastic stiffness
    elasticity.beInverseOf(compliance);

    // strain increment
    strainIncrement.beDifferenceOf(totalStrain, status->giveStrainVector());

    // --- loop over number of substeps ---
    // (in most cases, convergence is achieved with 1 substep only)
    nsubstep = status->giveNsubsteps(); // from previous iteration, or 1 in the first iteration
    do {
        // initialize the plastic strain and cumulative plastic strain
        // by values after the previous step
        tempPlasDef = * status->givePlasDef();
        tempKappa = status->giveKappa();
        // divide the strain increment into equal subincrements
        strainSubIncrement.beScaled( 1.0 / nsubstep, strainIncrement);
        // set the strain before the first substep
        strainAfterSubstep = status->giveStrainVector();

        isubstep = 1;
        convergence = true;
        // --- loop over substeps ---
        do {
            // apply the strain subincrement
            strainAfterSubstep.add(strainSubIncrement);
            // evaluate the trial stress
            tempStrainDiff.beDifferenceOf(strainAfterSubstep, tempPlasDef);
            trialEffectiveStress.beProductOf(elasticity, tempStrainDiff);
            tempEffectiveStress = trialEffectiveStress;
            // apply the iterative procedure that solves the system of nonlinear equations
            // consisting of the yield condition and discretized flow rule
            // and evaluates:
            // tempKappa ... cumulative plastic strain at the end of the substep
            // tempEffectiveStress ... effective stress at the end of the substep
            // tempPlasDef ... plastic strain at the end of the substep
            if ( printflag ) {
                printf( "substep %d/%d, strain %g %g\n", isubstep, nsubstep, strainAfterSubstep.at(1), strainAfterSubstep.at(3) );
            }

            convergence = projectOnYieldSurface(tempKappa, tempEffectiveStress, tempPlasDef, trialEffectiveStress, elasticity, compliance, status);
            isubstep++;
        } while ( convergence && isubstep <= nsubstep ); // end of loop over substeps

        nsubstep *= 2;
        if ( printflag ) {
            printf("\n");
        }
    } while ( !convergence && nsubstep <= max_num_substeps ); // end of loop over number of substeps

    if ( convergence ) {
        status->setTempPlasDef(tempPlasDef);
        status->setTempKappa(tempKappa);
        status->setTempEffectiveStress(tempEffectiveStress);
        status->setNsubsteps(nsubstep/2);
    } else {
        _error("No convergence of the stress return algorithm in TrabBone3D :: performPlasticityReturn\n");
    }
}


double
TrabBone3D :: computeDamageParam(double tempKappa, GaussPoint *gp)
{
    double tempDam;
    if ( tempKappa > 0. ) {
        tempDam = critDam * ( 1.0 - exp(-expDam * tempKappa) );
    } else {
        tempDam = 0.0;
    }

    return tempDam;
}


double
TrabBone3D :: computeDamage(GaussPoint *gp,  TimeStep *atTime)
{
    double tempKappa;

    computeCumPlastStrain(tempKappa, gp, atTime);

    double tempDam = computeDamageParam(tempKappa, gp);

    return tempDam;
}


void TrabBone3D :: computeCumPlastStrain(double &tempKappa, GaussPoint *gp, TimeStep *atTime)
{
    TrabBone3DStatus *status = ( TrabBone3DStatus * ) this->giveStatus(gp);
    tempKappa = status->giveTempKappa();
}


void TrabBone3D :: computeDensificationStress(FloatArray &answer, GaussPoint *gp, const FloatArray &totalStrain, TimeStep *atTime)
{
    double J;
    FloatArray Id(6);
    FloatMatrix C(3, 3), invC(3, 3), U(3, 3), invU(3, 3);

    // The following development is valid only with the nonlinear geometries (Green-Lagrange: E=0.5(C-I))
    C.at(1, 1) = 1. + 2. * totalStrain.at(1);
    C.at(2, 2) = 1. + 2. * totalStrain.at(2);
    C.at(3, 3) = 1. + 2. * totalStrain.at(3);
    C.at(1, 2) = C.at(2, 1) = totalStrain.at(6);
    C.at(1, 3) = C.at(3, 1) = totalStrain.at(5);
    C.at(2, 3) = C.at(3, 2) = totalStrain.at(4);
    invC.beInverseOf(C);
    J = sqrt( C.giveDeterminant() );

    if ( J < JCrit ) {
        Id.at(1) = invC.at(1, 1);
        Id.at(2) = invC.at(2, 2);
        Id.at(3) = invC.at(3, 3);
        Id.at(4) = 2 * invC.at(3, 2);
        Id.at(5) = 2 * invC.at(3, 1);
        Id.at(6) = 2 * invC.at(2, 1);
        answer.beScaled( gamDens * pow(J / JCrit - JCrit / J, tDens) * J, Id);
    } else {
        answer.resize(6);
    }
}


// BEGIN: SUBROUTINE FOR EVALUATION OF TOTAL STRESS
// returns real stress vector in 3d stress space of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
void
TrabBone3D :: giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                                   const FloatArray &totalStrain, TimeStep *atTime)
{
    double tempDam;
    FloatArray effStress, densStress;
    TrabBone3DStatus *status = ( TrabBone3DStatus * ) this->giveStatus(gp);
    this->initGpForNewStep(gp);

    // compute effective stress using the plasticity model
    performPlasticityReturn(gp, totalStrain);
    effStress =  * status->giveTempEffectiveStress();

    // evaluate damage variable
    tempDam = computeDamage(gp, atTime);

    // transform effective stress into nominal stress
    answer.beScaled( 1.0 - tempDam, effStress);

    computePlasStrainEnerDensity(gp, totalStrain, answer);

    // add densification stress, if the densificator is activated
    if ( JCrit != 0. ) {
        computeDensificationStress(densStress, gp, totalStrain, atTime);
        answer.add(densStress);
    }

    // store final damage, strain and stress in status
    status->setTempDam(tempDam);
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);
}

//
// END: SUBROUTINE FOR EVALUATION OF TOTAL STRESS
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: MATRIX DEFINITION
void
TrabBone3D :: constructAnisoComplTensor(FloatMatrix &answer)
{
    double m3 = 3. - m1 - m2;
    double eps0k = eps0 * pow(rho, expk);
    double mu0k = mu0 * pow(rho, expk);
    double m1l = pow(m1, expl);
    double m2l = pow(m2, expl);
    double m3l = pow(m3, expl);

    answer.resize(6, 6);

    answer.at(1, 1) = 1 / ( eps0k * m1l * m1l );
    answer.at(2, 2) = 1 / ( eps0k * m2l * m2l );
    answer.at(3, 3) = 1 / ( eps0k * m3l * m3l );
    answer.at(1, 2) = -nu0 / ( eps0k * m1l * m2l );
    answer.at(2, 1) = answer.at(1, 2);
    answer.at(1, 3) = -nu0 / ( eps0k * m1l * m3l );
    answer.at(3, 1) = answer.at(1, 3);
    answer.at(2, 3) = -nu0 / ( eps0k * m2l * m3l );
    answer.at(3, 2) = answer.at(2, 3);
    answer.at(4, 4) = 1 / ( mu0k * m2l * m3l );
    answer.at(5, 5) = 1 / ( mu0k * m3l * m1l );
    answer.at(6, 6) = 1 / ( mu0k * m1l * m2l );
}


void
TrabBone3D :: constructAnisoFabricTensor(FloatMatrix &answer, const int posSignFlag)
{
    double chi0, sig0;
    if ( posSignFlag ) {
        chi0 = chi0Pos;
        sig0 = sig0Pos;
    } else {
        chi0 = chi0Neg;
        sig0 = sig0Neg;
    }

    double m3 = 3. - m1 - m2;
    double sig0p = sig0 * sig0 * pow(rho, 2. * expp);
    double tau0p = tau0 * tau0 * pow(rho, 2. * expp);
    double m1q = pow(m1, 2. * expq);
    double m2q = pow(m2, 2. * expq);
    double m3q = pow(m3, 2. * expq);

    answer.resize(6, 6);

    answer.at(1, 1) = 1 / ( sig0p * m1q * m1q );
    answer.at(2, 2) = 1 / ( sig0p * m2q * m2q );
    answer.at(3, 3) = 1 / ( sig0p * m3q * m3q );
    answer.at(1, 2) = -chi0 / ( sig0p * m1q * m2q );
    answer.at(2, 1) = answer.at(1, 2);
    answer.at(1, 3) = -chi0 / ( sig0p * m1q * m3q );
    answer.at(3, 1) = answer.at(1, 3);
    answer.at(2, 3) = -chi0 / ( sig0p * m2q * m3q );
    answer.at(3, 2) = answer.at(2, 3);
    answer.at(4, 4) = 1 / ( tau0p * m2q * m3q );
    answer.at(5, 5) = 1 / ( tau0p * m3q * m1q );
    answer.at(6, 6) = 1 / ( tau0p * m1q * m2q );
}


void
TrabBone3D :: constructNormAdjustTensor(FloatMatrix &answer)
{
    int i;

    answer.resize(6, 6);
    answer.zero();

    for ( i = 1; i <= 3; i++ ) {
        answer.at(i, i) = 1.;
    }

    for ( i = 4; i <= 6; i++ ) {
        answer.at(i, i) = 0.5;
    }
}


IRResultType
TrabBone3D :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    // Mandatory parameters
    IR_GIVE_FIELD(ir, eps0, IFT_TrabBone3D_eps0, "eps0"); // Macro
    IR_GIVE_FIELD(ir, nu0, IFT_TrabBone3D_nu0, "nu0"); // Macro
    IR_GIVE_FIELD(ir, mu0, IFT_TrabBone3D_mu0, "mu0"); // Macro
    IR_GIVE_FIELD(ir, expk, IFT_TrabBone3D_expk, "expk"); // Macro
    IR_GIVE_FIELD(ir, expl, IFT_TrabBone3D_expl, "expl"); // Macro

    IR_GIVE_FIELD(ir, m1, IFT_TrabBone3D_m1, "m1"); // Macro
    IR_GIVE_FIELD(ir, m2, IFT_TrabBone3D_m2, "m2"); // Macro
    IR_GIVE_FIELD(ir, rho, IFT_TrabBone3D_rho, "rho"); // Macro

    IR_GIVE_FIELD(ir, sig0Pos, IFT_TrabBone3D_sig0Pos, "sig0pos"); // Macro
    IR_GIVE_FIELD(ir, sig0Neg, IFT_TrabBone3D_sig0Neg, "sig0neg"); // Macro
    IR_GIVE_FIELD(ir, chi0Pos, IFT_TrabBone3D_chi0Pos, "chi0pos"); // Macro
    IR_GIVE_FIELD(ir, tau0, IFT_TrabBone3D_tau0, "tau0"); // Macro
    IR_GIVE_FIELD(ir, expq, IFT_TrabBone3D_expq, "expq"); // Macro
    IR_GIVE_FIELD(ir, expp, IFT_TrabBone3D_expp, "expp"); // Macro

    IR_GIVE_FIELD(ir, plasHardFactor, IFT_TrabBone3D_plasHardFactor, "plashardfactor"); // Macro
    IR_GIVE_FIELD(ir, expPlasHard, IFT_TrabBone3D_expPlasHard, "expplashard"); // Macro

    IR_GIVE_FIELD(ir, expDam, IFT_TrabBone3D_expDam, "expdam"); // Macro
    IR_GIVE_FIELD(ir, critDam, IFT_TrabBone3D_critDam, "critdam"); // Macro

    // evaluation of dependent parameter
    chi0Neg = ( sig0Neg * sig0Neg ) / ( sig0Pos * sig0Pos ) * ( chi0Pos + 1 ) - 1;

    // optional densificator parameters
    gamDens = 0.;
    tDens = 1.;
    JCrit = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, gamDens, IFT_TrabBone3D_gamDens, "gamdens"); // Macro
    IR_GIVE_OPTIONAL_FIELD(ir, tDens, IFT_TrabBone3D_tDens, "tdens"); // Macro
    IR_GIVE_OPTIONAL_FIELD(ir, JCrit, IFT_TrabBone3D_JCrit, "jcrit"); // Macro

    // optional control parameters for printing and convergence
    printflag = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, printflag, IFT_TrabBone3D_JCrit, "printflag");
    max_num_iter = 20;
    IR_GIVE_OPTIONAL_FIELD(ir, max_num_iter, IFT_TrabBone3D_JCrit, "max_num_iter");
    max_num_substeps = 1024;
    IR_GIVE_OPTIONAL_FIELD(ir, max_num_substeps, IFT_TrabBone3D_JCrit, "max_num_substeps");
    rel_yield_tol = 1.e-10;
    IR_GIVE_OPTIONAL_FIELD(ir, rel_yield_tol, IFT_TrabBone3D_JCrit, "rel_yield_tol");
    strain_tol = 1.e-10;
    IR_GIVE_OPTIONAL_FIELD(ir, strain_tol, IFT_TrabBone3D_JCrit, "strain_tol");
    abaqus = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, abaqus, IFT_TrabBone3D_JCrit, "abaqus");

    return StructuralMaterial :: initializeFrom(ir);
}


int
TrabBone3D :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    TrabBone3DStatus *status = ( TrabBone3DStatus * ) this->giveStatus(aGaussPoint);
    if ( type == IST_DamageScalar ) {
        answer.resize(1);
        answer.at(1) = status->giveTempDam();
        return 1;
    } else if ( type == IST_PlasticStrainTensor ) {
        answer.resize(6);
        answer = * status->giveTempPlasDef();
        return 1;
    } else if ( type == IST_MaxEquivalentStrainLevel ) {
        answer.resize(1);
        answer.at(1) = status->giveTempKappa();
        return 1;
    } else if ( type == IST_BoneVolumeFraction ) {
        answer.resize(1);
        answer.at(1) = rho;
        return 1;
    } else if ( type == IST_PlasStrainEnerDens ) {
        answer.resize(1);
        answer.at(1) = status->giveTempPSED();
        return 1;
    } else if ( type == IST_ElasStrainEnerDens ) {
        answer.resize(1);
        answer.at(1) = status->giveTempTSED() - status->giveTempPSED();
        return 1;
    } else if ( type == IST_TotalStrainEnerDens ) {
        answer.resize(1);
        answer.at(1) = status->giveTempTSED();
        return 1;
    } else {
        return StructuralMaterial :: giveIPValue(answer, aGaussPoint, type, atTime);
    }
}

InternalStateValueType
TrabBone3D :: giveIPValueType(InternalStateType type)
{
    if ( type == IST_PlasticStrainTensor ) {
        return ISVT_TENSOR_S3;
    } else if ( ( type == IST_DamageScalar ) || ( type == IST_MaxEquivalentStrainLevel ) || ( type == IST_BoneVolumeFraction ) || ( type == IST_PlasStrainEnerDens ) || ( type == IST_ElasStrainEnerDens ) || ( type == IST_TotalStrainEnerDens ) ) {
        return ISVT_SCALAR;
    } else {
        return StructuralMaterial :: giveIPValueType(type);
    }
}

int
TrabBone3D :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode)
{
    if ( type == IST_DamageScalar ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;
    } else if ( type == IST_PlasticStrainTensor ) {
        answer.resize(6);
        answer.at(1) = 1;
        answer.at(2) = 2;
        answer.at(3) = 3;
        answer.at(4) = 4;
        answer.at(5) = 5;
        answer.at(6) = 6;
        return 1;
    } else if ( type == IST_MaxEquivalentStrainLevel ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;
    } else if ( type == IST_BoneVolumeFraction ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;
    } else if ( type == IST_PlasStrainEnerDens ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;
    } else if ( type == IST_ElasStrainEnerDens ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;
    } else if ( type == IST_TotalStrainEnerDens ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;
    } else {
        return StructuralMaterial :: giveIntVarCompFullIndx(answer, type, mmode);
    }
}

int
TrabBone3D :: giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint)
{
    if ( ( type == IST_DamageScalar ) || ( type == IST_MaxEquivalentStrainLevel ) || ( type == IST_BoneVolumeFraction ) || ( type == IST_ElasStrainEnerDens ) || ( type == IST_PlasStrainEnerDens ) || ( type == IST_TotalStrainEnerDens ) ) {
        return 1;
    } else if ( type == IST_PlasticStrainTensor ) {
        return 6;
    } else {
        return StructuralMaterial :: giveIPValueSize(type, aGaussPoint);
    }
}


TrabBone3DStatus :: TrabBone3DStatus(int n, Domain *d, GaussPoint *g) : StructuralMaterialStatus(n, d, g)
{
    kappa = 0.0;
    tempKappa = 0.0;
    tempDam = 0.0;
    tempPSED = 0.0;
    tsed = 0.0;
    tempTSED = 0.0;
    beta = 0.0;
    tempEffectiveStress.resize(6);
    plasDef.resize(6);
    tempPlasDef.resize(6);
    plasFlowDirec.resize(6);
    smtrx.resize(6, 6);
    tangentMatrix.resize(6, 6);
    SSaTensor.resize(6, 6);
    nss = 1;
}


TrabBone3DStatus :: ~TrabBone3DStatus()
{}


double
TrabBone3DStatus :: giveKappa()
{
    return kappa;
}

double
TrabBone3DStatus :: giveTempKappa()
{
    return tempKappa;
}

double
TrabBone3DStatus :: giveDam()
{
    return dam;
}

double
TrabBone3DStatus :: giveTempDam()
{
    return tempDam;
}

double
TrabBone3DStatus :: giveTempPSED()
{
    return tempPSED;
}

double
TrabBone3DStatus :: giveTSED()
{
    return tsed;
}

double
TrabBone3DStatus :: giveTempTSED()
{
    return tempTSED;
}

double
TrabBone3DStatus :: giveBeta()
{
    return beta;
}

const FloatArray *
TrabBone3DStatus :: giveTempEffectiveStress()
{
    return & tempEffectiveStress;
}

const FloatArray *
TrabBone3DStatus :: givePlasFlowDirec()
{
    return & plasFlowDirec;
}

const FloatArray *
TrabBone3DStatus :: givePlasDef()
{
    return & plasDef;
}

const FloatArray *
TrabBone3DStatus :: giveTempPlasDef()
{
    return & tempPlasDef;
}

const FloatMatrix *
TrabBone3DStatus :: giveSSaTensor()
{
    return & SSaTensor;
}


void
TrabBone3DStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { ");
    fprintf( file, "plastrains: %f  %f  %f  %f  %f  %f", this->plasDef.at(1), this->plasDef.at(2), this->plasDef.at(3), this->plasDef.at(4), this->plasDef.at(5), this->plasDef.at(6) );
    fprintf(file, " , kappa %f , dam %f , esed %f , psed %f , tsed %f ", this->tempKappa, this->tempDam, this->tempTSED - this->tempPSED, this->tempPSED, this->tempTSED);
    fprintf(file, "}\n");
}


void
TrabBone3DStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();
}


void
TrabBone3DStatus :: updateYourself(TimeStep *atTime)
{
    StructuralMaterialStatus :: updateYourself(atTime);
    this->kappa = this->tempKappa;
    this->dam = this->tempDam;
    this->tsed = this->tempTSED;
    this->plasDef = this->tempPlasDef;
    nss = 1;
}


contextIOResultType
TrabBone3DStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


contextIOResultType
TrabBone3DStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }


    return CIO_OK;
}


MaterialStatus *TrabBone3D :: CreateStatus(GaussPoint *gp) const
{
    TrabBone3DStatus *status =
        new  TrabBone3DStatus(1, StructuralMaterial :: giveDomain(), gp);
    return status;
}

}
