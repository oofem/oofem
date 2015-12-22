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

#include "trabbone3d.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "internalstatetype.h"
#include "contextioerr.h"
#include "intarray.h"
#include "datastream.h"
#include "classfactory.h"


namespace oofem {
REGISTER_Material(TrabBone3D);

TrabBone3D :: TrabBone3D(int n, Domain *d) : StructuralMaterial(n, d)
{ }


void TrabBone3D :: computePlasStrainEnerDensity(GaussPoint *gp, const FloatArray &totalStrain, const FloatArray &totalStress)
{
    TrabBone3DStatus *status = static_cast< TrabBone3DStatus * >( this->giveStatus(gp) );

    double tsed, tempPSED, tempTSED, tempESED;
    FloatArray oldTotalDef, tempPlasDef, oldStress;

    oldTotalDef = status->giveStrainVector();

    tsed = status->giveTSED();
    tempPlasDef = status->giveTempPlasDef();
    oldStress = status->giveTempStressVector();

    FloatArray elStrain, deltaStrain, tmpStress;

    elStrain.beDifferenceOf(totalStrain, tempPlasDef);
    deltaStrain.beDifferenceOf(totalStrain, oldTotalDef);

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
                                            MatResponseMode mode, GaussPoint *gp,
                                            TimeStep *tStep)
{
    double tempDam, beta, tempKappa, kappa;
    FloatArray tempEffectiveStress, tempTensor2, prodTensor, plasFlowDirec;
    FloatMatrix elasticity, compliance, SSaTensor, secondTerm, thirdTerm, tangentMatrix;
    TrabBone3DStatus *status = static_cast< TrabBone3DStatus * >( this->giveStatus(gp) );

    if ( mode == ElasticStiffness ) {
        this->constructAnisoComplTensor(compliance);
        //this->constructAnisoStiffnessTensor(elasticity);
        elasticity.beInverseOf(compliance);
        answer = elasticity;
    } else if ( mode == SecantStiffness ) {
        if ( printflag ) {
            printf("secant\n");
        }

        this->constructAnisoStiffnessTensor(elasticity);

        this->constructAnisoComplTensor(compliance);
        elasticity.beInverseOf(compliance);

        tempDam = status->giveTempDam();
        answer = elasticity;
        answer.times(1.0 - tempDam);
    } else if ( mode == TangentStiffness ) {
        kappa = status->giveKappa();
        tempKappa = status->giveTempKappa();
        tempDam = status->giveTempDam();
        if ( tempKappa > kappa ) {
            // plastic loading
            // Imports
            tempEffectiveStress = status->giveTempEffectiveStress();
            plasFlowDirec = status->givePlasFlowDirec();
            SSaTensor = status->giveSSaTensor();
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
            //answer.beTranspositionOf(tangentMatrix);
        } else {
            // elastic behavior with damage
            // Construction of the secant stiffness
            this->constructAnisoComplTensor(compliance);
            elasticity.beInverseOf(compliance);

            //this->constructAnisoStiffnessTensor(elasticity);
            answer = elasticity;
            answer.times(1.0 - tempDam);
        }

        double g = status->giveDensG();
        if ( g <= 0. ) {
            double factor = gammaL0 * pow(rho, rL) + gammaP0 *pow(rho, rP) * ( tDens - 1 ) * pow(g, tDens - 2);
            // printf("densification");
            tangentMatrix.resize(6, 6);
            tangentMatrix.zero();
            tangentMatrix.at(1, 1) = tangentMatrix.at(1, 2) = tangentMatrix.at(1, 3) = 1;
            tangentMatrix.at(2, 1) = tangentMatrix.at(2, 2) = tangentMatrix.at(2, 3) = 1;
            tangentMatrix.at(3, 1) = tangentMatrix.at(3, 2) = tangentMatrix.at(3, 3) = 1;
            tangentMatrix.times(factor);
            answer.add(tangentMatrix);
        }
    }

    status->setSmtrx(answer);
}


double
TrabBone3D :: evaluateCurrentYieldStress(const double kappa)
{
    if ( formulation == 0 ) {
        return hardFactor * ( 1. + plasHardFactor * ( 1.0 - exp(-kappa * expPlasHard) ) );
    }

#if 0
    else {
        //Hadi post-yield function
        double damage = this->computeDamageParam(kappa);
        if ( kappa < kappaMax ) {
            return ( yR + ( 1. - yR ) * ( 1 - pow( ( kappaMax - kappa ) / kappaMax, kappaSlope * kappaMax ) ) ) / ( 1. - damage );
        } else if ( kappa < ( kappaMin + kappaMax ) / 2. ) {
            return ( yR + ( 1. - yR ) * ( 1 - ( 1 - gMin ) / 2 * pow(2 * ( kappa - kappaMax ) / ( kappaMin - kappaMax ), N) ) ) / ( 1. - damage );
        } else if ( kappa < kappaMin ) {
            return ( yR + ( 1. - yR ) * ( gMin + ( 1 - gMin ) / 2. * pow(2 * ( kappaMin - kappa ) / ( kappaMin - kappaMax ), N) ) ) / ( 1. - damage );
        } else {
            return ( yR + ( 1. - yR ) * gMin ) / ( 1. - damage );
        }
    }
#endif
    else {
        double a = ( yR + ( pR - yR ) * ( 1. - exp(-expPlasHard * kappa) ) );
        return a;
    }
}

double
TrabBone3D :: evaluateCurrentPlasticModulus(const double kappa)
{
    if ( formulation == 0 ) {
        return hardFactor * ( plasHardFactor * expPlasHard * exp(-kappa * expPlasHard) );
    }

#if 0
    else {
        double damage = this->computeDamageParam(kappa);
        double gPrime;
        double damagePrime =  critDam * expDam * exp(-expDam * kappa);
        double g = evaluateCurrentYieldStress(kappa);

        if ( kappa < kappaMax ) {
            gPrime = kappaSlope * pow( ( kappaMax - kappa ) / kappaMax, kappaSlope * kappaMax - 1. );
        } else if ( kappa < ( ( kappaMin + kappaMax ) / 2. ) ) {
            gPrime   = ( ( gMin - 1. ) / ( kappaMin - kappaMax ) * N * pow(2. * ( kappa - kappaMax ) / ( kappaMin - kappaMax ), N - 1.) );
        } else if ( kappa < kappaMin ) {
            gPrime =  ( ( 1. - gMin ) / ( kappaMin - kappaMax ) * N * pow(2. * ( kappaMin - kappa ) / ( kappaMin - kappaMax ), N - 1.) );
        } else {
            gPrime = 0.;
        }

        return ( ( 1. - yR ) * gPrime / ( 1. - damage ) + damagePrime * ( yR + ( 1. - yR ) * g ) / ( 1. - damage ) / ( 1. - damage ) );
    }
#endif
    else {
        double a = ( ( pR - yR ) * expPlasHard * exp(-expPlasHard * kappa) );
        return a;
    }
}


double
TrabBone3D :: evaluateCurrentViscousStress(const double deltaKappa, TimeStep *tStep)
{
    //  double deltaT = 0.01;
    double deltaT =  tStep->giveTimeIncrement();
    double answer;
    //return answer;
    if ( deltaT == 0 ) {
        answer = 0;
    } else {
        answer = -viscosity * deltaKappa / deltaT;
    }

    return answer;
}

double
TrabBone3D :: evaluateCurrentViscousModulus(const double deltaKappa, TimeStep *tStep)
{
    double deltaT = tStep->giveTimeIncrement();
    double answer = -viscosity / deltaT;

    return answer;
}


void
TrabBone3D :: performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain, TimeStep *tStep)
{
    bool convergence;
    double tempKappa;
    FloatArray tempPlasDef, tempEffectiveStress, trialEffectiveStress;
    FloatMatrix elasticity, compliance;

    TrabBone3DStatus *status = static_cast< TrabBone3DStatus * >( this->giveStatus(gp) );

    // elastic compliance
    this->constructAnisoComplTensor(compliance);
    // elastic stiffness
    elasticity.beInverseOf(compliance);

    //  this->constructAnisoStiffnessTensor(elasticity);
    // initialize the plastic strain and cumulative plastic strain
    // by values after the previous step
    tempPlasDef = status->givePlasDef();
    tempKappa = status->giveKappa();
    // evaluate the trial stress
    trialEffectiveStress.beProductOf(elasticity, totalStrain - tempPlasDef);
    tempEffectiveStress = trialEffectiveStress;
    // apply the iterative procedure that solves the system of nonlinear equations
    // consisting of the yield condition and discretized flow rule
    // and evaluates:
    // tempKappa ... cumulative plastic strain at the end of the substep
    // tempEffectiveStress ... effective stress at the end of the substep
    // tempPlasDef ... plastic strain at the end of the substep
    convergence = projectOnYieldSurface(tempKappa, tempEffectiveStress, tempPlasDef, trialEffectiveStress, elasticity, compliance, status, tStep, gp, 0);
    if ( convergence ) {
        status->setTempPlasDef(tempPlasDef);
        status->setTempKappa(tempKappa);
        status->setTempEffectiveStress(tempEffectiveStress);
    } else {
        //printf("LineSearch \n");
        tempEffectiveStress = trialEffectiveStress;
        tempKappa = status->giveKappa();
        tempPlasDef = status->givePlasDef();
        convergence = projectOnYieldSurface(tempKappa, tempEffectiveStress, tempPlasDef, trialEffectiveStress, elasticity, compliance, status, tStep, gp, 1);
        if ( !convergence ) {
            //printf("No convergence %d",gp->giveNumber());
            //_error("No convergence of the stress return algorithm in TrabBone3D :: performPlasticityReturn\n");
        }

        status->setTempPlasDef(tempPlasDef);
        status->setTempKappa(tempKappa);
        status->setTempEffectiveStress(tempEffectiveStress);
    }
}


bool
TrabBone3D :: projectOnYieldSurface(double &tempKappa, FloatArray &tempEffectiveStress, FloatArray &tempPlasDef, const FloatArray &trialEffectiveStress, const FloatMatrix &elasticity, const FloatMatrix &compliance, TrabBone3DStatus *status, TimeStep *tStep, GaussPoint *gp, int lineSearchFlag)
{
    bool convergence;
    int flagLoop;
    double deltaKappa, incKappa,  beta, tempScalar, norm, FS, SFS;
    double plasCriterion, plasModulus, viscoModulus, toSolveScalar, errorF, errorR;
    FloatArray incTempEffectiveStress, errLoop;
    FloatArray toSolveTensor, plasFlowDirec, tempTensor2, tensorFF_S, F;
    FloatMatrix fabric, SSaTensor, tempTensor4, normAdjust, derivPlasFlowDirec;

    this->constructAnisoFabricTensor(fabric);
    this->constructAnisoFtensor(F);

    tempTensor2.beProductOf(fabric, trialEffectiveStress);
    FS = F.dotProduct(trialEffectiveStress);
    SFS = sqrt( trialEffectiveStress.dotProduct(tempTensor2) );
    plasCriterion = SFS + FS - this->evaluateCurrentYieldStress(tempKappa);

    if ( plasCriterion < rel_yield_tol ) {
        // trial stress in elastic domain
        convergence = true;
    } else {
        // return to the yield surface needed
        // Initial valuesr
        toSolveTensor.resize(6);
        toSolveTensor.zero();
        toSolveScalar = plasCriterion;
        errorF = plasCriterion;
        errorR = 0;
        SSaTensor = elasticity;
        this->constructPlasFlowDirec(plasFlowDirec, norm, fabric, F, tempEffectiveStress);
        tensorFF_S.beProductOf(fabric, tempEffectiveStress);
        deltaKappa = 0.;
        /***********************************************************************************************/
        double radialReturnFlag = lineSearchFlag;
        if ( radialReturnFlag == 2 ) {
            //printf("Radial return");
            double denom, dSYdA, dAlfa;
            double k = tempKappa;
            double f = plasCriterion;
            double alfa = 1;
            FloatArray helpArray, stress;
            stress.resize(6);
            stress.zero();
            double SFSr = sqrt( trialEffectiveStress.dotProduct(tempTensor2) );
            double FSr = F.dotProduct(trialEffectiveStress);
            tempTensor2.beProductOf(compliance, trialEffectiveStress);
            this->constructNormAdjustTensor(normAdjust);
            helpArray.beProductOf(normAdjust, tempTensor2);
            norm = sqrt( tempTensor2.dotProduct(helpArray) );
            while ( fabs(f) > 1.e-12 ) {
                dSYdA = norm * this->evaluateCurrentPlasticModulus(k);
                dSYdA *= -norm;
                denom = SFSr + FSr - dSYdA;
                dAlfa = -f / denom;
                alfa += dAlfa;
                stress = trialEffectiveStress;
                stress.times(alfa);
                k = k + ( 1 - alfa ) * norm;
                f =  evaluatePlasCriterion(fabric, F, stress, k, ( 1 - alfa ) * norm, tStep);
            }

            tempEffectiveStress  = stress;
            deltaKappa = k;
            toSolveScalar = 0;
            this->constructPlasFlowDirec(plasFlowDirec, norm, fabric, F, tempEffectiveStress);


            if ( tempEffectiveStress.giveSize() != trialEffectiveStress.giveSize() ) {
                printf( "tempS  %d \n", tempEffectiveStress.giveSize() );
                printf( "trial S %d \n", trialEffectiveStress.giveSize() );
            }

            tempTensor2.beProductOf( compliance, ( tempEffectiveStress - trialEffectiveStress ) );
            toSolveTensor = tempTensor2 + deltaKappa * plasFlowDirec;
            // Construction of the derivative of the plastic flow
            this->constructDerivativeOfPlasFlowDirec(derivPlasFlowDirec, fabric, F, tempEffectiveStress);
            // Construction of the gradient Nabla_S of R and SSa tensor
            tempTensor4 = derivPlasFlowDirec;
            tempTensor4.times(deltaKappa);
            tempTensor4.add(compliance);
            SSaTensor.beInverseOf(tempTensor4);
        }

        /***********************************************************************************************/
        // iteration loop - solution of a set of nonlinear equations
        flagLoop = 1;
        do {
            plasModulus = evaluateCurrentPlasticModulus(tempKappa + deltaKappa);
            viscoModulus = evaluateCurrentViscousModulus(deltaKappa, tStep);
            //*************************************
            //Evaluation of the Recursive Equations
            //*************************************
            tempTensor2.beProductOf(SSaTensor, plasFlowDirec);
            beta = plasFlowDirec.dotProduct(tempTensor2);
            beta += ( plasModulus - viscoModulus ) / norm; //+ viscoModulus;
            // Construction of the equation of Delta Kappa
            tempTensor2.beProductOf(SSaTensor, toSolveTensor);
            tempScalar = plasFlowDirec.dotProduct(tempTensor2);
            incKappa =  ( toSolveScalar / norm - tempScalar ) / beta;
            tempTensor2 = plasFlowDirec;
            tempTensor2.times(incKappa);
            tempTensor2 += toSolveTensor;
            incTempEffectiveStress.beProductOf(SSaTensor, tempTensor2);

            if ( lineSearchFlag == 1 ) {
                max_num_iter = 4000;
                double eta = 0.1;
                beta = 25.e-4;
                int jMax = 10;
                double alfa = 1;
                double M = ( errorR * errorR + errorF * errorF ) / 2;
                double dM = -2 * M;
                double newDeltaKappa;
                FloatArray tempStress;
                for ( int j = 0;; ++j ) {
                    tempStress = incTempEffectiveStress;
                    tempStress.times(-alfa);
                    tempStress.add(tempEffectiveStress);
                    newDeltaKappa = deltaKappa + alfa * incKappa;
                    //*************************************
                    // Evaluation of the f and R
                    //*************************************
                    // Construction of the derivative of the plastic flow
                    this->constructDerivativeOfPlasFlowDirec(derivPlasFlowDirec, fabric, F, tempStress);
                    // Construction of the gradient Nabla_S of R and SSa tensor
                    tempTensor4 = derivPlasFlowDirec;
                    tempTensor4.times(newDeltaKappa);
                    tempTensor4.add(compliance);
                    SSaTensor.beInverseOf(tempTensor4);
                    // Evaluation of R
                    this->constructPlasFlowDirec(plasFlowDirec, norm, fabric, F, tempStress);
                    tempTensor2.beProductOf( compliance, ( tempStress - trialEffectiveStress ) );
                    toSolveTensor = tempTensor2 + newDeltaKappa * plasFlowDirec;
                    // Evaluation of f
                    tempTensor2.beProductOf(fabric, tempStress);
                    SFS = sqrt( tempStress.dotProduct(tempTensor2) );
                    toSolveScalar = evaluatePlasCriterion(fabric, F, tempStress, tempKappa + newDeltaKappa, newDeltaKappa, tStep);
                    //*************************************
                    // Evaluation of the error
                    //*************************************
                    errLoop = toSolveTensor;
                    errorR = sqrt( errLoop.dotProduct(errLoop) );
                    errorF = fabs(toSolveScalar);
                    double newM = ( errorR * errorR + errorF * errorF ) / 2;
                    //check Goldstein's condition
                    //if(newM < M || j == jMax)
                    if ( newM < ( 1 - 2 * beta * alfa ) * M || j == jMax ) {
                        deltaKappa = newDeltaKappa;
                        tempEffectiveStress = tempStress;
                        break;
                    }

                    double alfa1 = eta * alfa;
                    double alfa2 = -alfa * alfa * dM / 2 / ( newM - M - alfa * dM );
                    alfa = alfa1;
                    if ( alfa1 < alfa2 ) {
                        alfa = alfa2;
                    }
                }
            } else {
                max_num_iter = 100;
                //////////////////////////////////////////////////////////
                deltaKappa += incKappa;
                tempEffectiveStress -= incTempEffectiveStress;
                //*************************************
                // Evaluation of the f and R
                //*************************************
                // Construction of the derivative of the plastic flow
                this->constructDerivativeOfPlasFlowDirec(derivPlasFlowDirec, fabric, F, tempEffectiveStress);
                // Construction of the gradient Nabla_S of R and SSa tensor
                tempTensor4 = derivPlasFlowDirec;
                tempTensor4.times(deltaKappa);
                tempTensor4.add(compliance);
                SSaTensor.beInverseOf(tempTensor4);
                // Evaluation of R
                this->constructPlasFlowDirec(plasFlowDirec, norm, fabric, F, tempEffectiveStress);
                tempTensor2.beProductOf( compliance, ( tempEffectiveStress - trialEffectiveStress ) );
                toSolveTensor = tempTensor2 + deltaKappa * plasFlowDirec;
                // Evaluation of f
                tempTensor2.beProductOf(fabric, tempEffectiveStress);
                SFS = sqrt( tempEffectiveStress.dotProduct(tempTensor2) );
                toSolveScalar = evaluatePlasCriterion(fabric, F, tempEffectiveStress, tempKappa + deltaKappa, deltaKappa, tStep);
                //*************************************
                // Evaluation of the error
                //*************************************
                errLoop = toSolveTensor;
                errorR = sqrt( errLoop.dotProduct(errLoop) );
                errorF = fabs(toSolveScalar);
            }

            if ( printflag ) {
                printf("   %d %g %g %g %g\n", flagLoop, tempEffectiveStress.at(1), tempEffectiveStress.at(3), incKappa, deltaKappa);
            }

            flagLoop++;
            convergence = ( fabs(errorF) < rel_yield_tol && errorR < strain_tol );
        } while ( flagLoop <= max_num_iter && !convergence );

        if ( convergence ) {
            plasModulus = evaluateCurrentPlasticModulus(tempKappa + deltaKappa);
            viscoModulus = evaluateCurrentViscousModulus(deltaKappa, tStep);
            tempTensor2.beProductOf(SSaTensor, plasFlowDirec);
            beta = plasFlowDirec.dotProduct(tempTensor2);
            beta += ( plasModulus - viscoModulus ) / norm;
            tempPlasDef += deltaKappa * plasFlowDirec;
            tempKappa += deltaKappa;
            status->setBeta(beta);
            status->setPlasFlowDirec(plasFlowDirec);
            status->setSSaTensor(SSaTensor);
        }
    }

    return convergence;
}

void
TrabBone3D :: constructPlasFlowDirec(FloatArray &answer, double  &norm, FloatMatrix &fabric, FloatArray &F, FloatArray &S)
{
    double SFS;
    FloatArray FFS, tempTensor2;
    FloatMatrix normAdjust;
    //////////////////////////////////////////////////
    FFS.beProductOf(fabric, S);
    //  FS = F.dotProduct(S);
    SFS = sqrt( S.dotProduct(FFS) );
    // scaling matrix
    this->constructNormAdjustTensor(normAdjust);
    //direction of Np
    answer = F;
    answer.add(1. / SFS, FFS);
    tempTensor2.beProductOf(normAdjust, answer);
    //norm Np
    norm = sqrt( answer.dotProduct(tempTensor2) );
    //plastic flow
    answer.times(1.0 / norm);
    //////////////////////////////////////////////////
}
void
TrabBone3D :: constructDerivativeOfPlasFlowDirec(FloatMatrix &answer, FloatMatrix &fabric, FloatArray &F, FloatArray &S)
{
    double SFS, norm, SGS, h;
    FloatArray FFS, FFSF, tempTensor2, tempTensor21, dNorm;
    FloatMatrix normAdjust, tempTensor4, dNp, FSFS;
    //////////////////////////////////////////////////
    FFS.beProductOf(fabric, S);
    SFS = sqrt( S.dotProduct(FFS) );

    SGS = pow(SFS, -3.);
    FloatArray gradientOfG, gradientOfH;
    FloatMatrix secondGradientOfG, helpMatrix;
    gradientOfG.zero();
    gradientOfG.add(FFS);
    gradientOfG.times(1. / SFS);
    gradientOfG.add(F);

    helpMatrix.zero();
    helpMatrix.add(fabric);
    helpMatrix.times(1. / SFS);
    secondGradientOfG.beDyadicProductOf(FFS, FFS);
    secondGradientOfG.times(-1. * SGS);
    secondGradientOfG.add(helpMatrix);

    h = gradientOfG.dotProduct(gradientOfG);
    h = sqrt(h);


    gradientOfH.beTProductOf(secondGradientOfG, gradientOfG);
    gradientOfH.times(1. / h);

    secondGradientOfG.times(h);
    FloatMatrix test;
    test.beDyadicProductOf(gradientOfG, gradientOfH);
    test.times(-1);
    test.add(secondGradientOfG);
    test.times(1. / h / h);
    //////////////////////////////////////////////////////////////////
    FFS.beProductOf(fabric, S);
    SFS = sqrt( S.dotProduct(FFS) );
    // scaling matrix
    this->constructNormAdjustTensor(normAdjust);
    //norm
    tempTensor2 = FFS;
    tempTensor2.times(1. / SFS);
    tempTensor2.add(F);
    tempTensor21.beProductOf(normAdjust, tempTensor2);
    //norm Np
    norm = sqrt( tempTensor2.dotProduct(tempTensor21) );
    ///////////////////////////////////////////////////////////////////
    FSFS.beDyadicProductOf(FFS, FFS);
    dNp.zero();
    dNp.add(FSFS);
    dNp.times(-1. / SFS / SFS);
    dNp.add(fabric);
    dNp.times(1. / SFS / norm);
    //////////////////////////////////////////////////////////////
    FFSF.zero();
    FFSF.add(FFS);
    FFSF.times(1. / SFS);
    FFSF.add(F);
    /////////////////////////////////////////////////////////////////
    tempTensor4.zero();
    tempTensor4.add(FSFS);
    tempTensor4.times(-1. / SFS / SFS);
    tempTensor4.add(fabric);
    tempTensor4.times(1. / SFS / norm);
    tempTensor2.beProductOf(normAdjust, FFSF);
    dNorm.beProductOf(tempTensor4, tempTensor2);
    ///////////////////////////////////////////////////////////////////
    tempTensor4.beDyadicProductOf(FFSF, dNorm);
    tempTensor4.times(-1. / norm / norm);
    ////////////////////////////////////////////////////////////////////
    answer.zero();
    answer.add(dNp);
    answer.add(tempTensor4);
}

double
TrabBone3D :: evaluatePlasCriterion(FloatMatrix &fabric, FloatArray &F, FloatArray &stress, double kappa, double deltaKappa, TimeStep *tStep)
{
    FloatArray FFS;
    double FS, SFS;
    FFS.beProductOf(fabric, stress);
    FS = F.dotProduct(stress);
    SFS =  sqrt( stress.dotProduct(FFS) );
    return SFS + FS - evaluateCurrentYieldStress(kappa) + this->evaluateCurrentViscousStress(deltaKappa, tStep);
}

double
TrabBone3D :: computeDamageParam(double tempKappa)
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
TrabBone3D :: computeDamageParamPrime(double tempKappa)
{
    double damagePrime;
    damagePrime = critDam * expDam * exp(-expDam * tempKappa);
    return damagePrime;
}
//
// END: FUNCTION FOR DAMAGE PARAMETER
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: FUNCTION FOR DAMAGE EVALUATION
//

double
TrabBone3D :: computeDamage(GaussPoint *gp,  TimeStep *tStep)
{
    double tempKappa;

    computeCumPlastStrain(tempKappa, gp, tStep);

    double tempDam = computeDamageParam(tempKappa);
    if ( tempDam < 0 ) {
        OOFEM_ERROR("negative damage");
    }

    return tempDam;
}


void TrabBone3D :: computeCumPlastStrain(double &tempKappa, GaussPoint *gp, TimeStep *tStep)
{
    TrabBone3DStatus *status = static_cast< TrabBone3DStatus * >( this->giveStatus(gp) );
    tempKappa = status->giveTempKappa();
}



void TrabBone3D :: computeDensificationStress(FloatArray &answer, GaussPoint *gp, const FloatArray &totalStrain, TimeStep *tStep)
{
    TrabBone3DStatus *status = static_cast< TrabBone3DStatus * >( this->giveStatus(gp) );
    answer.resize(6);
    answer.zero();
    double traceLnU =  ( totalStrain.at(1) + totalStrain.at(2) + totalStrain.at(3) );
    double g = traceLnU - densCrit;
    status->setDensG(g);
    if ( g <= 0 ) {
        answer.at(1) = answer.at(2) = answer.at(3) = 1;
        double factor = gammaL0 * pow(rho, rL) * g + gammaP0 *pow(rho, rP) * pow(g, tDens - 1);
        answer.times(factor);
    }
}


void
TrabBone3D :: giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp,
                                      const FloatArray &totalStrain, TimeStep *tStep)
{
    double tempDam;
    FloatArray effStress, densStress;
    TrabBone3DStatus *status = static_cast< TrabBone3DStatus * >( this->giveStatus(gp) );

    this->initTempStatus(gp);
    // compute effective stress using the plasticity model
    performPlasticityReturn(gp, totalStrain, tStep);
    effStress =  status->giveTempEffectiveStress();

    // evaluate damage variable
    tempDam = computeDamage(gp, tStep);

    // transform effective stress into nominal stress
    answer = ( 1 - tempDam ) * effStress;

    computePlasStrainEnerDensity(gp, totalStrain, answer);

    // add densification stress, if the densificator is activated
    if ( densCrit != 0 ) {
        computeDensificationStress(densStress, gp, totalStrain, tStep);
        answer.add(densStress);
    }

    // store final damage, strain and stress in status
    status->setTempDam(tempDam);
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);
}


void
TrabBone3D :: constructAnisoComplTensor(FloatMatrix &answer)
{
    FloatMatrix D(6, 6);
    FloatMatrix DT(6, 6);
    FloatMatrix T;

    this->constructFabricTransformationMatrix(T);

    double m3 = 3. - m1 - m2;
    double eps0k = eps0 * pow(rho, expk);
    double mu0k = mu0 * pow(rho, expk);
    double m1l = pow(m1, expl);
    double m2l = pow(m2, expl);
    double m3l = pow(m3, expl);

    answer.resize(6, 6);

    D.at(1, 1) = 1 / ( eps0k * m1l * m1l );
    D.at(2, 2) = 1 / ( eps0k * m2l * m2l );
    D.at(3, 3) = 1 / ( eps0k * m3l * m3l );
    D.at(1, 2) = -nu0 / ( eps0k * m1l * m2l );
    D.at(2, 1) = D.at(1, 2);
    D.at(1, 3) = -nu0 / ( eps0k * m1l * m3l );
    D.at(3, 1) = D.at(1, 3);
    D.at(2, 3) = -nu0 / ( eps0k * m2l * m3l );
    D.at(3, 2) = D.at(2, 3);
    D.at(4, 4) = 1 / ( mu0k * m2l * m3l );
    D.at(5, 5) = 1 / ( mu0k * m3l * m1l );
    D.at(6, 6) = 1 / ( mu0k * m1l * m2l );

    FloatMatrix stiffness, t;
    this->constructStiffnessTransformationMatrix(t);
    stiffness.beInverseOf(D);
    FloatMatrix ST, TST;
    ST.beProductOf(stiffness, t);
    TST.beTProductOf(t, ST);

    DT.beProductOf(D, T);
    answer.beTProductOf(T, DT);
}


void
TrabBone3D :: constructAnisoStiffnessTensor(FloatMatrix &answer)
{
    double m3 = 3. - m1 - m2;
    double eps0k = eps0 * pow(rho, expk);
    double mu0k = mu0 * pow(rho, expk);
    double m1l = pow(m1, expl);
    double m2l = pow(m2, expl);
    double m3l = pow(m3, expl);


    double eksi, n13, n23, n12, n31, n32, n21;

    double E1 = eps0k * m1l * m1l;
    double E2 = eps0k * m2l * m2l;
    double E3 = eps0k * m3l * m3l;

    double G23 = mu0k * m2l * m3l;
    double G13 = mu0k * m1l * m3l;
    double G12 = mu0k * m1l * m2l;

    n23 = nu0 * m2l / m3l;
    n12 = nu0 * m1l / m2l;
    n31 = nu0 * m3l / m1l;
    n32 = nu0 * m3l / m2l;
    n21 = nu0 * m2l / m1l;
    n13 = nu0 * m1l / m3l;

    eksi = 1. - ( n12 * n21 + n23 * n32 + n31 * n13 ) - ( n12 * n23 * n31 + n21 * n32 * n13 );

    //constitutiveMatrix = new FloatMatrix(6,6) ;
    FloatMatrix D(6, 6);
    D.zero();
    // switched letters from original oofem -> now produces same material stiffness matrix as Abaqus method
    D.at(1, 1) =  E1 * ( 1. - n23 * n32 ) / eksi;
    D.at(1, 2) =  E2 * ( n12 + n13 * n32 ) / eksi;
    D.at(1, 3) =  E3 * ( n13 + n23 * n12 ) / eksi;
    D.at(2, 2) =  E2 * ( 1. - n13 * n31 ) / eksi;
    D.at(2, 3) =  E3 * ( n23 + n21 * n13 ) / eksi;
    D.at(3, 3) =  E3 * ( 1. - n21 * n12 ) / eksi;

    // define the lower triangle
    for ( int i = 1; i < 4; i++ ) {
        for ( int j = 1; j < i; j++ ) {
            D.at(i, j) = D.at(j, i);
        }
    }

    D.at(4, 4) =  G23;
    D.at(5, 5) =  G13;
    D.at(6, 6) =  G12;

    FloatMatrix t;
    this->constructStiffnessTransformationMatrix(t);
    FloatMatrix DT;
    DT.beProductOf(D, t);
    answer.beTProductOf(t, DT);
}



void
TrabBone3D :: constructAnisoFabricTensor(FloatMatrix &answer)
{
    FloatMatrix T;
    double S0 = ( sig0Pos + sig0Neg ) / ( 2. * sig0Pos * sig0Neg );
    double rhoP = pow(rho, 2. * expp);
    double m3 = 3. - m1 - m2;
    double m1q = pow(m1, 2. * expq);
    double m2q = pow(m2, 2. * expq);
    double m3q = pow(m3, 2. * expq);


    this->constructFabricTransformationMatrix(T);

    FloatMatrix F(6, 6);
    FloatMatrix FT;

    F.at(1, 1) = S0 * S0 / ( rhoP * m1q * m1q );
    F.at(2, 2) = S0 * S0 / ( rhoP * m2q * m2q );
    F.at(3, 3) = S0 * S0 / ( rhoP * m3q * m3q );
    F.at(1, 2) = -chi0 * S0 * S0 / ( rhoP * m1q * m2q );
    F.at(2, 1) = F.at(1, 2);
    F.at(1, 3) = -chi0 * S0 * S0 / ( rhoP * m1q * m3q );
    F.at(3, 1) = F.at(1, 3);
    F.at(2, 3) = -chi0 * S0 * S0 / ( rhoP * m2q * m3q );
    F.at(3, 2) = F.at(2, 3);
    F.at(4, 4) = 1. / ( tau0 * tau0 * rhoP * m2q * m3q );
    F.at(5, 5) = 1. / ( tau0 * tau0 * rhoP * m1q * m3q );
    F.at(6, 6) = 1. / ( tau0 * tau0 * rhoP * m1q * m2q );


    FT.beProductOf(F, T);
    answer.beTProductOf(T, FT);
}


void
TrabBone3D :: constructAnisoFtensor(FloatArray &answer)
{
    FloatMatrix T;

    double rhoP = pow(rho, expp);
    double m3 = 3. - m1 - m2;
    double m1q = pow(m1, 2. * expq);
    double m2q = pow(m2, 2. * expq);
    double m3q = pow(m3, 2. * expq);


    this->constructFabricTransformationMatrix(T);
    FloatArray F(6);
    F.zero();

    F.at(1) = -( sig0Pos - sig0Neg ) / ( 2. * sig0Pos * sig0Neg * rhoP * m1q );
    F.at(2) = -( sig0Pos - sig0Neg ) / ( 2. * sig0Pos * sig0Neg * rhoP * m2q );
    F.at(3) = -( sig0Pos - sig0Neg ) / ( 2. * sig0Pos * sig0Neg * rhoP * m3q );

    answer.beTProductOf(T, F);
}


void
TrabBone3D :: constructStiffnessTransformationMatrix(FloatMatrix &answer)
{
    answer.resize(6, 6);

    answer.at(1, 1) = x1 * x1;
    answer.at(1, 2) = x2 * x2;
    answer.at(1, 3) = x3 * x3;
    answer.at(1, 4) = x2 * x3;
    answer.at(1, 5) = x1 * x3;
    answer.at(1, 6) = x1 * x2;
    //second row of pull back transformation matrix
    answer.at(2, 1) = y1 * y1;
    answer.at(2, 2) = y2 * y2;
    answer.at(2, 3) = y3 * y3;
    answer.at(2, 4) = y2 * y3;
    answer.at(2, 5) = y1 * y3;
    answer.at(2, 6) = y1 * y2;
    //third row of pull back transformation matrix
    answer.at(3, 1) = z1 * z1;
    answer.at(3, 2) = z2 * z2;
    answer.at(3, 3) = z3 * z3;
    answer.at(3, 4) = z2 * z3;
    answer.at(3, 5) = z1 * z3;
    answer.at(3, 6) = z1 * z2;
    //fourth row of pull back transformation matrix
    answer.at(4, 1) = 2. * y1 * z1;
    answer.at(4, 2) = 2. * y2 * z2;
    answer.at(4, 3) = 2. * y3 * z3;
    answer.at(4, 4) = ( y2 * z3 + y3 * z2 );
    answer.at(4, 5) = ( y1 * z3 + y3 * z1 );
    answer.at(4, 6) = ( y1 * z2 + y2 * z1 );
    //fifth row of pull back transformation matrix
    answer.at(5, 1) = 2. * x1 * z1;
    answer.at(5, 2) = 2. * x2 * z2;
    answer.at(5, 3) = 2. * x3 * z3;
    answer.at(5, 4) = ( x2 * z3 + x3 * z2 );
    answer.at(5, 5) = ( x1 * z3 + x3 * z1 );
    answer.at(5, 6) = ( x1 * z2 + x2 * z1 );
    //sixth row of pull back transformation matrix
    answer.at(6, 1) = 2. * x1 * y1;
    answer.at(6, 2) = 2. * x2 * y2;
    answer.at(6, 3) = 2. * x3 * y3;
    answer.at(6, 4) = ( x2 * y3 + x3 * y2 );
    answer.at(6, 5) = ( x1 * y3 + x3 * y1 );
    answer.at(6, 6) = ( x1 * y2 + x2 * y1 );
}


void
TrabBone3D :: constructNormAdjustTensor(FloatMatrix &answer)
{
    answer.resize(6, 6);
    answer.zero();

    for ( int i = 1; i <= 3; i++ ) {
        answer.at(i, i) = 1.;
    }

    for ( int i = 4; i <= 6; i++ ) {
        answer.at(i, i) = 0.5;
    }
}



void
TrabBone3D :: constructFabricTransformationMatrix(FloatMatrix &answer)
{
    answer.resize(6, 6);

    answer.at(1, 1) = x1 * x1;
    answer.at(1, 2) = x2 * x2;
    answer.at(1, 3) = x3 * x3;
    answer.at(1, 4) = 2 * x2 * x3; //2
    answer.at(1, 5) = 2 * x1 * x3; //2
    answer.at(1, 6) = 2 * x1 * x2; //2
    //second row of pull back transformation matrix
    answer.at(2, 1) = y1 * y1;
    answer.at(2, 2) = y2 * y2;
    answer.at(2, 3) = y3 * y3;
    answer.at(2, 4) = 2 * y2 * y3; //2
    answer.at(2, 5) = 2 * y1 * y3; //2
    answer.at(2, 6) = 2 * y1 * y2; //2
    //third row of pull back transformation matrix
    answer.at(3, 1) = z1 * z1;
    answer.at(3, 2) = z2 * z2;
    answer.at(3, 3) = z3 * z3;
    answer.at(3, 4) = 2 * z2 * z3; //2
    answer.at(3, 5) = 2 * z1 * z3; //2
    answer.at(3, 6) = 2 * z1 * z2; //2
    //fourth row of pull back transformation matrix
    answer.at(4, 1) = y1 * z1;
    answer.at(4, 2) = y2 * z2;
    answer.at(4, 3) = y3 * z3;
    answer.at(4, 4) = ( y2 * z3 + y3 * z2 );
    answer.at(4, 5) = ( y1 * z3 + y3 * z1 );
    answer.at(4, 6) = ( y1 * z2 + y2 * z1 );
    //fifth row of pull back transformation matrix
    answer.at(5, 1) = x1 * z1;
    answer.at(5, 2) = x2 * z2;
    answer.at(5, 3) = x3 * z3;
    answer.at(5, 4) = ( x2 * z3 + x3 * z2 );
    answer.at(5, 5) = ( x1 * z3 + x3 * z1 );
    answer.at(5, 6) = ( x1 * z2 + x2 * z1 );
    //sixth row of pull back transformation matrix
    answer.at(6, 1) = x1 * y1;
    answer.at(6, 2) = x2 * y2;
    answer.at(6, 3) = x3 * y3;
    answer.at(6, 4) = ( x2 * y3 + x3 * y2 );
    answer.at(6, 5) = ( x1 * y3 + x3 * y1 );
    answer.at(6, 6) = ( x1 * y2 + x2 * y1 );
}

IRResultType
TrabBone3D :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    // Mandatory parameters
    IR_GIVE_FIELD(ir, eps0, _IFT_TrabBone3D_eps0);
    IR_GIVE_FIELD(ir, nu0, _IFT_TrabBone3D_nu0);
    IR_GIVE_FIELD(ir, mu0, _IFT_TrabBone3D_mu0);
    IR_GIVE_FIELD(ir, expk, _IFT_TrabBone3D_expk);
    IR_GIVE_FIELD(ir, expl, _IFT_TrabBone3D_expl);

    eps0 = eps0 * 0.1069;
    mu0 = mu0 * 0.1069;


    IR_GIVE_FIELD(ir, m1, _IFT_TrabBone3D_m1);
    IR_GIVE_FIELD(ir, m2, _IFT_TrabBone3D_m2);
    IR_GIVE_FIELD(ir, rho, _IFT_TrabBone3D_rho);

    IR_GIVE_FIELD(ir, sig0Pos, _IFT_TrabBone3D_sig0Pos);
    IR_GIVE_FIELD(ir, sig0Neg, _IFT_TrabBone3D_sig0Neg);
    IR_GIVE_FIELD(ir, chi0Pos, _IFT_TrabBone3D_chi0Pos);
    IR_GIVE_FIELD(ir, tau0, _IFT_TrabBone3D_tau0);
    sig0Pos = sig0Pos * 0.15;
    sig0Neg = sig0Neg * 0.15;
    tau0 = tau0 * 0.15;


    IR_GIVE_FIELD(ir, expq, _IFT_TrabBone3D_expq);
    IR_GIVE_FIELD(ir, expp, _IFT_TrabBone3D_expp);



    IR_GIVE_FIELD(ir, plasHardFactor, _IFT_TrabBone3D_plasHardFactor);
    IR_GIVE_FIELD(ir, expPlasHard, _IFT_TrabBone3D_expPlasHard);

    IR_GIVE_FIELD(ir, expDam, _IFT_TrabBone3D_expDam);
    IR_GIVE_FIELD(ir, critDam, _IFT_TrabBone3D_critDam);

    // evaluation of dependent parameter
    chi0Neg = ( sig0Neg * sig0Neg ) / ( sig0Pos * sig0Pos ) * ( chi0Pos + 1 ) - 1;
    chi0 = chi0Pos;





    //local coordinate system
    //x'
    x1 = 1;
    x2 = 0;
    x3 = 0;
    //y'
    y1 = 0;
    y2 = 1;
    y3 = 0;




    ///@todo Why not use vectors?
    IR_GIVE_OPTIONAL_FIELD(ir, x1, _IFT_TrabBone3D_x1);
    IR_GIVE_OPTIONAL_FIELD(ir, x2, _IFT_TrabBone3D_x2);
    IR_GIVE_OPTIONAL_FIELD(ir, x3, _IFT_TrabBone3D_x3);

    IR_GIVE_OPTIONAL_FIELD(ir, y1, _IFT_TrabBone3D_y1);
    IR_GIVE_OPTIONAL_FIELD(ir, y2, _IFT_TrabBone3D_y2);
    IR_GIVE_OPTIONAL_FIELD(ir, y3, _IFT_TrabBone3D_y3);

    double normX = sqrt(x1 * x1 + x2 * x2 + x3 * x3);
    x1 = x1 / normX;
    x2 = x2 / normX;
    x3 = x3 / normX;

    //y'
    double normY = sqrt(y1 * y1 + y2 * y2 + y3 * y3);
    y1 = y1 / normY;
    y2 = y2 / normY;
    y3 = y3 / normY;
    /////////////////////////////////////////////
    /*
     * x1 = 1;
     * x2 = 0;
     * x3 = 0;
     * y1 = 0;
     * y2 = 1;
     * y3 = 0;
     */
    //////////////////////////////////////////////
    //z'
    z1 = x2 * y3 - x3 * y2;
    z2 = x3 * y1 - x1 * y3;
    z3 = x1 * y2 - x2 * y1;
    //viscosity parameter
    viscosity = 0.05;
    //    viscosity = 0.01;
    //    viscosity = 0.015;
    IR_GIVE_OPTIONAL_FIELD(ir, viscosity, _IFT_TrabBone3D_viscosity);
    // Hadi post-yield function parameters
    yR = 0.7;
    kappaMax = 0.01;
    kappaMin = 0.1;
    kappaSlope = 300;
    N = 2.;
    gMin = -2.;
    formulation  = 10;



    IR_GIVE_OPTIONAL_FIELD(ir, yR, _IFT_TrabBone3D_yR);
    IR_GIVE_OPTIONAL_FIELD(ir, kappaMax, _IFT_TrabBone3D_kappaMax);
    IR_GIVE_OPTIONAL_FIELD(ir, kappaMin, _IFT_TrabBone3D_kappaMin);
    IR_GIVE_OPTIONAL_FIELD(ir, kappaSlope, _IFT_TrabBone3D_kappaSlope);
    IR_GIVE_OPTIONAL_FIELD(ir, N, _IFT_TrabBone3D_N);
    IR_GIVE_OPTIONAL_FIELD(ir, gMin, _IFT_TrabBone3D_gMin);
    IR_GIVE_OPTIONAL_FIELD(ir, formulation, _IFT_TrabBone3D_formulation);


    // optional densificator parameters
    gammaL0 = 550;
    // gammaL0 = 1100;
    // gammaP0 = 1300;
    gammaP0 = 39.;
    tDens = 6.;
    // densCrit = -0.55;
    densCrit = -0.3;
    rL = 2.928;
    rP = 2.77;


    IR_GIVE_OPTIONAL_FIELD(ir, gammaL0, _IFT_TrabBone3D_gammaL);
    IR_GIVE_OPTIONAL_FIELD(ir, gammaP0, _IFT_TrabBone3D_gammaP);
    IR_GIVE_OPTIONAL_FIELD(ir, tDens, _IFT_TrabBone3D_tDens);
    IR_GIVE_OPTIONAL_FIELD(ir, densCrit, _IFT_TrabBone3D_densCrit);





    // optional control parameters for printing and convergence
    printflag = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, printflag, _IFT_TrabBone3D_printflag);
    max_num_iter = 100;
    IR_GIVE_OPTIONAL_FIELD(ir, max_num_iter, _IFT_TrabBone3D_max_num_iter);
    max_num_substeps = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, max_num_substeps, _IFT_TrabBone3D_max_num_substeps);
    rel_yield_tol = 1.e-9;
    IR_GIVE_OPTIONAL_FIELD(ir, rel_yield_tol, _IFT_TrabBone3D_rel_yield_tol);
    strain_tol = 1.e-9;
    IR_GIVE_OPTIONAL_FIELD(ir, strain_tol, _IFT_TrabBone3D_strain_tol);
    abaqus = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, abaqus, _IFT_TrabBone3D_abaqus);


    ////////////////////////////////
    double oM, oPM, s;
    s = expPlasHard;

    /*  P1 = exp(-expDam* kappaMax);
     * P2 = exp(-2*expDam*kappaMax);
     * P3 = exp(expPlasHard*kappaMax);
     * P4 = P3/(1-critDam*(1-P1));
     * citatel = (critDam*expDam*P1 + expPlasHard -expPlasHard*critDam + expPlasHard*critDam *P1);
     * jmenovatel = expPlasHard*(-1+2*critDam-2*critDam*P1-pow(critDam,2) + 2*pow(critDam,2)*P1-pow(critDam,2)*P2);
     */

    oM = this->computeDamageParam(kappaMax);
    oPM = this->computeDamageParamPrime(kappaMax);

    pR = ( s + oPM / ( 1 - oM ) ) / ( 1 - oM ) / s;
    //pR = pR;
    //hardFactor = 0.75;
    //    pR = -citatel/jmenovatel;
    yR = exp(expPlasHard * kappaMax) / ( 1 - oM ) - pR * ( exp(expPlasHard * kappaMax) - 1 );
    pR = 0.98 * pR;

    return StructuralMaterial :: initializeFrom(ir);
}



int
TrabBone3D :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    TrabBone3DStatus *status = static_cast< TrabBone3DStatus * >( this->giveStatus(gp) );
    if ( type == IST_DamageScalar ) {
        answer.resize(1);
        answer.at(1) = status->giveTempDam();
        return 1;
    } else if ( type == IST_PlasticStrainTensor ) {
        answer = status->giveTempPlasDef();
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
        return StructuralMaterial :: giveIPValue(answer, gp, type, tStep);
    }
}


double TrabBone3D :: predictRelativeComputationalCost(GaussPoint *gp)
{
    TrabBone3DStatus *status = static_cast< TrabBone3DStatus * >( this->giveStatus(gp) );

    if ( status->giveTempDam() > 0.0 ) {
        return 15.0;
    } else {
        return 1.0;
    }
}


double TrabBone3D :: predictRelativeRedistributionCost(GaussPoint *gp)
{
    return 1.0;
}



TrabBone3DStatus :: TrabBone3DStatus(int n, Domain *d, GaussPoint *g) : StructuralMaterialStatus(n, d, g)
{
    kappa = 0.0;
    tempKappa = 0.0;
    dam = 0.0;
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
    densG = 200;
}



TrabBone3DStatus :: ~TrabBone3DStatus()
{ }


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

const FloatArray &
TrabBone3DStatus :: giveTempEffectiveStress() const
{
    return tempEffectiveStress;
}

const FloatArray &
TrabBone3DStatus :: givePlasFlowDirec() const
{
    return plasFlowDirec;
}

const FloatArray &
TrabBone3DStatus :: givePlasDef() const
{
    return plasDef;
}

const FloatArray &
TrabBone3DStatus :: giveTempPlasDef() const
{
    return tempPlasDef;
}

const FloatMatrix &
TrabBone3DStatus :: giveSSaTensor() const
{
    return SSaTensor;
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
    densG = 200;
    this->tempKappa = this->kappa;
    this->tempDam = this->dam;
    this->tempTSED = this->tsed;
    this->tempPlasDef = this->plasDef;
    nss = 1;
}


void
TrabBone3DStatus :: updateYourself(TimeStep *tStep)
{
    StructuralMaterialStatus :: updateYourself(tStep);
    this->kappa = this->tempKappa;
    this->dam = this->tempDam;
    this->tsed = this->tempTSED;
    this->plasDef = this->tempPlasDef;
    nss = 1;
}


contextIOResultType
TrabBone3DStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = plasDef.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }


    if ( !stream.write(dam) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(beta) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( iores = effectiveStress.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = plasFlowDirec.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    /*
     * if ( ( iores = smtrx.storeYourself(stream) ) != CIO_OK ) {
     * THROW_CIOERR(iores);
     * }
     * if ( ( iores = tangentMatrix.storeYourself(stream) ) != CIO_OK ) {
     * THROW_CIOERR(iores);
     * }
     * if ( ( iores = SSaTensor.storeYourself(stream) ) != CIO_OK ) {
     * THROW_CIOERR(iores);
     * }
     *
     * if ( ( iores = tempStrain.storeYourself(stream) ) != CIO_OK ) {
     * THROW_CIOERR(iores);
     * }
     */
    return CIO_OK;
}


contextIOResultType
TrabBone3DStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
    if ( ( iores = plasDef.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }


    if ( !stream.read(dam) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(beta) ) {
        THROW_CIOERR(CIO_IOERR);
    }


    if ( ( iores = effectiveStress.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = plasFlowDirec.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    /*
     * if ( ( iores = smtrx.restoreYourself(stream) ) != CIO_OK ) {
     * THROW_CIOERR(iores);
     * }
     * if ( ( iores = tangentMatrix.restoreYourself(stream) ) != CIO_OK ) {
     * THROW_CIOERR(iores);
     * }
     * if ( ( iores = SSaTensor.restoreYourself(stream) ) != CIO_OK ) {
     * THROW_CIOERR(iores);
     * }
     * if ( ( iores = tempStrain.restoreYourself(stream) ) != CIO_OK ) {
     * THROW_CIOERR(iores);
     * }
     */



    return CIO_OK;
}


MaterialStatus *TrabBone3D :: CreateStatus(GaussPoint *gp) const
{
    return new TrabBone3DStatus(1, StructuralMaterial :: giveDomain(), gp);
}
} //end namespace oofem
