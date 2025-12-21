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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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


void TrabBone3D :: computePlasStrainEnerDensity(GaussPoint *gp, const FloatArrayF<6> &strain, const FloatArrayF<6> &stress) const
{
    auto status = static_cast< TrabBone3DStatus * >( this->giveStatus(gp) );

    auto tsed = status->giveTSED();
    auto &oldTotalDef = status->giveStrainVector();
    auto &tempPlasDef = status->giveTempPlasDef();
    auto &oldStress = status->giveTempStressVector();

    auto elStrain = strain - tempPlasDef;
    auto deltaStrain = strain - oldTotalDef;

    auto tmpStress = stress + oldStress;

    double tempESED = 0.5 * dot(elStrain, stress);
    double tempTSED = tsed + 0.5 * dot(deltaStrain, tmpStress);
    double tempPSED = tempTSED - tempESED;

    status->setTempTSED(tempTSED);
    status->setTempPSED(tempPSED);
}


FloatMatrixF<6,6>
TrabBone3D :: give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< TrabBone3DStatus * >( this->giveStatus(gp) );

    if ( mode == ElasticStiffness ) {
        auto compliance = this->constructAnisoComplTensor();
        auto elasticity = inv(compliance);
        return elasticity;
    } else if ( mode == SecantStiffness ) {
        if ( printflag ) {
            printf("secant\n");
        }

        auto compliance = this->constructAnisoComplTensor();
        auto elasticity = inv(compliance);
        auto tempDam = status->giveTempDam();
        return elasticity * (1.0 - tempDam);
    } else /*if ( mode == TangentStiffness )*/ {
        double kappa = status->giveKappa();
        double tempKappa = status->giveTempKappa();
        double tempDam = status->giveTempDam();

        FloatMatrixF<6,6> answer;
        if ( tempKappa > kappa ) {
            // plastic loading
            // Imports
            auto &tempEffectiveStress = status->giveTempEffectiveStress();
            auto &plasFlowDirec = status->givePlasFlowDirec();
            auto &SSaTensor = status->giveSSaTensor();
            auto beta = status->giveBeta();
            // Construction of the dyadic product tensor
            auto prodTensor = Tdot(SSaTensor, plasFlowDirec);
            // Construction of the tangent stiffness third term
            auto thirdTerm = dyad(tempEffectiveStress, prodTensor) * ((-expDam * critDam * exp(-expDam * tempKappa)) / beta);
            // Construction of the tangent stiffness second term
            auto tempTensor2 = dot(SSaTensor, plasFlowDirec);
            auto secondTerm = dyad(tempTensor2, prodTensor) * (( 1.0 - tempDam ) / beta);
            // Construction of the tangent stiffness
            auto tangentMatrix = SSaTensor * (1.0 - tempDam) + secondTerm + thirdTerm;
            answer = tangentMatrix;
        } else {
            // elastic behavior with damage
            // Construction of the secant stiffness
            auto compliance = this->constructAnisoComplTensor();
            auto elasticity = inv(compliance);
            answer = elasticity * (1.0 - tempDam);
        }

        double g = status->giveDensG();
        if ( g <= 0. ) {
            double factor = gammaL0 * pow(rho, rL) + gammaP0 *pow(rho, rP) * ( tDens - 1 ) * pow(g, tDens - 2);
            // printf("densification");
            answer += I6_I6 * factor;
        }
        return answer;
    }
}


double
TrabBone3D :: evaluateCurrentYieldStress(double kappa) const
{
    return  ( 1. + plasHardFactor * ( 1.0 - exp(-kappa * expPlasHard) ) );
}

double
TrabBone3D :: evaluateCurrentPlasticModulus(double kappa) const
{
    return ( plasHardFactor * expPlasHard * exp(-kappa * expPlasHard) );
}


double
TrabBone3D :: evaluateCurrentViscousStress(const double deltaKappa, TimeStep *tStep) const
{
    double deltaT = tStep->giveTimeIncrement();
    double answer;
    if ( deltaT == 0 ) {
        answer = 0;
    } else {
        answer = -viscosity * deltaKappa / deltaT;
    }

    return answer;
}

double
TrabBone3D :: evaluateCurrentViscousModulus(double deltaKappa, TimeStep *tStep) const
{
    double deltaT = tStep->giveTimeIncrement();
    double answer = -viscosity / deltaT;

    return answer;
}


void
TrabBone3D :: performPlasticityReturn(GaussPoint *gp, const FloatArrayF<6> &strain, TimeStep *tStep) const
{
    auto status = static_cast< TrabBone3DStatus * >( this->giveStatus(gp) );

    // elastic compliance
    auto compliance = this->constructAnisoComplTensor();
    // elastic stiffness
    auto elasticity = inv(compliance);

    //  this->constructAnisoStiffnessTensor(elasticity);
    // initialize the plastic strain and cumulative plastic strain
    // by values after the previous step
    auto tempPlasDef = status->givePlasDef();
    double tempKappa = status->giveKappa();
    // evaluate the trial stress
    auto trialEffectiveStress = dot(elasticity, strain - tempPlasDef);
    auto tempEffectiveStress = trialEffectiveStress;
    // apply the iterative procedure that solves the system of nonlinear equations
    // consisting of the yield condition and discretized flow rule
    // and evaluates:
    // tempKappa ... cumulative plastic strain at the end of the substep
    // tempEffectiveStress ... effective stress at the end of the substep
    // tempPlasDef ... plastic strain at the end of the substep
    bool convergence = projectOnYieldSurface(tempKappa, tempEffectiveStress, tempPlasDef, trialEffectiveStress, elasticity, compliance, status, tStep, gp, 0);
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
            OOFEM_ERROR("No convergence of the stress return algorithm in TrabBone3D :: performPlasticityReturn\n");
        }

        status->setTempPlasDef(tempPlasDef);
        status->setTempKappa(tempKappa);
        status->setTempEffectiveStress(tempEffectiveStress);
    }
}


bool
TrabBone3D :: projectOnYieldSurface(double &tempKappa, FloatArrayF<6> &tempEffectiveStress, FloatArrayF<6> &tempPlasDef, const FloatArrayF<6> &trialEffectiveStress, const FloatMatrixF<6,6> &elasticity, const FloatMatrixF<6,6> &compliance, TrabBone3DStatus *status, TimeStep *tStep, GaussPoint *gp, int lineSearchFlag) const
{
    bool convergence;

    auto fabric = this->constructAnisoFabricTensor();
    auto F = this->constructAnisoFtensor();

    auto tempTensor2 = dot(fabric, trialEffectiveStress);
    auto FS = dot(F, trialEffectiveStress);
    auto SFS = sqrt( dot(trialEffectiveStress, tempTensor2) );
    auto plasCriterion = SFS + FS - this->evaluateCurrentYieldStress(tempKappa);

    int current_max_num_iter = this->max_num_iter;
    if ( plasCriterion < rel_yield_tol ) {
        // trial stress in elastic domain
        convergence = true;
    } else {
        // return to the yield surface needed
        // Initial valuesr
        FloatArrayF<6> toSolveTensor;
        double toSolveScalar = plasCriterion;
        double errorF = plasCriterion;
        double errorR = 0;
        double deltaKappa = 0.;
        auto SSaTensor = elasticity;
        //auto [plastFloxDirec, norm] = this->constructPlasFlowDirec(fabric, F, tempEffectiveStress);
        auto tmp = this->constructPlasFlowDirec(fabric, F, tempEffectiveStress);
        auto plasFlowDirec = tmp.first;
        double norm = tmp.second;
        //auto tensorFF_S = dot(fabric, tempEffectiveStress); // FIXME: unusued in old code
        /***********************************************************************************************/
        double radialReturnFlag = lineSearchFlag;
        if ( radialReturnFlag == 2 ) {
            //printf("Radial return");
            double k = tempKappa;
            double f = plasCriterion;
            double SFSr = sqrt( dot(trialEffectiveStress, tempTensor2) );
            double FSr = dot(F, trialEffectiveStress);
            auto normAdjust = this->constructNormAdjustTensor();
            auto tempTensor2 = dot(compliance, trialEffectiveStress);
            auto helpArray = dot(normAdjust, tempTensor2);
            norm = sqrt( dot(tempTensor2, helpArray) );
            double alfa = 1;
            FloatArrayF<6> stress;
            while ( fabs(f) > 1.e-12 ) {
                double dSYdA = norm * this->evaluateCurrentPlasticModulus(k);
                dSYdA *= -norm;
                double denom = SFSr + FSr - dSYdA;
                double dAlfa = -f / denom;
                alfa += dAlfa;
                stress = trialEffectiveStress * alfa;
                k += ( 1 - alfa ) * norm;
                f = evaluatePlasCriterion(fabric, F, stress, k, ( 1 - alfa ) * norm, tStep);
            }

            tempEffectiveStress = stress;
            deltaKappa = k;
            toSolveScalar = 0;

            tmp = this->constructPlasFlowDirec(fabric, F, tempEffectiveStress);
            plasFlowDirec = tmp.first;
            norm = tmp.second;

            if ( tempEffectiveStress.giveSize() != trialEffectiveStress.giveSize() ) {
                printf( "tempS  %d \n", tempEffectiveStress.giveSize() );
                printf( "trial S %d \n", trialEffectiveStress.giveSize() );
            }

            toSolveTensor = dot( compliance, ( tempEffectiveStress - trialEffectiveStress ) ) + deltaKappa * plasFlowDirec;
            // Construction of the derivative of the plastic flow
            auto derivPlasFlowDirec = this->constructDerivativeOfPlasFlowDirec(fabric, F, tempEffectiveStress);
            // Construction of the gradient Nabla_S of R and SSa tensor
            SSaTensor = inv(derivPlasFlowDirec * deltaKappa + compliance);
        }

        /***********************************************************************************************/
        // iteration loop - solution of a set of nonlinear equations
        int flagLoop = 1;
        do {
            double plasModulus = evaluateCurrentPlasticModulus(tempKappa + deltaKappa);
            double viscoModulus = evaluateCurrentViscousModulus(deltaKappa, tStep);
            //*************************************
            //Evaluation of the Recursive Equations
            //*************************************
            double beta = dot(plasFlowDirec, dot(SSaTensor, plasFlowDirec)) + ( plasModulus - viscoModulus ) / norm; //+ viscoModulus;
            // Construction of the equation of Delta Kappa
            auto tempScalar = dot(plasFlowDirec, dot(SSaTensor, toSolveTensor));
            auto incKappa = ( toSolveScalar / norm - tempScalar ) / beta;
            auto incTempEffectiveStress = dot(SSaTensor, plasFlowDirec * incKappa + toSolveTensor);

            if ( lineSearchFlag == 1 ) {
                current_max_num_iter = 4000;
                double eta = 0.1;
                beta = 25.e-4;
                int jMax = 10;
                double alfa = 1;
                double M = ( errorR * errorR + errorF * errorF ) / 2;
                double dM = -2 * M;
                for ( int j = 0; ; ++j ) {
                    auto tempStress = incTempEffectiveStress * (-alfa) + tempEffectiveStress;
                    auto newDeltaKappa = deltaKappa + alfa * incKappa;
                    //*************************************
                    // Evaluation of the f and R
                    //*************************************
                    // Construction of the derivative of the plastic flow
                    auto derivPlasFlowDirec = this->constructDerivativeOfPlasFlowDirec(fabric, F, tempStress);

                    // Construction of the gradient Nabla_S of R and SSa tensor
                    SSaTensor = inv(derivPlasFlowDirec * newDeltaKappa + compliance);

                    // Evaluation of R
                    auto tmp = this->constructPlasFlowDirec(fabric, F, tempStress);
                    plasFlowDirec = tmp.first;
                    norm = tmp.second;  // FIXME: unusued in old code

                    toSolveTensor = dot( compliance, ( tempStress - trialEffectiveStress ) ) + newDeltaKappa * plasFlowDirec;
                    // Evaluation of f
                    //auto SFS = sqrt( dot(tempStress, dot(fabric, tempStress)) ); // FIXME: was unused in old code
                    toSolveScalar = evaluatePlasCriterion(fabric, F, tempStress, tempKappa + newDeltaKappa, newDeltaKappa, tStep);
                    //*************************************
                    // Evaluation of the error
                    //*************************************
                    errorR = sqrt(dot(toSolveTensor, toSolveTensor));
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
                //////////////////////////////////////////////////////////
                deltaKappa += incKappa;
                tempEffectiveStress -= incTempEffectiveStress;
                //*************************************
                // Evaluation of the f and R
                //*************************************
                // Construction of the derivative of the plastic flow
                //auto derivPlasFlowDirec = this->constructDerivativeOfPlasFlowDirec(fabric, F, tempEffectiveStress); // FIXME: unused in old code

                // Construction of the gradient Nabla_S of R and SSa tensor
                //auto SSaTensor = inv(derivPlasFlowDirec * deltaKappa + compliance); // FIXME: unused in old code

                // Evaluation of R
                auto tmp = this->constructPlasFlowDirec(fabric, F, tempEffectiveStress);
                plasFlowDirec = tmp.first;
                norm = tmp.second;
                toSolveTensor = dot( compliance, ( tempEffectiveStress - trialEffectiveStress ) ) + deltaKappa * plasFlowDirec;

                // Evaluation of f
                //auto SFS = sqrt( dot(tempEffectiveStress, dot(fabric, tempEffectiveStress)) );  // FIXME: was unsued in old code?
                toSolveScalar = evaluatePlasCriterion(fabric, F, tempEffectiveStress, tempKappa + deltaKappa, deltaKappa, tStep);
                //*************************************
                // Evaluation of the error
                //*************************************
                errorR = sqrt(dot(toSolveTensor, toSolveTensor));
                errorF = fabs(toSolveScalar);
            }

            if ( printflag ) {
                printf("   %d %g %g %g %g\n", flagLoop, tempEffectiveStress.at(1), tempEffectiveStress.at(3), incKappa, deltaKappa);
            }

            flagLoop++;
            convergence = ( fabs(errorF) < rel_yield_tol && errorR < strain_tol );
        } while ( flagLoop <= current_max_num_iter && !convergence );

        if ( convergence ) {
            auto plasModulus = evaluateCurrentPlasticModulus(tempKappa + deltaKappa);
            auto viscoModulus = evaluateCurrentViscousModulus(deltaKappa, tStep);
            auto tempTensor2 = dot(SSaTensor, plasFlowDirec);
            auto beta = dot(plasFlowDirec, tempTensor2) + ( plasModulus - viscoModulus ) / norm;
            tempPlasDef += deltaKappa * plasFlowDirec;
            tempKappa += deltaKappa;
            status->setBeta(beta);
            status->setPlasFlowDirec(plasFlowDirec);
            status->setSSaTensor(SSaTensor);
        }
    }

    return convergence;
}

std::pair<FloatArrayF<6>, double>
TrabBone3D :: constructPlasFlowDirec(FloatMatrixF<6,6> &fabric, const FloatArrayF<6> &F, const FloatArrayF<6> &S) const
{
    auto FFS = dot(fabric, S);
    //  FS = dot(F, S);
    double SFS = sqrt( dot(S, FFS) );
    // scaling matrix
    auto normAdjust = this->constructNormAdjustTensor();
    //direction of Np
    auto answer = F + FFS * (1. / SFS);
    auto tempTensor2 = dot(normAdjust, answer);
    //norm Np
    double norm = sqrt( dot(answer, tempTensor2) );
    //plastic flow
    answer *= 1.0 / norm;

    return {answer, norm};
}


FloatMatrixF<6,6>
TrabBone3D :: constructDerivativeOfPlasFlowDirec(const FloatMatrixF<6,6> &fabric, const FloatArrayF<6> &F, const FloatArrayF<6> &S) const
{
    // FIXME: unused in old code:
    //auto FFS = dot(fabric, S);
    //auto SFS = sqrt( dot(S, FFS) );
    //auto SGS = pow(SFS, -3.);
    //auto secondGradientOfG = dyad(FFS, FFS) * (-1. * SGS) + fabric * (1. / SFS);
    //to gradientOfG = FFS * (1. / SFS) + F;
    //auto h = norm(gradientOfG));
    //auto gradientOfH = Tdot(secondGradientOfG, gradientOfG) * (1. / h);
    //auto test = - dyad(gradientOfG, gradientOfH) + secondGradientOfG * (1./ h/.h);
    //secondGradientOfG *= h;
    //////////////////////////////////////////////////////////////////
    auto FFS = dot(fabric, S);
    auto SFS = sqrt( dot(S, FFS) );
    // scaling matrix
    auto normAdjust = this->constructNormAdjustTensor();
    //norm
    auto tempTensor2 = FFS * (1. / SFS) + F;
    auto tempTensor21 = dot(normAdjust, tempTensor2);
    //norm Np
    auto norm = sqrt( dot(tempTensor2, tempTensor21) );
    ///////////////////////////////////////////////////////////////////
    auto FSFS = dyad(FFS, FFS);
    auto dNp = (FSFS * (-1. / SFS / SFS) + fabric) * (1. / SFS / norm);
    //////////////////////////////////////////////////////////////
    auto FFSF = FFS * (1. / SFS) + F;
    /////////////////////////////////////////////////////////////////
    auto tempTensor4 = (FSFS * (-1. / SFS / SFS) + fabric) * (1. / SFS / norm);
    tempTensor2 = dot(normAdjust, FFSF);
    auto dNorm = dot(tempTensor4, tempTensor2);
    ///////////////////////////////////////////////////////////////////
    tempTensor4 = dyad(FFSF, dNorm) * (-1. / norm / norm);
    ////////////////////////////////////////////////////////////////////
    return dNp + tempTensor4;
}

double
TrabBone3D :: evaluatePlasCriterion(const FloatMatrixF<6,6> &fabric, const FloatArrayF<6> &F, const FloatArrayF<6> &stress, double kappa, double deltaKappa, TimeStep *tStep) const
{
    auto FFS = dot(fabric, stress);
    auto FS = dot(F, stress);
    auto SFS = sqrt( dot(stress, FFS) );
    return SFS + FS - evaluateCurrentYieldStress(kappa) + this->evaluateCurrentViscousStress(deltaKappa, tStep);
}

double
TrabBone3D :: computeDamageParam(double tempKappa) const
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
TrabBone3D :: computeDamageParamPrime(double tempKappa) const
{
    return critDam * expDam * exp(-expDam * tempKappa);
}
//
// END: FUNCTION FOR DAMAGE PARAMETER
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: FUNCTION FOR DAMAGE EVALUATION
//

double
TrabBone3D :: computeDamage(GaussPoint *gp, TimeStep *tStep) const
{
    double tempKappa = computeCumPlastStrain( gp, tStep);
    double tempDam = computeDamageParam(tempKappa);
    if ( tempDam < 0 ) {
        OOFEM_ERROR("negative damage");
    }

    return tempDam;
}


double TrabBone3D :: computeCumPlastStrain(GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< TrabBone3DStatus * >( this->giveStatus(gp) );
    return status->giveTempKappa();
}



FloatArrayF<6> TrabBone3D :: computeDensificationStress(GaussPoint *gp, const FloatArrayF<6> &strain, TimeStep *tStep) const
{
    auto status = static_cast< TrabBone3DStatus * >( this->giveStatus(gp) );

    FloatArrayF<6> answer;
    double traceLnU =  ( strain.at(1) + strain.at(2) + strain.at(3) );
    double g = traceLnU - densCrit;
    status->setDensG(g);
    if ( g <= 0 ) {
        double factor = gammaL0 * pow(rho, rL) * g + gammaP0 *pow(rho, rP) * pow(g, tDens - 1);
        answer.at(1) = answer.at(2) = answer.at(3) = factor;
    }
    return answer;
}


FloatArrayF<6>
TrabBone3D :: giveRealStressVector_3d(const FloatArrayF<6> &strain, GaussPoint *gp,
                                      TimeStep *tStep) const
{
    auto status = static_cast< TrabBone3DStatus * >( this->giveStatus(gp) );

    this->initTempStatus(gp);
    // compute effective stress using the plasticity model
    this->performPlasticityReturn(gp, strain, tStep);

    auto effStress = status->giveTempEffectiveStress();

    // evaluate damage variable
    double tempDam = computeDamage(gp, tStep);

    // transform effective stress into nominal stress
    auto stress = ( 1 - tempDam ) * effStress;

    computePlasStrainEnerDensity(gp, strain, stress);

    // add densification stress, if the densificator is activated
    if ( densCrit != 0 ) {
        stress += computeDensificationStress(gp, strain, tStep);
    }

    // store final damage, strain and stress in status
    status->setTempDam(tempDam);
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(stress);
    return stress;
}


FloatMatrixF<6,6>
TrabBone3D :: constructAnisoComplTensor() const
{
    double m3 = 3. - m1 - m2;
    double eps0k = eps0 * pow(rho, expk);
    double mu0k = mu0 * pow(rho, expk);
    double m1l = pow(m1, expl);
    double m2l = pow(m2, expl);
    double m3l = pow(m3, expl);

    FloatMatrixF<6,6> D;

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

#if 0
    /// Old code also computed this unused expression:
    auto t = this->constructStiffnessTransformationMatrix();
    auto stiffness = inv(D);
    auto TST = rotate(S, t);
#endif

    auto T = this->constructFabricTransformationMatrix();
    return rotate(D, T);
}


FloatMatrixF<6,6>
TrabBone3D :: constructAnisoStiffnessTensor() const
{
    double m3 = 3. - m1 - m2;
    double eps0k = eps0 * pow(rho, expk);
    double mu0k = mu0 * pow(rho, expk);
    double m1l = pow(m1, expl);
    double m2l = pow(m2, expl);
    double m3l = pow(m3, expl);

    double E1 = eps0k * m1l * m1l;
    double E2 = eps0k * m2l * m2l;
    double E3 = eps0k * m3l * m3l;

    double G23 = mu0k * m2l * m3l;
    double G13 = mu0k * m1l * m3l;
    double G12 = mu0k * m1l * m2l;

    double n23 = nu0 * m2l / m3l;
    double n12 = nu0 * m1l / m2l;
    double n31 = nu0 * m3l / m1l;
    double n32 = nu0 * m3l / m2l;
    double n21 = nu0 * m2l / m1l;
    double n13 = nu0 * m1l / m3l;

    double eksi = 1. - ( n12 * n21 + n23 * n32 + n31 * n13 ) - ( n12 * n23 * n31 + n21 * n32 * n13 );

    FloatMatrixF<6,6> D;
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

    D.at(4, 4) = G23;
    D.at(5, 5) = G13;
    D.at(6, 6) = G12;

    auto t = this->constructStiffnessTransformationMatrix();
    return rotate(D, t);
}


FloatMatrixF<6,6>
TrabBone3D :: constructAnisoFabricTensor() const
{
    double S0 = ( sig0Pos + sig0Neg ) / ( 2. * sig0Pos * sig0Neg );
    double rhoP = pow(rho, 2. * expp);
    double m3 = 3. - m1 - m2;
    double m1q = pow(m1, 2. * expq);
    double m2q = pow(m2, 2. * expq);
    double m3q = pow(m3, 2. * expq);

    FloatMatrixF<6,6> F;
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

    auto T = this->constructFabricTransformationMatrix();    
    return rotate(F, T);
}


FloatArrayF<6>
TrabBone3D :: constructAnisoFtensor() const
{
    double rhoP = pow(rho, expp);
    double m3 = 3. - m1 - m2;
    double m1q = pow(m1, 2. * expq);
    double m2q = pow(m2, 2. * expq);
    double m3q = pow(m3, 2. * expq);

    FloatArrayF<6> F;
    F.at(1) = -( sig0Pos - sig0Neg ) / ( 2. * sig0Pos * sig0Neg * rhoP * m1q );
    F.at(2) = -( sig0Pos - sig0Neg ) / ( 2. * sig0Pos * sig0Neg * rhoP * m2q );
    F.at(3) = -( sig0Pos - sig0Neg ) / ( 2. * sig0Pos * sig0Neg * rhoP * m3q );
    
    auto T = this->constructFabricTransformationMatrix();
    return Tdot(T, F);
}


FloatMatrixF<6,6>
TrabBone3D :: constructStiffnessTransformationMatrix() const
{
    FloatMatrixF<6,6> answer;

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
    
    return answer;
}


FloatMatrixF<6,6>
TrabBone3D :: constructNormAdjustTensor() const
{
    FloatMatrixF<6,6> answer;

    for ( int i = 1; i <= 3; i++ ) {
        answer.at(i, i) = 1.;
    }

    for ( int i = 4; i <= 6; i++ ) {
        answer.at(i, i) = 0.5;
    }
    return answer;
}



FloatMatrixF<6,6>
TrabBone3D :: constructFabricTransformationMatrix() const
{
    FloatMatrixF<6,6> answer;

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
    
    return answer;
}

void
TrabBone3D :: initializeFrom(InputRecord &ir)
{
    StructuralMaterial :: initializeFrom(ir);

    // Mandatory parameters
    IR_GIVE_FIELD(ir, eps0, _IFT_TrabBone3D_eps0);
    IR_GIVE_FIELD(ir, nu0, _IFT_TrabBone3D_nu0);
    IR_GIVE_FIELD(ir, mu0, _IFT_TrabBone3D_mu0);
    IR_GIVE_FIELD(ir, expk, _IFT_TrabBone3D_expk);
    IR_GIVE_FIELD(ir, expl, _IFT_TrabBone3D_expl);


    IR_GIVE_FIELD(ir, m1, _IFT_TrabBone3D_m1);
    IR_GIVE_FIELD(ir, m2, _IFT_TrabBone3D_m2);
    IR_GIVE_FIELD(ir, rho, _IFT_TrabBone3D_rho);

    IR_GIVE_FIELD(ir, sig0Pos, _IFT_TrabBone3D_sig0Pos);
    IR_GIVE_FIELD(ir, sig0Neg, _IFT_TrabBone3D_sig0Neg);
    IR_GIVE_FIELD(ir, chi0Pos, _IFT_TrabBone3D_chi0Pos);
    IR_GIVE_FIELD(ir, tau0, _IFT_TrabBone3D_tau0);



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
    //z'
    z1 = x2 * y3 - x3 * y2;
    z2 = x3 * y1 - x1 * y3;
    z3 = x1 * y2 - x2 * y1;
    //viscosity parameter
    viscosity = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, viscosity, _IFT_TrabBone3D_viscosity);

    // optional control parameters for printing and convergence
    printflag = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, printflag, _IFT_TrabBone3D_printflag);
    max_num_iter = 100;
    IR_GIVE_OPTIONAL_FIELD(ir, max_num_iter, _IFT_TrabBone3D_max_num_iter);
    rel_yield_tol = 1.e-9;
    IR_GIVE_OPTIONAL_FIELD(ir, rel_yield_tol, _IFT_TrabBone3D_rel_yield_tol);
    strain_tol = 1.e-9;
    IR_GIVE_OPTIONAL_FIELD(ir, strain_tol, _IFT_TrabBone3D_strain_tol);
}



int
TrabBone3D :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    auto status = static_cast< TrabBone3DStatus * >( this->giveStatus(gp) );
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
    auto status = static_cast< TrabBone3DStatus * >( this->giveStatus(gp) );

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



TrabBone3DStatus :: TrabBone3DStatus(GaussPoint *g) : StructuralMaterialStatus(g)
{
}


void
TrabBone3DStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
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
    densG = 1;
    this->tempKappa = this->kappa;
    this->tempDam = this->dam;
    this->tempTSED = this->tsed;
    this->tempPlasDef = this->plasDef;
}


void
TrabBone3DStatus :: updateYourself(TimeStep *tStep)
{
    StructuralMaterialStatus :: updateYourself(tStep);
    this->kappa = this->tempKappa;
    this->dam = this->tempDam;
    this->tsed = this->tempTSED;
    this->plasDef = this->tempPlasDef;
}


void
TrabBone3DStatus :: saveContext(DataStream &stream, ContextMode mode)
{
    contextIOResultType iores;

    // save parent class status
    StructuralMaterialStatus :: saveContext(stream, mode);

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
}


void
TrabBone3DStatus :: restoreContext(DataStream &stream, ContextMode mode)
{
    contextIOResultType iores;

    // read parent class status
    StructuralMaterialStatus :: restoreContext(stream, mode);

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
}


std::unique_ptr<MaterialStatus> TrabBone3D :: CreateStatus(GaussPoint *gp) const
{
    return std::make_unique<TrabBone3DStatus>(gp);
}
} //end namespace oofem
