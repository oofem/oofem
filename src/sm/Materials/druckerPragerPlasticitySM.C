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

#include "druckerPragerPlasticitySM.h"

#include "floatarray.h"
#include "floatmatrix.h"
#include "sm/Materials/structuralms.h"
#include "gausspoint.h"
#include "intarray.h"
#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/isolinearelasticmaterial.h"
#include "datastream.h"
#include "contextioerr.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(DruckerPragerPlasticitySM);

DruckerPragerPlasticitySMStatus :: DruckerPragerPlasticitySMStatus(GaussPoint *gp) :
    StructuralMaterialStatus(gp)
{
    stressVector.resize(6);
    strainVector.resize(6);
    tempStressVector = stressVector;
    tempStrainVector = strainVector;
}

void
DruckerPragerPlasticitySMStatus :: initTempStatus()
{
    // Call the function of the parent class to initialize the variables defined there.
    StructuralMaterialStatus :: initTempStatus();

    tempPlasticStrainDeviator = plasticStrainDeviator;
    tempVolumetricPlasticStrain = volumetricPlasticStrain;
    tempKappa = kappa;
    temp_state_flag = state_flag;
}

void
DruckerPragerPlasticitySMStatus :: updateYourself(TimeStep *tStep)
{
    // Call corresponding function of the parent class to update variables defined there.
    StructuralMaterialStatus :: updateYourself(tStep);

    // update variables defined in DruckerPragerPlasticitySMStatus
    plasticStrainDeviator = tempPlasticStrainDeviator;
    volumetricPlasticStrain = tempVolumetricPlasticStrain;
    kappa = tempKappa;
    state_flag = temp_state_flag;
}

void
DruckerPragerPlasticitySMStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
{
    // Call corresponding function of the parent class to print variables defined there.
    StructuralMaterialStatus :: printOutputAt(file, tStep);

    fprintf(file, "\tstatus { ");

    // print status flag
    switch ( state_flag ) {
    case DruckerPragerPlasticitySMStatus :: DP_Elastic:
        fprintf(file, " Elastic, ");
        break;
    case DruckerPragerPlasticitySMStatus :: DP_Yielding:
        fprintf(file, " Yielding, ");
        break;
    case DruckerPragerPlasticitySMStatus :: DP_Vertex:
        fprintf(file, " Vertex_return, ");
        break;
    case DruckerPragerPlasticitySMStatus :: DP_Unloading:
        fprintf(file, " Unloading, ");
        break;
    }

    // print plastic strain vector
    auto plasticStrainVector = givePlasticStrainVector();

    fprintf(file, "plasticStrains ");
    for ( auto &val : plasticStrainVector ) {
        fprintf( file, " %.4e", val );
    }

    fprintf(file, "}\n");

    fprintf(file, "\t\thardening_parameter ");
    // print hardening parameter
    fprintf(file, " %.4e\n", kappa);
}

void
DruckerPragerPlasticitySMStatus :: saveContext(DataStream &stream, ContextMode mode)
{
    StructuralMaterialStatus :: saveContext(stream, mode);

    contextIOResultType iores;
    if ( ( iores = plasticStrainDeviator.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream.write(volumetricPlasticStrain) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(temp_state_flag) ) {
        THROW_CIOERR(CIO_IOERR);
    }
}


void
DruckerPragerPlasticitySMStatus :: restoreContext(DataStream &stream, ContextMode mode)
{
    StructuralMaterialStatus :: restoreContext(stream, mode);

    contextIOResultType iores;
    if ( ( iores = plasticStrainDeviator.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream.read(volumetricPlasticStrain) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(temp_state_flag) ) {
        THROW_CIOERR(CIO_IOERR);
    }
}

//   *************************************************************
//   *** CLASS DRUCKER-PRAGER PLASTICITY STRUCTURAL MATERIAL   ***
//   *************************************************************


DruckerPragerPlasticitySM :: DruckerPragerPlasticitySM(int n, Domain *d) :
    StructuralMaterial(n, d),
    LEMaterial(n, d)
{}


void
DruckerPragerPlasticitySM :: initializeFrom(InputRecord &ir)
{
    // call the corresponding service of structural material
    StructuralMaterial :: initializeFrom(ir);
    // call the corresponding service for the linear elastic material
    LEMaterial.initializeFrom(ir);

    // initialize elastic constants
    //eM = LEMaterial.give(Ex,gp);
    //nu = LEMaterial.give(NYxz,gp);
    //gM = eM / ( 2. * ( 1. + nu ) );
    //kM = eM / ( 3. * ( 1. - 2. * nu ) );

    // instanciate the variables defined in DruckerPragerPlasticitySM
    IR_GIVE_FIELD(ir, initialYieldStress, _IFT_DruckerPragerPlasticitySM_iys); // initial yield stress under pure shear
    IR_GIVE_FIELD(ir, alpha, _IFT_DruckerPragerPlasticitySM_alpha); // friction coefficient
    IR_GIVE_FIELD(ir, alphaPsi, _IFT_DruckerPragerPlasticitySM_alphapsi); //dilatancy coefficient
    // this is valid for strain hardening/softening only (not for work hardening/softening)
    kFactor = sqrt(1. / 3. + 2. * alphaPsi * alphaPsi);

    IR_GIVE_FIELD(ir, hardeningType, _IFT_DruckerPragerPlasticitySM_ht);

    switch ( hardeningType ) {
    case 1:     // linear hardening/softening
        IR_GIVE_FIELD(ir, hardeningModulus, _IFT_DruckerPragerPlasticitySM_hm);
        break;
    case 2:     // exponential hardening/softening
        IR_GIVE_FIELD(ir, kappaC, _IFT_DruckerPragerPlasticitySM_kc);
        IR_GIVE_FIELD(ir, limitYieldStress, _IFT_DruckerPragerPlasticitySM_lys);
        break;
    default:
        throw ValueInputException(ir, _IFT_DruckerPragerPlasticitySM_ht,
                                  "must be 1 (linear hardening/softening) or 2 (exponential hardening/softening)");
    }

    yieldTol = 1.e-14;
    IR_GIVE_OPTIONAL_FIELD(ir, yieldTol, _IFT_DruckerPragerPlasticitySM_yieldtol);
    newtonIter = 30;
    IR_GIVE_OPTIONAL_FIELD(ir, newtonIter, _IFT_DruckerPragerPlasticitySM_newtoniter);
}

FloatArrayF<6>
DruckerPragerPlasticitySM :: giveRealStressVector_3d(const FloatArrayF<6> &strain, GaussPoint *gp,
                                                  TimeStep *tStep) const
{
    auto status = static_cast< DruckerPragerPlasticitySMStatus * >( this->giveStatus(gp) );

    // Initialize temp variables for this gauss point
    const_cast<DruckerPragerPlasticitySM*>(this)->initTempStatus(gp);

    // subtract stress independent part
    // note: eigenStrains (temperature) is not contained in mechanical strain stored in gp
    // therefore it is necessary to subtract always the total eigen strain value
    auto thermalStrain = this->computeStressIndependentStrainVector_3d(gp, tStep, VM_Total);
    auto strainVectorR = strain - thermalStrain;

    // perform the local stress return and update the history variables
    performLocalStressReturn(gp, strainVectorR);

    // copy total strain vector to the temp status
    status->letTempStrainVectorBe(strain);

    // give back correct form of stressVector to giveRealStressVector
    return status->giveTempStressVector();
}

void
DruckerPragerPlasticitySM :: performLocalStressReturn(GaussPoint *gp, const FloatArrayF<6> &strain) const
{
    auto status = static_cast< DruckerPragerPlasticitySMStatus * >( giveStatus(gp) );

    // split total strains in volumetric and deviatoric part
    // elastic constants
    double eM = LEMaterial.give(Ex, gp);
    double nu = LEMaterial.give(NYxz, gp);
    double gM = eM / ( 2. * ( 1. + nu ) );
    double kM = eM / ( 3. * ( 1. - 2. * nu ) );

    //auto [strainDeviator, volumetricStrain] = computeDeviatoricVolumetricSplit(strain); // c++17
    auto tmp = computeDeviatoricVolumetricSplit(strain);
    auto strainDeviator = tmp.first;
    auto volumetricStrain = tmp.second;

    // compute trial elastic strains
    double volumetricElasticTrialStrain = volumetricStrain - status->giveVolumetricPlasticStrain();
    auto plasticStrainDeviator = status->givePlasticStrainDeviator();
    auto elasticStrainDeviator = strainDeviator - plasticStrainDeviator;

    // compute trial stresses
    double volumetricStress = 3. * kM * volumetricElasticTrialStrain;
    auto stressDeviator = applyDeviatoricElasticStiffness(elasticStrainDeviator, gM);
    // norm of trial stress deviator
    double trialStressJTwo = computeSecondStressInvariant(stressDeviator);

    // initialize hardening parameter
    double kappa = status->giveKappa();
    double tempKappa = kappa;

    // choose and perform correct stress return and update state flag
    if ( computeYieldValue(volumetricStress, trialStressJTwo, tempKappa, eM) / eM
         > yieldTol ) {
      //printf("*");
        if ( checkForVertexCase(eM, gM, kM, trialStressJTwo, volumetricStress, tempKappa) ) {
            performVertexReturn(eM, gM, kM, trialStressJTwo, stressDeviator, volumetricStress, 
                                tempKappa, volumetricElasticTrialStrain, kappa);
            status->letTempStateFlagBe(DruckerPragerPlasticitySMStatus :: DP_Vertex);
        } else {
            performRegularReturn(eM, gM, kM, trialStressJTwo, stressDeviator, volumetricStress, tempKappa);
            status->letTempStateFlagBe(DruckerPragerPlasticitySMStatus :: DP_Yielding);
        }
    } else {
        const int state_flag = status->giveStateFlag();
        if ( state_flag == DruckerPragerPlasticitySMStatus :: DP_Elastic ) {
            status->letTempStateFlagBe(DruckerPragerPlasticitySMStatus :: DP_Elastic);
        } else {
            status->letTempStateFlagBe(DruckerPragerPlasticitySMStatus :: DP_Unloading);
        }
    }

    // update kappa
    status->letTempKappaBe(tempKappa);

    // compute full stresses from deviatoric and volumetric part and store them
    auto stress = computeDeviatoricVolumetricSum(stressDeviator, volumetricStress);
    status->letTempStressVectorBe(stress);

    // compute and update plastic strains, volumetric and deviatoric part
    elasticStrainDeviator = applyDeviatoricElasticCompliance(stressDeviator, gM);
    plasticStrainDeviator = strainDeviator - elasticStrainDeviator;
    status->letTempPlasticStrainDeviatorBe(plasticStrainDeviator);
    status->letTempVolumetricPlasticStrainBe(volumetricStrain - volumetricStress / 3. / kM);
}

bool
DruckerPragerPlasticitySM :: checkForVertexCase(double eM, double gM, double kM, double trialStressJTwo, double volumetricStress, double tempKappa) const
{
    // delta lambda max corresponds to the maximum value
    // of the rate of the plastic multiplier for regular plastic flow
    // and allows to distinguish between regular return and vertex case
    double deltaLambdaMax = sqrt(trialStressJTwo) / gM;

    // vertex case:
    // yield value positive under the assumption of maximum regular plastic flow
    double volConstant = 3. * kM * alphaPsi;
    double yieldValue =
        computeYieldValue(volumetricStress - volConstant * deltaLambdaMax,
                          0., tempKappa + kFactor * deltaLambdaMax, eM);
    if ( yieldValue / eM > -yieldTol ) {
        return true;
    }

    // regular case
    return false;
}

void
DruckerPragerPlasticitySM :: performRegularReturn(double eM, double gM, double kM, double trialStressJTwo, FloatArrayF<6> &stressDeviator, double &volumetricStress, double &tempKappa) const
{
    // delta lambda max controls the maximum plastic flow, see below
    double deltaLambdaMax = sqrt(trialStressJTwo) / gM;

    // declare some constants for faster use
    double volConstant = 3. * kM * alphaPsi;
    double devConstant = sqrt(2.) * gM;
    // yield value prime is derivative of yield value with respect to deltaLambda
    double yieldValuePrimeZero = -9. * alpha * alphaPsi * kM - gM;

    auto flowDir = stressDeviator * ( 1. / sqrt(2. * trialStressJTwo) );

    // initialize Newton iteration
    int iterationCount = 0;
    double deltaLambda = 0.;
    double deltaLambdaIncrement = 0.;
    double yieldValue = computeYieldValue(volumetricStress, trialStressJTwo, tempKappa, eM);
    double newtonError = fabs(yieldValue / eM);
    //printf("\nnewtonError = %e\n", newtonError) ;
    // Newton iteration to find deltaLambda
    while ( newtonError > yieldTol ) {
        if ( ++iterationCount > newtonIter ) {
            OOFEM_ERROR("Newton iteration for deltaLambda (regular stress return) did not converge after newtonIter iterations. You might want to try increasing the optional parameter newtoniter or yieldtol in the material record of your input file.");
        }

        double yieldValuePrime = yieldValuePrimeZero - kFactor *computeYieldStressPrime(tempKappa, eM);
        deltaLambdaIncrement = -yieldValue / yieldValuePrime;

        // deltaLambdaMax may be exceeded if the yield stress has almost vanished
        // If this happens, the stress deviator will evolve too much,
        // and will then be on the other side of the hydrostatic axis.
        // This causes the failure of the Newton-iteration and has to be avoided.
        if ( deltaLambda + deltaLambdaIncrement > deltaLambdaMax ) {
            OOFEM_LOG_DEBUG("Special case in Newton-iteration for regular return. This may cause loss of quadratic convergence.\n");
            deltaLambdaIncrement = deltaLambdaMax - deltaLambda;
        }

        deltaLambda += deltaLambdaIncrement;
        tempKappa += kFactor * deltaLambdaIncrement;
        volumetricStress -= volConstant * deltaLambdaIncrement;

        // auto plasticFlow = flowDir * (devConstant * deltaLambdaIncrement)
        //stressDeviator -= plasticFlow;

        stressDeviator += (-devConstant * deltaLambdaIncrement) * flowDir;
        double tempJTwo = computeSecondStressInvariant(stressDeviator);
        yieldValue = computeYieldValue(volumetricStress, tempJTwo, tempKappa, eM);
        newtonError = fabs(yieldValue / eM);
        //printf("newtonError = %e\n", newtonError) ;
    }

    OOFEM_LOG_DEBUG("IterationCount in regular return = %d\n", iterationCount);

    if ( deltaLambda < 0. ) {
        OOFEM_ERROR("Fatal error in the Newton iteration for regular stress return. deltaLambda is evaluated as negative, but should always be positive. This is most likely due to a softening law with local snapback, which is physically inadmissible.n");
    }
}

void
DruckerPragerPlasticitySM :: performVertexReturn(double eM, double gM, double kM, double trialStressJTwo, FloatArrayF<6> &stressDeviator, double &volumetricStress, double &tempKappa, double volumetricElasticTrialStrain, double kappa) const
{
    // declare some constants for faster use
    // yield value prime is derivative of yield value with respect to deltaLambda
    double yieldValuePrimeZero = 3. * alpha;
    // contribution of deviatoric strain to hardening
    double deviatorContribution = trialStressJTwo / 3. / gM / gM;
    // in the vertex case, deviatoric stresses are zero
    stressDeviator = zeros<6>();
    volumetricStress = 3. * kM * volumetricElasticTrialStrain;

    // initialize Newton iteration
    int iterationCount = 0;
    double deltaVolumetricStress = 0.;
    double deltaVolumetricStressIncrement = 0.;
    double deltaKappa = sqrt(deviatorContribution);
    tempKappa = kappa + deltaKappa;
    double yieldValue = computeYieldValue(volumetricStress, 0., tempKappa, eM);
    double newtonError = fabs(yieldValue / eM);

    // Newton iteration to find deltaLambda
    while ( newtonError > yieldTol ) {
        if ( ++iterationCount > newtonIter ) {
            OOFEM_ERROR("Newton iteration for deltaLambda (vertex stress return) did not converge after newtonIter iterations. You might want to try increasing the optional parameter newtoniter or yieldtol in the material record of your input file.");
        }

        // exclude division by zero
        double yieldValuePrime;
        if ( deltaKappa == 0. ) {
            yieldValuePrime = yieldValuePrimeZero
                              - sqrt(2.) / 3. / kM *computeYieldStressPrime(tempKappa, eM);
        } else {
            yieldValuePrime = yieldValuePrimeZero
                              - 2. / 9. / kM / kM *computeYieldStressPrime(tempKappa, eM)
            * deltaVolumetricStress / deltaKappa;
        }

        deltaVolumetricStressIncrement = -yieldValue / yieldValuePrime;
        deltaVolumetricStress += deltaVolumetricStressIncrement;
        volumetricStress += deltaVolumetricStressIncrement;
        deltaKappa = sqrt(2. / 9. / kM / kM * deltaVolumetricStress * deltaVolumetricStress
                          + deviatorContribution);
        tempKappa = kappa + deltaKappa;
        yieldValue = computeYieldValue(volumetricStress, 0., tempKappa, eM);
        newtonError = fabs(yieldValue / eM);
        OOFEM_LOG_DEBUG("NewtonError in iteration %d in vertex return = %e\n", iterationCount, newtonError);
    }

    OOFEM_LOG_DEBUG("Done iteration in vertex return, after %d\n", iterationCount);

    if ( deltaKappa < 0. ) {
        OOFEM_ERROR("Fatal error in the Newton iteration for vertex stress return. deltaKappa is evaluated as negative, but should always be positive. This is most likely due to a softening law with a local snapback, which is physically inadmissible.");
    }
}

double
DruckerPragerPlasticitySM :: computeYieldValue(double volumetricStress,
                                               double JTwo,
                                               double kappa,
                                               double eM) const
{
    return 3. * alpha * volumetricStress + sqrt(JTwo) - computeYieldStressInShear(kappa, eM);
}

double
DruckerPragerPlasticitySM :: computeYieldStressInShear(double kappa, double eM) const
{
    double yieldStress;
    switch ( hardeningType ) {
    case 1:     // linear hardening/softening.
        yieldStress = initialYieldStress + hardeningModulus * eM * kappa;
        if ( yieldStress < 0. ) {
            yieldStress = 0.;
            //printf("Yield stress zero reached in computeYieldStressInShear.\n") ;
        }
        return yieldStress;
    case 2:     // exponential hardening
        return limitYieldStress - ( limitYieldStress - initialYieldStress ) * exp(-kappa / kappaC);
    default:
        //StructuralMaterial :: OOFEM_ERROR( "Case failed: choose linear hardening/softening (1), exponential hardening/softening (2) in input file.") ;
        OOFEM_ERROR("Case failed: choose linear hardening/softening (1), exponential hardening/softening (2) in input file.");
        return 0.;
    }
}

double
DruckerPragerPlasticitySM :: computeYieldStressPrime(double kappa, double eM) const
{
    switch ( hardeningType ) {
    case 1:             // linear hardening/softening.
        // if the limit value for kappa is exceeded in the softening case, the derivative is zero
        if ( ( hardeningModulus < 0. ) && ( kappa >= -initialYieldStress / hardeningModulus / eM ) ) {
            return 0.0;
        } else {
            return eM * hardeningModulus;
        }
    case 2:             // exponential hardening
        return ( limitYieldStress - initialYieldStress ) / kappaC *exp(-kappa / kappaC);
    default:
        //StructuralMaterial :: OOFEM_ERROR( "Case failed: choose linear hardening/softening (1), exponential hardening/softening (2) in input file.") ;
        OOFEM_ERROR("Case failed: choose linear hardening/softening (1), exponential hardening/softening (2) in input file.");
        return 0.;
    }
}


FloatMatrixF<6,6>
DruckerPragerPlasticitySM :: give3dMaterialStiffnessMatrix(MatResponseMode mode,
                                                           GaussPoint *gp,
                                                           TimeStep *tStep) const
{
    switch ( mode ) {
    case ElasticStiffness:
        return LEMaterial.give3dMaterialStiffnessMatrix(mode, gp, tStep);

    case SecantStiffness:
        return LEMaterial.give3dMaterialStiffnessMatrix(mode, gp, tStep);

    case TangentStiffness:
        switch ( ( static_cast< DruckerPragerPlasticitySMStatus * >( this->giveStatus(gp) ) )
                ->giveTempStateFlag() ) {
        case DruckerPragerPlasticitySMStatus :: DP_Elastic:        // elastic stiffness
        case DruckerPragerPlasticitySMStatus :: DP_Unloading:        // elastic stiffness
            return LEMaterial.give3dMaterialStiffnessMatrix(mode, gp, tStep);
        case DruckerPragerPlasticitySMStatus :: DP_Yielding:
            // elasto-plastic stiffness for regular case
            //printf("\nAssembling regular algorithmic stiffness matrix.") ;
            return giveRegAlgorithmicStiffMatrix(mode, gp, tStep);
        case DruckerPragerPlasticitySMStatus :: DP_Vertex:
            // elasto-plastic stiffness for vertex case
            //printf("\nAssembling vertex case algorithmic stiffness matrix.") ;
            return giveVertexAlgorithmicStiffMatrix(mode, gp, tStep);
        default:
            OOFEM_ERROR("Case did not match.");
        }

        break;

    default:
        OOFEM_ERROR("Switch failed: Only elastic and tangent stiffness are supported.");
        break;
    }
    return FloatMatrixF<6,6>();
}

FloatMatrixF<6,6>
DruckerPragerPlasticitySM :: giveRegAlgorithmicStiffMatrix(MatResponseMode mode,
                                                           GaussPoint *gp,
                                                           TimeStep *tStep) const
{
    auto status = static_cast< DruckerPragerPlasticitySMStatus * >( this->giveStatus(gp) );

    const auto &stressVector = status->giveTempStressVector();

    auto deviatoricStress = computeDeviator(stressVector);
    
    double normOfStress = sqrt( 2. * computeSecondStressInvariant(deviatoricStress) );
    // elastic constants
    double eM = LEMaterial.give(Ex, gp);
    double nu = LEMaterial.give(NYxz, gp);
    double gM = eM / ( 2. * ( 1. + nu ) );
    double kM = eM / ( 3. * ( 1. - 2. * nu ) );

    auto flowDir = deviatoricStress / normOfStress;

    double kappa = status->giveKappa();
    double tempKappa = status->giveTempKappa();
    double deltaKappa = tempKappa - kappa;
    double deltaLambdaStar = sqrt(2.) * gM * deltaKappa / kFactor / normOfStress;
    double hStar = kFactor * computeYieldStressPrime(tempKappa, eM);

    //exclude division by zero
    if ( hStar == 0. ) {
        OOFEM_ERROR("computeYieldStressPrime is zero. This happens mainly due to excessive softening.");
    }

    double a_const = 1. + deltaLambdaStar;
    double b_const = 3. * alpha * alphaPsi * kM / hStar - deltaLambdaStar / 3.;
    double c_const = 3. * sqrt(2.) * alphaPsi * kM / 2. / hStar;
    double d_const = sqrt(2.) * alpha * gM / hStar;
    double e_const = gM / hStar - deltaLambdaStar;

    FloatMatrixF<6,6> A_Matrix;

    for ( int i = 1; i < 7; i++ ) {
        A_Matrix.at(i, i) = a_const;
    }

    for ( int i = 1; i < 4; i++ ) {
        for ( int j = 1; j < 4; j++ ) {
            A_Matrix.at(i, j) += b_const;
        }
    }

    for ( int i = 1; i < 4; i++ ) {
        for ( int j = 1; j < 4; j++ ) {
            A_Matrix.at(i, j) += c_const * flowDir.at(j);
        }

        for ( int j = 4; j < 7; j++ ) {
            A_Matrix.at(i, j) += 2. *c_const *flowDir.at(j);
        }
    }

    for ( int i = 1; i < 7; i++ ) {
        for ( int j = 1; j < 4; j++ ) {
            A_Matrix.at(i, j) += d_const * flowDir.at(i);
        }
    }

    for ( int i = 1; i < 7; i++ ) {
        for ( int j = 1; j < 4; j++ ) {
            A_Matrix.at(i, j) += e_const * flowDir.at(i) * flowDir.at(j);
        }

        // account for engineering notation
        for ( int j = 4; j < 7; j++ ) {
            A_Matrix.at(i, j) += 2. *e_const *flowDir.at(i) * flowDir.at(j);
        }
    }

    auto De = LEMaterial.give3dMaterialStiffnessMatrix(mode, gp, tStep);

    // answer is A_Matrix^-1 * De
    return dot(inv(A_Matrix), De);
}

FloatMatrixF<6,6>
DruckerPragerPlasticitySM :: giveVertexAlgorithmicStiffMatrix(MatResponseMode mode,
                                                              GaussPoint *gp,
                                                              TimeStep *tStep) const
{
    auto status = static_cast< DruckerPragerPlasticitySMStatus * >( this->giveStatus(gp) );

    double tempKappa = status->giveTempKappa();
    double deltaKappa = tempKappa - status->giveKappa();
    // elastic constants
    double eM = LEMaterial.give(Ex, gp);
    double nu = LEMaterial.give(NYxz, gp);
    //double gM = eM / ( 2. * ( 1. + nu ) );
    double kM = eM / ( 3. * ( 1. - 2. * nu ) );

    if ( deltaKappa <= 0. ) {
        // This case occurs in the first iteration of a step.
        // printf("deltaKappa<=0. for vertex case algorithmic stiffness, i.e. continuum tangent stiffness. Since the continuum tangent stiffness does not exist at the vertex, elastic stiffness is used instead. This will cause the loss of quadratic convergence.\n") ;
        return LEMaterial.give3dMaterialStiffnessMatrix( mode, gp, tStep);
        
    }

    double deltaVolumetricPlasticStrain =
        status->giveTempVolumetricPlasticStrain() - status->giveVolumetricPlasticStrain();
    double HBar = computeYieldStressPrime(tempKappa, eM);

    // compute elastic trial strain deviator of latest temp-state
    auto strainDeviator = computeDeviator(status->giveTempStrainVector());

    auto elasticStrainDeviator = strainDeviator - FloatArrayF<6>(status->givePlasticStrainDeviator());

    double a_const =
        kM * HBar / ( HBar * deltaVolumetricPlasticStrain + 9. / 2. * alpha * kM * deltaKappa );

    if ( ( HBar * deltaVolumetricPlasticStrain + 9. / 2. * alpha * kM * deltaKappa ) == 0. ) {
        OOFEM_ERROR("Tangent type is singular, material ID %d\n", this->giveNumber() );
    }
    // compute the algorithmic tangent stiffness


    FloatMatrixF<6,6> answer;
    for ( int i = 1; i < 4; i++ ) {
        for ( int j = 1; j < 4; j++ ) {
            answer.at(i, j) = deltaVolumetricPlasticStrain;
        }
    }

    for ( int i = 1; i < 4; i++ ) {
        for ( int j = 1; j < 4; j++ ) {
            answer.at(i, j) += elasticStrainDeviator.at(j);
        }

        for ( int j = 4; j < 7; j++ ) {
            answer.at(i, j) += .5 * elasticStrainDeviator.at(j);
        }
    }

    return a_const * answer;
}

int
DruckerPragerPlasticitySM :: giveIPValue(FloatArray &answer,
                                         GaussPoint *gp,
                                         InternalStateType type,
                                         TimeStep *tStep)
{
    const auto status = static_cast< DruckerPragerPlasticitySMStatus * >( giveStatus(gp) );

    switch ( type ) {
    case IST_PlasticStrainTensor:
        answer = status->givePlasticStrainVector();
        return 1;

    case IST_DamageTensor:
        answer.resize(6);
        answer.zero();
        answer.at(1) = answer.at(2) = answer.at(3) = status->giveKappa();
        return 1;

    default:
        return StructuralMaterial :: giveIPValue(answer, gp, type, tStep);
    }
}

MaterialStatus *
DruckerPragerPlasticitySM :: CreateStatus(GaussPoint *gp) const
{
    return new DruckerPragerPlasticitySMStatus(gp);
}


double
DruckerPragerPlasticitySM :: predictRelativeComputationalCost(GaussPoint *gp)
{
    auto status = static_cast< DruckerPragerPlasticitySMStatus * >( giveStatus(gp) );
    int state_flag = status->giveStateFlag();

    if ( ( state_flag == DruckerPragerPlasticitySMStatus :: DP_Vertex ) ||
        ( state_flag == DruckerPragerPlasticitySMStatus :: DP_Yielding ) ) {
        return 20.;
    } else {
        return 1.0;
    }
}

} // end namespace oofem

#undef PRINTFDP
