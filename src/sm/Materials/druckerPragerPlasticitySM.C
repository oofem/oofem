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
#include "../sm/Materials/structuralms.h"
#include "gausspoint.h"
#include "intarray.h"
#include "../sm/Materials/structuralmaterial.h"
#include "Materials/isolinearelasticmaterial.h"
#include "datastream.h"
#include "contextioerr.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(DruckerPragerPlasticitySM);

DruckerPragerPlasticitySMStatus :: DruckerPragerPlasticitySMStatus(int n, Domain *d, GaussPoint *gp) :
    StructuralMaterialStatus(n, d, gp),
    plasticStrainDeviator( 6 ),
    tempPlasticStrainDeviator( 6 )
{
    stressVector.resize(6);
    strainVector.resize(6);
    tempStressVector = stressVector;
    tempStrainVector = strainVector;

    kappa = tempKappa = 0.;
    state_flag = temp_state_flag = DruckerPragerPlasticitySMStatus :: DP_Elastic;
    volumetricPlasticStrain = tempVolumetricPlasticStrain = 0.;
}

DruckerPragerPlasticitySMStatus :: ~DruckerPragerPlasticitySMStatus()
{ }

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
DruckerPragerPlasticitySMStatus :: printOutputAt(FILE *file, TimeStep *tStep)
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
    FloatArray plasticStrainVector;
    givePlasticStrainVector(plasticStrainVector);

    fprintf(file, "plasticStrains ");
    for ( auto &val : plasticStrainVector ) {
        fprintf( file, " %.4e", val );
    }

    fprintf(file, "}\n");

    fprintf(file, "\t\thardening_parameter ");
    // print hardening parameter
    fprintf(file, " %.4e\n", kappa);
}

contextIOResultType
DruckerPragerPlasticitySMStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write raw data
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

    return CIO_OK;
}


contextIOResultType
DruckerPragerPlasticitySMStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
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

    return CIO_OK;
}

//   *************************************************************
//   *** CLASS DRUCKER-PRAGER PLASTICITY STRUCTURAL MATERIAL   ***
//   *************************************************************


DruckerPragerPlasticitySM :: DruckerPragerPlasticitySM(int n, Domain *d) :
    StructuralMaterial(n, d)
{
    LEMaterial = new IsotropicLinearElasticMaterial(n, d);
    kFactor = 0.;
    yieldTol = 0.;
    newtonIter = 0;
}

DruckerPragerPlasticitySM :: ~DruckerPragerPlasticitySM()
{
    delete LEMaterial;
}

IRResultType
DruckerPragerPlasticitySM :: initializeFrom(InputRecord *ir)
{
    // Required by IR_GIVE_FIELD macro
    IRResultType result;
    // call the corresponding service of structural material
    result = StructuralMaterial :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    // call the corresponding service for the linear elastic material
    result = LEMaterial->initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    // initialize elastic constants
    //eM = LEMaterial->give(Ex,gp);
    //nu = LEMaterial->give(NYxz,gp);
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
        OOFEM_WARNING("Choose hardeningType 1 (linear hardening/softening), 2 (exponential hardening/softening) in input file!");
        return IRRT_BAD_FORMAT;
    }

    yieldTol = 1.e-14;
    IR_GIVE_OPTIONAL_FIELD(ir, yieldTol, _IFT_DruckerPragerPlasticitySM_yieldtol);
    newtonIter = 30;
    IR_GIVE_OPTIONAL_FIELD(ir, newtonIter, _IFT_DruckerPragerPlasticitySM_newtoniter);

    return IRRT_OK;
}

void
DruckerPragerPlasticitySM :: giveRealStressVector_3d(FloatArray &answer,
                                                  GaussPoint *gp,
                                                  const FloatArray &totalStrain,
                                                  TimeStep *tStep)
{
    FloatArray strainVectorR;

    DruckerPragerPlasticitySMStatus *status =
        static_cast< DruckerPragerPlasticitySMStatus * >( this->giveStatus(gp) );

    // Initialize temp variables for this gauss point
    initTempStatus(gp);

    // subtract stress independent part
    // note: eigenStrains (temperature) is not contained in mechanical strain stored in gp
    // therefore it is necessary to subtract always the total eigen strain value
    this->giveStressDependentPartOfStrainVector(strainVectorR, gp, totalStrain,
                                                tStep, VM_Total);

    // perform the local stress return and update the history variables
    performLocalStressReturn(gp, strainVectorR);

    // copy total strain vector to the temp status
    status->letTempStrainVectorBe(totalStrain);

    // give back correct form of stressVector to giveRealStressVector
    answer = status->giveTempStressVector();
}

void
DruckerPragerPlasticitySM :: performLocalStressReturn(GaussPoint *gp,
                                                      const FloatArray &strain)
{
    DruckerPragerPlasticitySMStatus *status =
        static_cast< DruckerPragerPlasticitySMStatus * >( giveStatus(gp) );

    // split total strains in volumetric and deviatoric part
    FloatArray strainDeviator;
    double volumetricStrain;
    // elastic constants
    double eM = LEMaterial->give(Ex, gp);
    double nu = LEMaterial->give(NYxz, gp);
    double gM = eM / ( 2. * ( 1. + nu ) );
    double kM = eM / ( 3. * ( 1. - 2. * nu ) );

    volumetricStrain = computeDeviatoricVolumetricSplit(strainDeviator, strain);

    // compute trial elastic strains
    double volumetricElasticTrialStrain =
        volumetricStrain - status->giveVolumetricPlasticStrain();
    FloatArray plasticStrainDeviator = status->givePlasticStrainDeviator();
    FloatArray elasticStrainDeviator = strainDeviator;
    elasticStrainDeviator.subtract(plasticStrainDeviator);

    // compute trial stresses
    FloatArray stressDeviator;
    double volumetricStress = 3. * kM * volumetricElasticTrialStrain;
    applyDeviatoricElasticStiffness(stressDeviator, elasticStrainDeviator, gM);
    // norm of trial stress deviator
    double trialStressJTwo = computeSecondStressInvariant(stressDeviator);

    // initialize hardening parameter
    double kappa = status->giveKappa();
    double tempKappa = kappa;

    // choose and perform correct stress return and update state flag
    if ( computeYieldValue(volumetricStress, trialStressJTwo, tempKappa, eM) / eM
         > yieldTol ) {
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
    FloatArray stress;
    computeDeviatoricVolumetricSum(stress, stressDeviator, volumetricStress);
    status->letTempStressVectorBe(stress);

    // compute and update plastic strains, volumetric and deviatoric part
    applyDeviatoricElasticCompliance(elasticStrainDeviator, stressDeviator, gM);
    plasticStrainDeviator = strainDeviator;
    plasticStrainDeviator.subtract(elasticStrainDeviator);
    status->letTempPlasticStrainDeviatorBe(plasticStrainDeviator);
    status->letTempVolumetricPlasticStrainBe(volumetricStrain - volumetricStress / 3. / kM);
}

bool
DruckerPragerPlasticitySM :: checkForVertexCase(double eM, double gM, double kM, double trialStressJTwo, double volumetricStress, double tempKappa)
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
DruckerPragerPlasticitySM :: performRegularReturn(double eM, double gM, double kM, double trialStressJTwo, FloatArray &stressDeviator, double &volumetricStress, double &tempKappa)
{
    // delta lambda max controls the maximum plastic flow, see below
    double deltaLambdaMax = sqrt(trialStressJTwo) / gM;

    // declare some constants for faster use
    double volConstant = 3. * kM * alphaPsi;
    double devConstant = sqrt(2.) * gM;
    // yield value prime is derivative of yield value with respect to deltaLambda
    double yieldValuePrimeZero = -9. * alpha * alphaPsi * kM - gM;

    FloatArray flowDir = stressDeviator;
    flowDir.times( 1. / sqrt(2. * trialStressJTwo) );


    // some variables needed for iteration
    double yieldValuePrime;
    FloatArray plasticFlow;

    // initialize Newton iteration
    double tempJTwo;
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

        yieldValuePrime = yieldValuePrimeZero - kFactor *computeYieldStressPrime(tempKappa, eM);
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

        // plasticFlow = flowDir ;
        //plasticFlow.times( devConstant * deltaLambdaIncrement ) ;
        //stressDeviator.subtract( plasticFlow ) ;

        stressDeviator.add(-devConstant * deltaLambdaIncrement, flowDir);
        tempJTwo = computeSecondStressInvariant(stressDeviator);
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
DruckerPragerPlasticitySM :: performVertexReturn(double eM, double gM, double kM, double trialStressJTwo, FloatArray &stressDeviator, double &volumetricStress, double &tempKappa, double volumetricElasticTrialStrain, double kappa)
{
    // declare some constants for faster use
    // yield value prime is derivative of yield value with respect to deltaLambda
    double yieldValuePrimeZero = 3. * alpha;
    // contribution of deviatoric strain to hardening
    double deviatorContribution = trialStressJTwo / 3. / gM / gM;
    // in the vertex case, deviatoric stresses are zero
    stressDeviator.zero();
    volumetricStress = 3. * kM * volumetricElasticTrialStrain;

    // needed for iteration
    double yieldValuePrime;

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

        break;
    case 2:     // exponential hardening
        yieldStress =
            limitYieldStress - ( limitYieldStress - initialYieldStress ) * exp(-kappa / kappaC);
        break;
    default:
        //StructuralMaterial :: OOFEM_ERROR( "Case failed: choose linear hardening/softening (1), exponential hardening/softening (2) in input file.") ;
        OOFEM_ERROR("Case failed: choose linear hardening/softening (1), exponential hardening/softening (2) in input file.");
        return 0.;

        break;
    }

    return yieldStress;
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

        break;
    case 2:             // exponential hardening
        return ( limitYieldStress - initialYieldStress ) / kappaC *exp(-kappa / kappaC);

        break;
    default:
        //StructuralMaterial :: OOFEM_ERROR( "Case failed: choose linear hardening/softening (1), exponential hardening/softening (2) in input file.") ;
        OOFEM_ERROR("Case failed: choose linear hardening/softening (1), exponential hardening/softening (2) in input file.");
        return 0.;

        break;
    }
}

void
DruckerPragerPlasticitySM :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                                           MatResponseMode mode,
                                                           GaussPoint *gp,
                                                           TimeStep *tStep)
{
    switch ( mode ) {
    case ElasticStiffness:
        LEMaterial->give3dMaterialStiffnessMatrix(answer, mode, gp, tStep);
        break;

    case SecantStiffness:
        LEMaterial->give3dMaterialStiffnessMatrix(answer,  mode, gp, tStep);
        break;


    case TangentStiffness:
        switch ( ( static_cast< DruckerPragerPlasticitySMStatus * >( this->giveStatus(gp) ) )
                ->giveTempStateFlag() ) {
        case DruckerPragerPlasticitySMStatus :: DP_Elastic:        // elastic stiffness
        case DruckerPragerPlasticitySMStatus :: DP_Unloading:        // elastic stiffness
            LEMaterial->give3dMaterialStiffnessMatrix(answer, mode, gp, tStep);
            break;
        case DruckerPragerPlasticitySMStatus :: DP_Yielding:
            // elasto-plastic stiffness for regular case
            //printf("\nAssembling regular algorithmic stiffness matrix.") ;
            giveRegAlgorithmicStiffMatrix(answer, mode, gp, tStep);
            break;
        case DruckerPragerPlasticitySMStatus :: DP_Vertex:
            // elasto-plastic stiffness for vertex case
            //printf("\nAssembling vertex case algorithmic stiffness matrix.") ;
            giveVertexAlgorithmicStiffMatrix(answer, mode, gp, tStep);
            break;
        default:
            OOFEM_ERROR("Case did not match.");
            break;
        }

        break;

    default:
        OOFEM_ERROR("Switch failed: Only elastic and tangent stiffness are supported.");
        break;
    }
}

void
DruckerPragerPlasticitySM :: giveRegAlgorithmicStiffMatrix(FloatMatrix &answer,
                                                           MatResponseMode mode,
                                                           GaussPoint *gp,
                                                           TimeStep *tStep)
{
    int i, j;
    DruckerPragerPlasticitySMStatus *status =
        static_cast< DruckerPragerPlasticitySMStatus * >( this->giveStatus(gp) );

    const FloatArray &stressVector = status->giveTempStressVector();
    FloatArray deviatoricStress;
    computeDeviatoricVolumetricSplit(deviatoricStress, stressVector);
    double normOfStress = sqrt( 2. * computeSecondStressInvariant(deviatoricStress) );
    // elastic constants
    double eM = LEMaterial->give(Ex, gp);
    double nu = LEMaterial->give(NYxz, gp);
    double gM = eM / ( 2. * ( 1. + nu ) );
    double kM = eM / ( 3. * ( 1. - 2. * nu ) );



    FloatArray flowDir = deviatoricStress;
    flowDir.times(1. / normOfStress);

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

    FloatMatrix A_Matrix(6, 6);
    A_Matrix.zero();
    for ( i = 1; i < 7; i++ ) {
        A_Matrix.at(i, i) = a_const;
    }

    for ( i = 1; i < 4; i++ ) {
        for ( j = 1; j < 4; j++ ) {
            A_Matrix.at(i, j) += b_const;
        }
    }

    for ( i = 1; i < 4; i++ ) {
        for ( j = 1; j < 4; j++ ) {
            A_Matrix.at(i, j) += c_const * flowDir.at(j);
        }

        for ( j = 4; j < 7; j++ ) {
            A_Matrix.at(i, j) += 2. *c_const *flowDir.at(j);
        }
    }

    for ( i = 1; i < 7; i++ ) {
        for ( j = 1; j < 4; j++ ) {
            A_Matrix.at(i, j) += d_const * flowDir.at(i);
        }
    }

    for ( i = 1; i < 7; i++ ) {
        for ( j = 1; j < 4; j++ ) {
            A_Matrix.at(i, j) += e_const * flowDir.at(i) * flowDir.at(j);
        }

        // account for engineering notation
        for ( j = 4; j < 7; j++ ) {
            A_Matrix.at(i, j) += 2. *e_const *flowDir.at(i) * flowDir.at(j);
        }
    }

    FloatMatrix De;
    LEMaterial->give3dMaterialStiffnessMatrix(De, mode, gp, tStep);

    // answer is A_Matrix^-1 * De
    A_Matrix.solveForRhs(De, answer);
}

void
DruckerPragerPlasticitySM :: giveVertexAlgorithmicStiffMatrix(FloatMatrix &answer,
                                                              MatResponseMode mode,
                                                              GaussPoint *gp,
                                                              TimeStep *tStep)
{
    DruckerPragerPlasticitySMStatus *status =
        static_cast< DruckerPragerPlasticitySMStatus * >( this->giveStatus(gp) );

    double tempKappa = status->giveTempKappa();
    double deltaKappa = tempKappa - status->giveKappa();
    // elastic constants
    double eM = LEMaterial->give(Ex, gp);
    double nu = LEMaterial->give(NYxz, gp);
    //double gM = eM / ( 2. * ( 1. + nu ) );
    double kM = eM / ( 3. * ( 1. - 2. * nu ) );

    if ( deltaKappa <= 0. ) {
        // This case occurs in the first iteration of a step.
        // printf("deltaKappa<=0. for vertex case algorithmic stiffness, i.e. continuum tangent stiffness. Since the continuum tangent stiffness does not exist at the vertex, elastic stiffness is used instead. This will cause the loss of quadratic convergence.\n") ;
        LEMaterial->give3dMaterialStiffnessMatrix(answer, mode, gp, tStep);
    }

    double deltaVolumetricPlasticStrain =
        status->giveTempVolumetricPlasticStrain() - status->giveVolumetricPlasticStrain();
    double HBar = computeYieldStressPrime(tempKappa, eM);

    // compute elastic trial strain deviator of latest temp-state
    FloatArray strainDeviator;
    computeDeviatoricVolumetricSplit(strainDeviator, status->giveTempStrainVector());

    FloatArray elasticStrainDeviator = strainDeviator;
    elasticStrainDeviator.subtract( status->givePlasticStrainDeviator() );

    double a_const =
        kM * HBar / ( HBar * deltaVolumetricPlasticStrain + 9. / 2. * alpha * kM * deltaKappa );

    if ( ( HBar * deltaVolumetricPlasticStrain + 9. / 2. * alpha * kM * deltaKappa ) == 0. ) {
        OOFEM_ERROR("Tangent type is singular, material ID %d\n", this->giveNumber() );
    }
    // compute the algorithmic tangent stiffness

    answer.resize(6, 6);
    answer.zero();
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

    answer.times(a_const);
}

int
DruckerPragerPlasticitySM :: giveIPValue(FloatArray &answer,
                                         GaussPoint *gp,
                                         InternalStateType type,
                                         TimeStep *tStep)
{
    const DruckerPragerPlasticitySMStatus *status =
        static_cast< DruckerPragerPlasticitySMStatus * >( giveStatus(gp) );

    switch ( type ) {
    case IST_PlasticStrainTensor:
        status->givePlasticStrainVector(answer);
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
    return new  DruckerPragerPlasticitySMStatus(1, StructuralMaterial :: giveDomain(), gp);
}


double
DruckerPragerPlasticitySM :: predictRelativeComputationalCost(GaussPoint *gp)
{
    DruckerPragerPlasticitySMStatus *status =
        static_cast< DruckerPragerPlasticitySMStatus * >( giveStatus(gp) );
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
