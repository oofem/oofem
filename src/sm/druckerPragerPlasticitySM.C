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

#include "druckerPragerPlasticitySM.h"

#include "flotarry.h"
#include "flotmtrx.h"
#include "structuralms.h"
#include "gausspnt.h"
#include "intarray.h"
#include "structuralmaterial.h"
#include "isolinearelasticmaterial.h"
#include "structuralcrosssection.h"
#include "datastream.h"
#include "contextioerr.h"
#include "mathfem.h"

namespace oofem {
DruckerPragerPlasticitySMStatus :: DruckerPragerPlasticitySMStatus(int n, Domain *d, GaussPoint *gp) :
    StructuralMaterialStatus(n, d, gp),
    plasticStrainDeviator( gp->giveMaterialMode() ),
    tempPlasticStrainDeviator( gp->giveMaterialMode() )
{
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
DruckerPragerPlasticitySMStatus :: updateYourself(TimeStep *atTime)
{
    // Call corresponding function of the parent class to update variables defined there.
    StructuralMaterialStatus :: updateYourself(atTime);

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
    StrainVector plasticStrainVector( gp->giveMaterialMode() );
    giveFullPlasticStrainVector(plasticStrainVector);

    fprintf(file, "plasticStrains ");
    int n = plasticStrainVector.giveSize();
    for ( int i = 1; i <= n; i++ ) {
        fprintf( file, " % .4e", plasticStrainVector.at(i) );
    }

    fprintf(file, "}\n");

    fprintf(file, "\t\thardening_parameter ");
    // print hardening parameter
    fprintf(file, " % .4e\n", kappa);
}

contextIOResultType
DruckerPragerPlasticitySMStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write raw data
    if ( ( iores = plasticStrainDeviator.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream->write(& volumetricPlasticStrain, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& kappa, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& temp_state_flag, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}


contextIOResultType
DruckerPragerPlasticitySMStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
    if ( ( iores = plasticStrainDeviator.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream->read(& volumetricPlasticStrain, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& kappa, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& temp_state_flag, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}

//   *************************************************************
//   *** CLASS DRUCKER-PRAGER PLASTICITY STRUCTURAL MATERIAL   ***
//   *************************************************************


DruckerPragerPlasticitySM :: DruckerPragerPlasticitySM(int n, Domain *d) :
    StructuralMaterial(n, d),
    stressDeviator(_Unknown)
{
    LEMaterial = new IsotropicLinearElasticMaterial(n, d);
    volumetricElasticTrialStrain = 0.;
    volumetricStress = 0.;
    kappa = 0.;
    tempKappa = 0.;
    trialStressJTwo = 0.;
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
    const char *__proc = "initializeFrom";
    IRResultType result;
    // call the corresponding service of structural material
    StructuralMaterial :: initializeFrom(ir);

    // call the corresponding service for the linear elastic material
    LEMaterial->initializeFrom(ir);

    // initialize elastic constants
    //eM = LEMaterial->give(Ex,gp);
    //nu = LEMaterial->give(NYxz,gp);
    //gM = eM / ( 2. * ( 1. + nu ) );
    //kM = eM / ( 3. * ( 1. - 2. * nu ) );

    // instanciate the variables defined in DruckerPragerPlasticitySM
    IR_GIVE_FIELD(ir, initialYieldStress, IFT_DruckerPragerPlasticitySM_iys, "iys"); // initial yield stress under pure shear
    IR_GIVE_FIELD(ir, alpha, IFT_DruckerPragerPlasticitySM_alpha, "alpha"); // friction coefficient
    IR_GIVE_FIELD(ir, alphaPsi, IFT_DruckerPragerPlasticitySM_alphapsi, "alphapsi"); //dilatancy coefficient
    // this is valid for strain hardening/softening only (not for work hardening/softening)
    kFactor = sqrt(1. / 3. + 2. * alphaPsi * alphaPsi);

    IR_GIVE_FIELD(ir, hardeningType, IFT_DruckerPragerPlasticitySM_ht, "ht");

    switch ( hardeningType ) {
    case 1:     // linear hardening/softening
        IR_GIVE_FIELD(ir, hardeningModulus, IFT_DruckerPragerPlasticitySM_hm, "hm");
        break;
    case 2:     // exponential hardening/softening
        IR_GIVE_FIELD(ir, kappaC, IFT_DruckerPragerPlasticitySM_kc, "kc");
        IR_GIVE_FIELD(ir, limitYieldStress, IFT_DruckerPragerPlasticitySM_lys, "lys");
        break;
    default:
        _error("Choose hardeningType 1 (linear hardening/softening), 2 (exponential hardening/softening) in input file!");
        break;
    }

    yieldTol = 1.e-14;
    IR_GIVE_OPTIONAL_FIELD(ir, yieldTol, IFT_DruckerPragerPlasticitySM_yieldtol, "yieldtol");
    newtonIter = 30;
    IR_GIVE_OPTIONAL_FIELD(ir, newtonIter, IFT_DruckerPragerPlasticitySM_newtoniter, "newtoniter");

    return IRRT_OK;
}

int
DruckerPragerPlasticitySM :: hasMaterialModeCapability(MaterialMode mMode)
{
    if ( ( mMode == _3dMat ) || //mMode == _PlaneStress - this mode needs to be elaborated
        ( mMode == _PlaneStrain ) ||
        ( mMode == _3dRotContinuum ) ) {
        return 1;
    } else {
        return 0;
    }
}

void
DruckerPragerPlasticitySM :: giveRealStressVector(FloatArray &answer,
                                                  MatResponseForm form,
                                                  GaussPoint *gp,
                                                  const FloatArray &totalStrain,
                                                  TimeStep *atTime)
{
    FloatArray strainVectorR;

    if ( stressDeviator.giveStressStrainMode() == _Unknown ) {
        stressDeviator.letStressStrainModeBe( gp->giveMaterialMode() );
    }

    DruckerPragerPlasticitySMStatus *status =
        ( DruckerPragerPlasticitySMStatus * ) ( this->giveStatus(gp) );

    // Initialize temp variables for this gauss point
    initTempStatus(gp);

    // subtract stress independent part
    // note: eigenStrains (temperature) is not contained in mechanical strain stored in gp
    // therefore it is necessary to subtract always the total eigen strain value
    this->giveStressDependentPartOfStrainVector(strainVectorR, gp, totalStrain,
                                                atTime, VM_Total);

    // perform the local stress return and update the history variables
    StrainVector strain( strainVectorR, gp->giveMaterialMode() );
    performLocalStressReturn(gp, strain);

    // copy total strain vector to the temp status
    status->letTempStrainVectorBe(totalStrain);

    // give back correct form of stressVector to giveRealStressVector
    if ( form == ReducedForm ) {
        answer = status->giveTempStressVector();
    } else {
        ( ( StructuralCrossSection * ) ( gp->giveElement()->giveCrossSection() ) )
        ->giveFullCharacteristicVector( answer, gp, status->giveTempStressVector() );
    }
}

void
DruckerPragerPlasticitySM :: performLocalStressReturn(GaussPoint *gp,
                                                      const StrainVector &strain)
{
    DruckerPragerPlasticitySMStatus *status =
        ( DruckerPragerPlasticitySMStatus * ) ( giveStatus(gp) );

    // split total strains in volumetric and deviatoric part
    StrainVector strainDeviator( gp->giveMaterialMode() );
    double volumetricStrain;
    // elastic constants
    double eM = LEMaterial->give(Ex, gp);
    double nu = LEMaterial->give(NYxz, gp);
    double gM = eM / ( 2. * ( 1. + nu ) );
    double kM = eM / ( 3. * ( 1. - 2. * nu ) );

    strain.computeDeviatoricVolumetricSplit(strainDeviator, volumetricStrain);

    // compute trial elastic strains
    volumetricElasticTrialStrain =
        volumetricStrain - status->giveVolumetricPlasticStrain();
    StrainVector plasticStrainDeviator( gp->giveMaterialMode() );
    status->givePlasticStrainDeviator(plasticStrainDeviator);
    StrainVector elasticStrainDeviator = strainDeviator;
    elasticStrainDeviator.subtract(plasticStrainDeviator);

    // compute trial stresses
    volumetricStress = 3. * kM * volumetricElasticTrialStrain;
    elasticStrainDeviator.applyDeviatoricElasticStiffness(stressDeviator, gM);
    // norm of trial stress deviator
    trialStressJTwo = stressDeviator.computeSecondInvariant();

    // initialize hardening parameter
    kappa = status->giveKappa();
    tempKappa = kappa;

    // choose and perform correct stress return and update state flag
    if ( computeYieldValue(volumetricStress, trialStressJTwo, tempKappa, eM) / eM
         > yieldTol ) {
        if ( checkForVertexCase(eM, gM, kM) ) {
            performVertexReturn(eM, gM, kM);
            status->letTempStateFlagBe(DruckerPragerPlasticitySMStatus :: DP_Vertex);
        } else {
            performRegularReturn(eM, gM, kM);
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
    StressVector stress( gp->giveMaterialMode() );
    stressDeviator.computeDeviatoricVolumetricSum(stress, volumetricStress);
    status->letTempStressVectorBe(stress);

    // compute and update plastic strains, volumetric and deviatoric part
    stressDeviator.applyDeviatoricElasticCompliance(elasticStrainDeviator, gM);
    plasticStrainDeviator = strainDeviator;
    plasticStrainDeviator.subtract(elasticStrainDeviator);
    status->letTempPlasticStrainDeviatorBe(plasticStrainDeviator);
    status->letTempVolumetricPlasticStrainBe(volumetricStrain - volumetricStress / 3. / kM);
}

bool
DruckerPragerPlasticitySM :: checkForVertexCase(const double eM, const double gM, const double kM)
{
    // delta lambda max corresponds to the maximum value
    // of the rate of the plastic multiplier for regular plastic flow
    // and allows to distinguish between regular return and vertex case
    const double deltaLambdaMax = sqrt(trialStressJTwo) / gM;

    // vertex case:
    // yield value positive under the assumption of maximum regular plastic flow
    const double volConstant = 3. * kM * alphaPsi;
    const double yieldValue =
        computeYieldValue(volumetricStress - volConstant * deltaLambdaMax,
                          0., tempKappa + kFactor * deltaLambdaMax, eM);
    if ( yieldValue / eM > -yieldTol ) {
        return true;
    }

    // regular case
    return false;
}

void
DruckerPragerPlasticitySM :: performRegularReturn(double eM, double gM, double kM)
{
    // delta lambda max controls the maximum plastic flow, see below
    const double deltaLambdaMax = sqrt(trialStressJTwo) / gM;

    // declare some constants for faster use
    const double volConstant = 3. * kM * alphaPsi;
    const double devConstant = sqrt(2.) * gM;
    // yield value prime is derivative of yield value with respect to deltaLambda
    const double yieldValuePrimeZero = -9. * alpha * alphaPsi * kM - gM;

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
            _error("Newton iteration for deltaLambda (regular stress return) did not converge after newtonIter iterations. You might want to try increasing the optional parameter newtoniter or yieldtol in the material record of your input file.");
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
        tempJTwo = stressDeviator.computeSecondInvariant();
        yieldValue = computeYieldValue(volumetricStress, tempJTwo, tempKappa, eM);
        newtonError = fabs(yieldValue / eM);
        //printf("newtonError = %e\n", newtonError) ;
    }

    OOFEM_LOG_DEBUG("IterationCount in regular return = %d\n", iterationCount) ;

    if ( deltaLambda < 0. ) {
        _error("Fatal error in the Newton iteration for regular stress return. deltaLambda is evaluated as negative, but should always be positive. This is most likely due to a softening law with local snapback, which is physically inadmissible.n");
    }
}

void
DruckerPragerPlasticitySM :: performVertexReturn(double eM, double gM, double kM)
{
    // declare some constants for faster use
    // yield value prime is derivative of yield value with respect to deltaLambda
    const double yieldValuePrimeZero = 3. * alpha;
    // contribution of deviatoric strain to hardening
    const double deviatorContribution = trialStressJTwo / 3. / gM / gM;
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
            _error("Newton iteration for deltaLambda (vertex stress return) did not converge after newtonIter iterations. You might want to try increasing the optional parameter newtoniter or yieldtol in the material record of your input file.");
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
        OOFEM_LOG_DEBUG("NewtonError in iteration %d in vertex return = %e\n", iterationCount, newtonError) ;
    }

    OOFEM_LOG_DEBUG("Done iteration in vertex return, after %d\n", iterationCount);

    if ( deltaKappa < 0. ) {
        _error("Fatal error in the Newton iteration for vertex stress return. deltaKappa is evaluated as negative, but should always be positive. This is most likely due to a softening law with a local snapback, which is physically inadmissible.\n");
    }
}

double
DruckerPragerPlasticitySM :: computeYieldValue(const double volumetricStress,
                                               const double JTwo,
                                               const double kappa,
                                               const double eM) const
{
    return 3. * alpha * volumetricStress + sqrt(JTwo) - computeYieldStressInShear(kappa, eM);
}

double
DruckerPragerPlasticitySM :: computeYieldStressInShear(const double kappa, const double eM) const
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
        //StructuralMaterial :: _error( "Case failed: choose linear hardening/softening (1), exponential hardening/softening (2) in input file.") ;
        _error("Case failed: choose linear hardening/softening (1), exponential hardening/softening (2) in input file.");
        return 0.;

        break;
    }

    return yieldStress;
}

double
DruckerPragerPlasticitySM :: computeYieldStressPrime(const double kappa, const double eM) const
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
        //StructuralMaterial :: _error( "Case failed: choose linear hardening/softening (1), exponential hardening/softening (2) in input file.") ;
        _error("Case failed: choose linear hardening/softening (1), exponential hardening/softening (2) in input file.");
        return 0.;

        break;
    }
}

void
DruckerPragerPlasticitySM :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                                           MatResponseForm form,
                                                           MatResponseMode mode,
                                                           GaussPoint *gp,
                                                           TimeStep *atTime)
{
    switch ( mode ) {
    case ElasticStiffness:
        LEMaterial->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
        break;

    case TangentStiffness:
        switch ( ( ( DruckerPragerPlasticitySMStatus * ) ( this->giveStatus(gp) ) )
                ->giveTempStateFlag() ) {
        case DruckerPragerPlasticitySMStatus :: DP_Elastic:        // elastic stiffness
        case DruckerPragerPlasticitySMStatus :: DP_Unloading:        // elastic stiffness
            LEMaterial->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
            break;
        case DruckerPragerPlasticitySMStatus :: DP_Yielding:
            // elasto-plastic stiffness for regular case
            //printf("\nAssembling regular algorithmic stiffness matrix.") ;
            giveRegAlgorithmicStiffMatrix(answer, form, mode, gp, atTime);
            break;
        case DruckerPragerPlasticitySMStatus :: DP_Vertex:
            // elasto-plastic stiffness for vertex case
            //printf("\nAssembling vertex case algorithmic stiffness matrix.") ;
            giveVertexAlgorithmicStiffMatrix(answer, form, mode, gp, atTime);
            break;
        default:
            _error("Case did not match.\n");
            break;
        }

        break;

    default:
        _error("Switch failed: Only elastic and tangent stiffness are supported.\n");
        break;
    }
}

void
DruckerPragerPlasticitySM :: giveRegAlgorithmicStiffMatrix(FloatMatrix &answer,
                                                           MatResponseForm form,
                                                           MatResponseMode mode,
                                                           GaussPoint *gp,
                                                           TimeStep *atTime)
{
    int i, j;
    DruckerPragerPlasticitySMStatus *status =
        ( DruckerPragerPlasticitySMStatus * ) ( this->giveStatus(gp) );
    StructuralCrossSection *crossSection =
        ( StructuralCrossSection * ) ( gp->giveElement()->giveCrossSection() );

    const FloatArray stressVector = status->giveTempStressVector();
    FloatArray fullStressVector;
    crossSection->giveFullCharacteristicVector(fullStressVector, gp, stressVector);
    const StressVector stress(fullStressVector, _3dMat);
    StressVector deviatoricStress(_3dMat);
    double volumetricStress;
    stress.computeDeviatoricVolumetricSplit(deviatoricStress, volumetricStress);
    const double normOfStress = sqrt( 2. * deviatoricStress.computeSecondInvariant() );
    // elastic constants
    const double eM = LEMaterial->give(Ex, gp);
    const double nu = LEMaterial->give(NYxz, gp);
    const double gM = eM / ( 2. * ( 1. + nu ) );
    const double kM = eM / ( 3. * ( 1. - 2. * nu ) );



    FloatArray flowDir = deviatoricStress;
    flowDir.times(1. / normOfStress);

    const double kappa = status->giveKappa();
    const double tempKappa = status->giveTempKappa();
    const double deltaKappa = tempKappa - kappa;
    const double deltaLambdaStar = sqrt(2.) * gM * deltaKappa / kFactor / normOfStress;
    double hStar = kFactor * computeYieldStressPrime(tempKappa, eM);

    //exclude division by zero
    if(hStar == 0.){
        OOFEM_ERROR("DruckerPragerPlasticitySM :: computeYieldStressPrime is zero. This happens mainly due to excessive softening.");
    }

    const double a_const = 1. + deltaLambdaStar;
    const double b_const = 3. * alpha * alphaPsi * kM / hStar - deltaLambdaStar / 3.;
    const double c_const = 3. * sqrt(2.) * alphaPsi * kM / 2. / hStar;
    const double d_const = sqrt(2.) * alpha * gM / hStar;
    const double e_const = gM / hStar - deltaLambdaStar;

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
    LEMaterial->giveCharacteristicMatrix(De, form, mode, gp, atTime);

    // answer is A_Matrix^-1 * De
    A_Matrix.solveForRhs(De, answer);
}

void
DruckerPragerPlasticitySM :: giveVertexAlgorithmicStiffMatrix(FloatMatrix &answer,
                                                              MatResponseForm form,
                                                              MatResponseMode mode,
                                                              GaussPoint *gp,
                                                              TimeStep *atTime)
{
    int i, j;
    DruckerPragerPlasticitySMStatus *status =
        ( DruckerPragerPlasticitySMStatus * ) ( this->giveStatus(gp) );
    StructuralCrossSection *crossSection =
        ( StructuralCrossSection * ) ( gp->giveElement()->giveCrossSection() );

    const double tempKappa = status->giveTempKappa();
    const double deltaKappa = tempKappa - status->giveKappa();
    // elastic constants
    const double eM = LEMaterial->give(Ex, gp);
    const double nu = LEMaterial->give(NYxz, gp);
    //const double gM = eM / ( 2. * ( 1. + nu ) );
    const double kM = eM / ( 3. * ( 1. - 2. * nu ) );

    if ( deltaKappa <= 0. ) {
        // This case occurs in the first iteration of a step.
        // printf("deltaKappa<=0. for vertex case algorithmic stiffness, i.e. continuum tangent stiffness. Since the continuum tangent stiffness does not exist at the vertex, elastic stiffness is used instead. This will cause the loss of quadratic convergence.\n") ;
        LEMaterial->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
    }

    const double deltaVolumetricPlasticStrain =
        status->giveTempVolumetricPlasticStrain() - status->giveVolumetricPlasticStrain();
    const double HBar = computeYieldStressPrime(tempKappa, eM);

    // compute elastic trial strain deviator of latest temp-state
    FloatArray fullStrainVector;
    crossSection->giveFullCharacteristicVector( fullStrainVector, gp,
                                               status->giveTempStrainVector() );
    StrainVector strain(fullStrainVector, _3dMat);
    StrainVector strainDeviator(_3dMat);
    double volumetricStrain;
    strain.computeDeviatoricVolumetricSplit(strainDeviator, volumetricStrain);

    StrainVector plasticStrainDeviator(_3dMat);
    status->givePlasticStrainDeviator(plasticStrainDeviator);
    StrainVector elasticStrainDeviator = strainDeviator;
    elasticStrainDeviator.subtract(plasticStrainDeviator);

    const double a_const =
        kM * HBar / ( HBar * deltaVolumetricPlasticStrain + 9. / 2. * alpha * kM * deltaKappa );

    if ( ( HBar * deltaVolumetricPlasticStrain + 9. / 2. * alpha * kM * deltaKappa ) == 0.){
        OOFEM_ERROR2("DruckerPragerPlasticitySM :: giveVertexAlgorithmicStiffMatrix of tangent type is singular, material ID %d\n", this->giveNumber());
    }
    // compute the algorithmic tangent stiffness

    answer.resize(6, 6);
    answer.zero();
    for ( i = 1; i < 4; i++ ) {
        for ( j = 1; j < 4; j++ ) {
            answer.at(i, j) = deltaVolumetricPlasticStrain;
        }
    }

    for ( i = 1; i < 4; i++ ) {
        for ( j = 1; j < 4; j++ ) {
            answer.at(i, j) += elasticStrainDeviator.at(j);
        }

        for ( j = 4; j < 7; j++ ) {
            answer.at(i, j) += .5 * elasticStrainDeviator.at(j);
        }
    }

    answer.times(a_const);
}

int
DruckerPragerPlasticitySM :: giveIPValue(FloatArray &answer,
                                         GaussPoint *gp,
                                         InternalStateType type,
                                         TimeStep *atTime)
{
    const DruckerPragerPlasticitySMStatus *status =
        ( DruckerPragerPlasticitySMStatus * ) giveStatus(gp);
    StrainVector plasticStrainVector(_3dMat);

    switch ( type ) {
    case IST_PlasticStrainTensor:
        status->giveFullPlasticStrainVector(plasticStrainVector);
        answer = plasticStrainVector;
        if ( answer.giveSize() == 0 ) {
            answer.resize(6);
            answer.zero();
        }

        return 1;

    case IST_DamageTensor:
        answer.resize(1);
        answer.zero();
        answer.at(1) = status->giveKappa();
        return 1;

    default:
        return StructuralMaterial :: giveIPValue(answer, gp, type, atTime);

    }

    return 0;
}

int
DruckerPragerPlasticitySM :: giveIPValueSize(InternalStateType type,
                                             GaussPoint *gp)
{
    switch ( type ) {
    case IST_PlasticStrainTensor:
        return 6;

    case IST_DamageTensor:
        return 1;

    default:
        return StructuralMaterial :: giveIPValueSize(type, gp);

    }
}

int
DruckerPragerPlasticitySM :: giveIntVarCompFullIndx(IntArray &answer,
                                                    InternalStateType type,
                                                    MaterialMode mmode)
{
    switch ( type ) {
    case IST_PlasticStrainTensor:
        answer.resize(6);
        answer.zero();
        answer.at(1) = 1;
        answer.at(2) = 2;
        answer.at(3) = 3;
        answer.at(4) = 4;
        answer.at(5) = 5;
        answer.at(6) = 6;
        return 1;

    case IST_DamageTensor:
        answer.resize(1);
        answer.zero();
        answer.at(1) = 1;
        return 1;

    default:
        return StructuralMaterial :: giveIntVarCompFullIndx(answer, type, mmode);

    }
}

InternalStateValueType
DruckerPragerPlasticitySM :: giveIPValueType(InternalStateType type)
{
    switch ( type ) {
    case IST_PlasticStrainTensor:      // plastic strain tensor
        return ISVT_TENSOR_S3;

    case IST_DamageTensor:     // damage tensor used for internal variables
        return ISVT_TENSOR_G;

    default:
        return StructuralMaterial :: giveIPValueType(type);

    }
}

MaterialStatus *
DruckerPragerPlasticitySM :: CreateStatus(GaussPoint *gp) const
{
    DruckerPragerPlasticitySMStatus *status =
        new  DruckerPragerPlasticitySMStatus(1, StructuralMaterial :: giveDomain(), gp);
    return status;
}


#ifdef __PARALLEL_MODE
double
DruckerPragerPlasticitySM :: predictRelativeComputationalCost(GaussPoint *gp)
{
    //
    DruckerPragerPlasticitySMStatus *status =
        ( DruckerPragerPlasticitySMStatus * ) ( giveStatus(gp) );
    const int state_flag = status->giveStateFlag();

    if ( ( state_flag == DruckerPragerPlasticitySMStatus :: DP_Vertex ) ||
        ( state_flag == DruckerPragerPlasticitySMStatus :: DP_Yielding ) ) {
        return 20.;
    } else {
        return 1.0;
    }
}

#endif
} // end namespace oofem

#undef PRINTFDP