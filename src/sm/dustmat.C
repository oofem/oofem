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

#include "dustmat.h"

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
DustMaterialStatus :: DustMaterialStatus(int n, Domain *d, GaussPoint *gp) :
    StructuralMaterialStatus(n, d, gp),
    plasticStrain( gp->giveMaterialMode() ),
    tempPlasticStrain( gp->giveMaterialMode() )
{
    q = ( ( DustMaterial * ) gp->giveMaterial() )->giveQ0();
}

DustMaterialStatus :: ~DustMaterialStatus()
{ }

void
DustMaterialStatus :: initTempStatus()
{
    // Call the function of the parent class to initialize the variables defined there.
    StructuralMaterialStatus :: initTempStatus();
    tempPlasticStrain = plasticStrain;
    tempQ = q;
    tempStateFlag = stateFlag;
}

void
DustMaterialStatus :: updateYourself(TimeStep *atTime)
{
    // Call the corresponding function of the parent class to update variables defined there.
    StructuralMaterialStatus :: updateYourself(atTime);
    plasticStrain = tempPlasticStrain;
    q = tempQ;
    stateFlag = tempStateFlag;
}

void
DustMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    // Call the corresponding function of the parent class to print variables defined there.
    StructuralMaterialStatus :: printOutputAt(file, tStep);

    fprintf(file, "\tstatus { ");

    // print status flag
    switch ( stateFlag ) {
    case DustMaterialStatus :: DM_Elastic:
        fprintf(file, " Elastic");
        break;
    case DustMaterialStatus :: DM_Yielding1:
        fprintf(file, " Yielding1");
        break;
    case DustMaterialStatus :: DM_Yielding2:
        fprintf(file, " Yielding2");
        break;
    case DustMaterialStatus :: DM_Yielding3:
        fprintf(file, " Yielding3");
        break;
    case DustMaterialStatus :: DM_Unloading:
        fprintf(file, " Unloading");
        break;
    }

    // print plastic strain vector
    StrainVector plasticStrain( gp->giveMaterialMode() );
    givePlasticStrain(plasticStrain);

    fprintf(file, ", plasticStrains ");
    int n = plasticStrain.giveSize();
    for ( int i = 1; i <= n; i++ ) {
        fprintf( file, " % .4e", plasticStrain.at(i) );
    }

    fprintf(file, ", q  % .4e", q);

    fprintf(file, "}\n");
}

contextIOResultType
DustMaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


contextIOResultType
DustMaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}



//   *************************************************************
//   *** CLASS DUST MATERIAL   ***
//   *************************************************************


DustMaterial :: DustMaterial(int n, Domain *d) : StructuralMaterial(n, d)
{
    LEMaterial = new IsotropicLinearElasticMaterial(n, d);
}

DustMaterial :: ~DustMaterial()
{
    delete LEMaterial;
}

IRResultType
DustMaterial :: initializeFrom(InputRecord *ir)
{
    // Required by IR_GIVE_FIELD macro
    const char *__proc = "initializeFrom";
    IRResultType result;
    // call the corresponding service of structural material
    StructuralMaterial :: initializeFrom(ir);

    // call the corresponding service for the linear elastic material
    this->LEMaterial->initializeFrom(ir);

    // instanciate the variables defined in DustMaterial
    ft = 3e6;
    x0 = 150e6;
    rEll = .5;
    beta = .008e-6;
    theta = .32;
    alpha = 44e6;
    lambda = 34e6;
    wHard = .1;
    dHard = .0003e-6;
    mStiff = 0.;
    newtonTol = 1e-8;
    newtonIter = 200;

    IR_GIVE_OPTIONAL_FIELD(ir, alpha, IFT_DustMaterial_alpha, "alpha");
    IR_GIVE_OPTIONAL_FIELD(ir, beta, IFT_DustMaterial_beta, "beta");
    IR_GIVE_OPTIONAL_FIELD(ir, lambda, IFT_DustMaterial_lambda, "lambda");
    IR_GIVE_OPTIONAL_FIELD(ir, theta, IFT_DustMaterial_theta, "theta");
    IR_GIVE_OPTIONAL_FIELD(ir, ft, IFT_DustMaterial_ft, "ft");
    IR_GIVE_OPTIONAL_FIELD(ir, rEll, IFT_DustMaterial_rEll, "rell");
    IR_GIVE_OPTIONAL_FIELD(ir, x0, IFT_DustMaterial_x0, "x0");
    IR_GIVE_OPTIONAL_FIELD(ir, wHard, IFT_DustMaterial_wHard, "whard");
    IR_GIVE_OPTIONAL_FIELD(ir, dHard, IFT_DustMaterial_dHard, "dhard");
    IR_GIVE_OPTIONAL_FIELD(ir, mStiff, IFT_DustMaterial_mStiff, "mstiff");
    IR_GIVE_OPTIONAL_FIELD(ir, mStiff, IFT_DustMaterial_newtonTol, "newtontol");
    IR_GIVE_OPTIONAL_FIELD(ir, newtonIter, IFT_DustMaterial_newtonIter, "newtoniter");

    // check parameters admissibility
    if ( ft < 0 ) {
        OOFEM_ERROR("parameter 'ft' must be positive")
    }

    ;
    if ( x0 < 0 ) {
        OOFEM_ERROR("parameter 'x0' must be positive")
    }

    ;
    if ( rEll < 0 ) {
        OOFEM_ERROR("parameter 'rEll' must be positive")
    }

    ;
    if ( theta < 0 ) {
        OOFEM_ERROR("parameter 'theta' must be positive")
    }

    ;
    if ( beta < 0 ) {
        OOFEM_ERROR("parameter 'beta' must be positive")
    }

    ;
    if ( lambda < 0 ) {
        OOFEM_ERROR("parameter 'lambda' must be positive")
    }

    ;
    if ( alpha < lambda ) {
        OOFEM_ERROR("parameter 'alpha' must be greater than parameter 'lambda'")
    }

    ;
    x0 = -x0; // compressive strength is negative, although on input it is a positive number

    hardeningType = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, hardeningType, IFT_DustMaterial_hardeningType, "ht");

    q0 = x0;
    solveQ0(q0);

    return IRRT_OK;
}

int
DustMaterial :: hasMaterialModeCapability(MaterialMode mMode)
{
    if ( ( mMode == _3dMat ) ||
        ( mMode == _PlaneStrain ) ||
        ( mMode == _3dRotContinuum ) ) {
        return 1;
    } else {
        return 0;
    }
}

void
DustMaterial :: giveRealStressVector(FloatArray &answer,
                                     MatResponseForm form,
                                     GaussPoint *gp,
                                     const FloatArray &totalStrain,
                                     TimeStep *atTime)
{
    FloatArray strainVectorR;

    DustMaterialStatus *status = ( DustMaterialStatus * )( this->giveStatus(gp) );

    // Initialize temp variables for this Gauss point
    this->initTempStatus(gp);

    // subtract stress-independent part of strain
    this->giveStressDependentPartOfStrainVector(strainVectorR, gp, totalStrain, atTime, VM_Total);

    // perform the local stress return and update the history variables
    StrainVector strain( strainVectorR, gp->giveMaterialMode() );
    performStressReturn(gp, strain);

    // copy total strain vector to the temp status
    status->letTempStrainVectorBe(totalStrain);

    // pass the correct form of stressVector to giveRealStressVector
    if ( form == ReducedForm ) {
        answer = status->giveTempStressVector();
    } else {
        ( ( StructuralCrossSection * )( gp->giveElement()->giveCrossSection() ) )
        ->giveFullCharacteristicVector( answer, gp, status->giveTempStressVector() );
    }
}

void
DustMaterial :: performStressReturn(GaussPoint *gp, StrainVector strain)
{
    DustMaterialStatus *status = ( DustMaterialStatus * ) giveStatus(gp);
    MaterialMode mode = gp->giveMaterialMode();

    // compute total strain components
    StrainVector strainDeviator(mode);
    double volumetricStrain;
    strain.computeDeviatoricVolumetricSplit(strainDeviator, volumetricStrain);

    // compute trial elastic strains
    StrainVector plasticStrain(mode);
    double volumetricPlasticStrain;
    StrainVector plasticStrainDeviator(mode);
    status->givePlasticStrain(plasticStrain);
    plasticStrain.computeDeviatoricVolumetricSplit(plasticStrainDeviator, volumetricPlasticStrain);
    double volumetricElasticStrain = volumetricStrain - volumetricPlasticStrain;
    StrainVector elasticStrainDeviator = strainDeviator;
    elasticStrainDeviator.subtract(plasticStrainDeviator);

    // compute trial stresses
    double bulkModulus, shearModulus;
    computeAndSetBulkAndShearModuli(bulkModulus, shearModulus, gp);
    double volumetricStress = 3. * bulkModulus * volumetricElasticStrain;
    StressVector stressDeviator(mode);
    elasticStrainDeviator.applyDeviatoricElasticStiffness(stressDeviator, shearModulus);

    // norm of trial stress deviator
    double rho = stressDeviator.computeSecondCoordinate();
    double i1 = 3 * volumetricStress;
    double f1, f2, f3, q, tempQ;
    q = tempQ = status->giveQ();
    f1 = yieldFunction1(rho, i1);
    f2 = yieldFunction2(rho, i1, q);
    f3 = yieldFunction3(i1);

    // actual stress return
    double lambda = 0.;
    StrainVector m(mode);
    m.zero();
    double feft = functionFe(ft);
    double auxModulus = 2 * shearModulus / ( 9 * bulkModulus );
    double temp = feft - auxModulus * ( i1 - ft ) / functionFeDI1(ft);

    if ( f1 > 0 && i1 >= q && rho >= temp ) { // yield function 1
        status->letTempStateFlagBe(DustMaterialStatus :: DM_Yielding1);
        performF1return(i1, rho, gp);
        tempQ = status->giveTempQ();
        lambda = ( rho - functionFe( functionI1(q, tempQ, i1, bulkModulus) ) ) / ( 2 * shearModulus );
        computePlastStrainDirM1(m, stressDeviator, rho, i1, q);
    } else if ( f2 > 0 && i1 < q ) { // yield function 2
        status->letTempStateFlagBe(DustMaterialStatus :: DM_Yielding2);
        performF2return(i1, rho, gp);
        tempQ = status->giveTempQ();
        lambda = computeDeltaGamma2(tempQ, q, i1, bulkModulus);
        computePlastStrainDirM2(m, stressDeviator, rho, i1, tempQ);
    } else if ( f3 > 0 && rho < temp ) { // yield function 3
        status->letTempStateFlagBe(DustMaterialStatus :: DM_Yielding3);
        double fFeFt = functionFe(ft);
        double fFeDqFt = functionFeDI1(ft);
        double b = fFeFt - 2 * shearModulus / ( 9 * bulkModulus ) * ( i1 - ft ) - ( 1 + fFeDqFt ) * rho;
        double c = -fFeFt / ( 9 * bulkModulus ) * ( i1 - ft );
        lambda = 1 / ( 4 * shearModulus ) * ( -b + sqrt(b * b - 8 * shearModulus * c) );
        double deltaVolumetricPlasticStrain = 3 * ( 1 - ( 1 + fFeDqFt ) * rho / fFeFt ) * lambda;
        if ( volumetricPlasticStrain + deltaVolumetricPlasticStrain < -.5 * wHard ) {
            deltaVolumetricPlasticStrain = -.5 * wHard - volumetricPlasticStrain;
        }

        if ( q <= 0.0 ) {
            computeQFromPlastVolEps(tempQ, q, deltaVolumetricPlasticStrain);
        }

        status->letTempQBe(tempQ);
        computePlastStrainDirM3(m, stressDeviator, rho, i1, tempQ);
    } else { // elastic case
        int stateFlag = status->giveStateFlag();
        if ( ( stateFlag == DustMaterialStatus :: DM_Unloading ) || ( stateFlag == DustMaterialStatus :: DM_Elastic ) ) {
            status->letTempStateFlagBe(DustMaterialStatus :: DM_Elastic);
        } else {
            status->letTempStateFlagBe(DustMaterialStatus :: DM_Unloading);
        }
    }

    if ( lambda < 0 ) {
        OOFEM_ERROR("TODO");
    }

    // compute correct stress
    m.times(lambda);
    plasticStrain.add(m);
    double mVol;
    StrainVector mDeviator(mode);
    m.computeDeviatoricVolumetricSplit(mDeviator, mVol);
    i1 -= 3 * bulkModulus * mVol;
    volumetricStress = i1 / 3.;
    mDeviator.times(-2 * shearModulus);
    stressDeviator.add(mDeviator);

    // compute full stresses from deviatoric and volumetric part and store them
    StressVector stress(mode);
    stressDeviator.computeDeviatoricVolumetricSum(stress, volumetricStress);
    status->letTempStressVectorBe(stress);

    // compute and update plastic strain and q
    status->letTempPlasticStrainBe(plasticStrain);
}

void
DustMaterial :: performF1return(double i1, double rho, GaussPoint *gp)
{
    DustMaterialStatus *status = ( DustMaterialStatus * ) giveStatus(gp);
    double bulkModulus = status->giveBulkModulus();
    double shearModulus = status->giveShearModulus();
    double q = status->giveQ();
    double tempQ = status->giveTempQ();
    double fx, dfx;
    int i;
    double m = 9 * bulkModulus / ( 2 * shearModulus );
    double vfI1, vfI1DQ, a, b, c, d, da, db, dc;
    int positiveFlag = 0;

    for ( i = 0; i < newtonIter; i++ ) {
        vfI1 = functionI1(q, tempQ, i1, bulkModulus);
        vfI1DQ = functionI1DQ(tempQ, bulkModulus);
        a = ( vfI1 - tempQ ) / ( ft - tempQ );
        b = functionFeDI1(vfI1);
        c = rho - functionFe(vfI1);
        da = ( ( vfI1DQ - 1 ) * ( ft - tempQ ) + ( vfI1 - tempQ ) ) / ( ( ft - tempQ ) * ( ft - tempQ ) );
        db = functionFeDI1DI1(vfI1) * vfI1DQ;
        dc = -functionFeDI1(vfI1) * vfI1DQ;
        d = da * b * c + a * db * c + a * b * dc;
        fx  = -3 *bulkModulus *functionH(q, tempQ) - m * a * b * c;
        dfx = -3 *bulkModulus *functionHDQ(tempQ) - m * d;
        tempQ -= fx / dfx;
        if ( tempQ >= 0 ) {
            if ( positiveFlag >= 1 ) {
                status->letTempQBe(0.0);
                return;
            }

            tempQ = 0;
            positiveFlag += 1;
        }

        if ( fabs(fx / dfx / tempQ) < newtonTol ) {
            status->letTempQBe(tempQ);
            return;
        }
    }

    OOFEM_LOG_DEBUG("  i1 %e rho %e  bulkM %e  shearM %e\n", i1, rho, bulkModulus, shearModulus);
    OOFEM_ERROR("performF1return: Newton's method did not converge\n");
}

void
DustMaterial :: performF2return(double i1, double rho, GaussPoint *gp)
{
    DustMaterialStatus *status = ( DustMaterialStatus * ) giveStatus(gp);
    double bulkModulus = status->giveBulkModulus();
    double shearModulus = status->giveShearModulus();
    double q = status->giveQ();
    double qRight = q;
    double qLeft = q;
    double tempQ = .5 * ( qLeft + qRight );
    double fqRight, fqLeft, fq;
    int i;
    double fx, dfx;
    for ( i = 0; i < newtonIter; i++ ) {
        fx = i1 - 3 *bulkModulus *functionH(q, qLeft) - qLeft;
        dfx =   -3 *bulkModulus *functionHDQ(qLeft) - 1;
        qLeft -= fx / dfx;
        if (  fabs(fx / dfx / q0) < newtonTol ) {
            break;
        }
    }

    for ( i = 0; i < newtonIter; i++ ) {
        fq = fTempR2(tempQ, q, i1, rho, bulkModulus, shearModulus);
        if ( fabs( ( qRight - qLeft ) / qRight ) < newtonTol ) {
            status->letTempQBe(tempQ);
            return;
        }

        fqRight = fTempR2(qRight, q, i1, rho, bulkModulus, shearModulus);
        fqLeft = fTempR2(qLeft, q, i1, rho, bulkModulus, shearModulus);
        if ( fq > 0 ) {
            qRight = tempQ;
        } else {
            qLeft = tempQ;
        }

        tempQ = .5 * ( qLeft + qRight );
    }

    OOFEM_ERROR("performF2return: bisection method did not converge\n");
}

void
DustMaterial :: computeQFromPlastVolEps(double &answer, double q, double deltaVolumetricPlasticStrain)
{
    if ( q >= 0. ) {
        answer = 0.;
        return;
    }

    double fx, dfx;
    int i;
    for ( i = 0; i <= newtonIter; i++ ) {
        fx = functionH(q, answer) - deltaVolumetricPlasticStrain;
        dfx = functionHDQ(answer);
        answer -= fx / dfx;
        if (  fabs(fx / dfx / answer) < newtonTol ) {
            if ( answer > 0 ) {
                answer = 0.;
            }

            return;
        }
    }

    OOFEM_LOG_DEBUG("  dVolEpsPl: %e\n", deltaVolumetricPlasticStrain);
    OOFEM_ERROR("computeQFromPlastVolEps: Newton's method did not converge\n");
}

void
DustMaterial :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                              MatResponseForm form,
                                              MatResponseMode mode,
                                              GaussPoint *gp,
                                              TimeStep *atTime)
{
    DustMaterialStatus *status = ( DustMaterialStatus * ) giveStatus(gp);
    double ym0 = LEMaterial->giveYoungsModulus();
    double ym = status->giveYoungsModulus();
    double coeff = status->giveVolumetricPlasticStrain() < 0 ? ym / ym0 : 1.0;
    if ( mode == ElasticStiffness ) {
        LEMaterial->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
        answer.times(coeff);
    } else if ( mode == SecantStiffness || mode == TangentStiffness ) {
        LEMaterial->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
        answer.times(coeff);
    } else {
        _error("Unsupported MatResponseMode\n");
    }
}

int
DustMaterial :: setIPValue(const FloatArray value, GaussPoint *gp, InternalStateType type)
{
    DustMaterialStatus *status = ( DustMaterialStatus * ) giveStatus(gp);
    if ( type == IST_PlasticStrainTensor ) {
        StrainVector plasticStrain = StrainVector( value, gp->giveMaterialMode() );
        status->letPlasticStrainBe(plasticStrain);
        return 1;
    } else if ( type == IST_StressCapPos ) {
        status->letQBe( value.at(1) );
        return 1;
    } else {
        return StructuralMaterial :: setIPValue(value, gp, type);
    }
}

int
DustMaterial :: giveIPValue(FloatArray &answer,
                            GaussPoint *gp,
                            InternalStateType type,
                            TimeStep *atTime)
{
    const DustMaterialStatus *status = ( DustMaterialStatus * ) giveStatus(gp);
    if ( type == IST_PlasticStrainTensor ||
         type == IST_VolumetricPlasticStrain ||
         type == IST_PrincipalPlasticStrainTensor ) {
        StrainVector plasticStrain( gp->giveMaterialMode() );
        status->givePlasticStrain(plasticStrain);
        if ( type == IST_PlasticStrainTensor ) {
            answer = plasticStrain;
            return 1;
        }

        if ( type == IST_PrincipalPlasticStrainTensor ) {
            plasticStrain.computePrincipalValues(answer);
            return 1;
        }

        StrainVector plasticStrainDeviator( gp->giveMaterialMode() );
        double volumetricPlasticStrain;
        plasticStrain.computeDeviatoricVolumetricSplit(plasticStrainDeviator, volumetricPlasticStrain);
        if ( type == IST_VolumetricPlasticStrain ) {
            answer.resize(1);
            answer.at(1) = volumetricPlasticStrain;
            return 1;
        }
    } else if ( type == IST_StressCapPos ) {
        answer.resize(1);
        answer.at(1) = status->giveQ();
        return 1;
    } else {
        return StructuralMaterial :: giveIPValue(answer, gp, type, atTime);
    }

    return 0;
}

int
DustMaterial :: giveIPValueSize(InternalStateType type,
                                GaussPoint *gp)
{
    if ( type == IST_PlasticStrainTensor ) {
        return this->giveSizeOfReducedStressStrainVector( gp->giveMaterialMode() );
    } else if ( type == IST_PrincipalPlasticStrainTensor ) {
        return 3;
    } else if ( type == IST_VolumetricPlasticStrain || type == IST_StressCapPos ) {
        return 1;
    } else {
        return StructuralMaterial :: giveIPValueSize(type, gp);
    }
}

int
DustMaterial :: giveIntVarCompFullIndx(IntArray &answer,
                                       InternalStateType type,
                                       MaterialMode mmode)
{
    switch ( type ) {
    default:
        return StructuralMaterial :: giveIntVarCompFullIndx(answer, type, mmode);
    }
}

InternalStateValueType
DustMaterial :: giveIPValueType(InternalStateType type)
{
    if ( type == IST_PlasticStrainTensor ) {
        return ISVT_TENSOR_S3E;
    } else if ( type == IST_PrincipalPlasticStrainTensor ) {
        return ISVT_VECTOR;
    } else if ( ( type == IST_StressCapPos ) || ( type == IST_VolumetricPlasticStrain ) ) {
        return ISVT_SCALAR;
    } else {
        return StructuralMaterial :: giveIPValueType(type);
    }
}

MaterialStatus *
DustMaterial :: CreateStatus(GaussPoint *gp) const
{
    DustMaterialStatus *status =
        new  DustMaterialStatus(1, StructuralMaterial :: giveDomain(), gp);
    return status;
}

double
DustMaterial :: functionFe(double i1)
{
    return alpha - lambda *exp(beta *i1) - theta * i1;
}

double
DustMaterial :: functionFeDI1(double i1)
{
    return -lambda *beta *exp(beta *i1) - theta;
}

double
DustMaterial :: functionFeDI1DI1(double i1)
{
    return -lambda *beta *beta *exp(beta *i1);
}

double
DustMaterial :: functionFc(double rho, double i1, double q)
{
    return sqrt( rho * rho + 1 / rEll / rEll * ( q - i1 ) * ( q - i1 ) );
}

double
DustMaterial :: yieldFunction1(double rho, double i1) {
    return rho - functionFe(i1);
}

double
DustMaterial :: yieldFunction2(double rho, double i1, double q)
{
    return functionFc(rho, i1, q) - functionFe(q);
}

double
DustMaterial :: yieldFunction3(double i1)
{
    return i1 - ft;
}

double
DustMaterial :: functionX(double q)
{
    return q - rEll * ( alpha - lambda * exp(beta * q) - theta * q );
}

double
DustMaterial :: functionXDQ(double q)
{
    return 1 - rEll * ( -lambda * beta * exp(beta * q) - theta );
}

void
DustMaterial :: solveQ0(double &answer)
{
    double fx, dfx;
    int i;
    for ( i = 0; i < newtonIter; i++ ) {
        fx = -x0 + answer - rEll * ( alpha -      lambda * exp(beta * answer) - theta * answer );
        dfx =       1 - rEll * (      -beta * lambda * exp(beta * answer) - theta );
        answer -= fx / dfx;
        if (  fabs(fx / dfx / answer) < newtonTol ) {
            if ( answer >= 0 ) {
                OOFEM_ERROR("internal parameter q has to be negative\n");
            }

            return;
        }
    }

    OOFEM_ERROR("solveQ0: Newton's method did not converge\n");
}

void
DustMaterial :: computeAndSetBulkAndShearModuli(double &bulkModulus, double &shearModulus, GaussPoint *gp)
{
    double ym = LEMaterial->giveYoungsModulus();
    double nu = LEMaterial->givePoissonsRatio();
    StrainVector plasticStrain( gp->giveMaterialMode() );
    DustMaterialStatus *status =  ( DustMaterialStatus * ) giveStatus(gp);
    status->givePlasticStrain(plasticStrain);
    double volumetricPlasticStrain = plasticStrain.computeVolumetricPart();
    if ( volumetricPlasticStrain < 0. ) {
        ym -= mStiff * volumetricPlasticStrain;
    }

    bulkModulus = IsotropicLinearElasticMaterial :: computeBulkModulusFromYoungAndPoisson(ym, nu);
    shearModulus = IsotropicLinearElasticMaterial :: computeShearModulusFromYoungAndPoisson(ym, nu);
    status->setBulkModulus(bulkModulus);
    status->setShearModulus(shearModulus);
    status->setYoungsModulus(ym);
}

void
DustMaterial :: computePlastStrainDirM1(StrainVector &answer, const StressVector &stressDeviator, double rho, double i1, double q)
{
    if ( ( answer.giveStressStrainMode() == _3dMat ) && ( stressDeviator.giveStressStrainMode() == _3dMat ) ) {
        for ( int i = 1; i <= 6; i++ ) {
            answer.at(i) = stressDeviator.at(i) / rho;
        }

        double temp = ( lambda * beta * exp(beta * i1) + theta ) * ( i1 - q ) / ( ft - q );
        answer.at(1) += temp;
        answer.at(2) += temp;
        answer.at(3) += temp;
    } else   {
        OOFEM_ERROR("Incorrect mode of stressstrainvector in DustMaterial :: computePlastStrainDirM1\n");
    }
}

void
DustMaterial :: computePlastStrainDirM2(StrainVector &answer, const StressVector &stressDeviator, double rho, double i1, double q)
{
    if ( ( answer.giveStressStrainMode() == _3dMat ) && ( stressDeviator.giveStressStrainMode() == _3dMat ) ) {
        double fc = functionFc(rho, i1, q);
        for ( int i = 1; i <= 6; i++ ) {
            answer.at(i) = stressDeviator.at(i) / fc;
        }

        double temp = ( q - i1 ) / ( rEll * rEll * fc );
        answer.at(1) -= temp;
        answer.at(2) -= temp;
        answer.at(3) -= temp;
    } else   {
        OOFEM_ERROR("Incorrect mode of stressstrainvector in DustMaterial :: computePlastStrainDirM2\n");
    }
}

void
DustMaterial :: computePlastStrainDirM3(StrainVector &answer, const StressVector &stressDeviator, double rho, double i1, double q)
{
    if ( ( answer.giveStressStrainMode() == _3dMat ) && ( stressDeviator.giveStressStrainMode() == _3dMat ) ) {
        double feft = functionFe(ft);
        for ( int i = 1; i <= 6; i++ ) {
            answer.at(i) = stressDeviator.at(i) / feft;
        }

        double dfeft = functionFeDI1(ft);
        double temp = 1 - ( 1 + dfeft ) * rho / feft;
        answer.at(1) += temp;
        answer.at(2) += temp;
        answer.at(3) += temp;
    } else   {
        OOFEM_ERROR("Incorrect mode of stressstrainvector in DustMaterial :: computePlastStrainDirM3\n");
    }
}

double
DustMaterial :: functionH(double q, double tempQ)
{
    double xq = functionX(q);
    double xtq = functionX(tempQ);
    switch ( hardeningType ) {
    case 1:
        return wHard * ( exp(dHard * xtq) - exp(dHard * xq) );

    default:     // 0
        return dHard * wHard * ( xtq / ( 1 - dHard * xtq ) - xq / ( 1 - dHard * xq ) );
    }
}

double
DustMaterial :: functionHDQ(double tempQ)
{
    double xtq = functionX(tempQ);
    double dxtq = functionXDQ(tempQ);
    switch ( hardeningType ) {
    case 1:
        return wHard * dHard * exp(dHard * xtq) * dHard * dxtq;

    default:     // 0
        return dHard * wHard * ( dxtq * ( 1 - dHard * xtq ) - xtq * ( -dHard * dxtq ) ) / ( 1 - dHard * xtq ) / ( 1 - dHard * xtq );
    }
}

double
DustMaterial :: functionI1(double q, double tempQ, double i1, double bulkModulus)
{
    return i1 - 3 *bulkModulus *functionH(q, tempQ);
}

double
DustMaterial :: functionI1DQ(double tempQ, double bulkModulus)
{
    return -3 *bulkModulus *functionHDQ(tempQ);
}

double
DustMaterial :: computeDeltaGamma2(double tempQ, double q, double i1, double bulkModulus)
{
    double vfH = functionH(q, tempQ);
    return rEll * rEll * functionFe(tempQ) * vfH / ( 3 * ( i1 - 3 * bulkModulus * vfH - tempQ ) );
}

double
DustMaterial :: computeDeltaGamma2DQ(double tempQ, double q, double i1, double bulkModulus)
{
    double vfH = functionH(q, tempQ);
    double vdfH = functionHDQ(tempQ);
    double vfEq = functionFe(tempQ);
    double vdfEq = functionFeDI1(tempQ);
    double frac1 = rEll * rEll * vfEq * vfH;
    double dfrac1 = rEll * rEll * ( vdfEq * vfH + vfEq * vdfH );
    double frac2 = 3 * ( i1 - 3 * bulkModulus * vfH - tempQ );
    double dfrac2 = 3 * ( -3 * bulkModulus * vdfH - 1 );
    return ( dfrac1 * frac2 - frac1 * dfrac2 ) / frac2 / frac2;
}

double
DustMaterial :: fTempR2(double tempQ, double q, double i1, double rho, double bulkModulus, double shearModulus)
{
    double vfEq = functionFe(tempQ);
    double dgq = computeDeltaGamma2(tempQ, q, i1, bulkModulus);
    double frac = ( vfEq + 2 * shearModulus * dgq );
    double a = rho * vfEq / frac;
    frac = ( rEll * rEll * vfEq + 9 * bulkModulus * dgq );
    double b = ( tempQ - i1 ) * rEll * vfEq / frac;
    return a * a + b * b - vfEq * vfEq;
}
} // end namespace oofem

