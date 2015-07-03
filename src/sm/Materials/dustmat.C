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

#include "dustmat.h"

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
REGISTER_Material(DustMaterial);

DustMaterialStatus :: DustMaterialStatus(int n, Domain *d, GaussPoint *gp, double q0) :
    StructuralMaterialStatus(n, d, gp),
    plasticStrain( 6 ),
    tempPlasticStrain( 6 )
{
    stressVector.resize(6);
    strainVector.resize(6);
    tempStressVector.resize(6);
    tempStrainVector.resize(6);
    q = q0;
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
DustMaterialStatus :: updateYourself(TimeStep *tStep)
{
    // Call the corresponding function of the parent class to update variables defined there.
    StructuralMaterialStatus :: updateYourself(tStep);
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

    fprintf(file, ", plasticStrains ");
    for ( auto &val : this->givePlasticStrain() ) {
        fprintf( file, " %.4e", val );
    }

    fprintf(file, ", q  %.4e", q);

    fprintf(file, "}\n");
}

contextIOResultType
DustMaterialStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


contextIOResultType
DustMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
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
    IRResultType result;
    // call the corresponding service of structural material
    result = StructuralMaterial :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    // call the corresponding service for the linear elastic material
    result = this->LEMaterial->initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

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

    IR_GIVE_OPTIONAL_FIELD(ir, alpha, _IFT_DustMaterial_alpha);
    IR_GIVE_OPTIONAL_FIELD(ir, beta, _IFT_DustMaterial_beta);
    IR_GIVE_OPTIONAL_FIELD(ir, lambda, _IFT_DustMaterial_lambda);
    IR_GIVE_OPTIONAL_FIELD(ir, theta, _IFT_DustMaterial_theta);
    IR_GIVE_OPTIONAL_FIELD(ir, ft, _IFT_DustMaterial_ft);
    IR_GIVE_OPTIONAL_FIELD(ir, rEll, _IFT_DustMaterial_rEll);
    IR_GIVE_OPTIONAL_FIELD(ir, x0, _IFT_DustMaterial_x0);
    IR_GIVE_OPTIONAL_FIELD(ir, wHard, _IFT_DustMaterial_wHard);
    IR_GIVE_OPTIONAL_FIELD(ir, dHard, _IFT_DustMaterial_dHard);
    IR_GIVE_OPTIONAL_FIELD(ir, mStiff, _IFT_DustMaterial_mStiff);
    IR_GIVE_OPTIONAL_FIELD(ir, mStiff, _IFT_DustMaterial_newtonTol);
    IR_GIVE_OPTIONAL_FIELD(ir, newtonIter, _IFT_DustMaterial_newtonIter);

    // check parameters admissibility
    if ( ft < 0 ) {
        OOFEM_WARNING("parameter 'ft' must be positive");
        return IRRT_BAD_FORMAT;
    }

    if ( x0 < 0 ) {
        OOFEM_WARNING("parameter 'x0' must be positive");
        return IRRT_BAD_FORMAT;
    }

    if ( rEll < 0 ) {
        OOFEM_WARNING("parameter 'rEll' must be positive");
        return IRRT_BAD_FORMAT;
    }

    if ( theta < 0 ) {
        OOFEM_WARNING("parameter 'theta' must be positive");
        return IRRT_BAD_FORMAT;
    }

    if ( beta < 0 ) {
        OOFEM_WARNING("parameter 'beta' must be positive");
        return IRRT_BAD_FORMAT;
    }

    if ( lambda < 0 ) {
        OOFEM_WARNING("parameter 'lambda' must be positive");
        return IRRT_BAD_FORMAT;
    }

    if ( alpha < lambda ) {
        OOFEM_WARNING("parameter 'alpha' must be greater than parameter 'lambda'");
        return IRRT_BAD_FORMAT;
    }

    x0 = -x0; // compressive strength is negative, although on input it is a positive number

    hardeningType = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, hardeningType, _IFT_DustMaterial_hardeningType);

    q0 = x0;
    solveQ0(q0);

    return IRRT_OK;
}

void
DustMaterial :: giveRealStressVector_3d(FloatArray &answer,
                                     GaussPoint *gp,
                                     const FloatArray &totalStrain,
                                     TimeStep *tStep)
{
    FloatArray strainVectorR;

    DustMaterialStatus *status = static_cast< DustMaterialStatus * >( this->giveStatus(gp) );

    // Initialize temp variables for this Gauss point
    this->initTempStatus(gp);

    // subtract stress-independent part of strain
    this->giveStressDependentPartOfStrainVector(strainVectorR, gp, totalStrain, tStep, VM_Total);

    // perform the local stress return and update the history variables
    performStressReturn(gp, strainVectorR);

    // copy total strain vector to the temp status
    status->letTempStrainVectorBe(totalStrain);

    // pass the correct form of stressVector to giveRealStressVector
    answer = status->giveTempStressVector();
}

void
DustMaterial :: performStressReturn(GaussPoint *gp, const FloatArray &strain)
{
    DustMaterialStatus *status = static_cast< DustMaterialStatus * >( giveStatus(gp) );

    // compute total strain components
    FloatArray strainDeviator;
    double volumetricStrain;
    volumetricStrain = computeDeviatoricVolumetricSplit(strainDeviator, strain);

    // compute trial elastic strains
    FloatArray plasticStrain = status->givePlasticStrain();
    double volumetricPlasticStrain;
    FloatArray plasticStrainDeviator;
    volumetricPlasticStrain = computeDeviatoricVolumetricSplit(plasticStrainDeviator, plasticStrain);
    double volumetricElasticStrain = volumetricStrain - volumetricPlasticStrain;
    FloatArray elasticStrainDeviator = strainDeviator;
    elasticStrainDeviator.subtract(plasticStrainDeviator);

    // compute trial stresses
    double bulkModulus, shearModulus;
    computeAndSetBulkAndShearModuli(bulkModulus, shearModulus, gp);
    double volumetricStress = 3. * bulkModulus * volumetricElasticStrain;
    FloatArray stressDeviator = {2 * elasticStrainDeviator[0], 2 * elasticStrainDeviator[1], 2 * elasticStrainDeviator[2], 
                                elasticStrainDeviator[3], elasticStrainDeviator[4], elasticStrainDeviator[5]};
    stressDeviator.times(shearModulus);

    // norm of trial stress deviator
    double rho = computeSecondCoordinate(stressDeviator);
    double i1 = 3 * volumetricStress;
    double f1, f2, f3, q, tempQ;
    q = tempQ = status->giveQ();
    f1 = yieldFunction1(rho, i1);
    f2 = yieldFunction2(rho, i1, q);
    f3 = yieldFunction3(i1);

    // actual stress return
    double lambda = 0.;
    FloatArray m(6);
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
    FloatArray mDeviator;
    mVol = computeDeviatoricVolumetricSplit(mDeviator, m);
    i1 -= 3 * bulkModulus * mVol;
    volumetricStress = i1 / 3.;
    mDeviator.times(-2 * shearModulus);
    stressDeviator.add(mDeviator);

    // compute full stresses from deviatoric and volumetric part and store them
    FloatArray stress = stressDeviator;
    stress.at(1) += volumetricStress;
    stress.at(2) += volumetricStress;
    stress.at(3) += volumetricStress;
    status->letTempStressVectorBe(stress);

    // compute and update plastic strain and q
    status->letTempPlasticStrainBe(plasticStrain);
}

void
DustMaterial :: performF1return(double i1, double rho, GaussPoint *gp)
{
    DustMaterialStatus *status = static_cast< DustMaterialStatus * >( giveStatus(gp) );
    double bulkModulus = status->giveBulkModulus();
    double shearModulus = status->giveShearModulus();
    double q = status->giveQ();
    double tempQ = status->giveTempQ();
    double fx, dfx;
    double m = 9 * bulkModulus / ( 2 * shearModulus );
    double vfI1, vfI1DQ, a, b, c, d, da, db, dc;
    int positiveFlag = 0;

    for ( int i = 0; i < newtonIter; i++ ) {
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
    OOFEM_ERROR("Newton's method did not converge");
}

void
DustMaterial :: performF2return(double i1, double rho, GaussPoint *gp)
{
    DustMaterialStatus *status = static_cast< DustMaterialStatus * >( giveStatus(gp) );
    double bulkModulus = status->giveBulkModulus();
    double shearModulus = status->giveShearModulus();
    double q = status->giveQ();
    double qRight = q;
    double qLeft = q;
    double tempQ = .5 * ( qLeft + qRight );
    double fq;
    double fx, dfx;
    for ( int i = 0; i < newtonIter; i++ ) {
        fx = i1 - 3 *bulkModulus *functionH(q, qLeft) - qLeft;
        dfx =   -3 *bulkModulus *functionHDQ(qLeft) - 1;
        qLeft -= fx / dfx;
        if (  fabs(fx / dfx / q0) < newtonTol ) {
            break;
        }
    }

    for ( int i = 0; i < newtonIter; i++ ) {
        fq = fTempR2(tempQ, q, i1, rho, bulkModulus, shearModulus);
        if ( fabs( ( qRight - qLeft ) / qRight ) < newtonTol ) {
            status->letTempQBe(tempQ);
            return;
        }

        if ( fq > 0 ) {
            qRight = tempQ;
        } else {
            qLeft = tempQ;
        }

        tempQ = .5 * ( qLeft + qRight );
    }

    OOFEM_ERROR("bisection method did not converge");
}

void
DustMaterial :: computeQFromPlastVolEps(double &answer, double q, double deltaVolumetricPlasticStrain)
{
    if ( q >= 0. ) {
        answer = 0.;
        return;
    }

    double fx, dfx;
    for ( int i = 0; i <= newtonIter; i++ ) {
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
    OOFEM_ERROR("Newton's method did not converge");
}

void
DustMaterial :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                              MatResponseMode mode,
                                              GaussPoint *gp,
                                              TimeStep *tStep)
{
    DustMaterialStatus *status = static_cast< DustMaterialStatus * >( giveStatus(gp) );
    double ym0 = LEMaterial->giveYoungsModulus();
    double ym = status->giveYoungsModulus();
    double coeff = status->giveVolumetricPlasticStrain() < 0 ? ym / ym0 : 1.0;
    if ( mode == ElasticStiffness ) {
        LEMaterial->give3dMaterialStiffnessMatrix(answer, mode, gp, tStep);
    } else if ( mode == SecantStiffness || mode == TangentStiffness ) {
        LEMaterial->give3dMaterialStiffnessMatrix(answer, mode, gp, tStep);
        answer.times(coeff);
    } else {
        OOFEM_ERROR("Unsupported MatResponseMode");
    }
}

int
DustMaterial :: setIPValue(const FloatArray &value, GaussPoint *gp, InternalStateType type)
{
    DustMaterialStatus *status = static_cast< DustMaterialStatus * >( giveStatus(gp) );
    if ( type == IST_PlasticStrainTensor ) {
        status->letPlasticStrainBe(value);
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
                            TimeStep *tStep)
{
    const DustMaterialStatus *status = static_cast< DustMaterialStatus * >( giveStatus(gp) );
    if ( type == IST_PlasticStrainTensor ) {
        answer = status->givePlasticStrain();
        return 1;
    } else if ( type == IST_PrincipalPlasticStrainTensor ) {
        computePrincipalValues(answer, status->givePlasticStrain(), principal_strain);
        return 1;
    } else if ( type == IST_VolumetricPlasticStrain ) {
        answer.resize(1);
        answer.at(1) = status->giveVolumetricPlasticStrain(); ///@todo This is actually the mean, not the volumetric part. / Mikael
        return 1;
    } else if ( type == IST_StressCapPos ) {
        answer.resize(1);
        answer.at(1) = status->giveQ();
        return 1;
    } else {
        return StructuralMaterial :: giveIPValue(answer, gp, type, tStep);
    }

    return 0;
}

MaterialStatus *
DustMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new DustMaterialStatus( 1, StructuralMaterial :: giveDomain(), gp, this->giveQ0() );
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
DustMaterial :: yieldFunction1(double rho, double i1)
{
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
    for ( int i = 0; i < newtonIter; i++ ) {
        fx = -x0 + answer - rEll * ( alpha -      lambda * exp(beta * answer) - theta * answer );
        dfx =       1 - rEll * (      -beta * lambda * exp(beta * answer) - theta );
        answer -= fx / dfx;
        if (  fabs(fx / dfx / answer) < newtonTol ) {
            if ( answer >= 0 ) {
                OOFEM_ERROR("internal parameter q has to be negative");
            }

            return;
        }
    }

    OOFEM_ERROR("Newton's method did not converge");
}

void
DustMaterial :: computeAndSetBulkAndShearModuli(double &bulkModulus, double &shearModulus, GaussPoint *gp)
{
    DustMaterialStatus *status = static_cast< DustMaterialStatus * >( giveStatus(gp) );
    double ym = LEMaterial->giveYoungsModulus();
    double nu = LEMaterial->givePoissonsRatio();
    double volumetricPlasticStrain = status->giveVolumetricPlasticStrain();
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
DustMaterial :: computePlastStrainDirM1(FloatArray &answer, const FloatArray &stressDeviator, double rho, double i1, double q)
{
    answer.beScaled(1./rho, stressDeviator);

    double temp = ( lambda * beta * exp(beta * i1) + theta ) * ( i1 - q ) / ( ft - q );
    answer.at(1) += temp;
    answer.at(2) += temp;
    answer.at(3) += temp;
}

void
DustMaterial :: computePlastStrainDirM2(FloatArray &answer, const FloatArray &stressDeviator, double rho, double i1, double q)
{
    double fc = functionFc(rho, i1, q);
    answer.beScaled(1./fc, stressDeviator);

    double temp = ( q - i1 ) / ( rEll * rEll * fc );
    answer.at(1) -= temp;
    answer.at(2) -= temp;
    answer.at(3) -= temp;
}

void
DustMaterial :: computePlastStrainDirM3(FloatArray &answer, const FloatArray &stressDeviator, double rho, double i1, double q)
{
    double feft = functionFe(ft);
    answer.beScaled(1./feft, stressDeviator);

    double dfeft = functionFeDI1(ft);
    double temp = 1 - ( 1 + dfeft ) * rho / feft;
    answer.at(1) += temp;
    answer.at(2) += temp;
    answer.at(3) += temp;
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
        return wHard *dHard *exp(dHard *xtq) * dHard * dxtq;

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
    return rEll *rEll *functionFe(tempQ) * vfH / ( 3 * ( i1 - 3 * bulkModulus * vfH - tempQ ) );
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
