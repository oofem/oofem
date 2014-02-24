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

#include "simpleinterfacemat.h"
#include "structuralelement.h"
#include "interfaceelement1d.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "contextioerr.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"

namespace oofem {
REGISTER_Material(SimpleInterfaceMaterial);

SimpleInterfaceMaterial :: SimpleInterfaceMaterial(int n, Domain *d) : StructuralMaterial(n, d)
    //
    // constructor
    //
{ }


SimpleInterfaceMaterial :: ~SimpleInterfaceMaterial()
//
// destructor
//
{ }

int
SimpleInterfaceMaterial :: hasMaterialModeCapability(MaterialMode mode)
//
// returns whether receiver supports given mode
//
{
    return mode == _1dInterface || mode == _2dInterface || mode == _3dInterface;
}


void
SimpleInterfaceMaterial :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                                         MatResponseMode mode,
                                                         GaussPoint *gp,
                                                         TimeStep *tStep)
//
// computes full constitutive matrix for case of gp stress-strain state.
//
{
    OOFEM_ERROR("not implemented");
}


void
SimpleInterfaceMaterial :: giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                                                //const FloatArray &totalStrain,// @todo temporary -should not be here /JB
                                                const FloatArray &strainVector,
                                                TimeStep *tStep)
//
// returns real stress vector in 3d stress space of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
{
    SimpleInterfaceMaterialStatus *status = static_cast< SimpleInterfaceMaterialStatus * >( this->giveStatus(gp) );
    this->initGpForNewStep(gp);
    FloatArray shearStrain(2), shearStress; //, strainVector;
    StructuralElement *el = static_cast< StructuralElement * >( gp->giveElement() );
    //el->computeStrainVector(strainVector, gp, tStep);

    FloatArray tempShearStressShift = status->giveTempShearStressShift();
    const double normalStrain = strainVector.at(1);
    double normalStress, maxShearStress, dp;
    double shift = -this->kn * this->stiffCoeff * normalClearance;

    MaterialMode mMode = el->giveMaterialMode();
    //answer.resize(giveSizeOfReducedStressStrainVector(mMode));
    answer.zero();
    if ( normalStrain + normalClearance <= 0. ) {
        normalStress = this->kn * ( normalStrain + normalClearance ) + shift; //in compression and after the clearance gap closed
        maxShearStress = fabs(normalStress) * this->frictCoeff;
    } else {
        normalStress = this->kn * this->stiffCoeff * ( normalStrain + normalClearance ) + shift;
        maxShearStress = 0.;
    }

    switch ( mMode ) {
    case _1dInterface:
        answer.resize(1);
        break;
    case _2dInterface:
        answer.resize(2);
        shearStrain.at(1) = strainVector.at(2);
        shearStress.beScaled(this->kn, shearStrain);
        shearStress.subtract(tempShearStressShift);
        dp = shearStress.dotProduct(shearStress, 1);
        if ( dp > maxShearStress * maxShearStress ) {
            shearStress.times( maxShearStress / sqrt(dp) );
        }

        tempShearStressShift.beScaled(this->kn, shearStrain);
        tempShearStressShift.subtract(shearStress);
        answer.at(2) = shearStress.at(1);
        break;
    case _3dInterface:
    case _3dMat: //JB
        answer.resize(3);
        shearStrain.at(1) = strainVector.at(2);
        shearStrain.at(2) = strainVector.at(3);
        shearStress.beScaled(this->kn, shearStrain);
        shearStress.subtract(tempShearStressShift);
        dp = shearStress.dotProduct(shearStress, 2);
        if ( dp > maxShearStress * maxShearStress ) {
            shearStress.times( maxShearStress / sqrt(dp) );
        }

        tempShearStressShift.beScaled(this->kn, shearStrain);
        tempShearStressShift.subtract(shearStress);
        answer.at(2) = shearStress.at(1);
        answer.at(3) = shearStress.at(2);
        break;
    default:
        OOFEM_ERROR("Unsupported interface mode");
    }

    double lim = 1.e+50;
    answer.at(1) = min(normalStress, lim);  //threshold on maximum
    answer.at(1) = max(answer.at(1), -lim);  //threshold on minimum
    //answer.at(1) = normalStress > lim ? lim : normalStress < -lim ? -lim : normalStress;
    // update gp
    status->setTempShearStressShift(tempShearStressShift);
    status->letTempStrainVectorBe(strainVector);
    status->letTempStressVectorBe(answer);
}


void
SimpleInterfaceMaterial :: giveStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode rMode,
                                               GaussPoint *gp, TimeStep *tStep)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    MaterialMode mMode = gp->giveElement()->giveMaterialMode();

    FloatArray strainVector;
    StructuralElement *el = static_cast< StructuralElement * >( gp->giveElement() );
    double normalStrain;

    el->computeStrainVector(strainVector, gp, tStep);
    normalStrain = strainVector.at(1);
    answer.zero();
    switch ( mMode ) {
    case _1dInterface:
        answer.resize(1, 1);
        if ( rMode == SecantStiffness || rMode == TangentStiffness ) {
            if ( normalStrain + normalClearance <= 0 ) {
                answer.at(1, 1) = this->kn; //in compression and after the clearance gap closed
            } else {
                answer.at(1, 1) = this->kn * this->stiffCoeff;
            }
        } else {
            if ( rMode == ElasticStiffness ) {
                answer.at(1, 1) = this->kn;
            } else {
                OOFEM_ERROR("unknown MatResponseMode (%s)", __MatResponseModeToString(rMode) );
            }
        }

        return;

    case _2dInterface:
        answer.resize(2, 2);
        if ( rMode == SecantStiffness || rMode == TangentStiffness ) {
            if ( normalStrain + normalClearance <= 0. ) {
                answer.at(1, 1) = answer.at(2, 2) = this->kn; //in compression and after the clearance gap closed
            } else {
                answer.at(1, 1) = answer.at(2, 2) = this->kn * this->stiffCoeff;
            }
        } else {
            if ( rMode == ElasticStiffness ) {
                answer.at(1, 1) = answer.at(2, 2) = this->kn;
            } else {
                OOFEM_ERROR("unknown MatResponseMode (%s)", __MatResponseModeToString(rMode) );
            }
        }

        return;

    case _3dInterface:
        answer.resize(3, 3);
        if ( rMode == SecantStiffness || rMode == TangentStiffness ) {
            if ( normalStrain + normalClearance <= 0. ) {
                answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = this->kn; //in compression and after the clearance gap closed
            } else {
                answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = this->kn * this->stiffCoeff;
            }
        } else {
            if ( rMode == ElasticStiffness ) {
                answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = this->kn;
            } else {
                OOFEM_ERROR("unknown MatResponseMode (%s)", __MatResponseModeToString(rMode) );
            }
        }

        return;

    default:
        StructuralMaterial :: giveStiffnessMatrix(answer, rMode, gp, tStep);
    }
}


int
SimpleInterfaceMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    return StructuralMaterial :: giveIPValue(answer, gp, type, tStep);
}


void
SimpleInterfaceMaterial :: giveThermalDilatationVector(FloatArray &answer,
                                                       GaussPoint *gp,  TimeStep *tStep)
//
// returns a strain vector
// eps_0 = {exx_0, eyy_0, ezz_0, gyz_0, gxz_0, gxy_0}^T
// caused by unit temperature in direction of
// gp (element) local axes
//
{
    answer.resize(1);
    answer.at(1) = 0.0;
}


IRResultType
SimpleInterfaceMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    frictCoeff = 0.;
    stiffCoeff = 0.;
    normalClearance = 0.;
    IR_GIVE_FIELD(ir, kn, _IFT_SimpleInterfaceMaterial_kn);
    IR_GIVE_OPTIONAL_FIELD(ir, frictCoeff, _IFT_SimpleInterfaceMaterial_frictCoeff);
    IR_GIVE_OPTIONAL_FIELD(ir, stiffCoeff, _IFT_SimpleInterfaceMaterial_stiffCoeff);
    IR_GIVE_OPTIONAL_FIELD(ir, normalClearance, _IFT_SimpleInterfaceMaterial_normalClearance);

    return StructuralMaterial :: initializeFrom(ir);
}


void
SimpleInterfaceMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralMaterial :: giveInputRecord(input);
    input.setField(this->kn, _IFT_SimpleInterfaceMaterial_kn);
    input.setField(this->frictCoeff, _IFT_SimpleInterfaceMaterial_frictCoeff);
    input.setField(this->stiffCoeff, _IFT_SimpleInterfaceMaterial_stiffCoeff);
    input.setField(this->normalClearance, _IFT_SimpleInterfaceMaterial_normalClearance);
}


SimpleInterfaceMaterialStatus :: SimpleInterfaceMaterialStatus(int n, Domain *d, GaussPoint *g) : StructuralMaterialStatus(n, d, g)
{
    shearStressShift.resize(2);
    tempShearStressShift.resize(2);
    shearStressShift.zero();
    tempShearStressShift.zero();
}


SimpleInterfaceMaterialStatus :: ~SimpleInterfaceMaterialStatus()
{ }


void
SimpleInterfaceMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { ");
    fprintf( file, "shearStressShift (%f, %f)", this->shearStressShift.at(1), this->shearStressShift.at(2) );
    fprintf(file, "}\n");
}


void
SimpleInterfaceMaterialStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();
    tempShearStressShift = shearStressShift;
}


void
SimpleInterfaceMaterialStatus :: updateYourself(TimeStep *tStep)
{
    StructuralMaterialStatus :: updateYourself(tStep);
    shearStressShift = tempShearStressShift;
}


FloatArray
SimpleInterfaceMaterialStatus :: giveShearStressShift()
{
    FloatArray answer = shearStressShift;
    return answer;
}


FloatArray
SimpleInterfaceMaterialStatus :: giveTempShearStressShift()
{
    FloatArray answer = tempShearStressShift;
    return answer;
}


contextIOResultType
SimpleInterfaceMaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write a raw data
    //if ( !stream->write(& kappa, 1) ) {
    //THROW_CIOERR(CIO_IOERR);
    //}

    return CIO_OK;
}


contextIOResultType
SimpleInterfaceMaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
    //if ( !stream->read(& kappa, 1) ) {
    //THROW_CIOERR(CIO_IOERR);
    //}

    return CIO_OK;
}
} // end namespace oofem
