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

#include "simpleinterfacemat.h"
#include "structuralelement.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "mathfem.h"
#include "contextioerr.h"

namespace oofem {
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
    if ( mode == _1dInterface || mode == _2dInterface || mode == _3dInterface ) {
        return 1;
    }

    return 0;
}


void
SimpleInterfaceMaterial :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                                         MatResponseForm form,
                                                         MatResponseMode mode,
                                                         GaussPoint *gp,
                                                         TimeStep *atTime)
//
// computes full constitutive matrix for case of gp stress-strain state.
//
{
    _error("give3dMaterialStiffnessMatrix: not implemented");
}


void
SimpleInterfaceMaterial :: giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                                                const FloatArray &totalStrain,
                                                TimeStep *atTime)
//
// returns real stress vector in 3d stress space of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
{
    SimpleInterfaceMaterialStatus *status = ( SimpleInterfaceMaterialStatus * ) this->giveStatus(gp);
    this->initGpForNewStep(gp);
    FloatArray shearStrain(2), shearStress, strainVector;
    StructuralElement *el = ( StructuralElement * ) gp->giveElement();
    el->computeStrainVector(strainVector, gp, atTime);
    FloatArray tempShearStressShift = status->giveTempShearStressShift();
    const double normalStrain = strainVector.at(1);
    double normalStress, maxShearStress, dp;

    MaterialMode mMode = el->giveMaterialMode();
    //answer.resize(giveSizeOfReducedStressStrainVector(mMode));
    answer.zero();
    if ( normalStrain <= 0. ) {
        normalStress = this->kn * normalStrain;
        maxShearStress = fabs(normalStress) * this->frictCoeff;
    } else {
        normalStress = this->kn * this->stiffCoeff * normalStrain;
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
        _error("giveMaterialMode: Unsupported coord mode");
    }

    double lim = 1e50;
    answer.at(1) = normalStress > lim ? lim : normalStress < -lim ? -lim : normalStress;

    // update gp
    status->setTempShearStressShift(tempShearStressShift);
    status->letTempStrainVectorBe(strainVector);
    status->letTempStressVectorBe(answer);
}


void
SimpleInterfaceMaterial :: giveCharacteristicMatrix(FloatMatrix &answer,
                                                    MatResponseForm form, MatResponseMode rMode,
                                                    GaussPoint *gp, TimeStep *atTime)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    //MaterialMode mMode = gp->giveMaterialMode();
    MaterialMode mMode = gp->giveElement()->giveMaterialMode();

    //SimpleInterfaceMaterialStatus *status = ( SimpleInterfaceMaterialStatus * ) this->giveStatus(gp);
    FloatArray strainVector;
    StructuralElement *el = ( StructuralElement * ) gp->giveElement();
    el->computeStrainVector(strainVector, gp, atTime);
    double normalStrain = strainVector.at(1);
    answer.zero();
    switch ( mMode ) {
    case _1dInterface:
        answer.resize(1, 1);
        if ( rMode == SecantStiffness || rMode == TangentStiffness ) {
            if ( normalStrain <= 0 ) {
                answer.at(1, 1) = this->kn;
            } else {
                answer.at(1, 1) = this->kn * this->frictCoeff;
            }
        } else {
            if ( rMode == ElasticStiffness ) {
                answer.at(1, 1) = this->kn;
            } else {
                _error2( "give2dInterfaceMaterialStiffnessMatrix: unknown MatResponseMode (%s)", __MatResponseModeToString(rMode) );
            }
        }

        return;

    case _2dInterface:
        answer.resize(2, 2);
        if ( rMode == SecantStiffness || rMode == TangentStiffness ) {
            if ( normalStrain <= 0. ) {
                answer.at(1, 1) = answer.at(2, 2) = this->kn;
            } else {
                answer.at(1, 1) = answer.at(2, 2) = this->kn * this->frictCoeff;
            }
        } else {
            if ( rMode == ElasticStiffness ) {
                answer.at(1, 1) = answer.at(2, 2) = this->kn;
            } else {
                _error2( "give2dInterfaceMaterialStiffnessMatrix: unknown MatResponseMode (%s)", __MatResponseModeToString(rMode) );
            }
        }

        return;

    case _3dInterface:
        answer.resize(3, 3);
        if ( rMode == SecantStiffness || rMode == TangentStiffness ) {
            if ( normalStrain <= 0. ) {
                answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = this->kn;
            } else {
                answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = this->kn * this->frictCoeff;
            }
        } else {
            if ( rMode == ElasticStiffness ) {
                answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = this->kn;
            } else {
                _error2( "give2dInterfaceMaterialStiffnessMatrix: unknown MatResponseMode (%s)", __MatResponseModeToString(rMode) );
            }
        }

        return;

    default:
        StructuralMaterial :: giveCharacteristicMatrix(answer, form, rMode, gp, atTime);
    }
}


int
SimpleInterfaceMaterial :: giveSizeOfReducedStressStrainVector(MaterialMode mode)
//
// returns the size of reduced stress-strain vector
// according to mode given by gp.
//
{
    switch ( mode ) {
    case _1dInterface: return 1;

    case _2dInterface: return 2;

    case _3dInterface: return 3;

    default:
        return StructuralMaterial :: giveSizeOfReducedStressStrainVector(mode);
    }
}


int
SimpleInterfaceMaterial :: giveStressStrainComponentIndOf(MatResponseForm form, MaterialMode mMode, int ind)
//
// this function returns index of reduced(if form == ReducedForm)
// or Full(if form==FullForm) stressStrain component in Full or reduced
// stressStrainVector acording to stressStrain mode of given gp.
//
{
    if ( mMode == _1dInterface ) {
        return ind;
    } else {
        return StructuralMaterial :: giveStressStrainComponentIndOf(form, mMode, ind);
    }
}

void
SimpleInterfaceMaterial :: giveStressStrainMask(IntArray &answer, MatResponseForm form,
                                                MaterialMode mMode) const
//
// this function returns mask of reduced(if form == ReducedForm)
// or Full(if form==FullForm) stressStrain vector in full or
// reduced StressStrainVector
// acording to stressStrain mode of given gp.
//
//
// mask has size of reduced or full StressStrain Vector and  i-th component
// is index to full or reduced StressStrainVector where corresponding
// stressStrain resides.
//
// Reduced form is sub-vector (of stress or strain components),
// where components corresponding to imposed zero stress (plane stress,...)
// are not included. On the other hand, if zero strain component is imposed
// (Plane strain, ..) this condition must be taken into account in geometrical
// relations, and corresponding component is included in reduced vector.
//
{
    if ( mMode == _1dInterface ) {
        answer.resize(1);
        answer.at(1) = 1;
    } else {
        StructuralMaterial :: giveStressStrainMask(answer, form, mMode);
    }
}


void
SimpleInterfaceMaterial :: giveReducedCharacteristicVector(FloatArray &answer, GaussPoint *gp,
                                                           const FloatArray &charVector3d)
//
// returns reduced stressVector or strainVector from full 3d vector reduced
// to vector required by gp->giveStressStrainMode()
//
{
    MaterialMode mode = gp->giveMaterialMode();

    if ( mode == _1dInterface ) {
        answer = charVector3d;
        return;
    } else {
        StructuralMaterial :: giveReducedCharacteristicVector(answer, gp, charVector3d);
    }
}


void
SimpleInterfaceMaterial :: giveFullCharacteristicVector(FloatArray &answer,
                                                        GaussPoint *gp,
                                                        const FloatArray &strainVector)
//
// returns full 3d general strain vector from strainVector in reducedMode
// based on StressStrainMode in gp. Included are strains which
// perform nonzero work.
// General strain vector has one of the following forms:
// 1) strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
// 2) strainVectorShell {eps_x,eps_y,gamma_xy, kappa_x, kappa_y, kappa_xy, gamma_zx, gamma_zy}
//
// you must assigng your stress strain mode to one of the folloving modes (or add new)
// FullForm of MaterialStiffnessMatrix must have the same form.
//
{
    MaterialMode mode = gp->giveMaterialMode();
    if ( mode == _1dInterface ) {
        answer = strainVector;
        return;
    } else {
        StructuralMaterial :: giveFullCharacteristicVector(answer, gp, strainVector);
    }
}


int
SimpleInterfaceMaterial :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    return StructuralMaterial :: giveIPValue(answer, aGaussPoint, type, atTime);
}


InternalStateValueType
SimpleInterfaceMaterial :: giveIPValueType(InternalStateType type)
{
    return StructuralMaterial :: giveIPValueType(type);
}


int
SimpleInterfaceMaterial :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mMode)
{
    return StructuralMaterial :: giveIntVarCompFullIndx(answer, type, mMode);
}


int
SimpleInterfaceMaterial :: giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint)
{
    return StructuralMaterial :: giveIPValueSize(type, aGaussPoint);
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
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    frictCoeff = 0.;
    stiffCoeff = 0.;
    IR_GIVE_FIELD(ir, kn, IFT_SimpleInterfaceMaterial_kn, "kn");
    IR_GIVE_OPTIONAL_FIELD(ir, frictCoeff, IFT_SimpleInterfaceMaterial_frictCoeff, "fc");
    IR_GIVE_OPTIONAL_FIELD(ir, stiffCoeff, IFT_SimpleInterfaceMaterial_frictCoeff, "stiffcoeff");

    return StructuralMaterial :: initializeFrom(ir);
}


int
SimpleInterfaceMaterial :: giveInputRecordString(std :: string &str, bool keyword)
{
    char buff [ 1024 ];

    StructuralMaterial :: giveInputRecordString(str, keyword);

    sprintf(buff, " kn %e", kn);
    sprintf(buff, " frictCoeff %e", frictCoeff);
    sprintf(buff, " stiffCoeff %e", stiffCoeff);
    str += buff;

    return 1;
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
SimpleInterfaceMaterialStatus :: updateYourself(TimeStep *atTime)
{
    StructuralMaterialStatus :: updateYourself(atTime);
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
