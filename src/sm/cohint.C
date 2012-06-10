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

#include "cohint.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "contextioerr.h"

namespace oofem {
//---------------------------------------------------------------------------------------------------
// c l a s s   CohesiveInterfaceMaterial
//---------------------------------------------------------------------------------------------------

CohesiveInterfaceMaterial :: CohesiveInterfaceMaterial(int n, Domain *d) : StructuralMaterial(n, d)
//
// constructor
//
{}

void
CohesiveInterfaceMaterial :: giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                                                  const FloatArray &totalStrain,
                                                  TimeStep *atTime)
//
// returns real stress vector in 3d stress space of receiver according to
// the previous level of stress and current strain increment
//
{
    CohesiveInterfaceMaterialStatus *status = ( CohesiveInterfaceMaterialStatus * ) this->giveStatus(gp);

    // initialize
    this->initGpForNewStep(gp);
    answer.resize(3);

    // normal part of elastic stress-strain law
    answer.at(1) = kn * totalStrain.at(1);

    // shear part of elastic stress-strain law
    answer.at(2) = ks * totalStrain.at(2);
    answer.at(3) = ks * totalStrain.at(3);

    // update gp
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);
}

void
CohesiveInterfaceMaterial :: giveCharacteristicMatrix(FloatMatrix &answer,
                                                      MatResponseForm form, MatResponseMode rMode,
                                                      GaussPoint *gp, TimeStep *atTime)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _3dInterface:
        give3dInterfaceMaterialStiffnessMatrix(answer, form, rMode, gp, atTime);
        break;
    default:
        StructuralMaterial :: giveCharacteristicMatrix(answer, form, rMode, gp, atTime);
    }
}


void
CohesiveInterfaceMaterial :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
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

int
CohesiveInterfaceMaterial :: giveSizeOfReducedStressStrainVector(MaterialMode mode)
//
// returns the size of reduced stress-strain vector
// acording to mode given by gp.
//
{
    switch ( mode ) {
    case _3dInterface:
        return 3;

    default:
        return StructuralMaterial :: giveSizeOfReducedStressStrainVector(mode);
    }
}

int
CohesiveInterfaceMaterial :: giveStressStrainComponentIndOf(MatResponseForm form, MaterialMode mmode, int ind)
//
// this function returns index of reduced(if form == ReducedForm)
// or Full(if form==FullForm) stressStrain component in Full or reduced
// stressStrainVector acording to stressStrain mode of given gp.
//
{
    if ( mmode == _3dInterface ) {
        return ind;
    } else {
        return StructuralMaterial :: giveStressStrainComponentIndOf(form, mmode, ind);
    }
}

void
CohesiveInterfaceMaterial :: giveStressStrainMask(IntArray &answer, MatResponseForm form,
                                                  MaterialMode mmode) const
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
    int i;

    if ( mmode == _3dInterface ) {
        answer.resize(3);
        for ( i = 1; i <= 3; i++ ) {
            answer.at(i) = i;
        }
    } else {
        StructuralMaterial :: giveStressStrainMask(answer, form, mmode);
    }
}

void
CohesiveInterfaceMaterial :: giveReducedCharacteristicVector(FloatArray &answer, GaussPoint *gp,
                                                             const FloatArray &charVector3d)
//
// returns reduced stressVector or strainVector from full 3d vector reduced
// to vector required by gp->giveStressStrainMode()
//
{
    MaterialMode mode = gp->giveMaterialMode();

    if ( mode == _3dInterface ) {
        answer = charVector3d;
        return;
    } else {
        StructuralMaterial :: giveReducedCharacteristicVector(answer, gp, charVector3d);
    }
}

void
CohesiveInterfaceMaterial :: giveFullCharacteristicVector(FloatArray &answer,
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
    if ( mode == _3dInterface ) {
        answer = strainVector;
        return;
    } else                                                             {
        StructuralMaterial :: giveFullCharacteristicVector(answer, gp, strainVector);
    }
}

void
CohesiveInterfaceMaterial :: give3dInterfaceMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode rMode,
                                                                    GaussPoint *gp, TimeStep *atTime)
{
    // assemble elastic stiffness
    answer.resize(3, 3);
    answer.zero();
    answer.at(1, 1) = kn;
    answer.at(2, 2) = answer.at(3, 3) = ks;
}

int
CohesiveInterfaceMaterial :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    return StructuralMaterial :: giveIPValue(answer, aGaussPoint, type, atTime);
}

InternalStateValueType
CohesiveInterfaceMaterial :: giveIPValueType(InternalStateType type)
{
    return StructuralMaterial :: giveIPValueType(type);
}

int
CohesiveInterfaceMaterial :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode)
{
    return StructuralMaterial :: giveIntVarCompFullIndx(answer, type, mmode);
}

int
CohesiveInterfaceMaterial :: giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint)
{
    return StructuralMaterial :: giveIPValueSize(type, aGaussPoint);
}

void
CohesiveInterfaceMaterial :: giveThermalDilatationVector(FloatArray &answer,
                                                         GaussPoint *gp,  TimeStep *tStep)
//
// returns a FloatArray(6) of initial strain vector
// eps_0 = {exx_0, eyy_0, ezz_0, gyz_0, gxz_0, gxy_0}^T
// caused by unit temperature in direction of
// gp (element) local axes
//
{
    answer.resize(6);
    answer.zero();
    answer.at(1) = this->tempDillatCoeff;
    answer.at(2) = this->tempDillatCoeff;
    answer.at(3) = this->tempDillatCoeff;
}

IRResultType
CohesiveInterfaceMaterial :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    // elastic parameters
    IR_GIVE_FIELD(ir, kn, IFT_IsoInterfaceDamageMaterial_kn, "kn");
    IR_GIVE_FIELD(ir, ks, IFT_IsoInterfaceDamageMaterial_ks, "ks");

    // thermal coefficient
    tempDillatCoeff = 0.0;

    return StructuralMaterial :: initializeFrom(ir);
}

int
CohesiveInterfaceMaterial :: giveInputRecordString(std :: string &str, bool keyword)
{
    char buff [ 1024 ];

    StructuralMaterial :: giveInputRecordString(str, keyword);

    sprintf(buff, " talpha %e kn %e ks %e", this->tempDillatCoeff, kn, ks);
    str += buff;

    return 1;
}

//---------------------------------------------------------------------------------------------------
// c l a s s   CohesiveInterfaceMaterialStatus
//---------------------------------------------------------------------------------------------------

CohesiveInterfaceMaterialStatus :: CohesiveInterfaceMaterialStatus(int n, Domain *d, GaussPoint *g) : StructuralMaterialStatus(n, d, g)
{}

CohesiveInterfaceMaterialStatus :: ~CohesiveInterfaceMaterialStatus()
{}

void
CohesiveInterfaceMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { ");
    fprintf(file, "}\n");
}

void
CohesiveInterfaceMaterialStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();
}

void
CohesiveInterfaceMaterialStatus :: updateYourself(TimeStep *atTime)
{
    StructuralMaterialStatus :: updateYourself(atTime);
}

contextIOResultType
CohesiveInterfaceMaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}

contextIOResultType
CohesiveInterfaceMaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}
} // namespace oofem

