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

#include "cohint.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "contextioerr.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"

namespace oofem {
//---------------------------------------------------------------------------------------------------
// c l a s s   CohesiveInterfaceMaterial
//---------------------------------------------------------------------------------------------------

REGISTER_Material( CohesiveInterfaceMaterial );

CohesiveInterfaceMaterial :: CohesiveInterfaceMaterial(int n, Domain *d) : StructuralMaterial(n, d)
//
// constructor
//
{}

void
CohesiveInterfaceMaterial :: giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                                                  const FloatArray &totalStrain,
                                                  TimeStep *atTime)
//
// returns real stress vector in 3d stress space of receiver according to
// the previous level of stress and current strain increment
//
{
    CohesiveInterfaceMaterialStatus *status = static_cast< CohesiveInterfaceMaterialStatus * >( this->giveStatus(gp) );

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
CohesiveInterfaceMaterial :: giveStiffnessMatrix(FloatMatrix &answer,
                                                      MatResponseMode rMode,
                                                      GaussPoint *gp, TimeStep *atTime)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _3dInterface:
        give3dInterfaceMaterialStiffnessMatrix(answer, rMode, gp, atTime);
        break;
    default:
        StructuralMaterial :: giveStiffnessMatrix(answer, rMode, gp, atTime);
    }
}


void
CohesiveInterfaceMaterial :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
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
CohesiveInterfaceMaterial :: give3dInterfaceMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
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
    IR_GIVE_FIELD(ir, kn, _IFT_IsoInterfaceDamageMaterial_kn);
    IR_GIVE_FIELD(ir, ks, _IFT_IsoInterfaceDamageMaterial_ks);

    // thermal coefficient
    tempDillatCoeff = 0.0;

    return StructuralMaterial :: initializeFrom(ir);
}

void
CohesiveInterfaceMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralMaterial :: giveInputRecord(input);
    //input.setField(this->tempDillatCoeff, _IFT_CohesiveInterfaceMaterial_temperatureDillationCoeff);
    input.setField(this->kn, _IFT_IsoInterfaceDamageMaterial_kn);
    input.setField(this->ks, _IFT_IsoInterfaceDamageMaterial_ks);
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

