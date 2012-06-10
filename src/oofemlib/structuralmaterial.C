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

#include "structuralmaterial.h"
#include "structuralcrosssection.h"
#include "domain.h"
#include "verbose.h"
#include "structuralms.h"
#include "structuralelement.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "mathfem.h"
#include "engngm.h"
#include "fieldmanager.h"

namespace oofem {
int
StructuralMaterial :: hasMaterialModeCapability(MaterialMode mode)
//
// returns whether receiver supports given mode
//
{
    if ( ( mode == _3dMat ) || ( mode == _PlaneStress ) ||
        ( mode == _PlaneStrain ) || ( mode == _1dMat ) ||
        ( mode == _2dPlateLayer ) || ( mode == _2dBeamLayer ) ||
        ( mode == _3dShellLayer ) || ( mode == _1dFiber ) ) {
        return 1;
    }

    return 0;
}

void
StructuralMaterial :: giveCharacteristicMatrix(FloatMatrix &answer,
                                               MatResponseForm form, MatResponseMode rMode,
                                               GaussPoint *gp, TimeStep *atTime)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _3dMat:
    case _3dMatGrad:
    case _3dMat_F: // even if material uses deformation gradient, stiffness is computed in the usual way
        this->give3dMaterialStiffnessMatrix(answer, form, rMode, gp, atTime);
        break;
    case _PlaneStress:
    case _PlaneStressGrad:
        this->givePlaneStressStiffMtrx(answer, form, rMode, gp, atTime);
        break;
    case _PlaneStrain:
    case _PlaneStrainGrad:
        this->givePlaneStrainStiffMtrx(answer, form, rMode, gp, atTime);
        break;
    case _1dMat:
        this->give1dStressStiffMtrx(answer, form, rMode, gp, atTime);
        break;
    case _2dPlateLayer:
        this->give2dPlateLayerStiffMtrx(answer, form, rMode, gp, atTime);
        break;
    case _3dShellLayer:
        this->give3dShellLayerStiffMtrx(answer, form, rMode, gp, atTime);
        break;
    case _2dBeamLayer:
        this->give2dBeamLayerStiffMtrx(answer, form, rMode, gp, atTime);
        break;
    case _1dFiber:
        this->give1dFiberStiffMtrx(answer, form, rMode, gp, atTime);
        break;
    default:
        OOFEM_ERROR2("StructuralMaterial :: giveCharacteristicMatrix : unknown mode (%s)", __MaterialModeToString(mMode) );
        return;
    }
}


void
StructuralMaterial :: giveCharacteristicComplianceMatrix(FloatMatrix &answer,
                                                         MatResponseForm form, MatResponseMode rMode,
                                                         GaussPoint *gp,
                                                         TimeStep *atTime)
//
// Returns characteristic material compliance matrix of the receiver
// works for positive definite associated stiffnesses only
//
{
    FloatMatrix redInvAnswer, redAnswer;
    IntArray mask;

    this->giveCharacteristicMatrix(redInvAnswer, ReducedForm, rMode, gp, atTime);
    redAnswer.beInverseOf(redInvAnswer);

    if ( form == FullForm ) {
        this->giveStressStrainMask( mask, FullForm, gp->giveMaterialMode() );
        answer.beSubMatrixOf(redAnswer, mask);
    } else if ( form == ReducedForm ) {
        answer = redAnswer;
    } else {
        OOFEM_ERROR("StructuralMaterial :: giveCharacteristicComplianceMatrix - unsupported form mode");
    }
}


void
StructuralMaterial ::  reduceStiffMtrx3d(FloatMatrix &answer, MatResponseForm form, GaussPoint *gp,
                                         FloatMatrix &stiffMtrx3d) const
//
// Returns characteristic material stiffness matrix of the receiver
// reduced to corresponding mode obtained from gp.
{
    MaterialMode mode = gp->giveMaterialMode();
    switch ( mode ) {
    case _3dMat:
        answer = stiffMtrx3d;
        break;
    case _PlaneStress:
        this->reduceToPlaneStressStiffMtrx(answer, form, gp, stiffMtrx3d);
        break;
    case _PlaneStrain:
        this->reduceToPlaneStrainStiffMtrx(answer, form, gp, stiffMtrx3d);
        break;
    case _1dMat:
        this->reduceTo1dStressStiffMtrx(answer, form, gp, stiffMtrx3d);
        break;
    case _2dPlateLayer:
        this->reduceTo2dPlateLayerStiffMtrx(answer, form, gp, stiffMtrx3d);
        break;
    case _3dShellLayer:
        this->reduceTo3dShellLayerStiffMtrx(answer, form, gp, stiffMtrx3d);
        break;
    case _2dBeamLayer:
        this->reduceTo2dBeamLayerStiffMtrx(answer, form, gp, stiffMtrx3d);
        break;
    case _1dFiber:
        this->reduceTo1dFiberStiffMtrx(answer, form, gp, stiffMtrx3d);
        break;
    default:
        OOFEM_ERROR2("StructuralMaterial :: reduceStiffMtrx3d : unknown mode (%s)", __MaterialModeToString(mode) );
        return;
    }
}


void
StructuralMaterial :: reduceComplMtrx3d(FloatMatrix &answer, MatResponseForm form, GaussPoint *gp,
                                        FloatMatrix &complMtrx3d) const
//
// Returns characteristic material compliance matrix of the receiver
// reduced to corresponding mode obtained from gp.
{
    MaterialMode mode = gp->giveMaterialMode();
    switch ( mode ) {
    case _3dMat:
        answer = complMtrx3d;
        break;
    case _PlaneStress:
        this->reduceToPlaneStressComplMtrx(answer, form, gp, complMtrx3d);
        break;
    case _PlaneStrain:
        this->reduceToPlaneStrainComplMtrx(answer, form, gp, complMtrx3d);
        break;
    case _1dMat:
        this->reduceTo1dStressComplMtrx(answer, form, gp, complMtrx3d);
        break;
    case _2dPlateLayer:
        this->reduceTo2dPlateLayerComplMtrx(answer, form, gp, complMtrx3d);
        break;
    case _3dShellLayer:
        this->reduceTo3dShellLayerComplMtrx(answer, form, gp, complMtrx3d);
        break;
    case _2dBeamLayer:
        this->reduceTo2dBeamLayerComplMtrx(answer, form, gp, complMtrx3d);
        break;
    case _1dFiber:
        this->reduceTo1dFiberComplMtrx(answer, form, gp, complMtrx3d);
        break;
    default:
        OOFEM_ERROR2("StructuralMaterial :: reduceComplMtrx3d : unknown mode (%s)", __MaterialModeToString(mode) );
        return;
    }
}


void
StructuralMaterial :: giveStressDependentPartOfStrainVector(FloatArray &answer, GaussPoint *gp,
                                                            const FloatArray &reducedStrainVector,
                                                            TimeStep *stepN, ValueModeType mode)
{
    /*
     * This functions subtract from reducedStrainVector its stress independent part
     * caused by temperature, shrinkage and possibly by other phenomenas.
     */
    FloatArray epsilonTemperature;
    StructuralCrossSection *cs = ( StructuralCrossSection * ) gp->giveElement()->giveCrossSection();

    answer = reducedStrainVector;
    cs->computeStressIndependentStrainVector(epsilonTemperature, gp, stepN, mode);
    if ( epsilonTemperature.giveSize() ) {
        answer.subtract(epsilonTemperature);
    }
}


int
StructuralMaterial :: giveSizeOfReducedStressStrainVector(MaterialMode mode)
//
// returns the size of reduced stress-strain vector
// according to mode given by gp.
//
{
    switch ( mode ) {
    case _3dMat:
    case _3dMat_F:
        return 6;

    case _PlaneStress:
        return 3;

    case _PlaneStrain:
        return 4;

    case _PlaneStressRot:
        return 4;

    case _1dMat:
        return 1;

    case _2dPlateLayer:
    case _3dShellLayer:
        return 5;

    case _2dBeamLayer:
        return 2;

    case _2dPlate:
        return 5;

    case _2dBeam:
        return 3;

    case _3dShell:
        return 8;

    case _3dBeam:
        return 6;

    case _1dFiber:
        return 3;

    case _1dMatGrad:
        return 2;

    case _PlaneStressGrad:
        return 4;

    case _PlaneStrainGrad:
        return 5;

    case _3dMatGrad:
        return 7;

    default:
        OOFEM_ERROR2("StructuralMaterial :: giveSizeOfReducedStressStrainVector : unknown mode (%s)", __MaterialModeToString(mode) );
    }

    return 0;
}


int
StructuralMaterial :: giveStressStrainComponentIndOf(MatResponseForm form, MaterialMode mmode, int ind)
//
// this function returns index of reduced(if form == ReducedForm)
// or Full(if form==FullForm) stressStrain component in Full or reduced
// stressStrainVector acording to stressStrain mode of given gp.
//
{
    //MaterialMode mode  = gp -> giveMaterialMode ();

    if ( form == ReducedForm ) {
        switch ( mmode ) {
        case _3dMat:
            return ind;

        case _PlaneStress:
            if ( ind == 1 ) {
                return 1;
            } else if ( ind == 2 ) {
                return 2;
            } else if ( ind == 3 ) {
                return 6;
            }

            break;
        case _PlaneStrain:
            if ( ind == 1 ) {
                return 1;
            } else if ( ind == 2 ) {
                return 2;
            } else if ( ind == 3 ) {
                return 3;
            } else if ( ind == 4 ) {
                return 6;
            }

            break;
        case _1dMat:
            if ( ind == 1 ) {
                return 1;
            }

            break;
        case _2dPlateLayer:
        case _3dShellLayer:
            if ( ind == 1 ) {
                return 1;
            } else if ( ind == 2 ) {
                return 2;
            } else if ( ind == 3 ) {
                return 4;
            } else if ( ind == 4 ) {
                return 5;
            } else if ( ind == 5 ) {
                return 6;
            }

            break;
        case _2dBeamLayer:
            if ( ind == 1 ) {
                return 1;
            } else if ( ind == 2 ) {
                return 5;
            }

            break;
        case _2dPlate:
            if ( ind == 1 ) {
                return 7;
            } else if ( ind == 2 ) {
                return 8;
            } else if ( ind == 3 ) {
                return 12;
            } else if ( ind == 4 ) {
                return 5;
            } else if ( ind == 5 ) {
                return 4;
            }

            break;
        case _2dBeam:
            if ( ind == 1 ) {
                return 1;
            } else if ( ind == 2 ) {
                return 8;
            } else if ( ind == 3 ) {
                return 5;
            }

            break;
        case _1dFiber:
            if ( ind == 1 ) {
                return 1;
            } else if ( ind == 2 ) {
                return 5;
            } else if ( ind == 3 ) {
                return 6;
            }

            break;
        case _3dShell:
            if ( ind == 1 ) {
                return 1;
            } else if ( ind == 2 ) {
                return 2;
            } else if ( ind == 3 ) {
                return 6;
            } else if ( ind == 4 ) {
                return 7;
            } else if ( ind == 5 ) {
                return 8;
            } else if ( ind == 6 ) {
                return 12;
            } else if ( ind == 7 ) {
                return 5;
            } else if ( ind == 8 ) {
                return 4;
            }

            return ind;

        case _3dMatGrad:
            return ind;

        case _PlaneStressGrad:
            if ( ind == 1 ) {
                return 1;
            } else if ( ind == 2 ) {
                return 2;
            } else if ( ind == 3 ) {
                return 6;
            } else if ( ind == 4 ) {
                return 7;
            }

            break;

        case _PlaneStrainGrad:
            if ( ind == 1 ) {
                return 1;
            } else if ( ind == 2 ) {
                return 2;
            } else if ( ind == 3 ) {
                return 3;
            } else if ( ind == 4 ) {
                return 6;
            } else if ( ind == 5 ) {
                return 7;
            }

            break;

        case _1dMatGrad:
            if ( ind == 1 ) {
                return 1;
            } else if ( ind == 2 ) {
                return 7;
            }

            break;

        default:
            OOFEM_ERROR2("StructuralMaterial :: giveStressStrainComponentIndIn : unknown mode (%s)", __MaterialModeToString(mmode) );
        }

        return 0;
    } else if ( form == FullForm ) {
        switch ( mmode ) {
        case _3dMat:
            return ind;

        case _PlaneStress:
            if ( ind == 1 ) {
                return 1;
            } else if ( ind == 2 ) {
                return 2;
            } else if ( ind == 6 ) {
                return 3;
            }

            break;
        case _PlaneStrain:
            if ( ind == 1 ) {
                return 1;
            } else if ( ind == 2 ) {
                return 2;
            } else if ( ind == 3 ) {
                return 3;
            } else if ( ind == 6 ) {
                return 4;
            }

            break;
        case _1dMat:
            if ( ind == 1 ) {
                return 1;
            }

            break;
        case _2dPlateLayer:
        case _3dShellLayer:
            if ( ind == 1 ) {
                return 1;
            } else if ( ind == 2 ) {
                return 2;
            } else if ( ind == 4 ) {
                return 3;
            } else if ( ind == 5 ) {
                return 4;
            } else if ( ind == 6 ) {
                return 5;
            }

            break;
        case _2dBeamLayer:
            if ( ind == 1 ) {
                return 1;
            } else if ( ind == 5 ) {
                return 2;
            }

            break;
        case _2dPlate:
            if ( ind == 7 ) {
                return 1;
            } else if ( ind == 8 ) {
                return 2;
            } else if ( ind == 12 ) {
                return 3;
            } else if ( ind == 5 ) {
                return 4;
            } else if ( ind == 4 ) {
                return 5;
            }

            break;
        case _2dBeam:
            if ( ind == 1 ) {
                return 1;
            } else if ( ind == 8 ) {
                return 2;
            } else if ( ind == 5 ) {
                return 3;
            }

            break;
        case _1dFiber:
            if ( ind == 1 ) {
                return 1;
            } else if ( ind == 5 ) {
                return 2;
            } else if ( ind == 6 ) {
                return 3;
            }

            break;
        case _3dShell:
            if ( ind == 1 ) {
                return 1;
            } else if ( ind == 2 ) {
                return 2;
            } else if ( ind == 6 ) {
                return 3;
            } else if ( ind == 7 ) {
                return 4;
            } else if ( ind == 8 ) {
                return 5;
            } else if ( ind == 12 ) {
                return 6;
            } else if ( ind == 5 ) {
                return 7;
            } else if ( ind == 4 ) {
                return 8;
            }

            break;

        case _3dMatGrad:
            return ind;

        case _PlaneStressGrad:
            if ( ind == 1 ) {
                return 1;
            } else if ( ind == 2 ) {
                return 2;
            } else if ( ind == 6 ) {
                return 3;
            } else if ( ind == 7 ) {
                return 4;
            }

            break;

        case _PlaneStrainGrad:
            if ( ind == 1 ) {
                return 1;
            } else if ( ind == 2 ) {
                return 2;
            } else if ( ind == 3 ) {
                return 3;
            } else if ( ind == 6 ) {
                return 4;
            } else if ( ind == 7 ) {
                return 5;
            }

            break;

        case _1dMatGrad:
            if ( ind == 1 ) {
                return 1;
            } else if ( ind == 7 ) {
                return 2;
            }

            break;


        default:
            OOFEM_ERROR2("StructuralMaterial :: giveStressStrainComponentIndIn : unknown mode (%s)", __MaterialModeToString(mmode) );
        }

        return 0;
    } else {
        OOFEM_ERROR("StructuralMaterial :: giveStressStrainComponentIndIn : unknown form mode");
    }

    return 0;
}


void
StructuralMaterial :: giveStressStrainMask(IntArray &answer, MatResponseForm form,
                                           MaterialMode mmode) const
//
// this function returns mask of reduced(if form == ReducedForm)
// or Full(if form==FullForm) stressStrain vector in full or
// reduced StressStrainVector
// according to stressStrain mode of given gp.
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

    if ( form == ReducedForm ) {
        switch ( mmode ) {
        case _3dMat:
        case _3dMat_F:
            answer.resize(6);
            for ( i = 1; i <= 6; i++ ) {
                answer.at(i) = i;
            }

            break;
        case _PlaneStress:
            answer.resize(3);
            answer.at(1) = 1;
            answer.at(2) = 2;
            answer.at(3) = 6;
            break;
        case _PlaneStrain:
            answer.resize(4);
            answer.at(1) = 1;
            answer.at(2) = 2;
            answer.at(3) = 3;
            answer.at(4) = 6;
            break;
        case _1dMat:
            answer.resize(1);
            answer.at(1) = 1;
            break;
        case _2dPlateLayer:
        case _3dShellLayer:
            answer.resize(5);
            answer.at(1) = 1;
            answer.at(2) = 2;
            answer.at(3) = 4;
            answer.at(4) = 5;
            answer.at(5) = 6;
            break;
        case _2dBeamLayer:
            answer.resize(2);
            answer.at(1) = 1;
            answer.at(2) = 5;
            break;
        case _2dPlate:
            answer.resize(5);
            answer.at(1) = 7;
            answer.at(2) = 8;
            answer.at(3) = 12;
            answer.at(4) = 5;
            answer.at(5) = 4;
            break;
        case _2dBeam:
            answer.resize(3);
            answer.at(1) = 1;
            answer.at(2) = 8;
            answer.at(3) = 5;
            break;
        case _3dBeam:
            answer.resize(6);
            answer.at(1) = 1;
            answer.at(2) = 5;
            answer.at(3) = 6;
            answer.at(4) = 7;
            answer.at(5) = 8;
            answer.at(6) = 9;
            break;
        case _1dFiber:
            answer.resize(3);
            answer.at(1) = 1;
            answer.at(2) = 5;
            answer.at(3) = 6;
            break;
        case _3dShell:
            answer.resize(8);
            answer.at(1) = 1;
            answer.at(2) = 2;
            answer.at(3) = 6;
            answer.at(4) = 7;
            answer.at(5) = 8;
            answer.at(6) = 12;
            answer.at(7) = 5;
            answer.at(8) = 4;
            break;
        case _3dMatGrad:
            answer.resize(7);
            for ( i = 1; i <= 7; i++ ) {
                answer.at(i) = i;
            }

            break;
        case _PlaneStressGrad:
            answer.resize(4);
            answer.at(1) = 1;
            answer.at(2) = 2;
            answer.at(3) = 6;
            answer.at(4) = 7;
            break;
        case _PlaneStrainGrad:
            answer.resize(5);
            answer.at(1) = 1;
            answer.at(2) = 2;
            answer.at(3) = 3;
            answer.at(4) = 6;
            answer.at(5) = 7;
            break;
        case _1dMatGrad:
            answer.resize(2);
            answer.at(1) = 1;
            answer.at(2) = 7;
            break;
        default:
            OOFEM_ERROR2("StructuralMaterial :: giveStressStrainMask : unknown mode (%s)", __MaterialModeToString(mmode) );
        }
    } else if ( form == FullForm ) {
        switch ( mmode ) {
        case _3dMat:
            answer.resize(6);
            answer.zero();
            for ( i = 1; i <= 6; i++ ) {
                answer.at(i) = i;
            }

            break;
        case _PlaneStress:
            answer.resize(6);
            answer.zero();
            answer.at(1) = 1;
            answer.at(2) = 2;
            answer.at(6) = 3;
            break;
        case _PlaneStrain:
            answer.resize(6);
            answer.zero();
            answer.at(1) = 1;
            answer.at(2) = 2;
            answer.at(3) = 3;
            answer.at(6) = 4;
            break;
        case _1dMat:
            answer.resize(6);
            answer.zero();
            answer.at(1) = 1;
            break;
        case _2dPlateLayer:
        case _3dShellLayer:
            answer.resize(6);
            answer.zero();
            answer.at(1) = 1;
            answer.at(2) = 2;
            answer.at(4) = 3;
            answer.at(5) = 4;
            answer.at(6) = 5;
            break;
        case _2dBeamLayer:
            answer.resize(6);
            answer.zero();
            answer.at(1) = 1;
            answer.at(5) = 2;
            break;
        case _2dPlate:
            answer.resize(12);
            answer.zero();
            answer.at(7) = 1;
            answer.at(8) = 2;
            answer.at(12) = 3;
            answer.at(5) = 4;
            answer.at(4) = 5;
            break;
        case _2dBeam:
            answer.resize(12);
            answer.zero();
            answer.at(1) = 1;
            answer.at(8) = 2;
            answer.at(5) = 3;
            break;
        case _3dBeam:
            answer.resize(6);
            answer.at(1) = 1;
            answer.at(5) = 2;
            answer.at(6) = 3;
            answer.at(7) = 4;
            answer.at(8) = 5;
            answer.at(9) = 6;
            break;
        case _1dFiber:
            answer.resize(6);
            answer.zero();
            answer.at(1) = 1;
            answer.at(5) = 2;
            answer.at(6) = 3;
            break;
        case _3dShell:
            answer.resize(12);
            answer.zero();
            answer.at(1) = 1;
            answer.at(2) = 2;
            answer.at(6) = 3;
            answer.at(7) = 4;
            answer.at(8) = 5;
            answer.at(12) = 6;
            answer.at(5) = 7;
            answer.at(4) = 8;
            break;
        case _3dMatGrad:
            answer.resize(7);
            answer.zero();
            for ( i = 1; i <= 7; i++ ) {
                answer.at(i) = i;
            }

            break;
        case _PlaneStressGrad:
            answer.resize(7);
            answer.zero();
            answer.at(1) = 1;
            answer.at(2) = 2;
            answer.at(6) = 3;
            answer.at(7) = 4;
            break;
        case _PlaneStrainGrad:
            answer.resize(7);
            answer.zero();
            answer.at(1) = 1;
            answer.at(2) = 2;
            answer.at(3) = 3;
            answer.at(6) = 4;
            answer.at(7) = 5;
            break;
        case _1dMatGrad:
            answer.resize(7);
            answer.zero();
            answer.at(1) = 1;
            answer.at(7) = 2;
            break;

        default:
            OOFEM_ERROR2("StructuralMaterial :: giveStressStrainMask : unknown mode (%s)", __MaterialModeToString(mmode) );
        }
    } else {
        OOFEM_ERROR("StructuralMaterial :: giveStressStrainMask : unknown form mode");
    }
}


void
StructuralMaterial :: reduceToPlaneStressStiffMtrx(FloatMatrix &answer,
                                                   MatResponseForm form, GaussPoint *gp,
                                                   FloatMatrix &stiffMtrx3d) const
//
// returns receiver's 2dPlaneStressMtrx constructed from
// stiffMtrx3d (general 3dMatrialStiffnessMatrix)
// (2dPlaneStres ==> sigma_z = tau_xz = tau_yz = 0.)
// This method works for general 3d stiff matrix
//
{
    FloatMatrix inv3d, invAnswer, reducedAnswer;

    // check if stiffMtrx is proper
    if ( ( stiffMtrx3d.isSquare() ) && ( stiffMtrx3d.giveNumberOfRows() == 6 ) ) {
        inv3d.beInverseOf(stiffMtrx3d);

        invAnswer.resize(3, 3);

        invAnswer.at(1, 1) = inv3d.at(1, 1);
        invAnswer.at(1, 2) = inv3d.at(1, 2);
        invAnswer.at(1, 3) = inv3d.at(1, 6);

        invAnswer.at(2, 1) = inv3d.at(2, 1);
        invAnswer.at(2, 2) = inv3d.at(2, 2);
        invAnswer.at(2, 3) = inv3d.at(2, 6);

        invAnswer.at(3, 1) = inv3d.at(6, 1);
        invAnswer.at(3, 2) = inv3d.at(6, 2);
        invAnswer.at(3, 3) = inv3d.at(6, 6);

        reducedAnswer.beInverseOf(invAnswer);

        if ( form == ReducedForm ) {
            answer = reducedAnswer;
        } else {
            answer.resize(6, 6);
            answer.zero();

            answer.at(1, 1) = reducedAnswer.at(1, 1);
            answer.at(1, 2) = reducedAnswer.at(1, 2);
            answer.at(1, 6) = reducedAnswer.at(1, 3);
            answer.at(2, 1) = reducedAnswer.at(2, 1);
            answer.at(2, 2) = reducedAnswer.at(2, 2);
            answer.at(2, 6) = reducedAnswer.at(2, 3);
            answer.at(6, 6) = reducedAnswer.at(3, 3);
            answer.at(6, 1) = reducedAnswer.at(3, 1);
            answer.at(6, 2) = reducedAnswer.at(3, 2);
        }

        return;
    } else {
        OOFEM_ERROR("StructuralMaterial :: reduceToPlaneStressStiffMtrx : stiffMtrx size mismatch");
    }
}


void
StructuralMaterial :: reduceToPlaneStrainStiffMtrx(FloatMatrix &answer,
                                                   MatResponseForm form, GaussPoint *gp,
                                                   FloatMatrix &stiffMtrx3d) const
//
// returns receiver's 2dPlaneStrainMtrx constructed from
// general 3dMatrialStiffnessMatrix
// (2dPlaneStrain ==> eps_z = gamma_xz = gamma_yz = 0.)
// but we take ez, SigmaZ into account.
//
{
    // check if stiffMtrx is proper
    if ( ( stiffMtrx3d.isSquare() ) && ( stiffMtrx3d.giveNumberOfRows() == 6 ) ) {
        if ( form == ReducedForm ) {
            answer.resize(4, 4);
            answer.zero();

            answer.at(1, 1) = stiffMtrx3d.at(1, 1);
            answer.at(1, 2) = stiffMtrx3d.at(1, 2);
            answer.at(1, 4) = stiffMtrx3d.at(1, 6);

            answer.at(2, 1) = stiffMtrx3d.at(2, 1);
            answer.at(2, 2) = stiffMtrx3d.at(2, 2);
            answer.at(2, 4) = stiffMtrx3d.at(2, 6);

            answer.at(3, 1) = stiffMtrx3d.at(3, 1);
            answer.at(3, 2) = stiffMtrx3d.at(3, 2);
            answer.at(3, 4) = stiffMtrx3d.at(3, 6);

            answer.at(4, 1) = stiffMtrx3d.at(6, 1);
            answer.at(4, 2) = stiffMtrx3d.at(6, 2);
            answer.at(4, 4) = stiffMtrx3d.at(6, 6);
        } else {
            answer.resize(6, 6);
            answer.zero();

            answer.at(1, 1) = stiffMtrx3d.at(1, 1);
            answer.at(1, 2) = stiffMtrx3d.at(1, 2);
            answer.at(1, 6) = stiffMtrx3d.at(1, 6);

            answer.at(2, 1) = stiffMtrx3d.at(2, 1);
            answer.at(2, 2) = stiffMtrx3d.at(2, 2);
            answer.at(2, 6) = stiffMtrx3d.at(2, 6);

            answer.at(3, 1) = stiffMtrx3d.at(3, 1);
            answer.at(3, 2) = stiffMtrx3d.at(3, 2);
            answer.at(3, 6) = stiffMtrx3d.at(3, 6);

            answer.at(6, 6) = stiffMtrx3d.at(3, 3);
            answer.at(6, 1) = stiffMtrx3d.at(6, 1);
            answer.at(6, 2) = stiffMtrx3d.at(6, 2);
        }

        return;
    } else {
        OOFEM_ERROR("StructuralMaterial :: reduceToPlaneStrainStiffMtrx :: stiffMtrx size mismatch");
    }
}


void
StructuralMaterial :: reduceTo1dStressStiffMtrx(FloatMatrix &answer,
                                                MatResponseForm form, GaussPoint *gp,
                                                FloatMatrix &stiffMtrx3d) const
//
//
// returns receiver's 1dMaterialStiffnessMAtrix constructed from
// general 3dMatrialStiffnessMatrix
// (1d case ==> sigma_y = sigma_z = tau_yz = tau_zx = tau_xy  = 0.)
// This method works only if 3dMateriallStiffnessMatrix
// has two 3x3 independent blocks
{
    FloatMatrix m3d11, inv3d;
    double val11;

    // check if stiffMtrx is proper
    if ( ( stiffMtrx3d.isSquare() ) && ( stiffMtrx3d.giveNumberOfRows() == 6 ) ) {
        inv3d.beInverseOf(stiffMtrx3d);
        val11 = inv3d.at(1, 1);

        if ( form == ReducedForm ) {
            answer.resize(1, 1);

            answer.at(1, 1) = 1. / val11;
        } else {
            answer.resize(6, 6);
            answer.zero();

            answer.at(1, 1) = 1. / val11;
        }

        return;
    } else {
        OOFEM_ERROR("StructuralMaterial :: reduceTo1dStressStiffMtrx:: stiffMtrx3d size mismatch");
    }
}


void
StructuralMaterial :: reduceTo2dPlateLayerStiffMtrx(FloatMatrix &answer,
                                                    MatResponseForm form,
                                                    GaussPoint *gp,
                                                    FloatMatrix &stiffMtrx3d) const
//
// return material stiffness matrix for derived types of stressStreinState
// assumption sigma_z = 0.
//
// General strain vector has one of the following forms:
// 1) strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
//
{
    MaterialMode mode = gp->giveMaterialMode();
    FloatMatrix invMat3d, invMatLayer(5, 5), matLayer;

    int i, j;

    if ( !( ( mode == _2dPlateLayer ) || ( mode == _3dShellLayer ) ) ) {
        _error("ReduceTo2dPlateLayerStiffMtrx : unsupported mode");
    }


    // check if stiffMtrx is proper
    if ( ( stiffMtrx3d.isSquare() ) && ( stiffMtrx3d.giveNumberOfRows() == 6 ) ) {
        invMat3d.beInverseOf(stiffMtrx3d);



        for ( i = 1; i <= 2; i++ ) {
            for ( j = 1; j <= 2; j++ ) {
                invMatLayer.at(i, j) = invMat3d.at(i, j);
            }
        }

        for ( i = 4; i <= 6; i++ ) {
            for ( j = 4; j <= 6; j++ ) {
                invMatLayer.at(i - 1, j - 1) = invMat3d.at(i, j);
            }
        }

        for ( i = 1; i <= 2; i++ ) {
            for ( j = 4; j <= 6; j++ ) {
                invMatLayer.at(i, j - 1) = invMat3d.at(i, j);
                invMatLayer.at(j - 1, i) = invMat3d.at(j, i);
            }
        }

        matLayer.beInverseOf(invMatLayer);

        if ( form == ReducedForm ) {
            answer = matLayer;
        } else {
            answer.resize(6, 6);
            answer.zero();

            for ( i = 1; i <= 2; i++ ) {
                for ( j = 1; j <= 2; j++ ) {
                    answer.at(i, j) = matLayer.at(i, j);
                }
            }

            for ( i = 4; i <= 6; i++ ) {
                for ( j = 4; j <= 6; j++ ) {
                    answer.at(i, j) = matLayer.at(i - 1, j - 1);
                }
            }

            for ( i = 1; i <= 2; i++ ) {
                for ( j = 4; j <= 6; j++ ) {
                    answer.at(i, j) = matLayer.at(i, j - 1);
                    answer.at(j, i) = matLayer.at(j - 1, i);
                }
            }
        }

        return;
    } else {
        OOFEM_ERROR("StructuralMaterial :: reduceTo2dPlateLayerStiffMtrx : stiffMtrx size mismatch");
    }
}


void
StructuralMaterial :: reduceTo3dShellLayerStiffMtrx(FloatMatrix &answer,
                                                    MatResponseForm form,
                                                    GaussPoint *gp,
                                                    FloatMatrix &stiffMtrx3d) const
//
// return material stiffness matrix for derived types of stressStreinState
// assumption sigma_z = 0.
//
// General strain vector has one of the following forms:
// 1) strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}

{
    this->reduceTo2dPlateLayerStiffMtrx(answer, form, gp, stiffMtrx3d);
}


void
StructuralMaterial :: reduceTo2dBeamLayerStiffMtrx(FloatMatrix &answer,
                                                   MatResponseForm form,
                                                   GaussPoint *gp,
                                                   FloatMatrix &stiffMtrx3d) const
//
// return material stiffness matrix for derived types of stressStreinState
// assumption sigma_z = 0.
//
// General strain vector has one of the following forms:
// 1) strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
//
{
    MaterialMode mode = gp->giveMaterialMode();
    FloatMatrix invMat3d, invMatLayer(2, 2), matLayer;

    if ( mode != _2dBeamLayer ) {
        OOFEM_ERROR("StructuralMaterial :: ReduceTo2dBeamLayerStiffMtrx : unsupported mode");
    }

    if ( ( stiffMtrx3d.isSquare() ) && ( stiffMtrx3d.giveNumberOfRows() == 6 ) ) {
        invMat3d.beInverseOf(stiffMtrx3d);

        invMatLayer.at(1, 1) = invMat3d.at(1, 1);
        invMatLayer.at(1, 2) = invMat3d.at(1, 5);
        invMatLayer.at(2, 1) = invMat3d.at(5, 1);
        invMatLayer.at(2, 2) = invMat3d.at(5, 5);

        matLayer.beInverseOf(invMatLayer);

        if ( form == ReducedForm ) {
            answer = matLayer;
        } else {
            answer.resize(6, 6);
            answer.zero();

            answer.at(1, 1) = matLayer.at(1, 1);
            answer.at(1, 5) = matLayer.at(1, 2);
            answer.at(5, 1) = matLayer.at(2, 1);
            answer.at(5, 5) = matLayer.at(2, 2);
        }

        return;
    } else {
        OOFEM_ERROR("StructuralMaterial :: reduceTo2dBeamLayerStiffMtrx: stiffMtrx3d size mismatch");
    }
}


void
StructuralMaterial :: reduceTo1dFiberStiffMtrx(FloatMatrix &answer,
                                               MatResponseForm form,
                                               GaussPoint *gp,
                                               FloatMatrix &stiffMtrx3d) const
//
// return material stiffness matrix for derived types of stressStreinState
// assumption sigma_z = 0.
//
// General strain vector has one of the following forms:
// 1) strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
//
{
    MaterialMode mode = gp->giveMaterialMode();
    FloatMatrix invMat3d, invMatLayer(3, 3), matLayer;

    if ( mode != _1dFiber ) {
        OOFEM_ERROR("StructuralMaterial :: reduceTo1dFiberStiffMtrx : unsupported mode");
    }

    if ( ( stiffMtrx3d.isSquare() ) && ( stiffMtrx3d.giveNumberOfRows() == 6 ) ) {
        invMat3d.beInverseOf(stiffMtrx3d);

        invMatLayer.at(1, 1) = invMat3d.at(1, 1);
        invMatLayer.at(1, 2) = invMat3d.at(1, 5);
        invMatLayer.at(1, 3) = invMat3d.at(1, 6);
        invMatLayer.at(2, 1) = invMat3d.at(5, 1);
        invMatLayer.at(2, 2) = invMat3d.at(5, 5);
        invMatLayer.at(2, 3) = invMat3d.at(5, 6);
        invMatLayer.at(3, 1) = invMat3d.at(6, 1);
        invMatLayer.at(3, 2) = invMat3d.at(6, 5);
        invMatLayer.at(3, 3) = invMat3d.at(6, 6);

        matLayer.beInverseOf(invMatLayer);

        if ( form == ReducedForm ) {
            answer = matLayer;
        } else {
            answer.resize(6, 6);
            answer.zero();

            answer.at(1, 1) = matLayer.at(1, 1);
            answer.at(1, 5) = matLayer.at(1, 2);
            answer.at(1, 6) = matLayer.at(1, 3);
            answer.at(5, 1) = matLayer.at(2, 1);
            answer.at(5, 5) = matLayer.at(2, 2);
            answer.at(5, 6) = matLayer.at(2, 3);
            answer.at(6, 1) = matLayer.at(3, 1);
            answer.at(6, 5) = matLayer.at(3, 2);
            answer.at(6, 6) = matLayer.at(3, 3);
        }

        return;
    } else {
        OOFEM_ERROR("StructuralMaterial :: reduceTo1dFiberStiffMtrx: stiffMtrx3d size mismatch");
    }
}


//
// Compliance reduction functions
//
void
StructuralMaterial :: reduceToPlaneStressComplMtrx(FloatMatrix &answer,
                                                   MatResponseForm form, GaussPoint *gp,
                                                   FloatMatrix &complMtrx3d) const
//
// returns receiver's 2dPlaneComplMtrx constructed from
// complMtrx3d (general 3dMatrialComplianceMatrix)
// (2dPlaneStres ==> sigma_z = tau_xz = tau_yz = 0.)
// This method works for general 3d compl matrix
//
{
    // check if complMtrx is proper
    if ( ( complMtrx3d.isSquare() ) && ( complMtrx3d.giveNumberOfRows() == 6 ) ) {
        if ( form == ReducedForm ) {
            answer.resize(3, 3);
            answer.zero();

            answer.at(1, 1) = complMtrx3d.at(1, 1);
            answer.at(1, 2) = complMtrx3d.at(1, 2);
            answer.at(1, 3) = complMtrx3d.at(1, 6);
            answer.at(2, 1) = complMtrx3d.at(2, 1);
            answer.at(2, 2) = complMtrx3d.at(2, 2);
            answer.at(2, 3) = complMtrx3d.at(2, 6);
            answer.at(3, 3) = complMtrx3d.at(6, 6);
            answer.at(3, 1) = complMtrx3d.at(6, 1);
            answer.at(3, 2) = complMtrx3d.at(6, 2);
        } else {
            answer.resize(6, 6);
            answer.zero();

            answer.at(1, 1) = complMtrx3d.at(1, 1);
            answer.at(1, 2) = complMtrx3d.at(1, 2);
            answer.at(1, 6) = complMtrx3d.at(1, 6);
            answer.at(2, 1) = complMtrx3d.at(2, 1);
            answer.at(2, 2) = complMtrx3d.at(2, 2);
            answer.at(2, 6) = complMtrx3d.at(2, 6);
            answer.at(6, 6) = complMtrx3d.at(6, 6);
            answer.at(6, 1) = complMtrx3d.at(6, 1);
            answer.at(6, 2) = complMtrx3d.at(6, 2);
        }

        return;
    } else {
        OOFEM_ERROR("StructuralMaterial :: reduceToPlaneStressComplMtrx : complMtrx size mismatch");
    }
}


void
StructuralMaterial :: reduceToPlaneStrainComplMtrx(FloatMatrix &answer,
                                                   MatResponseForm form, GaussPoint *gp,
                                                   FloatMatrix &complMtrx3d) const
//
// returns receiver's 2dPlaneStrainMtrx constructed from
// general 3dMatrialComplianceMatrix
// (2dPlaneStrain ==> eps_z = gamma_xz = gamma_yz = 0.)
//
{
    int i, j;
    FloatMatrix inv3d, invAnswer(3, 3), reducedAnswer;

    // check if complMtrx is proper
    if ( ( complMtrx3d.isSquare() ) && ( complMtrx3d.giveNumberOfRows() == 6 ) ) {
        inv3d.beInverseOf(complMtrx3d);

        invAnswer.at(1, 1) = inv3d.at(1, 1);
        invAnswer.at(1, 2) = inv3d.at(1, 2);
        invAnswer.at(1, 3) = inv3d.at(1, 6);

        invAnswer.at(2, 1) = inv3d.at(2, 1);
        invAnswer.at(2, 2) = inv3d.at(2, 2);
        invAnswer.at(2, 3) = inv3d.at(2, 6);

        invAnswer.at(3, 1) = inv3d.at(6, 1);
        invAnswer.at(3, 2) = inv3d.at(6, 2);
        invAnswer.at(3, 3) = inv3d.at(6, 6);

        reducedAnswer.beInverseOf(invAnswer);

        if ( form == ReducedForm ) {
            answer.resize(4, 4);
            answer.zero();

            for ( i = 1; i <= 2; i++ ) {
                for ( j = 1; j <= 2; j++ ) {
                    answer.at(i, j) = reducedAnswer.at(i, j);
                }
            }

            answer.at(1, 4) = reducedAnswer.at(1, 3);
            answer.at(2, 4) = reducedAnswer.at(2, 3);
            answer.at(4, 1) = reducedAnswer.at(3, 1);
            answer.at(4, 2) = reducedAnswer.at(3, 2);
            answer.at(4, 4) = reducedAnswer.at(3, 3);
        } else {
            answer.resize(6, 6);
            answer.zero();

            answer.at(1, 1) = reducedAnswer.at(1, 1);
            answer.at(1, 2) = reducedAnswer.at(1, 2);
            answer.at(1, 6) = reducedAnswer.at(1, 3);
            answer.at(2, 1) = reducedAnswer.at(2, 1);
            answer.at(2, 2) = reducedAnswer.at(2, 2);
            answer.at(2, 6) = reducedAnswer.at(2, 3);
            answer.at(6, 6) = reducedAnswer.at(3, 3);
            answer.at(6, 1) = reducedAnswer.at(3, 1);
            answer.at(6, 2) = reducedAnswer.at(3, 2);
        }

        return;
    } else {
        OOFEM_ERROR("StructuralMaterial :: reduceToPlaneStrainComplMtrx :: complMtrx size mismatch");
    }
}


void
StructuralMaterial :: reduceTo1dStressComplMtrx(FloatMatrix &answer,
                                                MatResponseForm form, GaussPoint *gp,
                                                FloatMatrix &complMtrx3d) const
//
//
// returns receiver's 1dMaterialComplianceMAtrix constructed from
// general 3dMatrialComplianceMatrix
// (1d case ==> sigma_y = sigma_z = tau_yz = tau_zx = tau_xy  = 0.)
{
    // check if complMtrx is proper
    if ( ( complMtrx3d.isSquare() ) && ( complMtrx3d.giveNumberOfRows() == 6 ) ) {
        if ( form == ReducedForm ) {
            answer.resize(1, 1);

            answer.at(1, 1) = complMtrx3d.at(1, 1);
        } else {
            answer.resize(6, 6);
            answer.zero();

            answer.at(1, 1) = complMtrx3d.at(1, 1);
        }

        return;
    } else {
        OOFEM_ERROR("StructuralMaterial :: reduceTo1dStressComplMtrx:: complMtrx3d size mismatch");
    }
}



void
StructuralMaterial :: reduceTo2dPlateLayerComplMtrx(FloatMatrix &answer,
                                                    MatResponseForm form,
                                                    GaussPoint *gp,
                                                    FloatMatrix &complMtrx3d) const
//
// return material compliance matrix for derived types of stressStreinState
// assumption sigma_z = 0.
//
// General strain vector has one of the following forms:
// 1) strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
//
{
    MaterialMode mode = gp->giveMaterialMode();
    int i, j;

    if ( !( ( mode == _2dPlateLayer ) || ( mode == _3dShellLayer ) ) ) {
        OOFEM_ERROR("StructuralMaterial :: reduceTo2dPlateLayerComplMtrx : unsupported mode");
    }


    // check if complMtrx is proper
    if ( ( complMtrx3d.isSquare() ) && ( complMtrx3d.giveNumberOfRows() == 6 ) ) {
        if ( form == ReducedForm ) {
            answer.resize(5, 5);
            answer.zero();

            for ( i = 1; i <= 2; i++ ) {
                for ( j = 1; j <= 2; j++ ) {
                    answer.at(i, j) = complMtrx3d.at(i, j);
                }
            }

            for ( i = 4; i <= 6; i++ ) {
                for ( j = 4; j <= 6; j++ ) {
                    answer.at(i - 1, j - 1) = complMtrx3d.at(i, j);
                }
            }

            for ( i = 1; i <= 2; i++ ) {
                for ( j = 4; j <= 6; j++ ) {
                    answer.at(i, j - 1) = complMtrx3d.at(i, j);
                    answer.at(j - 1, i) = complMtrx3d.at(j, i);
                }
            }
        } else {
            answer.resize(6, 6);
            answer.zero();

            for ( i = 1; i <= 2; i++ ) {
                for ( j = 1; j <= 2; j++ ) {
                    answer.at(i, j) = complMtrx3d.at(i, j);
                }
            }

            for ( i = 4; i <= 6; i++ ) {
                for ( j = 4; j <= 6; j++ ) {
                    answer.at(i, j) = complMtrx3d.at(i, j);
                }
            }

            for ( i = 1; i <= 2; i++ ) {
                for ( j = 4; j <= 6; j++ ) {
                    answer.at(i, j) = complMtrx3d.at(i, j);
                    answer.at(j, i) = complMtrx3d.at(j, i);
                }
            }
        }

        return;
    } else {
        OOFEM_ERROR("StructuralMaterial :: reduceTo2dPlateLayerComplMtrx : stiffMtrx size mismatch");
    }
}


void
StructuralMaterial :: reduceTo3dShellLayerComplMtrx(FloatMatrix &answer,
                                                    MatResponseForm form,
                                                    GaussPoint *gp,
                                                    FloatMatrix &complMtrx3d) const
//
// return material compliance matrix for derived types of stressStreinState
// assumption sigma_z = 0.
//
// General strain vector has one of the following forms:
// 1) strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}

{
    this->reduceTo2dPlateLayerComplMtrx(answer, form, gp, complMtrx3d);
}



void
StructuralMaterial :: reduceTo2dBeamLayerComplMtrx(FloatMatrix &answer,
                                                   MatResponseForm form,
                                                   GaussPoint *gp,
                                                   FloatMatrix &complMtrx3d) const
//
// return material compliance matrix for derived types of stressStreinState
// assumption sigma_z = 0.
//
// General strain vector has one of the following forms:
// 1) strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
//
{
    MaterialMode mode = gp->giveMaterialMode();

    if ( mode != _2dBeamLayer ) {
        OOFEM_ERROR("StructuralMaterial :: reduceTo2dBeamLayerComplMtrx : unsupported mode");
    }

    if ( ( complMtrx3d.isSquare() ) && ( complMtrx3d.giveNumberOfRows() == 6 ) ) {
        if ( form == ReducedForm ) {
            answer.resize(2, 2);

            answer.at(1, 1) = complMtrx3d.at(1, 1);
            answer.at(1, 2) = complMtrx3d.at(1, 5);
            answer.at(2, 1) = complMtrx3d.at(5, 1);
            answer.at(2, 2) = complMtrx3d.at(5, 5);
        } else {
            answer.resize(6, 6);
            answer.zero();

            answer.at(1, 1) = complMtrx3d.at(1, 1);
            answer.at(1, 5) = complMtrx3d.at(1, 5);
            answer.at(5, 1) = complMtrx3d.at(5, 1);
            answer.at(5, 5) = complMtrx3d.at(5, 5);
        }

        return;
    } else {
        OOFEM_ERROR("StructuralMaterial :: reduceTo2dBeamLayerStiffMtrx: stiffMtrx3d size mismatch");
    }
}


void
StructuralMaterial :: reduceTo1dFiberComplMtrx(FloatMatrix &answer,
                                               MatResponseForm form,
                                               GaussPoint *gp,
                                               FloatMatrix &complMtrx3d) const
//
// return material compliance matrix for derived types of stressStreinState
// assumption sigma_z = 0.
//
// General strain vector has one of the following forms:
// 1) strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
//
{
    MaterialMode mode = gp->giveMaterialMode();

    if ( mode != _1dFiber ) {
        OOFEM_ERROR("StructuralMaterial :: reduceTo1dFiberComplMtrx : unsupported mode");
    }

    if ( ( complMtrx3d.isSquare() ) && ( complMtrx3d.giveNumberOfRows() == 6 ) ) {
        if ( form == ReducedForm ) {
            answer.resize(3, 3);

            answer.at(1, 1) = complMtrx3d.at(1, 1);
            answer.at(1, 2) = complMtrx3d.at(1, 5);
            answer.at(1, 3) = complMtrx3d.at(1, 6);
            answer.at(2, 1) = complMtrx3d.at(5, 1);
            answer.at(2, 2) = complMtrx3d.at(5, 5);
            answer.at(2, 3) = complMtrx3d.at(5, 6);
            answer.at(3, 1) = complMtrx3d.at(6, 1);
            answer.at(3, 2) = complMtrx3d.at(6, 5);
            answer.at(3, 3) = complMtrx3d.at(6, 6);
        } else {
            answer.resize(6, 6);
            answer.zero();

            answer.at(1, 1) = complMtrx3d.at(1, 1);
            answer.at(1, 5) = complMtrx3d.at(1, 5);
            answer.at(1, 6) = complMtrx3d.at(1, 6);
            answer.at(5, 1) = complMtrx3d.at(5, 1);
            answer.at(5, 5) = complMtrx3d.at(5, 5);
            answer.at(5, 6) = complMtrx3d.at(5, 6);
            answer.at(6, 1) = complMtrx3d.at(6, 1);
            answer.at(6, 5) = complMtrx3d.at(6, 5);
            answer.at(6, 6) = complMtrx3d.at(6, 6);
        }

        return;
    } else {
        OOFEM_ERROR("StructuralMaterial :: reduceTo1dFiberComplMtrx: stiffMtrx3d size mismatch");
    }
}

//
//
//


void
StructuralMaterial :: givePlaneStressStiffMtrx(FloatMatrix &answer,
                                               MatResponseForm form, MatResponseMode mode,
                                               GaussPoint *gp,
                                               TimeStep *atTime)
//
// returns Mat stiffness for PlaneStress
//
{
    FloatMatrix m3d;

    this->give3dMaterialStiffnessMatrix(m3d, FullForm, mode, gp, atTime);
    this->reduceToPlaneStressStiffMtrx(answer, form, gp, m3d);
}

void
StructuralMaterial :: givePlaneStrainStiffMtrx(FloatMatrix &answer,
                                               MatResponseForm form, MatResponseMode mode,
                                               GaussPoint *gp,
                                               TimeStep *atTime)
//
// return material stiffness matrix for PlaneStrain mode
//
{
    FloatMatrix m3d;

    this->give3dMaterialStiffnessMatrix(m3d, FullForm, mode, gp, atTime);
    this->reduceToPlaneStrainStiffMtrx(answer, form, gp, m3d);
}

void
StructuralMaterial :: give1dStressStiffMtrx(FloatMatrix &answer,
                                            MatResponseForm form, MatResponseMode mode,
                                            GaussPoint *gp,
                                            TimeStep *atTime)
//
// return material stiffness matrix for 1d stress strain mode
//
{
    FloatMatrix m3d;

    this->give3dMaterialStiffnessMatrix(m3d, FullForm, mode, gp, atTime);
    this->reduceTo1dStressStiffMtrx(answer, form, gp, m3d);
}


void
StructuralMaterial :: give2dBeamLayerStiffMtrx(FloatMatrix &answer,
                                               MatResponseForm form, MatResponseMode mode,
                                               GaussPoint *gp,
                                               TimeStep *atTime)
//
// return material stiffness matrix for2dBeamLayer mode
//
{
    FloatMatrix m3d;

    this->give3dMaterialStiffnessMatrix(m3d, FullForm, mode, gp, atTime);
    this->reduceTo2dBeamLayerStiffMtrx(answer, form, gp, m3d);
}


void
StructuralMaterial :: give2dPlateLayerStiffMtrx(FloatMatrix &answer,
                                                MatResponseForm form, MatResponseMode mode,
                                                GaussPoint *gp,
                                                TimeStep *atTime)
//
// return material stiffness matrix for 2dPlateLayer
//
{
    FloatMatrix m3d;

    this->give3dMaterialStiffnessMatrix(m3d, FullForm, mode, gp, atTime);
    this->reduceTo2dPlateLayerStiffMtrx(answer, form, gp, m3d);
}

void
StructuralMaterial :: give1dFiberStiffMtrx(FloatMatrix &answer,
                                           MatResponseForm form, MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *atTime)
//
// return material stiffness matrix for 2dPlateLayer
//
{
    FloatMatrix m3d;

    this->give3dMaterialStiffnessMatrix(m3d, FullForm, mode, gp, atTime);
    this->reduceTo1dFiberStiffMtrx(answer, form, gp, m3d);
}


void
StructuralMaterial :: give3dShellLayerStiffMtrx(FloatMatrix &answer,
                                                MatResponseForm form, MatResponseMode mode,
                                                GaussPoint *gp,
                                                TimeStep *atTime)
//
// returns material stiffness matrix for 3dShellLayer
//
{
    FloatMatrix m3d;

    this->give3dMaterialStiffnessMatrix(m3d, FullForm, mode, gp, atTime);
    this->reduceTo3dShellLayerStiffMtrx(answer, form, gp, m3d);
}


void
StructuralMaterial :: computePrincipalValues(FloatArray &answer, const FloatArray &s, stressStrainPrincMode mode)
//
// This function computes Principal values of strains or stresses.
// strains/stresses are stored in vector form in array s.
// Engineering notation is used.
//
// Problem size (3D/2D) is recognized automatically according to
// vector size.
// If size = 6 -> 3D problem, then array s contains:
//                            {Sxx,Syy,Szz,Syz,Szx,Sxy} if mode = stress_mode
//                            {Exx,Eyy,Ezz,GMyz,GMzx,GMxy} if mode = strain mode
// if size = 2 -> 2D problem, then array s contains:
//                            {Sxx,Syy,Sxy} if mode = stress_mode
//                            {Exx,Eyy,GMxy} if mode = strain mode
//
// mode      - principal strains
//           - principal stress
//           - principal deviatoric stress
//           - ..
//
// Return Values:
//
//    array sp -> principal strains or stresses
//
{
    int size = s.giveSize();
    double swap;
    int nonzeroFlag = 0;
    if ( !( ( size == 3 ) || ( size == 6 ) ) ) {
        OOFEM_ERROR("StructuralMaterial :: ComputePrincipalValues: Vector size mismatch");
    }

    if ( s.giveSize() == 3 ) {
        // 2D problem
        double ast, dst, D = 0.0;
        answer.resize(2);

        for ( int i = 1; i <= size; i++ ) {
            if ( fabs( s.at(i) ) > 1.e-20 ) {
                nonzeroFlag = 1;
            }
        }

        if ( nonzeroFlag == 0 ) {
            answer.zero();
            return;
        }

        ast = s.at(1) + s.at(2);
        dst = s.at(1) - s.at(2);
        if ( mode == principal_strain ) {
            D = dst * dst + s.at(3) * s.at(3);
        } else if ( mode == principal_stress ) {
            D = dst * dst + 4.0 * s.at(3) * s.at(3);
        } else {
            OOFEM_ERROR("StructuralMaterial :: ComputePrincipalValues: not supported");
        }

        if ( D < 0. ) {
            OOFEM_ERROR("StructuralMaterial :: ComputePrincipalValues: Imaginar roots ");
        }

        D = sqrt(D);
        answer.at(1) = 0.5 * ( ast - D );
        answer.at(2) = 0.5 * ( ast + D );

        // sort result
        if ( answer.at(1) > answer.at(2) ) {
            return;
        } else {
            swap = answer.at(1);
            answer.at(1) = answer.at(2);
            answer.at(2) = swap;
            return;
        }
    } else {
        // 3D problem
        double I1 = 0.0, I2 = 0.0, I3 = 0.0, help, s1, s2, s3;
        int i, j;

        for ( i = 1; i <= size; i++ ) {
            if ( fabs( s.at(i) ) > 1.e-20 ) {
                nonzeroFlag = 1;
            }
        }

        answer.resize(3);
        answer.zero();
        if ( nonzeroFlag == 0 ) {
            return;
        }



        if ( mode == principal_stress ) {
            I1 = s.at(1) + s.at(2) + s.at(3);
            I2 = s.at(1) * s.at(2) + s.at(2) * s.at(3) + s.at(3) * s.at(1) -
                 ( s.at(4) * s.at(4) + s.at(5) * s.at(5) + s.at(6) * s.at(6) );
            I3 = s.at(1) * s.at(2) * s.at(3) + 2. * s.at(4) * s.at(5) * s.at(6) -
                 ( s.at(1) * s.at(4) * s.at(4) + s.at(2) * s.at(5) * s.at(5) +
                  s.at(3) * s.at(6) * s.at(6) );
        } else if ( mode == principal_deviatoricstress ) {
            help = ( s.at(1) + s.at(2) + s.at(3) ) / 3.0;
            I1 = 0.;
            I2 = -( 1. / 6. ) * ( ( s.at(1) - s.at(2) ) * ( s.at(1) - s.at(2) ) + ( s.at(2) - s.at(3) ) * ( s.at(2) - s.at(3) ) +
                                 ( s.at(3) - s.at(1) ) * ( s.at(3) - s.at(1) ) ) - s.at(4) * s.at(4) - s.at(5) * s.at(5) -
                 s.at(6) * s.at(6);
            I3 = ( s.at(1) - help ) * ( s.at(2) - help ) * ( s.at(3) - help ) + 2. * s.at(4) * s.at(5) * s.at(6) -
                 s.at(5) * s.at(5) * ( s.at(2) - help ) - s.at(4) * s.at(4) * ( s.at(1) - help ) -
                 s.at(6) * s.at(6) * ( s.at(3) - help );
        } else if ( mode == principal_strain ) {
            I1 = s.at(1) + s.at(2) + s.at(3);
            I2 = s.at(1) * s.at(2) + s.at(2) * s.at(3) + s.at(3) * s.at(1) -
                 0.25 * ( s.at(4) * s.at(4) + s.at(5) * s.at(5) + s.at(6) * s.at(6) );
            I3 = s.at(1) * s.at(2) * s.at(3) +
                 0.25 * ( s.at(4) * s.at(5) * s.at(6) - s.at(1) * s.at(4) * s.at(4) -
                         s.at(2) * s.at(5) * s.at(5) - s.at(3) * s.at(6) * s.at(6) );
        } else {
            OOFEM_ERROR("StructuralMaterial :: ComputePrincipalValues: not supported");
        }

        /*
         * Call cubic3r to ensure, that all three real eigenvalues will be found, because we have symmetric tensor.
         * This aloows to overcome various rounding errors when solving general cubic equation.
         */
        cubic3r( ( double ) -1., I1, -I2, I3, & s1, & s2, & s3, & i );

        if ( i > 0 ) {
            answer.at(1) = s1;
        }

        if ( i > 1 ) {
            answer.at(2) = s2;
        }

        if ( i > 2 ) {
            answer.at(3) = s3;
        }

        // sort results
        for ( i = 1; i < 3; i++ ) {
            for ( j = 1; j < 3; j++ ) {
                if ( answer.at(j + 1) > answer.at(j) ) {
                    swap = answer.at(j + 1);
                    answer.at(j + 1) = answer.at(j);
                    answer.at(j) = swap;
                }
            }
        }

        return;
    }
}

void
StructuralMaterial :: computePrincipalValDir(FloatArray &answer, FloatMatrix &dir,
                                             const FloatArray &s,
                                             stressStrainPrincMode mode)
//
// This function cumputes Principal values & directions corresponding to principal values
// of strains or streses.
// strains/streses are stored in vector form in array s.
// Engineering notation is used.
//
// Problem size (3D/2D) is recognized automatically according to
// vector size.
// If size = 6 -> 3D problem, then array s contains:
//                            {Sxx,Syy,Szz,Syz,Szx,Sxy} if mode = stress_mode
//                            {Exx,Eyy,Ezz,GMyz,GMzx,GMxy} if mode = strain mode
// if size = 2 -> 2D problem, then array s contains:
//                            {Sxx,Syy,Sxy} if mode = stress_mode
//                            {Exx,Eyy,GMxy} if mode = strain mode
//
// mode      - principal strains
//           - principal stress
//           - principal deviatoric stress
//           - ..
// Input Values:
// mode
// s
//
// Return Values:
//
// matrix dir -> principal directions of strains or stresses
// array sp -> principal strains or stresses
//
{
    FloatMatrix ss;
    FloatArray sp;
    double swap;
    int i, ii, jj, kk, nval, size = s.giveSize();
    int nonzeroFlag = 0;

    // printf ("size is %d\n",size);
    if ( !( ( size == 3 ) || ( size == 6 ) ) ) {
        OOFEM_ERROR("StructuralMaterial :: computePrincipalValDir: Vector size mismatch");
    }

    if ( s.giveSize() == 3 ) {
        // 2D problem
        ss.resize(2, 2);
        answer.resize(2);

        for ( i = 1; i <= size; i++ ) {
            if ( fabs( s.at(i) ) > 1.e-20 ) {
                nonzeroFlag = 1;
            }
        }

        if ( nonzeroFlag == 0 ) {
            answer.zero();
            ss.zero();
            return;
        }

        ss.at(1, 1) = s.at(1);
        ss.at(2, 2) = s.at(2);

        if ( mode == principal_strain ) {
            ss.at(1, 2) = ss.at(2, 1) = 0.5 * s.at(3);
        } else if ( mode == principal_stress ) {
            ss.at(1, 2) = ss.at(2, 1) = s.at(3);
        } else {
            OOFEM_ERROR("StructuralMaterial :: computePrincipalValDir: not supported");
        }
    } else {
        // 3D problem
        double help;
        ss.resize(3, 3);
        answer.resize(3);

        for ( i = 1; i <= size; i++ ) {
            if ( fabs( s.at(i) ) > 1.e-20 ) {
                nonzeroFlag = 1;
            }
        }

        if ( nonzeroFlag == 0 ) {
            answer.zero();
            ss.zero();
            return;
        }

        if ( mode == principal_stress ) {
            ss.at(1, 1) = s.at(1);
            ss.at(2, 2) = s.at(2);
            ss.at(3, 3) = s.at(3);
            ss.at(1, 2) = ss.at(2, 1) = s.at(6);
            ss.at(1, 3) = ss.at(3, 1) = s.at(5);
            ss.at(2, 3) = ss.at(3, 2) = s.at(4);
        } else if ( mode == principal_deviatoricstress ) {
            help = ( s.at(1) + s.at(2) + s.at(3) ) / 3.0;
            ss.at(1, 1) = s.at(1) - help;
            ss.at(2, 2) = s.at(2) - help;
            ss.at(3, 3) = s.at(3) - help;
            ss.at(1, 2) = ss.at(2, 1) = s.at(6);
            ss.at(1, 3) = ss.at(3, 1) = s.at(5);
            ss.at(2, 3) = ss.at(3, 2) = s.at(4);
        } else if ( mode == principal_strain ) {
            ss.at(1, 1) = s.at(1);
            ss.at(2, 2) = s.at(2);
            ss.at(3, 3) = s.at(3);
            ss.at(1, 2) = ss.at(2, 1) = 0.5 * s.at(6);
            ss.at(1, 3) = ss.at(3, 1) = 0.5 * s.at(5);
            ss.at(2, 3) = ss.at(3, 2) = 0.5 * s.at(4);
        } else {
            OOFEM_ERROR("StructuralMaterial :: computePrincipalDirection: not supported");
        }
    }

#if 0
    ss.Jacobi(& answer, & dir, & i);
#else
    ss.jaco_(answer, dir, 10);
#endif
    // sort results
    nval = 2;
    if ( size == 6 ) {
        nval = 3;
    }

    for ( ii = 1; ii < nval; ii++ ) {
        for ( jj = 1; jj < nval; jj++ ) {
            if ( answer.at(jj + 1) > answer.at(jj) ) {
                // swap eigen values and eigen vectors
                swap = answer.at(jj + 1);
                answer.at(jj + 1) = answer.at(jj);
                answer.at(jj) = swap;
                for ( kk = 1; kk <= nval; kk++ ) {
                    swap = dir.at(kk, jj + 1);
                    dir.at(kk, jj + 1) = dir.at(kk, jj);
                    dir.at(kk, jj) = swap;
                }
            }
        }
    }
}


double
StructuralMaterial :: computeVonMisesStress(const FloatArray *currentStress) {
    double J2;
    double v1, v2, v3;

    if ( currentStress == NULL || currentStress->giveSize() != 6) {
        return 0.0;
    }

    v1 = ( ( currentStress->at(1) - currentStress->at(2) ) * ( currentStress->at(1) - currentStress->at(2) ) );
    v2 = ( ( currentStress->at(2) - currentStress->at(3) ) * ( currentStress->at(2) - currentStress->at(3) ) );
    v3 = ( ( currentStress->at(3) - currentStress->at(1) ) * ( currentStress->at(3) - currentStress->at(1) ) );

    J2 = ( 1. / 6. ) * ( v1 + v2 + v3 ) + currentStress->at(4) * currentStress->at(4) +
             currentStress->at(5) * currentStress->at(5) + currentStress->at(6) * currentStress->at(6);

    return sqrt(3*J2);
}


void
StructuralMaterial :: giveStrainVectorTranformationMtrx(FloatMatrix &answer,
                                                        const FloatMatrix &base,
                                                        bool transpose) const
//
// returns transformation matrix for 3d - strains to another system of axes,
// given by base.
// In base (FloatMatrix[3,3]) there are on each column stored vectors of
// coordinate system to which we do transformation.
//
// If transpose == 1 we transpose base matrix before transforming
//
{
    FloatMatrix t;
    answer.resize(6, 6);
    answer.zero();

    if ( transpose ) {
        t.beTranspositionOf(base);
    } else {
        t = base;
    }

    answer.at(1, 1) = t.at(1, 1) * t.at(1, 1);
    answer.at(1, 2) = t.at(2, 1) * t.at(2, 1);
    answer.at(1, 3) = t.at(3, 1) * t.at(3, 1);
    answer.at(1, 4) = t.at(2, 1) * t.at(3, 1);
    answer.at(1, 5) = t.at(1, 1) * t.at(3, 1);
    answer.at(1, 6) = t.at(1, 1) * t.at(2, 1);

    answer.at(2, 1) = t.at(1, 2) * t.at(1, 2);
    answer.at(2, 2) = t.at(2, 2) * t.at(2, 2);
    answer.at(2, 3) = t.at(3, 2) * t.at(3, 2);
    answer.at(2, 4) = t.at(2, 2) * t.at(3, 2);
    answer.at(2, 5) = t.at(1, 2) * t.at(3, 2);
    answer.at(2, 6) = t.at(1, 2) * t.at(2, 2);

    answer.at(3, 1) = t.at(1, 3) * t.at(1, 3);
    answer.at(3, 2) = t.at(2, 3) * t.at(2, 3);
    answer.at(3, 3) = t.at(3, 3) * t.at(3, 3);
    answer.at(3, 4) = t.at(2, 3) * t.at(3, 3);
    answer.at(3, 5) = t.at(1, 3) * t.at(3, 3);
    answer.at(3, 6) = t.at(1, 3) * t.at(2, 3);

    answer.at(4, 1) = 2.0 * t.at(1, 2) * t.at(1, 3);
    answer.at(4, 2) = 2.0 * t.at(2, 2) * t.at(2, 3);
    answer.at(4, 3) = 2.0 * t.at(3, 2) * t.at(3, 3);
    answer.at(4, 4) = ( t.at(2, 2) * t.at(3, 3) + t.at(3, 2) * t.at(2, 3) );
    answer.at(4, 5) = ( t.at(1, 2) * t.at(3, 3) + t.at(3, 2) * t.at(1, 3) );
    answer.at(4, 6) = ( t.at(1, 2) * t.at(2, 3) + t.at(2, 2) * t.at(1, 3) );

    answer.at(5, 1) = 2.0 * t.at(1, 1) * t.at(1, 3);
    answer.at(5, 2) = 2.0 * t.at(2, 1) * t.at(2, 3);
    answer.at(5, 3) = 2.0 * t.at(3, 1) * t.at(3, 3);
    answer.at(5, 4) = ( t.at(2, 1) * t.at(3, 3) + t.at(3, 1) * t.at(2, 3) );
    answer.at(5, 5) = ( t.at(1, 1) * t.at(3, 3) + t.at(3, 1) * t.at(1, 3) );
    answer.at(5, 6) = ( t.at(1, 1) * t.at(2, 3) + t.at(2, 1) * t.at(1, 3) );

    answer.at(6, 1) = 2.0 * t.at(1, 1) * t.at(1, 2);
    answer.at(6, 2) = 2.0 * t.at(2, 1) * t.at(2, 2);
    answer.at(6, 3) = 2.0 * t.at(3, 1) * t.at(3, 2);
    answer.at(6, 4) = ( t.at(2, 1) * t.at(3, 2) + t.at(3, 1) * t.at(2, 2) );
    answer.at(6, 5) = ( t.at(1, 1) * t.at(3, 2) + t.at(3, 1) * t.at(1, 2) );
    answer.at(6, 6) = ( t.at(1, 1) * t.at(2, 2) + t.at(2, 1) * t.at(1, 2) );
}


void
StructuralMaterial :: giveStressVectorTranformationMtrx(FloatMatrix &answer,
                                                        const FloatMatrix &base,
                                                        bool transpose) const
//
// returns transformation matrix for 3d - stress to another system of axes,
// given by base.
// In base (FloatMatrix[3,3]) there are on each column stored vectors of
// coordinate system to which we do transformation.
//
// If transpose == 1 we transpose base matrix before transforming
//
{
    FloatMatrix t;
    answer.resize(6, 6);
    answer.zero();

    if ( transpose ) {
        t.beTranspositionOf(base);
    } else {
        t = base;
    }

    answer.at(1, 1) = t.at(1, 1) * t.at(1, 1);
    answer.at(1, 2) = t.at(2, 1) * t.at(2, 1);
    answer.at(1, 3) = t.at(3, 1) * t.at(3, 1);
    answer.at(1, 4) = 2.0 * t.at(2, 1) * t.at(3, 1);
    answer.at(1, 5) = 2.0 * t.at(1, 1) * t.at(3, 1);
    answer.at(1, 6) = 2.0 * t.at(1, 1) * t.at(2, 1);

    answer.at(2, 1) = t.at(1, 2) * t.at(1, 2);
    answer.at(2, 2) = t.at(2, 2) * t.at(2, 2);
    answer.at(2, 3) = t.at(3, 2) * t.at(3, 2);
    answer.at(2, 4) = 2.0 * t.at(2, 2) * t.at(3, 2);
    answer.at(2, 5) = 2.0 * t.at(1, 2) * t.at(3, 2);
    answer.at(2, 6) = 2.0 * t.at(1, 2) * t.at(2, 2);

    answer.at(3, 1) = t.at(1, 3) * t.at(1, 3);
    answer.at(3, 2) = t.at(2, 3) * t.at(2, 3);
    answer.at(3, 3) = t.at(3, 3) * t.at(3, 3);
    answer.at(3, 4) = 2.0 * t.at(2, 3) * t.at(3, 3);
    answer.at(3, 5) = 2.0 * t.at(1, 3) * t.at(3, 3);
    answer.at(3, 6) = 2.0 * t.at(1, 3) * t.at(2, 3);

    answer.at(4, 1) = t.at(1, 2) * t.at(1, 3);
    answer.at(4, 2) = t.at(2, 2) * t.at(2, 3);
    answer.at(4, 3) = t.at(3, 2) * t.at(3, 3);
    answer.at(4, 4) = ( t.at(2, 2) * t.at(3, 3) + t.at(3, 2) * t.at(2, 3) );
    answer.at(4, 5) = ( t.at(1, 2) * t.at(3, 3) + t.at(3, 2) * t.at(1, 3) );
    answer.at(4, 6) = ( t.at(1, 2) * t.at(2, 3) + t.at(2, 2) * t.at(1, 3) );

    answer.at(5, 1) = t.at(1, 1) * t.at(1, 3);
    answer.at(5, 2) = t.at(2, 1) * t.at(2, 3);
    answer.at(5, 3) = t.at(3, 1) * t.at(3, 3);
    answer.at(5, 4) = ( t.at(2, 1) * t.at(3, 3) + t.at(3, 1) * t.at(2, 3) );
    answer.at(5, 5) = ( t.at(1, 1) * t.at(3, 3) + t.at(3, 1) * t.at(1, 3) );
    answer.at(5, 6) = ( t.at(1, 1) * t.at(2, 3) + t.at(2, 1) * t.at(1, 3) );

    answer.at(6, 1) = t.at(1, 1) * t.at(1, 2);
    answer.at(6, 2) = t.at(2, 1) * t.at(2, 2);
    answer.at(6, 3) = t.at(3, 1) * t.at(3, 2);
    answer.at(6, 4) = ( t.at(2, 1) * t.at(3, 2) + t.at(3, 1) * t.at(2, 2) );
    answer.at(6, 5) = ( t.at(1, 1) * t.at(3, 2) + t.at(3, 1) * t.at(1, 2) );
    answer.at(6, 6) = ( t.at(1, 1) * t.at(2, 2) + t.at(2, 1) * t.at(1, 2) );
}


void
StructuralMaterial :: givePlaneStressVectorTranformationMtrx(FloatMatrix &answer,
                                                             const FloatMatrix &base,
                                                             bool transpose) const
//
// returns transformation matrix for 2d - stress to another system of axes,
// given by base.
// In base (FloatMatrix[2,2]) there are on each column stored vectors of
// coordinate system to which we do transformation.
//
// If transpose == 1 we transpose base matrix before transforming
//
{
    FloatMatrix t;
    answer.resize(3, 3);
    answer.zero();

    if ( transpose ) {
        t.beTranspositionOf(base);
    } else {
        t = base;
    }

    answer.at(1, 1) = t.at(1, 1) * t.at(1, 1);
    answer.at(1, 2) = t.at(2, 1) * t.at(2, 1);
    answer.at(1, 3) = 2.0 * t.at(1, 1) * t.at(2, 1);

    answer.at(2, 1) = t.at(1, 2) * t.at(1, 2);
    answer.at(2, 2) = t.at(2, 2) * t.at(2, 2);
    answer.at(2, 3) = 2.0 * t.at(1, 2) * t.at(2, 2);

    answer.at(3, 1) = t.at(1, 1) * t.at(1, 2);
    answer.at(3, 2) = t.at(2, 1) * t.at(2, 2);
    answer.at(3, 3) = t.at(1, 1) * t.at(2, 2) + t.at(2, 1) * t.at(1, 2);
}


void
StructuralMaterial :: transformStrainVectorTo(FloatArray &answer, const FloatMatrix &base,
                                              const FloatArray &strainVector, bool transpose) const
//
// performs transformation of 3d-strain vector to another system of axes,
// given by base.
// In base (FloatMatrix[3,3]) there are on each column stored vectors of
// coordinate system to which we do transformation. These vectors must
// be expressed in the same coordinate system as strainVector
//
// If transpose == 1 we transpose base matrix before transforming
{
    FloatMatrix tt;

    this->giveStrainVectorTranformationMtrx(tt, base, transpose);
    answer.beProductOf(tt, strainVector);
}


void
StructuralMaterial :: transformStressVectorTo(FloatArray &answer, const FloatMatrix &base,
                                              const FloatArray &stressVector, bool transpose) const
//
//
// performs transformation of 3d-stress vector to another system of axes,
// given by base.
// In base (FloatMatrix[3,3]) there are on each column stored vectors of
// coordinate system to which we do transformation. These vectors must
// be expressed in the same coordinate system as strainVector
// If transpose == 1 we transpose base matrix before transforming
//

{
    FloatMatrix tt;

    this->giveStressVectorTranformationMtrx(tt, base, transpose);
    answer.beProductOf(tt, stressVector);
}


void
StructuralMaterial :: sortPrincDirAndValCloseTo(FloatArray *pVal, FloatMatrix *pDir,
                                                FloatMatrix *toPDir)
//
// this method sorts newly computed principal values (pVal) and
// corresponding principal directions (pDir) to be closed to
// some (often previous) principal directions (toPDir).
//
// remark : pDir and toPDir should have eigen vectors stored in columns
// and normalized.
//
{
    int i, j, k, maxJ = 0, size;
    double cosine, maxCosine, swap;

#ifdef DEBUG
    if ( ( !pDir->isSquare() ) || ( !toPDir->isSquare() ) ) {
        OOFEM_ERROR("StructuralMaterial :: sortPrincDirandValCloseTo - Not square matrix");
    }

    if ( pDir->giveNumberOfRows() != toPDir->giveNumberOfRows() ) {
        OOFEM_ERROR("StructuralMaterial :: sortPrincDirandValCloseTo - Incompatible matrices");
    }

    if ( pDir->giveNumberOfRows() != pVal->giveSize() ) {
        OOFEM_ERROR("StructuralMaterial :: sortPrincDirandValCloseTo - Incompatible pVal Array size");
    }

#endif

    //
    // compute cosine matrix, where member i,j is cosine of angle
    // between toPDir i th eigen vector and j th pDir eigen vector
    //
    // sort pVal and pDir
    size = pDir->giveNumberOfRows();
    for ( i = 1; i <= size - 1; i++ ) {
        // find closest pDir vector to toPDir i-th vector
        maxCosine = 0.0;
        for ( j = i; j <= size; j++ ) {
            for ( k = 1, cosine = 0.; k <= size; k++ ) {
                cosine += toPDir->at(k, i) * pDir->at(k, j);
            }

            cosine = fabs(cosine);
            if ( cosine > maxCosine ) {
                maxJ = j;
                maxCosine = cosine;
            }
        }

        // swap entries
        if ( maxJ != i ) {
            // swap eigenVectors and values
            swap = pVal->at(maxJ);
            pVal->at(maxJ) = pVal->at(i);
            pVal->at(i) = swap;
            for ( k = 1; k <= size; k++ ) {
                swap = pDir->at(k, maxJ);
                pDir->at(k, maxJ) = pDir->at(k, i);
                pDir->at(k, i) = swap;
            }
        }
    }
}


int
StructuralMaterial :: setIPValue(const FloatArray value, GaussPoint *aGaussPoint, InternalStateType type)
{
    StructuralMaterialStatus *status = ( StructuralMaterialStatus * ) this->giveStatus(aGaussPoint);
    if ( type == IST_StressTensor ) {
        status->letStressVectorBe(value);
        return 1;
    } else if ( type == IST_StrainTensor ) {
        status->letStrainVectorBe(value);
        return 1;
    } else if ( type == IST_StressTensorTemp ) {
        status->letTempStressVectorBe(value);
        return 1;
    } else if ( type == IST_StrainTensorTemp ) {
        status->letTempStrainVectorBe(value);
        return 1;
    } else {
        return 0;
    }
}

int
StructuralMaterial :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    StructuralMaterialStatus *status = ( StructuralMaterialStatus * ) this->giveStatus(aGaussPoint);
    if ( type == IST_StressTensor ) {
        answer = status->giveStressVector();
        return 1;
    } else if ( type == IST_vonMisesStress ) {
        answer.resize(1);
        answer.at(1) = this->computeVonMisesStress(&status->giveStressVector());
        return 1;
    } else if ( type == IST_StrainTensor ) {
        answer = status->giveStrainVector();
        return 1;
    } else if ( type == IST_StressTensorTemp ) {
        answer = status->giveTempStressVector();
        return 1;
    } else if ( type == IST_StrainTensorTemp ) {
        answer = status->giveTempStrainVector();
        return 1;
    } else if ( ( type == IST_PrincipalStressTensor ) || ( type == IST_PrincipalStressTempTensor ) ) {
        int indx;
        FloatArray st(6), s;

        if ( type == IST_PrincipalStressTensor ) {
            s = status->giveStressVector();
        } else {
            s = status->giveTempStressVector();
        }

        for ( int i = 1; i <= s.giveSize(); i++ ) {
            indx = this->giveStressStrainComponentIndOf(ReducedForm, aGaussPoint->giveMaterialMode(), i);
            if ( indx ) {
                st.at(indx) = s.at(i);
            }
        }

        this->computePrincipalValues(answer, st, principal_stress);
        return 1;
    } else if ( ( type == IST_PrincipalStrainTensor ) || ( type == IST_PrincipalStrainTempTensor ) ) {
        int indx;
        FloatArray st(6), s;

        if ( type == IST_PrincipalStrainTensor ) {
            s = status->giveStrainVector();
        } else {
            s = status->giveTempStrainVector();
        }

        for ( int i = 1; i <= s.giveSize(); i++ ) {
            indx = this->giveStressStrainComponentIndOf(ReducedForm, aGaussPoint->giveMaterialMode(), i);
            if ( indx ) {
                st.at(indx) = s.at(i);
            }
        }

        this->computePrincipalValues(answer, st, principal_strain);
        return 1;
    } else if ( type == IST_Temperature ) {
        /* add external source, if provided, such as staggered analysis */
        FieldManager *fm = domain->giveEngngModel()->giveContext()->giveFieldManager();
        Field *tf;
        int err;
        if ( ( tf = fm->giveField(FT_Temperature) ) ) {
            // temperature field registered
            FloatArray gcoords, et2;
            ( ( StructuralElement * ) aGaussPoint->giveElement() )->computeGlobalCoordinates( gcoords, * aGaussPoint->giveCoordinates() );
            if ( ( err = tf->evaluateAt(answer, gcoords, VM_Total, atTime) ) ) {
                OOFEM_ERROR3("StructuralMaterial :: giveIPValue: tf->evaluateAt failed, element %d, error code %d", aGaussPoint->giveElement()->giveNumber(), err);
            }
        } else {
            answer.resize(1);
            answer.zero();
        }

        return 1;
    } else if ( ( type == IST_CylindricalStressTensor ) || ( type == IST_CylindricalStrainTensor ) ) {
        FloatArray gc, val = status->giveStressVector();
        FloatMatrix base(3, 3);
        ( ( StructuralElement * ) aGaussPoint->giveElement() )->computeGlobalCoordinates( gc, * aGaussPoint->giveCoordinates() );
        double l = sqrt( gc.at(1) * gc.at(1) + gc.at(2) * gc.at(2) );
        if ( l > 1.e-4 ) {
            base.at(1, 1) = gc.at(1) / l;
            base.at(2, 1) = gc.at(2) / l;
            base.at(3, 1) = 0.0;

            base.at(1, 2) = -1.0 * base.at(2, 1);
            base.at(2, 2) = base.at(1, 1);
            base.at(3, 2) = 0.0;

            base.at(1, 3) = 0.0;
            base.at(2, 3) = 0.0;
            base.at(3, 3) = 1.0;

            if ( type == IST_CylindricalStressTensor ) {
                transformStressVectorTo(answer, base, val, 0);
            } else {
                transformStrainVectorTo(answer, base, val, 0);
            }
        } else {
            answer = val;
        }

        return 1;
    } else {
        return Material :: giveIPValue(answer, aGaussPoint, type, atTime);
    }
}


InternalStateValueType
StructuralMaterial :: giveIPValueType(InternalStateType type)
{
    if ( ( type == IST_StressTensor ) || ( type == IST_StressTensorTemp ) ||
        ( type == IST_PrincipalStressTensor ) || ( type == IST_PrincipalStrainTensor ) ||
        ( type == IST_PrincipalStressTempTensor ) || ( type == IST_PrincipalStrainTempTensor ) ||
        ( type == IST_CylindricalStressTensor ) ) {
        return ISVT_TENSOR_S3;
    }
    // strains components packed in engineering notation
    else if ( ( type == IST_StrainTensor ) || ( type == IST_StrainTensorTemp ) || ( type == IST_CylindricalStrainTensor ) ) {
        return ISVT_TENSOR_S3E;
    } else if ( ( type == IST_Temperature ) || ( type == IST_vonMisesStress ) ) {
        return ISVT_SCALAR;
    } else {
        return Material :: giveIPValueType(type);
    }
}


int
StructuralMaterial :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode)
{
    if ( ( type == IST_StressTensor ) || ( type == IST_StrainTensor ) ||
        ( type == IST_StressTensorTemp ) || ( type == IST_StrainTensorTemp ) ||
        ( type == IST_CylindricalStressTensor ) || ( type == IST_CylindricalStrainTensor ) ||
        ( type == IST_ShellForceMomentumTensor ) ) {
        this->giveStressStrainMask(answer, FullForm, mmode);
        return 1;
    } else if ( ( type == IST_PrincipalStressTensor ) || ( type == IST_PrincipalStrainTensor ) ||
               ( type == IST_PrincipalStressTempTensor ) || ( type == IST_PrincipalStrainTempTensor ) ) {
        answer.resize(6);
        answer.at(1) = 1;
        answer.at(2) = 2;
        answer.at(3) = 3;
        return 1;
    } else if ( type == IST_Temperature ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;
    } else {
        return Material :: giveIntVarCompFullIndx(answer, type, mmode);
    }
}


int
StructuralMaterial :: giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint)
{
    if ( ( type == IST_StressTensor ) || ( type == IST_StrainTensor ) ||
        ( type == IST_StressTensorTemp ) || ( type == IST_StrainTensorTemp ) ||
        ( type == IST_CylindricalStressTensor ) || ( type == IST_CylindricalStrainTensor ) ||
        ( type == IST_ShellForceMomentumTensor ) ) {
        return this->giveSizeOfReducedStressStrainVector( aGaussPoint->giveMaterialMode() );
    } else if ( ( type == IST_PrincipalStressTensor ) || ( type == IST_PrincipalStrainTensor ) ||
               ( type == IST_PrincipalStressTempTensor ) || ( type == IST_PrincipalStrainTempTensor ) ) {
        return 3;
    } else if ( ( type == IST_Temperature ) || ( type == IST_vonMisesStress ) ) {
        return 1;
    } else {
        return Material :: giveIPValueSize(type, aGaussPoint);
    }
}


void
StructuralMaterial :: computeStressIndependentStrainVector(FloatArray &answer,
                                                           GaussPoint *gp, TimeStep *stepN, ValueModeType mode)
{
    FloatArray fullAnswer, et, e0, eigenstrain, answerTemper, answerEigenstrain;
    FloatMatrix GCS;
    MaterialMode matmode = gp->giveMaterialMode();
    StructuralCrossSection *crossSection =  dynamic_cast< StructuralCrossSection * >( gp->giveCrossSection() );
    Element *elem = gp->giveElement();
    StructuralElement *selem = dynamic_cast< StructuralElement * >( gp->giveElement() );


    answer.resize(0);
    answerTemper.resize(0);
    answerEigenstrain.resize(0);

    if ( stepN->giveIntrinsicTime() < this->castingTime ) {
        answer.zero();
        return;
    }

    //sum up all prescribed temperatures over an element
    //elem->computeResultingIPTemperatureAt(et, stepN, gp, mode);
    if ( selem ) {
        selem->computeResultingIPTemperatureAt(et, stepN, gp, mode);        // HUHU
    }

    //sum up all prescribed eigenstrain over an element
    if ( selem ) {
        selem->computeResultingIPEigenstrainAt(eigenstrain, stepN, gp, mode);
    }

    if ( eigenstrain.giveSize() != 0 && eigenstrain.giveSize() != giveSizeOfReducedStressStrainVector(matmode) ) {
        OOFEM_ERROR5("StructuralMaterial :: Number of given eigenstrain components %d is different than required %d by material mode %s, element %d", eigenstrain.giveSize(), giveSizeOfReducedStressStrainVector(matmode), __MaterialModeToString(matmode), elem->giveNumber() );
    }

    /* add external source, if provided */
    FieldManager *fm = domain->giveEngngModel()->giveContext()->giveFieldManager();
    Field *tf;
    if ( ( tf = fm->giveField(FT_Temperature) ) ) {
        // temperature field registered
        FloatArray gcoords, et2;
        int err;
        elem->computeGlobalCoordinates( gcoords, * gp->giveCoordinates() );
        if ( ( err = tf->evaluateAt(et2, gcoords, mode, stepN) ) ) {
            OOFEM_ERROR3("StructuralMaterial :: computeStressIndependentStrainVector: tf->evaluateAt failed, element %d, error code %d", elem->giveNumber(), err);
        }

        if ( et2.isNotEmpty() ) {
            if ( et.isEmpty() ) {
                et = et2;
            } else {
                et.at(1) += et2.at(1);
            }
        }
    }


    if ( et.giveSize() ) { //found temperature boundary conditions or prescribed field
        double thick, width;

        this->giveThermalDilatationVector(e0, gp, stepN);

        switch ( matmode ) {
        case _2dBeam:
            thick = crossSection->give(CS_Thickness);
            answerTemper.resize(3);
            answerTemper.zero();
            answerTemper.at(1) = e0.at(1) * ( et.at(1) - this->giveReferenceTemperature() );
            if ( et.giveSize() > 1 ) {
                answerTemper.at(2) = e0.at(1) * et.at(2) / thick;     // kappa_x
            }

            break;
        case _3dBeam:
            thick = crossSection->give(CS_Thickness);
            width = crossSection->give(CS_Width);
            answerTemper.resize(6);
            answerTemper.zero();
            answerTemper.at(1) = e0.at(1) * ( et.at(1) - this->giveReferenceTemperature() );
            if ( et.giveSize() > 1 ) {
                answerTemper.at(5) = e0.at(1) * et.at(2) / thick;     // kappa_y
                if ( et.giveSize() > 2 ) {
                    answerTemper.at(6) = e0.at(1) * et.at(3) / width;     // kappa_z
                }
            }

            break;
        case _2dPlate:
            thick = crossSection->give(CS_Thickness);
            if ( et.giveSize() > 1 ) {
                answerTemper.resize(5);
                answerTemper.zero();

                if ( et.giveSize() > 1 ) {
                    answerTemper.at(1) = e0.at(1) * et.at(2) / thick;     // kappa_x
                    answerTemper.at(2) = e0.at(2) * et.at(2) / thick;     // kappa_y
                }
            }

            break;
        case _3dShell:
            thick = crossSection->give(CS_Thickness);
            answerTemper.resize(8);
            answerTemper.zero();

            answerTemper.at(1) = e0.at(1) * ( et.at(1) - this->giveReferenceTemperature() );
            answerTemper.at(2) = e0.at(2) * ( et.at(1) - this->giveReferenceTemperature() );
            if ( et.giveSize() > 1 ) {
                answerTemper.at(4) = e0.at(1) * et.at(2) / thick;     // kappa_x
                answerTemper.at(5) = e0.at(2) * et.at(2) / thick;     // kappa_y
            }

            break;
        default:
            if ( e0.giveSize() ) {
                fullAnswer = e0;
                if ( mode == VM_Total ) {
                    fullAnswer.times(et.at(1) - this->referenceTemperature);
                } else {
                    fullAnswer.times( et.at(1) );
                }

                this->giveReducedCharacteristicVector(answerTemper, gp, fullAnswer);
            }
        }
    }

    if ( eigenstrain.giveSize() ) { //found prescribed eigenstrain
        switch ( matmode ) {
        case _1dMat:
        case _3dMat:
        case _3dMat_F:
        case _PlaneStress:
        case _PlaneStrain:
        case _3dRotContinuum:
        case _3dMicroplane:
            fullAnswer = eigenstrain;
            break;

        default:
            OOFEM_ERROR2("StructuralMaterial :: Material mode %s for eigenstrains not supported", __MaterialModeToString(matmode) );
        }

        answerEigenstrain = fullAnswer;
        //this->giveReducedCharacteristicVector(answerEigenstrain, gp, fullAnswer);
    }

    //join temperature and eigenstrain vectors, compare vector sizes
    if ( answerTemper.giveSize() ) {
        answer = answerTemper;
        if ( answerEigenstrain.giveSize() ) {
            if ( answerTemper.giveSize() != answerEigenstrain.giveSize() ) {
                OOFEM_ERROR4("StructuralMaterial :: Vector of temperature strains has the size %d which is different with the size of eigenstrain vector %d, element %d", answerTemper.giveSize(), answerEigenstrain.giveSize(), elem->giveNumber() );
            }

            answer.add(answerEigenstrain);
        }
    } else {
        if ( answerEigenstrain.giveSize() ) {
            answer = answerEigenstrain;
        }
    }
}


void
StructuralMaterial :: giveFullCharacteristicVector(FloatArray &answer,
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
    IntArray indx;
    int i, j, answerSize = 6;

    if ( this->hasMaterialModeCapability(mode) ) {
        if ( mode == _3dMat ) {
            answer = strainVector;
            return;
        }

        answer.resize(answerSize);
        answer.zero();

        this->giveStressStrainMask( indx, ReducedForm, gp->giveMaterialMode() );
        for ( i = 1; i <= indx.giveSize(); i++ ) {
            if ( ( j = indx.at(i) ) ) {
                answer.at(j) = strainVector.at(i);
            }
        }

        return;
    } else {
        OOFEM_ERROR("StructuralMaterial :: giveFullCharacteristicVector - invalid mode");
    }
}


void
StructuralMaterial :: giveReducedCharacteristicVector(FloatArray &answer, GaussPoint *gp,
                                                      const FloatArray &charVector3d)
//
// returns reduced stressVector or strainVector from full 3d vector reduced
// to vector required by gp->giveStressStrainMode()
//
{
    MaterialMode mode = gp->giveMaterialMode();
    IntArray indx;
    int size = charVector3d.giveSize();
    int i, j;

    if ( this->hasMaterialModeCapability(mode) ) {
        if ( ( mode == _3dMat ) || ( mode == _3dMicroplane ) ) {
            if ( size != 6 ) {
                OOFEM_ERROR("StructuralMaterial :: giveReducedCharacteristicVector - charVector3d size mismatch");
            }

            answer = charVector3d;
            return;
        }

        this->giveStressStrainMask( indx, ReducedForm, gp->giveMaterialMode() );
        answer.resize( indx.giveSize() );
        answer.zero();

        for ( i = 1; i <= indx.giveSize(); i++ ) {
            if ( ( j = indx.at(i) ) ) {
                answer.at(i) = charVector3d.at(j);
            }
        }

        return;
    } else {
        OOFEM_ERROR("StructuralMaterial :: giveFullCharacteristicVector - invalid mode");
    }
}


IRResultType
StructuralMaterial :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

#  ifdef VERBOSE
    // VERBOSE_PRINT1 ("Instanciating material ",this->giveNumber())
#  endif
    this->Material :: initializeFrom(ir);

    referenceTemperature = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, referenceTemperature, IFT_StructuralMaterial_referencetemperature, "referencetemperature"); // Macro

    return IRRT_OK;
}


int
StructuralMaterial :: giveInputRecordString(std :: string &str, bool keyword)
{
    char buff [ 1024 ];

    Material :: giveInputRecordString(str, keyword);
    sprintf(buff, " referencetemperature %e", this->referenceTemperature);
    str += buff;

    return 1;
}

} // end namespace oofem
