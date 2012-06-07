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

#ifndef materialmode_h
#define materialmode_h

#include "enumitem.h"

namespace oofem {

#define MaterialMode_DEF \
    ENUM_ITEM(_Unknown)   \
    ENUM_ITEM(_3dMat) \
    ENUM_ITEM(_3dMat_F)  /* 3d deformation gradient */ \
    ENUM_ITEM(_3dMatGrad)  /* 3d model with gradient of internal variable */ \
    ENUM_ITEM(_1dMatGrad)  /* 1d model with gradient of internal variable */ \
    ENUM_ITEM(_PlaneStressGrad)  /* plane stress with gradient of internal variable */ \
    ENUM_ITEM(_PlaneStrainGrad)  /* plane strain with gradient of internal variable */ \
    ENUM_ITEM(_PlaneStress) \
    ENUM_ITEM(_PlaneStrain) \
    ENUM_ITEM(_2dPlate) \
    ENUM_ITEM(_1dMat) \
    ENUM_ITEM(_2dBeam) \
    ENUM_ITEM(_3dBeam) \
    ENUM_ITEM(_3dShell) \
    ENUM_ITEM(_3dRotContinuum) /* axisymmetry */ \
  \
    ENUM_ITEM(_2dPlateLayer) \
    ENUM_ITEM(_2dBeamLayer) \
    ENUM_ITEM(_3dShellLayer) \
    ENUM_ITEM(_PlaneStressRot) \
  \
    ENUM_ITEM(_1dFiber) \
    ENUM_ITEM(_3dMicroplane) \
    ENUM_ITEM(_3dInterface) \
    ENUM_ITEM(_2dInterface) \
    ENUM_ITEM(_1dInterface) \
  \
    ENUM_ITEM(_1dHeat) /* 1d heat */ \
    ENUM_ITEM(_1dHeMo) /* 1d heat and mass (one component) transfer */ \
    ENUM_ITEM(_2dHeat) /* 2d heat */ \
    ENUM_ITEM(_2dHeMo) /* 2d heat and mass (one component) transfer */ \
    ENUM_ITEM(_3dHeat) \
    ENUM_ITEM(_3dHeMo) \
  \
    ENUM_ITEM(_2dFlow) \
    ENUM_ITEM(_2dAxiFlow) \
    ENUM_ITEM(_3dFlow) \
\
    ENUM_ITEM(_2dLattice) \
/**
 * Type representing material mode of integration point.
 */
enum MaterialMode {
    MaterialMode_DEF
};

#undef ENUM_ITEM
#undef ENUM_ITEM_WITH_VALUE
#undef enumitem_h


const char *__MaterialModeToString(MaterialMode _value);
} // end namespace oofem
#endif // materialmode_h
