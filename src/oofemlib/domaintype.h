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

#ifndef domaintype_h
#define domaintype_h

#include "enumitem.h"

namespace oofem {

#define domainType_DEF \
    ENUM_ITEM(_unknownMode) \
    ENUM_ITEM(_2dPlaneStressMode) \
    ENUM_ITEM(_PlaneStrainMode) \
    ENUM_ITEM(_2dPlaneStressRotMode) \
    ENUM_ITEM(_3dMode) \
    ENUM_ITEM(_3dAxisymmMode) \
    ENUM_ITEM(_2dMindlinPlateMode) \
    ENUM_ITEM(_3dShellMode) \
    ENUM_ITEM(_2dTrussMode) \
    ENUM_ITEM(_1dTrussMode) \
    ENUM_ITEM(_2dBeamMode) \
    ENUM_ITEM(_HeatTransferMode) \
    ENUM_ITEM(_HeatMass1Mode) /* Coupled heat and mass (1 matter) transfer */   \
    ENUM_ITEM(_2dIncompressibleFlow) /* 2d Incompressible flow, no energy eq */ \
    ENUM_ITEM(_3dIncompressibleFlow) /* 3d Incompressible flow, no energy eq */ \
    ENUM_ITEM(_2dLatticeMode)\
    ENUM_ITEM(_3dDirShellMode) /* 7 parameter shell based on director fields */ \
/**
 * Type representing type of domain.
 * Domain type (the member value of Domain class) is used to determine the default
 * number of DOFs per node and side and to determine their corresponding physical meaning.
 */
enum domainType {
    domainType_DEF
};

#undef ENUM_ITEM
#undef ENUM_ITEM_WITH_VALUE
#undef enumitem_h


const char *__domainTypeToString(domainType _value);
} // end namespace oofem
#endif // domaintype_h
