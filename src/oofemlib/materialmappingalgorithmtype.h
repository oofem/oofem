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

#ifndef materialmappingalgorithmtype_h
#define materialmappingalgorithmtype_h

#include "enumitem.h"

namespace oofem {
#define MaterialMappingAlgorithmType_DEF \
    ENUM_ITEM(MMA_ClosestPoint) \
    ENUM_ITEM(MMA_LeastSquareProjection) \
    ENUM_ITEM(MMA_ShapeFunctionProjection)

/**
 * Enumerative type used to classify supported
 * MaterialMappingAlgorithms
 */
enum MaterialMappingAlgorithmType {
    MaterialMappingAlgorithmType_DEF
};

#undef ENUM_ITEM
#undef ENUM_ITEM_WITH_VALUE
#undef enumitem_h


const char *__MaterialMappingAlgorithmTypeToString(MaterialMappingAlgorithmType _value);
} // end namespace oofem
#endif // materialmappingalgorithmtype_h
