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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#ifndef matesponsemode_h
#define matesponsemode_h

#include "enum.h"

namespace oofem {
#define ENUM_TYPE MatResponseMode
#define ENUM_DEF \
    ENUM_ITEM_WITH_VALUE(TangentStiffness, 0) \
    ENUM_ITEM_WITH_VALUE(SecantStiffness, 1) \
    ENUM_ITEM_WITH_VALUE(ElasticStiffness, 2) \
    ENUM_ITEM_WITH_VALUE(Stress, 3) \
    ENUM_ITEM_WITH_VALUE(Conductivity, 4)  /* element level conductivity matrix */ \
    ENUM_ITEM_WITH_VALUE(Conductivity_ww, 5) /* material level conductivity submatrix */ \
    ENUM_ITEM_WITH_VALUE(Conductivity_hh, 6) /* material level conductivity submatrix */ \
    ENUM_ITEM_WITH_VALUE(Conductivity_hw, 7) /* material level conductivity submatrix */ \
    ENUM_ITEM_WITH_VALUE(Conductivity_wh, 8) /* material level conductivity submatrix */ \
    ENUM_ITEM_WITH_VALUE(Capacity, 9)                                                   \
    ENUM_ITEM_WITH_VALUE(Capacity_ww, 10) /* material level capacity submatrix */ \
    ENUM_ITEM_WITH_VALUE(Capacity_hh, 11) /* material level capacity submatrix */ \
    ENUM_ITEM_WITH_VALUE(Capacity_hw, 12) /* material level capacity submatrix */ \
    ENUM_ITEM_WITH_VALUE(Capacity_wh, 13) /* material level capacity submatrix */ \
    ENUM_ITEM_WITH_VALUE(IntSource, 14)                                                  \
    ENUM_ITEM_WITH_VALUE(IntSource_ww, 15) /* material level internal source submatrix - water source */ \
    ENUM_ITEM_WITH_VALUE(IntSource_hh, 16) /*  - heat source */ \
    ENUM_ITEM_WITH_VALUE(IntSource_hw, 17) /*  - heat source dependency on water content change */ \
    ENUM_ITEM_WITH_VALUE(IntSource_wh, 18) /*  - water source dependency on temperature change */ \
    ENUM_ITEM_WITH_VALUE(Permeability, 19) \
    ENUM_ITEM_WITH_VALUE(FluidMassBalancePressureContribution, 20) \
    ENUM_ITEM_WITH_VALUE(BiotConstant, 21) \
    ENUM_ITEM_WITH_VALUE(CompressibilityCoefficient, 22) \
    ENUM_ITEM_WITH_VALUE(FluidViscosity, 23) \
    ENUM_ITEM_WITH_VALUE(Flux, 24) \
    ENUM_ITEM_WITH_VALUE(DSigmaDT, 25) \
    ENUM_ITEM_WITH_VALUE(ElasticBulkModulus, 26) \
    ENUM_ITEM_WITH_VALUE(ElasticBulkModulusInverse, 27) \
    ENUM_ITEM_WITH_VALUE(MRM_ScalarOne, 28) \
    ENUM_ITEM_WITH_VALUE(DeviatoricStiffness, 29) \
    ENUM_ITEM_WITH_VALUE(DeviatoricStress, 30) 

#include "enum-impl.h"

} // end namespace oofem
#endif // matesponsemode_h
