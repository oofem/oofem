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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//
// FILE: interfacetype.h
//

#ifndef interfacetype_h
#define interfacetype_h

/**
 * Enumerative type, used to identify interface type.
 * @see Interface class for details
 */
enum InterfaceType {
    UnknownInterfaceType,

    LayeredCrossSectionInterfaceType,
    FiberedCrossSectionInterfaceType,

    ZZNodalRecoveryModelInterfaceType,
    NodalAveragingRecoveryModelInterfaceType,
    SPRNodalRecoveryModelInterfaceType,

    ZZErrorEstimatorInterfaceType,
    HuertaErrorEstimatorInterfaceType,
    Huerta1dErrorEstimatorInterfaceType, // experimental

    DirectErrorIndicatorRCInterfaceType,
    ZZRemeshingCriteriaInterfaceType,
    HuertaRemeshingCriteriaInterfaceType,

    SpatialLocalizerInterfaceType,

    EIPrimaryUnknownMapperInterfaceType,
    MMAShapeFunctProjectionInterfaceType,
    EIPrimaryFieldInterfaceType,

    NonlocalMaterialStatusExtensionInterfaceType,

    NonlocalMaterialExtensionInterfaceType,
    NonlocalMaterialStiffnessInterfaceType,
    MaterialModelMapperInterfaceType,
    RandomMaterialStatusExtensionInterfaceType,

    HydrationModelInterfaceType,
    HydrationModelStatusInterfaceType,

    LEPlicElementInterfaceType,
    LevelSetPCSElementInterfaceType,

    XfemElementInterfaceType
};

#endif // interfacetype_h
