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

#ifndef interfacetype_h
#define interfacetype_h

namespace oofem {
/**
 * Enumerative type, used to identify interface type.
 * @see Interface More details.
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

    SpatialLocalizerInterfaceType,

    EIPrimaryUnknownMapperInterfaceType,
    EIPrimaryFieldInterfaceType,

    NonlocalMaterialStatusExtensionInterfaceType,
    GradDpMaterialExtensionInterfaceType,
    GradDpMaterialStatusExtensionInterfaceType,

    NonlocalMaterialExtensionInterfaceType,
    NonlocalMaterialStiffnessInterfaceType,
    MaterialModelMapperInterfaceType,
    RandomMaterialStatusExtensionInterfaceType,

    HydrationModelInterfaceType,
    HydrationModelStatusInterfaceType,

    LEPlicElementInterfaceType,
    LevelSetPCSElementInterfaceType,

    XfemElementInterfaceType,
    VTKXMLExportModuleElementInterfaceType,
    FailureModuleElementInterfaceType
};
} // end namespace oofem
#endif // interfacetype_h
