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

#ifndef drawmode_h
#define drawmode_h

namespace oofem {
enum DrawMode {
    unknown,

    rawGeometry,
    deformedGeometry,
    eigenVectorGeometry,
    nodeAnnotation,
    appliedPrimaryBc,

    internalStateBegin,
    mxForce,
    myForce,
    mzForce,
    myzForce,
    mzxForce,
    mxyForce,
    //qxzForce,
    //qyzForce,
    sxForce,
    syForce,
    szForce,
    syzForce,
    szxForce,
    sxyForce,
    yieldState,
    crackedState,
    stressErrorState,
    requiredAdaptiveMeshSizeState,
    damageLevel,
    errorIndicatorLevel,
    relativeMeshSizeDensity,
    temperatureField,
    massConcentration1Field,
    velocityField,
    pressureField,
    vofField,
    densityField,

    hydrationDegreeState,
    humidityState,

    internalStateEnd
};
} // end namespace oofem
#endif // drawmode_h
