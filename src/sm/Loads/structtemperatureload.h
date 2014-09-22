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

#ifndef structtemperatureload_h
#define structtemperatureload_h

#include "load.h"

#define _IFT_StructuralTemperatureLoad_Name "structtemperatureload"

namespace oofem {
class Element;
class TimeStep;

/**
 * This class implements temperature (constant) load over the element
 * component array contains one or two numbers;
 * componentArray->at(1) contains increment of temperature in mid-surface
 * componentArray->at(2) contains increment of gradient of temperature
 * over the thickness of element (optional)
 */
class StructuralTemperatureLoad : public Load
{
public:
    StructuralTemperatureLoad(int n, Domain * d) : Load(n, d) { }

    virtual void computeValueAt(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, ValueModeType mode);

    virtual bcValType giveBCValType() const { return TemperatureBVT; }
    virtual bcGeomType giveBCGeoType() const { return BodyLoadBGT; }
    virtual const char *giveInputRecordName() const { return _IFT_StructuralTemperatureLoad_Name; }
    virtual const char *giveClassName() const { return "StructuralTemperatureLoad"; }
};
} // end namespace oofem
#endif // structtemperatureload_h
