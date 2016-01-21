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

#ifndef emptymaterial_h
#define emptymaterial_h

#include "material.h"

#define _IFT_DummyMaterial_Name "dummymat"

namespace oofem {
class GaussPoint;
class Domain;
class InputRecord;

/**
 * Dummy material model, no functionality. Convenient for special-purpose elements not
 * requiring real material.
 */
class OOFEM_EXPORT DummyMaterial : public Material
{
public:
    DummyMaterial(int n, Domain * d) : Material(n, d) { }
    virtual int hasMaterialModeCapability(MaterialMode mode) { return 0; }

    virtual const char *giveClassName() const { return "DummyMaterial"; }
    virtual const char *giveInputRecordName() const { return _IFT_DummyMaterial_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir) { return IRRT_OK; }
};
} // end namespace oofem
#endif // material_h
