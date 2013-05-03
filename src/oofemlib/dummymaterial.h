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
class DummyMaterial : public Material
{
public:
    DummyMaterial (int n, Domain* d) : Material(n,d) {};
    virtual int hasMaterialModeCapability(MaterialMode mode) { return 0;}

    virtual const char *giveClassName() const { return "DummyMaterial"; }
    virtual classType giveClassID() const { return DummyMaterialClass; }
    virtual IRResultType initializeFrom(InputRecord *ir) { return IRRT_OK; }
};
} // end namespace oofem
#endif // material_h
