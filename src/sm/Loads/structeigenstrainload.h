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

#ifndef structeigenstrainload_h
#define structeigenstrainload_h

#include "load.h"

#define _IFT_StructuralEigenstrainLoad_Name "structeigenstrainload"

namespace oofem {
class Element;
class TimeStep;

/**
 * This class implements prescribed eigenstrain (stress-free strain). It reads six
 * strain components (xx, yy, zz, yz, zx, xy) in the global coordinate system. 2D, 1D?
 *
 * @author Vit Smilauer
 */
class StructuralEigenstrainLoad : public Load
{
public:
    StructuralEigenstrainLoad(int i, Domain * d) : Load(i, d) { }

    /**
     * Computes components values of eigenstrain field at given point (coordinates given in Global c.s.).
     * taking into account corresponding load time function value while respecting load response mode.
     * @param answer Component values at given point and time.
     * @param tStep Time step representing time.
     * @param coords Integration point global coordinates, which are used to evaluate components values.
     * @param mode Determines response mode.
     */
    virtual void computeValueAt(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, ValueModeType mode);

    virtual const char *giveInputRecordName() const { return _IFT_StructuralEigenstrainLoad_Name; }
    virtual const char *giveClassName() const { return "StructuralEigenstrainLoad"; }
    virtual bcValType giveBCValType() const { return EigenstrainBVT; }
    virtual bcGeomType giveBCGeoType() const { return BodyLoadBGT; }
};
} // end namespace oofem
#endif // structeigenstrainload_h
