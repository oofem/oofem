/* v 1.9 2010/01/05 vs */
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
 *               Copyright (C) 2009   Vit Smilauer
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

#ifndef structeigenstrainload_h
#define structeigenstrainload_h


#include "load.h"

namespace oofem {
class Element;
class TimeStep;


class StructuralEigenstrainLoad : public Load
{
    /* This class implements prescribed eigenstrain (stress-free strain). It reads six
     * strain components (xx, yy, zz, yz, zx, xy) in the global coordinate system. 2D, 1D?
     */

public:
    StructuralEigenstrainLoad(int i, Domain *d) : Load(i, d) { }  // constructor

    /**
     * Computes components values of eigenstrain field at given point (coordinates given in Global c.s.).
     * taking into account corresponding load time function value while respecting load response mode.
     * @param answer component values at given point and time
     * @param stepN time step representing time
     * @param coords gp global coordinates, which are used to evaluate components values
     * @param mode determines response mode
     */
    virtual void         computeValueAt(FloatArray &answer, TimeStep *tStep, FloatArray &coords, ValueModeType mode)
    { this->computeComponentArrayAt(answer, tStep, mode); }

    ///Returns type of class
    classType    giveClassID() const { return StructuralEigenstrainLoadClass; }
    ///Returns name of class
    const char *giveClassName() const { return "StructuralEigenstrainLoad"; }
    ///Returns type of boundary condition
    bcValType    giveBCValType() const { return EigenstrainBVT; }
    ///Returns type of load
    bcGeomType   giveBCGeoType() const { return BodyLoadBGT; }
};
} // end namespace oofem
#endif // structeigenstrainload_h
