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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#ifndef tf1_h
#define tf1_h

#include "structtemperatureload.h"
#include "domain.h"

namespace oofem {
/**
 * Class representing user defined temperature field.
 * No user input. The expression is hard - coded in the class body
 * as a function of global x,y and z coordinates and time t.
 *
 * The load time function is not used here, the function provided is
 * supposed to be function of time and coordinates.
 */
class TF1 : public StructuralTemperatureLoad
{
public:
    /**
     * Constructor. Creates temperature load function with given number, belonging to given domain.
     * @param n Load time function number.
     * @param d Domain to which new object will belongs.
     */
    TF1(int n, Domain *d) : StructuralTemperatureLoad(n, d) { }
    /// Destructor
    virtual ~TF1()  { }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void computeValueAt(FloatArray &answer, TimeStep *tStep, FloatArray &coords, ValueModeType mode);

    virtual classType giveClassID() const { return TF1Class; }
    virtual const char *giveClassName() const { return "TF1"; }

};
} // end namespace oofem
#endif // tf1_h
