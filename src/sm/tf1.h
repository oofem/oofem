/* $Header: /home/cvs/bp/oofem/sm/src/tf1.h,v 1.4 2003/04/06 14:08:31 bp Exp $ */
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

//   ************************************************
//   *** CLASS USER DEFINED TEMPERATURE FIELD Tf1 ***
//   ************************************************


#ifndef tf1_h
#define tf1_h

#include "structtemperatureload.h"
#include "domain.h"

namespace oofem {
/**
 * Class representing usr defined temparature field.
 * No user input. The expression is hard - coded in the class body
 * as a function of global x,y and z coordinates and time t.
 *
 * THE LOAD TIME FUNCTION is not used here, the function provided is
 * supposed to be function of time and coordinates.
 *
 */
class TF1 : public StructuralTemperatureLoad
{
private:
public:

    /**
     * Constructor. Creates temperature load function with given number, belonging to given domain.
     * @param n load time function number
     * @param d domain to which new object will belongs.
     */
    TF1(int i, Domain *d) : StructuralTemperatureLoad(i, d) { }
    /// Destructor
    virtual ~TF1()  { }

    // computations
    /**
     * Computes components values of temperature field at given point (coordinates given in Global c.s.).
     * taking into account corresponding load time function value respecting load response mode.
     * @param answer component values at given point and time
     * @param stepN time step representing time
     * @param coords gp global coordinates, which are used to evaluate components values.
     * @param mode determines response mode.
     */
    virtual void         computeValueAt(FloatArray &answer, TimeStep *tStep, FloatArray &coords, ValueModeType mode);

    /// Returns classType id of receiver.
    classType   giveClassID() const { return UserDefinedTemperatureFieldClass; }
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "TF1"; }
    /**
     * Initializes receiver acording to object description stored in input record.
     * Must be implemented in derived classes
     */
    IRResultType initializeFrom(InputRecord *ir);
};
} // end namespace oofem
#endif // tf1_h
