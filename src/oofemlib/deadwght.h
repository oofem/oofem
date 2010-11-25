/* $Header: /home/cvs/bp/oofem/oofemlib/src/deadwght.h,v 1.8 2003/04/06 14:08:23 bp Exp $ */
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

//   *************************
//   *** CLASS DEAD WEIGHT ***
//   *************************


#ifndef deadwght_h
#define deadwght_h

#include "bodyload.h"
#include "bcgeomtype.h"
#include "valuemodetype.h"

namespace oofem {
/**
 * This class implements a gravity-like load, or internal source (heat etc.) for transport problems.
 * The inheritted attribute 'componentArray' contains the components of an
 * loading prescribed per unit volume.
 */
class DeadWeight : public BodyLoad
{
    /*
     * This class implements a gravity-like load.
     * DESCRIPTION
     * The attribute 'componentArray' contains the components of an acceleration
     * 'a', expected (but not required) to be downwards vertical.
     * TASK
     * returning the body force  rho*a  acting on a given element.
     */
public:
    /// Constructor
    DeadWeight(int i, Domain *d) : BodyLoad(i, d) { }         // constructor
    /**
     * Computes components values of deadweight field at given point (coordinates given in Global c.s.).
     * taking into account corresponding load time function value while respecting load response mode.
     * @param answer component values at given point and time
     * @param stepN time step representing time
     * @param coords gp global coordinates, which are used to evaluate components values
     * @param mode determines response mode
     */
    void         computeValueAt(FloatArray &answer, TimeStep *atTime, FloatArray &coords, ValueModeType mode)
    { computeComponentArrayAt(answer, atTime, mode); }
    ///Returns type of class
    classType    giveClassID() const { return DeadWeightClass; }
    ///Returns name of class
    const char *giveClassName() const { return "DeadWeight"; }
    ///Returns input record name of the receiver.
    const char *giveInputRecordName() const { return "DeadWeight"; }
    ///Returns type of boundary condition
    bcValType    giveBCValType() const { return ForceLoadBVT; }
    ///Returns type of load
    bcGeomType    giveBCGeoType() const { return BodyLoadBGT; }
};
} // end namespace oofem
#endif // deadwght_h






