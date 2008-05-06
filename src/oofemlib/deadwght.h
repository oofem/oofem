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

/**
 * This class implements a gravity-like load.
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
     * Returns receiver load type.
     * @return StructuralLoadLT;
     */
    bcGeomType    giveBCGeoType() const { return BodyLoadBGT; }
    void         computeValueAt(FloatArray &answer, TimeStep *atTime, FloatArray &coords, ValueModeType mode)
    { computeComponentArrayAt(answer, atTime, mode); }

    /// Returns input record name of the receiver.
    const char *giveInputRecordName() const { return "DeadWeight"; }
};


#endif // deadwght_h






