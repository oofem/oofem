/* $Header: /home/cvs/bp/oofem/oofemlib/src/bodyload.h,v 1.8 2003/04/06 14:08:23 bp Exp $ */
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


//   ***********************
//   *** CLASS BODY LOAD ***
//   ***********************

#ifndef bodyload_h
#define bodyload_h

#include "load.h"

#include "classtype.h"
#include "bcgeomtype.h"

namespace oofem {
class Element;
class TimeStep;

/**
 * Class implementing element body load, acting over whole element volume (e.g., the dead weight).
 * Body load is usually attribute of one or more elements.
 *
 * This base body load class only defines the common services common to all derived classes.
 * Derived classes need to implement services declared by base Load class.
 */
class BodyLoad : public Load
{
    /*
     * This abstract class is the superclass of all loads (e.g., the dead weight)
     * that act on the whole volume of finite elements. A body load is usually
     * attribute of many elements.
     */
public:
    /**
     * Constructor. Creates Body Load object with given number, belonging to given domain.
     * @param n load number
     * @param d domain to which new object will belongs.
     */
    BodyLoad(int i, Domain *d) : Load(i, d) { }     // constructor


    classType    giveClassID() const { return BodyLoadClass; }
    /**
     * Returns receiver's load gemetry type.
     * @return returns BodyLoadGT value.
     */
    bcGeomType giveBCGeoType() const { return BodyLoadBGT; }
    const char *giveClassName() const { return "BodyLoad"; }
};
} // end namespace oofem
#endif // bodyload_h








