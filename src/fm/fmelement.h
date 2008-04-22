/* $Header: /home/cvs/bp/oofem/tm/src/transportelement.h,v 1.3 2003/04/23 14:22:15 bp Exp $ */
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

//   *********************************
//   *** CLASS GENERAL CFD ELEMENT ***
//   *********************************

#ifndef fmelement_h
#define fmelement_h


#include "element.h"
#include "femcmpnn.h"
#include "domain.h"
#include "flotmtrx.h"
#include "cltypes.h"
#include "primaryfield.h"


// declaration of basic boundary codes
#define FMElement_PrescribedTractionBC ( 1 << 0 )
#define FMElement_PrescribedUnBC       ( 1 << 1 )
#define FMElement_PrescribedUsBC       ( 1 << 2 )
#define FMElement_PrescribedPressureBC ( 1 << 3 )

class TimeStep;
class Node;
class Material;
class GaussPoint;
class FloatMatrix;
class FloatArray;
class IntArray;

/**
 * This abstract class represent a general base element class for
 * fluid dynamic problems.
 */
class FMElement : public Element
{
public:
protected:
    /// Array of Boundary sides
    IntArray boundarySides;
    /// boundary sides codes
    IntArray boundaryCodes;

public:
    // constructor
    FMElement(int, Domain *);
    ~FMElement();                        // destructor

    ///Initializes receiver acording to object description stored in input record.
    IRResultType initializeFrom(InputRecord *ir);

    // definition
    const char *giveClassName() const { return "FMElement"; }
    classType                giveClassID() const { return FMElementClass; }
};

#endif // fmelement_h







