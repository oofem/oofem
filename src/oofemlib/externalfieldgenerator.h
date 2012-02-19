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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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

#ifndef externalfieldgenerator_h
#define externalfieldgenerator_h

#include "randomfieldgenerator.h"
#include "gausspnt.h"

#define NAME_MAX_LENGTH 200

namespace oofem {
/**
 * This class implements a randomfieldgenerator which reads 
 * an externally generated field interpolates to determine 
 * at Gausspoints.
 */
class ExternalFieldGenerator : public RandomFieldGenerator
{
protected:
    Domain *domain;
    FloatArray field;
    IntArray numberReal;
public:
    /// Constructor. Creates empty RandomFieldGenerator
    ExternalFieldGenerator(int num, Domain *d);
    /// Destructor
    virtual ~ExternalFieldGenerator();

    virtual void generateRandomValue(double &value, FloatArray *position);

    virtual IRResultType initializeFrom(InputRecord *ir);
};
} // end namespace oofem

#endif
