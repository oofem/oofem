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

#ifndef eleminterpunknownmapper_h
#define eleminterpunknownmapper_h

#include "primaryunknownmapper.h"
#include "compiler.h"

namespace oofem {
class Domain;
class Element;
class TimeStep;

/**
 * The class implementing the primary unknown mapper using element interpolation functions.
 * The basic task is to map the primary unknowns from one (old) mesh to the new one.
 * This task requires the special element algorithms, these are to be included using interface concept.
 * Requires the element support via EIPrimaryUnknownMapperInterface.
 */
class EIPrimaryUnknownMapper : public PrimaryUnknownMapper
{
public:
    /// Constructor
    EIPrimaryUnknownMapper();
    /// Destructor
    virtual ~EIPrimaryUnknownMapper() { }

    virtual int mapAndUpdate(FloatArray &answer, ValueModeType mode, EquationID ut,
                             Domain *oldd, Domain *newd,  TimeStep *tStep);
    virtual int evaluateAt(FloatArray &answer, IntArray &dofMask, EquationID ut, ValueModeType mode,
                           Domain *oldd, FloatArray &coords, IntArray &regList, TimeStep *tStep);
};
} // end namespace oofem
#endif // eleminterpunknownmapper_h
