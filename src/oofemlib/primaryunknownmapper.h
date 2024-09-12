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

#ifndef primaryunknownmapper_h
#define primaryunknownmapper_h

#include "oofemenv.h"
#include "interface.h"
#include "valuemodetype.h"

#include <string>

namespace oofem {
class Domain;
class Element;
class TimeStep;
class FloatArray;
class IntArray;

/**
 * The base class for all primary unknowns mappers.
 * The basic task is to map the primary unknowns from one (old) mesh to the new one.
 * If this task requires the special element algorithms, these should be included using interface concept.
 *
 */
class OOFEM_EXPORT PrimaryUnknownMapper
{
public:
    /// Constructor
    PrimaryUnknownMapper() { }
    /// Destructor
    virtual ~PrimaryUnknownMapper() { }
    /**
     * Maps and updates the vector(s) of primary unknowns from old mesh oldd to new mesh newd.
     * The result is stored in answer array.
     * The interpolation of the primary unknowns is determined by element interpolation.
     * The physical meaning of primary unknowns is determined by DofManagers.
     * The ordering of unknowns in answer is determined by code numbers of
     * new mesh dofmanagers.
     * @param answer Resulting array with primary unknowns.
     * @param mode Determines the mode of unknown.
     * @param oldd Old mesh reference.
     * @param newd New mesh reference.
     * @param tStep Time step.
     * @return Nonzero if o.k.
     */
    virtual int mapAndUpdate(FloatArray &answer, ValueModeType mode,
                             Domain *oldd, Domain *newd,  TimeStep *tStep) = 0;
    /**
     * Evaluates the vector of primary unknowns, determined by domain, at given point.
     * The physical meaning of primary unknowns mapped is determined by background.
     * element containing given point.
     * @param answer Contains evaluated unknown vector.
     * @param dofMask Parameter containing dofIDs of mapped values.
     * @param mode Determines the type of mode of unknown.
     * @param oldd Old mesh reference (mesh with unknown field).
     * @param coords Coordinates of point of interest.
     * @param regList List of regions where to search, if empty all region search performed.
     * @param tStep Solution step.
     * @return Nonzero if o.k.
     */
    virtual int evaluateAt(FloatArray &answer, IntArray &dofMask, ValueModeType mode,
                           Domain *oldd, const FloatArray &coords, IntArray &regList, TimeStep *tStep) = 0;
protected:
    /// Returns string for prepending output (used by error reporting macros).
    std :: string errorInfo(const char *func) const;
};
} // end namespace oofem
#endif // primaryunknownmapper_h
