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

#ifndef eleminterpunknownmapper_h
#define eleminterpunknownmapper_h

#include "primaryunknownmapper.h"

namespace oofem {
class Domain;
class Element;
class TimeStep;

/**
 * The class implementing the primary unknown mapper using element interpolation functions.
 * The basic task is to map the primary unknowns from one (old) mesh to the new one.
 * This task requires the special element algorithms, these are to be included using interface concept.
 */
class OOFEM_EXPORT EIPrimaryUnknownMapper : public PrimaryUnknownMapper
{
public:
    /// Constructor
    EIPrimaryUnknownMapper();
    /// Destructor
    virtual ~EIPrimaryUnknownMapper() { }

    virtual int mapAndUpdate(FloatArray &answer, ValueModeType mode,
                             Domain *oldd, Domain *newd,  TimeStep *tStep);
    virtual int evaluateAt(FloatArray &answer, IntArray &dofMask, ValueModeType mode,
                           Domain *oldd, FloatArray &coords, IntArray &regList, TimeStep *tStep);
};
} // end namespace oofem
#endif // eleminterpunknownmapper_h
