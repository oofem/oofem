/* $Header: /home/cvs/bp/oofem/sm/src/eleminterpunknownmapper.h,v 1.3 2003/04/06 14:08:30 bp Exp $ */
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

//   **********************************************************
//   *** CLASS ELEMENT INTERPOLATION PRIMARY UNKNOWN MAPPER ***
//   **********************************************************

#ifndef eleminterpunknownmapper_h
#define eleminterpunknownmapper_h

#include "primaryunknownmapper.h"
#include "compiler.h"

namespace oofem {

class Domain;
class Element;
class TimeStep;


/**
 * The class implementing the primary unknown mapper using element intrpolation functions.
 * The basic task is to map the primary unknowns from one (old) mesh to the new one.
 * This task requires the special element algorithms, these are to be included using interface concept.
 * Requires the element support via EIPrimaryUnknownMapperInterface.
 */
class EIPrimaryUnknownMapper : public PrimaryUnknownMapper
{
protected:

public:
    /// Constructor
    EIPrimaryUnknownMapper();
    /// Destructor
    ~EIPrimaryUnknownMapper() { }
    /** Maps and updates the vector(s) of primary unknowns from old mesh oldd to new mesh newd.
     * The result is stored in answer array.
     * The interpolation of the primary unknowns is determined by element interpolation.
     * The physical meaning of primary unknowns is determined by DofManagers.
     * The ordering of unknowns in answer is determined by code numbers of
     * new mesh dofmanagers.
     * @param answer contains results
     * @param mode determines the mode of unknown
     * @param ut   determines unknown type
     * @param oldd old mesh reference
     * @param newd new mesh reference
     * @param tStep time step
     * @return nonzero if o.k.
     */
    virtual int mapAndUpdate(FloatArray &answer, ValueModeType mode, EquationID ut,
                             Domain *oldd, Domain *newd,  TimeStep *tStep);
    /**
     * Evaluates the vector of primary unknowns, determined by domain, at given point.
     * The physical meaning of primary unknowns mapped is determined by bacground
     * element containing given point.
     * @param answer contains evaluated unknown vector
     * @param output parameter comtaining dofIDs of mapped values
     * @param mode determines the type of unknown
     * @param oldd old mesh reference (mesh with unknown field)
     * @param coords coordinates of point of interest
     * @param regList - list of regions where to search, if empty all region serach performed.
     * @param tStep solution step
     * @return nonzero if o.k.
     */
    virtual int evaluateAt(FloatArray &answer, IntArray &dofMAsk, EquationID ut, ValueModeType mode,
                           Domain *oldd, FloatArray &coords, IntArray &regList, TimeStep *tStep);

protected:
};

} // end namespace oofem
#endif // eleminterpunknownmapper_h
