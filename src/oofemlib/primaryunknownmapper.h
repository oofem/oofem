/* $Header: /home/cvs/bp/oofem/oofemlib/src/primaryunknownmapper.h,v 1.5.4.1 2004/04/05 15:19:43 bp Exp $ */
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


//   ************************************
//   *** CLASS PRIMARY UNKNOWN MAPPER ***
//   ************************************

#ifndef primaryunknownmapper_h
#define primaryunknownmapper_h

#include "compiler.h"
#include "cltypes.h"
#include "interface.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#endif

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
class PrimaryUnknownMapper
{
protected:

public:
    /// Constructor
    PrimaryUnknownMapper() { }
    /// Destructor
    virtual ~PrimaryUnknownMapper() { }
    /** Maps and updates the vector(s) of primary unknowns from old mesh oldd to new mesh newd.
     * The physical meaning of primary unknowns is determined by corresponding DofManagers.
     * The result is stored in answer array.
     * The ordering of unknowns in answer is determined by code numbers of
     * new mesh dofmanagers.
     * @param answer contains results
     * @param type determines the physical meaning of mapped unknown
     * @param mode determines the mode of unknown
     * @param ut   determines unknown type
     * @param oldd old mesh reference
     * @param newd new mesh reference
     * @param tStep time step
     * @return nonzero if o.k.
     */
    virtual int mapAndUpdate(FloatArray &answer, ValueModeType mode, EquationID ut,
                             Domain *oldd, Domain *newd,  TimeStep *tStep) = 0;
    /**
     * Evaluates the vector of primary unknowns, determined by domain, at given point.
     * The physical meaning of primary unknowns mapped is determined by bacground
     * element containing given point.
     * @param answer contains evaluated unknown vector
     * @param output parameter comtaining dofIDs of mapped values
     * @param ut   determines unknown type
     * @param mode determines the type of unknown
     * @param oldd old mesh reference (mesh with unknown field)
     * @param coords coordinates of point of interest
     * @param regList - list of regions where to search, if empty all region serach performed.
     * @param tStep solution step
     * @return nonzero if o.k.
     */
    virtual int evaluateAt(FloatArray &answer, IntArray &dofMAsk, EquationID ut, ValueModeType mode,
                           Domain *oldd, FloatArray &coords, IntArray &regList, TimeStep *tStep) = 0;
protected:
    /// prints error message and exits
    void error(const char *file, int line, const char *format, ...) const;
    /// prints warning message
    void warning(const char *file, int line, const char *format, ...) const;
};

#endif // primaryunknownmapper_h






