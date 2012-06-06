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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#ifndef appordering_h
#define appordering_h

#include "equationid.h"

namespace oofem {
class EngngModel;
class IntArray;

/**
 * Base class for ordering of equation for parallelization.
 * @todo Document more
 */
class ApplicationOrdering
{
public:
    enum EquationType { et_standard, et_prescribed };

    ApplicationOrdering() { }
    virtual ~ApplicationOrdering() { }

    /**
     * Initiates the receiver.
     * @param em Engineering model to determine general information about the problem.
     * @param ut Equation to solve.
     * @param di Domain index.
     * @param et Equation type.
     */
    virtual void init(EngngModel *em, EquationID ut, int di, EquationType et = et_standard) = 0;

    /**
     * Returns number of local eqs; ie. those that belong to receiver processor;
     * Note that some eqs may be owned by remote processors (some shared nodes,...).
     * The sum of local eqs for all processors should give total number of eqs.
     * @return Numbering of local equations.
     */
    virtual int giveNumberOfLocalEqs() { return 0; }
    /**
     * @return the total number of equations of the problem.
     */
    virtual int giveNumberOfGlobalEqs() { return 0; }

    /**
     * Finds the global equation from a local equation.
     * @param leq Local equation number.
     * @return Global equation number.
     */
    virtual int giveNewEq(int leq) = 0;
    /**
     * Finds the local equation number from a global equation.
     * @param eq Global equation number.
     * @return Local equation number.
     */
    virtual int giveOldEq(int eq) = 0;

    virtual void map2New(IntArray &answer, const IntArray &src, int baseOffset = 0) = 0;
    virtual void map2Old(IntArray &answer, const IntArray &src, int baseOffset = 0) = 0;
};
} // end namespace oofem
#endif // appordering_h
