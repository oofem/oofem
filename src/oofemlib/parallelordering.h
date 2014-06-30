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

#ifndef parallelordering_h
#define parallelordering_h

#include "oofemcfg.h"
#include "intarray.h"
#include "dofmanager.h"

#include <map>

namespace oofem {
class EngngModel;
class IntArray;
class UnknownNumberingScheme;

class OOFEM_EXPORT ParallelOrdering
{
public:
    ParallelOrdering() { }
    virtual ~ParallelOrdering() { }

    /// Returns true if given DofManager is local (ie maintained by the receiver processor).
    bool isLocal(DofManager *dman);
    /// Returns true if given DofManager is shared between partitions.
    bool isShared(DofManager *dman);

    /**
     * Initiates the receiver.
     * @param em Engineering model to determine general information about the problem.
     * @param di Domain index.
     */
    virtual void init(EngngModel *em, int di, const UnknownNumberingScheme &u) = 0;

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

/**
 * Ordering from oofem natural ordering (includes all local and shared eqs)
 * to global ordering.
 */
class OOFEM_EXPORT Natural2GlobalOrdering : public ParallelOrdering
{
protected:
    /// Old to new mapping; uses 0-based global eq ordering; 1-based local ordering.
    IntArray locGlobMap;
    /// New to old mapping.
    std :: map< int, int >globLocMap;

    /// Number of local and global eqs.
    int l_neqs, g_neqs;

public:
    Natural2GlobalOrdering();
    virtual ~Natural2GlobalOrdering() { }

    virtual void init(EngngModel *, int di, const UnknownNumberingScheme &n);

    virtual int giveNewEq(int leq);
    virtual int giveOldEq(int eq);

    virtual void map2New(IntArray &answer, const IntArray &src, int baseOffset = 0);
    virtual void map2Old(IntArray &answer, const IntArray &src, int baseOffset = 0);

    int giveNumberOfLocalEqs() { return l_neqs; }
    int giveNumberOfGlobalEqs() { return g_neqs; }

    IntArray *giveN2Gmap() { return & locGlobMap; }
};

/**
 * Ordering from oofem natural ordering (includes all local and shared eqs)
 * to local ordering, where only locally maintained eqs are considered.
 */
class OOFEM_EXPORT Natural2LocalOrdering : public ParallelOrdering
{
protected:
    /// Natural to local
    IntArray n2l;

public:
    Natural2LocalOrdering();
    virtual ~Natural2LocalOrdering() { }

    virtual void init(EngngModel *, int di, const UnknownNumberingScheme &n);

    virtual int giveNewEq(int leq);
    virtual int giveOldEq(int eq);

    virtual void map2New(IntArray &answer, const IntArray &src, int baseOffset = 0);
    virtual void map2Old(IntArray &answer, const IntArray &src, int baseOffset = 0);

    IntArray *giveN2Lmap() { return & n2l; }
};
} // end namespace oofem

#endif // parallelordering_h
