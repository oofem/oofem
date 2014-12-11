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

#ifndef parallelcontext_h
#define parallelcontext_h

#include "oofemcfg.h"
#ifdef __PARALLEL_MODE
 #include "parallelordering.h"
#endif

namespace oofem {
class EngngModel;
class FloatArray;
class DofManager;

/**
 * This class provides an communication context for distributed memory parallelism.
 * Tasks:
 * - Keeping track of the parallel communicator.
 * - Determining owner for shared dof managers.
 */
class OOFEM_EXPORT ParallelContext
{
protected:
    int di;
    EngngModel *emodel;
#ifdef __PARALLEL_MODE
    Natural2GlobalOrdering n2g;
    Natural2LocalOrdering n2l;
#endif

public:
    /**
     * Creates a context belonging to a system of equations in a given engineering model.
     * @param e Engineering model to work with.
     */
    ParallelContext(EngngModel * e);
    ~ParallelContext();

    /**
     * Initiates the mapping for given domain.
     * @param di Domain index.
     */
    void init(int newDi);

    int giveNumberOfLocalEqs();
    int giveNumberOfGlobalEqs();
    int giveNumberOfNaturalEqs();

    ///@name Convenience functions for working with distributed arrays.
    //@{
    bool isLocal(DofManager *dman);
    /**
     * Norm for a locally distributed array.
     * Common for convergence criterion and such.
     */
    double localNorm(const FloatArray &src);

    /**
     * Dot product for a distributed array.
     * Common for convergence criterion and such.
     */
    double dotProduct(const FloatArray &a, const FloatArray &b);
    /**
     * Dot product for a locally distributed array.
     * Common for convergence criterion and such.
     */
    double localDotProduct(const FloatArray &a, const FloatArray &b);
    /**
     * Accumulates the global value.
     */
    double accumulate(double local);
    /**
     * Accumulates the global value.
     */
    void accumulate(const FloatArray &local, FloatArray &global);
    //@}

#ifdef __PARALLEL_MODE
    Natural2GlobalOrdering *giveN2Gmap() { return & n2g; }
    Natural2LocalOrdering *giveN2Lmap() { return & n2l; }
#endif
};
} // end namespace oofem

#endif // parallelcontext_h
