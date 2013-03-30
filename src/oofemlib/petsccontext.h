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

#ifndef petscscontext_h
#define petscscontext_h

#ifdef __PETSC_MODULE

 #ifndef __MAKE_DEPEND
  #include <petscvec.h>
 #endif

 #include "petscordering.h"
 #include "equationid.h"

namespace oofem {
class EngngModel;
class FloatArray;
class DofManager;

/**
 * This class provides an communication context to PETSc library.
 * Tasks:
 * - Keeping track of the parallel communicator.
 * - Creating suitable parallel objects.
 * - Determining owner for shared dof managers.
 */
class PetscContext
{
protected:
    int di;
    EngngModel *emodel;
    EquationID ut;
    VecScatter n2gvecscat;
    VecScatter l2gvecscat;
    /// True if vectors are assumed to be natural distribution.
    bool naturalVectors;

    /// Communicator used for parallel objects.
    MPI_Comm comm;

#ifdef __PARALLEL_MODE
    PetscNatural2GlobalOrdering n2g;
    PetscNatural2LocalOrdering n2l;

    PetscNatural2GlobalOrdering n2g_prescribed;
    PetscNatural2LocalOrdering n2l_prescribed;
 #endif

public:
    /**
     * Creates a context belonging to a system of equations in a given engineering model.
     * @param e Engineering model to work with.
     * @param eid Equation ID to work with.
     * @param naturalVectors Should be true if shared dofs only contain the local contributions.
     * Some engineering models manually scatter local vectors to their global value, in which case this would be false.
     */
    PetscContext(EngngModel *e, EquationID eid, bool naturalVectors = true);
    ~PetscContext();

    /**
     * Initiates the mapping for given domain.
     * @param di Domain index.
     */
    void init(int di);

    /**
     * Gives the communicator for parallel (if active).
     * @return Parallel communicator object (typically PETSC_COMM_WORLD). Gives the self communicator if engineering problem isn't parallel.
     */
    MPI_Comm giveComm() { return comm; };

    int giveNumberOfLocalEqs();
    int giveNumberOfGlobalEqs();
    int giveNumberOfNaturalEqs();

    /// Scatters global vector to natural one.
    int scatterG2N(Vec src, Vec dest, InsertMode mode);
    /// Scatters global vector to natural one.
    int scatterG2N(Vec src, FloatArray *dest, InsertMode mode);

    /// Scatters vectors (natural or local) to a global one. Uses the naturalVectors variable to determine.
    int scatter2G(const FloatArray *src, Vec dest, InsertMode mode);

    /// Scatters and gathers vector in natural ordering (sequential) to global (parallel) one.
    int scatterN2G(Vec src, Vec dest, InsertMode mode);
    /// Scatters and gathers vector in natural ordering (sequential) to global (parallel) one.
    int scatterN2G(const FloatArray *src, Vec dest, InsertMode mode);

    /**
     * Scatters and gathers vector in natural ordering to global (parallel) one,
     * but only local entries are processed.
     */
    int scatterL2G(const FloatArray *src, Vec dest, InsertMode mode);

    ///@name Convenience functions for working with distributed arrays.
    //@{
    bool isLocal(DofManager *dman);
    /**
     * Norm for a distributed array.
     * Common for convergence criterion and such.
     */
    double norm(const FloatArray &src);
    /**
     * Norm for a naturally distributed array.
     * Common for convergence criterion and such.
     */
    double naturalNorm(const FloatArray &src);
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
     * Dot product for a naturally distributed array.
     * Common for convergence criterion and such.
     */
    double naturalDotProduct(const FloatArray &a, const FloatArray &b);

    /**
     * Accumulates the global value.
     */
    double accumulate(double local);
    /**
     * Accumulates the global value.
     */
    void accumulate(const FloatArray &local, FloatArray &global);
    //@}

    void createVecGlobal(Vec *answer);

 #ifdef __PARALLEL_MODE
    PetscNatural2GlobalOrdering *giveN2Gmap() { return & n2g; }
    PetscNatural2LocalOrdering *giveN2Lmap() { return & n2l; }

    PetscNatural2GlobalOrdering *giveN2GPrescribedmap() { return & n2g_prescribed; }
    PetscNatural2LocalOrdering *giveN2LPrescribedmap() { return & n2l_prescribed; }

 #endif
};
} // end namespace oofem
#endif
#endif // petscscontext_h
