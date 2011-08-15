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

#ifndef petscscontext_h
#define petscscontext_h

#ifdef __PETSC_MODULE

 #ifndef __MAKEDPEND
  #include "petscksp.h"
  #include "petscordering.h"
 #endif

 #include "equationid.h"

namespace oofem {
class EngngModel;
class FloatArray;

/**
 * This class provides an communication context to PETSc library.
 */
class PetscContext
{
protected:
    EngngModel *emodel;
    EquationID ut;
    VecScatter n2gvecscat;
    VecScatter l2gvecscat;
 #ifdef __PARALLEL_MODE
    PetscNatural2GlobalOrdering n2g;
    PetscNatural2LocalOrdering n2l;

    PetscNatural2GlobalOrdering n2g_prescribed;
    PetscNatural2LocalOrdering n2l_prescribed;
 #endif

public:
    PetscContext(EngngModel *e, EquationID ut);
    ~PetscContext();

    void init(int di);

    int giveNumberOfLocalEqs();
    int giveNumberOfGlobalEqs();
    int giveNumberOfNaturalEqs();


    /// Scatters global vector to natural one.
    int scatterG2N(Vec src, Vec dest, InsertMode mode);
    int scatterG2N(Vec src, FloatArray *dest, InsertMode mode);
    /// Scatters and gathers vector in natural ordering (sequential) to global (parallel) one.
    int scatterN2G(Vec src, Vec dest, InsertMode mode);
    int scatterN2G(FloatArray *src, Vec dest, InsertMode mode);
    /**
     * Scatters and gathers vector in natural ordering to global (parallel) one,
     * but only local entries are processed.
     */
    int scatterL2G(FloatArray *src, Vec dest, InsertMode mode);

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
