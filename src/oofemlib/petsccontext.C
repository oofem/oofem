/* $Header: /home/cvs/bp/oofem/oofemlib/src/Attic/petsccontext.C,v 1.1.2.1 2004/04/05 15:19:43 bp Exp $ */
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
#ifdef __PETSC_MODULE
#include "petsccontext.h"
#include "engngm.h"

//#define PetscContext_debug_print

PetscContext :: PetscContext(EngngModel *e, EquationID ut)
#ifdef __PARALLEL_MODE
    : n2g(), n2l()
#endif
{
    this->emodel = e;
    this->ut = ut;
    n2gvecscat = NULL;
}

PetscContext :: ~PetscContext()
{
    if ( n2gvecscat ) {
        VecScatterDestroy(n2gvecscat);
    }
}

void
PetscContext :: init(int di)
{
#ifdef __PARALLEL_MODE
    if ( emodel->isParallel() ) {
        n2g.init(emodel, ut, di);
        n2l.init(emodel, ut, di);

        n2g_prescribed.init(emodel, ut, di, ApplicationOrdering :: et_prescribed);
        n2l_prescribed.init(emodel, ut, di, ApplicationOrdering :: et_prescribed);

#ifdef  PetscContext_debug_print
        fprintf( stderr, "[%d] Petsccontext::init leq:%d, neq:%d, geq:%d\n", emodel->giveRank(), giveNumberOfLocalEqs(), giveNumberOfNaturalEqs(), giveNumberOfGlobalEqs() );
#endif
    }

#endif
    if ( n2gvecscat ) {
        VecScatterDestroy(n2gvecscat);
        n2gvecscat = NULL;
    }
}



int
PetscContext :: giveNumberOfLocalEqs()
{
#ifdef __PARALLEL_MODE
    if ( emodel->isParallel() ) {
        return n2g.giveNumberOfLocalEqs();
    } else {
#endif
    return emodel->giveNumberOfEquations(ut);

#ifdef __PARALLEL_MODE
}
#endif
}


int
PetscContext :: giveNumberOfGlobalEqs()
{
#ifdef __PARALLEL_MODE
    if ( emodel->isParallel() ) {
        return n2g.giveNumberOfGlobalEqs();
    } else {
#endif
    return emodel->giveNumberOfEquations(ut);

#ifdef __PARALLEL_MODE
}
#endif
}

int
PetscContext :: giveNumberOfNaturalEqs()
{
    return emodel->giveNumberOfEquations(ut);
}


int
PetscContext :: scatterG2N(Vec src, Vec dest, InsertMode mode)
{
#ifdef __PARALLEL_MODE
    if ( emodel->isParallel() ) {
        int neqs = giveNumberOfNaturalEqs();

        if ( n2gvecscat == NULL ) {
            //
            IS naturalIS, globalIS;
            ISCreateGeneral(PETSC_COMM_WORLD, neqs, this->giveN2Gmap()->giveN2Gmap()->givePointer(), & globalIS);
            ISCreateStride(PETSC_COMM_WORLD, neqs, 0, 1, & naturalIS);
            VecScatterCreate(dest, naturalIS, src, globalIS, & n2gvecscat);
        }

        VecScatterBegin(n2gvecscat, src, dest, mode, SCATTER_REVERSE); //
        return VecScatterEnd(n2gvecscat, src, dest, mode, SCATTER_REVERSE); //
    } else {
#endif
    return VecCopy(src, dest);

#ifdef __PARALLEL_MODE
}
#endif
}


int
PetscContext :: scatterG2N(Vec src, FloatArray *dest, InsertMode mode)
{
    PetscScalar *ptr;
    int i;

#ifdef __PARALLEL_MODE
    if ( emodel->isParallel() ) {
        int neqs = giveNumberOfNaturalEqs();
        Vec natVec;
        VecCreateSeq(PETSC_COMM_SELF, neqs, & natVec);

        this->scatterG2N(src, natVec, mode);

        dest->resize( giveNumberOfNaturalEqs() );
        VecGetArray(natVec, & ptr);
        for ( i = 0; i < neqs; i++ ) {
            dest->at(i + 1) = ptr [ i ];
        }

        VecRestoreArray(natVec, & ptr);
        VecDestroy(natVec);
    } else {
#endif
    int neqs = giveNumberOfNaturalEqs();
    dest->resize(neqs);
    VecGetArray(src, & ptr);
    for ( i = 0; i < neqs; i++ ) {
        dest->at(i + 1) = ptr [ i ];
    }

    VecRestoreArray(src, & ptr);
#ifdef __PARALLEL_MODE
}
#endif
    return 1;
}

int
PetscContext :: scatterN2G(Vec src, Vec dest, InsertMode mode)
{
#ifdef __PARALLEL_MODE
    if ( emodel->isParallel() ) {
        int neqs = giveNumberOfNaturalEqs();
        if ( n2gvecscat == NULL ) {
            //
            IS naturalIS, globalIS;
            ISCreateGeneral(PETSC_COMM_WORLD, neqs, this->giveN2Gmap()->giveN2Gmap()->givePointer(), & globalIS);
            ISCreateStride(PETSC_COMM_WORLD, neqs, 0, 1, & naturalIS);
            VecScatterCreate(src, naturalIS, dest, globalIS, & n2gvecscat);
        }

        VecScatterBegin(n2gvecscat, src, dest, mode, SCATTER_FORWARD); //
        VecScatterEnd(n2gvecscat, src, dest, mode, SCATTER_FORWARD); //
    } else {
#endif
    VecCopy(src, dest);
#ifdef __PARALLEL_MODE
}
#endif
    return 1;
}


int
PetscContext :: scatterN2G(FloatArray *src, Vec dest, InsertMode mode)
{
    PetscScalar *ptr;
    int i;

#ifdef __PARALLEL_MODE
    if ( emodel->isParallel() ) {
        int size = src->giveSize();
        Vec natVec;
        VecCreateSeq(PETSC_COMM_SELF, giveNumberOfNaturalEqs(), & natVec);

        ptr = src->givePointer();
        for ( i = 0; i < size; i++ ) {
            VecSetValues(natVec, 1, & i, ptr + i, ADD_VALUES);
        }

        VecAssemblyBegin(natVec);
        VecAssemblyEnd(natVec);
        //VecView(natVec,PETSC_VIEWER_STDOUT_SELF);

        this->scatterN2G(natVec, dest, mode);
        VecDestroy(natVec);
    } else {
#endif
    int size = src->giveSize();
    ptr = src->givePointer();
    for ( i = 0; i < size; i++ ) {
        VecSetValues(dest, 1, & i, ptr + i, ADD_VALUES);
    }

    VecAssemblyBegin(dest);
    VecAssemblyEnd(dest);
#ifdef __PARALLEL_MODE
}
#endif
    return 1;
}


int
PetscContext :: scatterL2G(FloatArray *src, Vec dest, InsertMode mode)
{
    PetscScalar *ptr;
    int i;

#ifdef __PARALLEL_MODE
    if ( emodel->isParallel() ) {
        int eqg, size = src->giveSize();
        ptr = src->givePointer();

        PetscNatural2LocalOrdering *n2l = this->giveN2Lmap();
        PetscNatural2GlobalOrdering *n2g = this->giveN2Gmap();
        for ( i = 0; i < size; i++ ) {
            if ( n2l->giveNewEq(i + 1) ) {
                eqg = n2g->giveNewEq(i + 1);
                VecSetValues(dest, 1, & eqg, ptr + i, ADD_VALUES);
            }
        }

        VecAssemblyBegin(dest);
        VecAssemblyEnd(dest);
        //VecView(natVec,PETSC_VIEWER_STDOUT_SELF);
    } else {
#endif

    int size = src->giveSize();
    ptr = src->givePointer();
    for ( i = 0; i < size; i++ ) {
        VecSetValues(dest, 1, & i, ptr + i, ADD_VALUES);
    }

    VecAssemblyBegin(dest);
    VecAssemblyEnd(dest);
#ifdef __PARALLEL_MODE
}
#endif
    return 1;
}



void
PetscContext :: createVecGlobal(Vec *answer)
{
#ifdef __PARALLEL_MODE
    if ( emodel->isParallel() ) {
        VecCreate(PETSC_COMM_WORLD, answer);
        VecSetSizes( * answer, giveNumberOfLocalEqs(), giveNumberOfGlobalEqs() );
        VecSetFromOptions(* answer);
    } else {
#endif
    VecCreateSeq(PETSC_COMM_SELF, giveNumberOfNaturalEqs(), answer);
#ifdef __PARALLEL_MODE
}
#endif
}

#endif
