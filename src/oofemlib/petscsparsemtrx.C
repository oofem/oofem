/* $Header: /home/cvs/bp/oofem/oofemlib/src/Attic/petscsparsemtrx.C,v 1.1.2.1 2004/04/05 15:19:43 bp Exp $ */
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

#include "parallel.h"
#include "petscsparsemtrx.h"
#include "petscordering.h"
#include "engngm.h"

#ifndef __MAKEDEPEND
#include <stdio.h>
#include <set>
#include <vector>
#include "petscksp.h"
#endif

namespace oofem {

SparseMtrx *
PetscSparseMtrx :: GiveCopy() const
{
    OOFEM_ERROR("PetscSparseMtrx :: GiveCopy - Not implemented");
    return NULL;
}

void
PetscSparseMtrx :: times(const FloatArray &x, FloatArray &answer) const
{
    OOFEM_ERROR("PetscSparseMtrx :: times - Not implemented");
}


void
PetscSparseMtrx :: times(double x)
{
    OOFEM_ERROR("PetscSparseMtrx::times(double x) - unsupported");
}

int
PetscSparseMtrx :: buildInternalStructure(EngngModel *eModel, int di, EquationID ut, const UnknownNumberingScheme&s)
{
    IntArray loc;
    Domain *domain = eModel->giveDomain(di);
    int nelem;

    if ( mtrx ) {
        MatDestroy(mtrx);
    }

    this->ut = ut;
    this->emodel = eModel;
    this->di = di;


#ifdef __PARALLEL_MODE
    int rank;
    PetscNatural2GlobalOrdering *n2g = NULL;
    rank = eModel->giveRank();
    if ( eModel->isParallel() ) {
        n2g = eModel->givePetscContext(di, ut)->giveN2Gmap();
    }

#ifdef __VERBOSE_PARALLEL
    VERBOSEPARALLEL_PRINT("PetscSparseMtrx:: buildInternalStructure", "", rank);
#endif


    // initialize n2l map
    PetscNatural2LocalOrdering n2l;
    n2l.init(eModel, ut, di);

    // initialize n2g map
    // n2g.init(eModel, di);


    leqs = n2g->giveNumberOfLocalEqs();
    geqs = n2g->giveNumberOfGlobalEqs();

#ifdef __VERBOSE_PARALLEL
    OOFEM_LOG_INFO( "[%d]PetscSparseMtrx:: buildInternalStructure: l_eqs = %d, g_eqs = %d, n_eqs = %d\n", rank, leqs, geqs, eModel->giveNumberOfEquations(ut) );
#endif

#else
    leqs = geqs = eModel->giveNumberOfEquations(ut);
    //map = new LocalizationMap();
#endif


    //determine nonzero structure of a "local (maintained)" part of matrix

#ifdef __PARALLEL_MODE
    if ( eModel->isParallel() ) {
        int i, ii, j, jj, n;
        Element *elem;
        // allocation map
        std :: vector< std :: set< int > > rows(leqs); // ??

        IntArray d_nnz(leqs), lloc;

        //fprintf (stderr,"[%d] n2l map: ",rank);
        //for (n=1; n<=n2l.giveN2Lmap()->giveSize(); n++) fprintf (stderr, "%d ", n2l.giveN2Lmap()->at(n));

        nelem = domain->giveNumberOfElements();
        for ( n = 1; n <= nelem; n++ ) {
            //fprintf (stderr, "(elem %d) ", n);
            elem = domain->giveElement(n);
            elem->giveLocationArray(loc, ut,s);
            n2l.map2New(lloc, loc, 0); // translate natural->local numbering
            for ( i = 1; i <= lloc.giveSize(); i++ ) {
                if ( ( ii = lloc.at(i) ) ) {
                    for ( j = 1; j <= lloc.giveSize(); j++ ) {
                        if ( ( jj = lloc.at(j) ) ) {
                            //fprintf (stderr, "{[%d] %d-%d %d-%d} ", rank, loc.at(i),ii-1,loc.at(j),jj-1);
                            rows [ ii - 1 ].insert(jj - 1);
                        }
                    }
                }
            }

            //fprintf (stderr, "\n");
        }

        for ( i = 0; i < leqs; i++ ) {
            d_nnz(i) = rows [ i ].size();
        }

        //fprintf (stderr,"\n[%d]PetscSparseMtrx: Profile ...",rank);
        //for (i=0; i<leqs; i++) fprintf(stderr, "%d ", d_nnz(i));



        //fprintf (stderr,"\n[%d]PetscSparseMtrx: Creating MPIAIJ Matrix ...\n",rank);

        // create PETSc mat
        /*
         * MatCreateMPIAIJ(PETSC_COMM_WORLD,leqs,leqs,geqs,geqs,
         *              leqs,d_nnz.givePointer(), // counted
         *              6,PETSC_NULL, // guess
         *              &mtrx);
         */
        MatCreate(PETSC_COMM_WORLD, & mtrx);
        MatSetSizes(mtrx, leqs, leqs, geqs, geqs);
        MatSetType(mtrx, MATMPIAIJ);
        MatSetFromOptions(mtrx);
        MatMPIAIJSetPreallocation(mtrx, 1, d_nnz.givePointer(), 6, PETSC_NULL);

#ifdef __VERBOSE_PARALLEL
        VERBOSEPARALLEL_PRINT("PetscSparseMtrx:: buildInternalStructure", "done", rank);
#endif
    } else {
#endif


    int i, ii, j, jj, n;
    Element *elem;
    // allocation map
    std :: vector< std :: set< int > > rows(leqs); // ??
    IntArray d_nnz(leqs);

    nelem = domain->giveNumberOfElements();
    for ( n = 1; n <= nelem; n++ ) {
        elem = domain->giveElement(n);
        elem->giveLocationArray(loc, ut,s);
        for ( i = 1; i <= loc.giveSize(); i++ ) {
            if ( ( ii = loc.at(i) ) ) {
                for ( j = 1; j <= loc.giveSize(); j++ ) {
                    jj = loc.at(j);
                    if ( ( jj = loc.at(j) ) ) {
                        rows [ ii - 1 ].insert(jj - 1);
                    }
                }
            }
        }
    }

    for ( i = 0; i < leqs; i++ ) {
        d_nnz(i) = rows [ i ].size();
    }

    // create PETSc mat
    /*
     * MatCreateSeqAIJ(PETSC_COMM_SELF,leqs,leqs,0,d_nnz.givePointer(),&mtrx);
     */
    MatCreate(PETSC_COMM_WORLD, & mtrx);
    MatSetSizes(mtrx, leqs, leqs, leqs, leqs);
    MatSetType(mtrx, MATSEQAIJ);
    MatSetFromOptions(mtrx);
    MatSeqAIJSetPreallocation( mtrx, 0, d_nnz.givePointer() );

#ifdef __PARALLEL_MODE
}
#endif

    nRows = nColumns = geqs;
    return TRUE;
}

int
PetscSparseMtrx :: assemble(const IntArray &loc, const FloatMatrix &mat)
{
    int ndofe = mat.giveNumberOfRows();

#ifdef __PARALLEL_MODE
    if ( emodel->isParallel() ) {
        // translate local code numbers to global ones
        IntArray gloc(ndofe);
        emodel->givePetscContext(this->di, ut)->giveN2Gmap()->map2New(gloc, loc, 0);

        //fprintf (stderr, "[?] gloc=");
        //for (int i=1; i<=ndofe; i++) fprintf (stderr, "%d ", gloc.at(i));

        // To allow the insertion of values using MatSetValues in column major order
        MatSetOption(mtrx, MAT_ROW_ORIENTED, PETSC_FALSE);
        MatSetValues(this->mtrx, ndofe, gloc.givePointer(), ndofe, gloc.givePointer(), mat.givePointer(), ADD_VALUES);
    } else {
#endif
    IntArray gloc(ndofe);
    for ( int i = 1; i <= ndofe; i++ ) {
        gloc.at(i) = loc.at(i) - 1;
    }

    // To allow the insertion of values using MatSetValues in column major order
    MatSetOption(mtrx, MAT_ROW_ORIENTED, PETSC_FALSE);
    MatSetValues(this->mtrx, ndofe, gloc.givePointer(), ndofe, gloc.givePointer(), mat.givePointer(), ADD_VALUES);

    //mat.printYourself();
    //loc.printYourself();
    //MatView(this->mtrx,PETSC_VIEWER_STDOUT_SELF);
#ifdef __PARALLEL_MODE
}
#endif
    this->version++;
    return 1;
}

int
PetscSparseMtrx :: assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat)
{
#ifdef __PARALLEL_MODE
    if ( emodel->isParallel() ) {
        // translate eq numbers
        IntArray grloc( rloc.giveSize() ), gcloc( cloc.giveSize() );
        emodel->givePetscContext(this->di, ut)->giveN2Gmap()->map2New(grloc, rloc, 0);
        emodel->givePetscContext(this->di, ut)->giveN2Gmap()->map2New(gcloc, cloc, 0);

        // To allow the insertion of values using MatSetValues in column major order
        MatSetOption(mtrx, MAT_ROW_ORIENTED, PETSC_FALSE);
        MatSetValues(this->mtrx, grloc.giveSize(), grloc.givePointer(),
                     gcloc.giveSize(), gcloc.givePointer(), mat.givePointer(), ADD_VALUES);
    } else {
#endif
    int i, rsize = rloc.giveSize(), csize = cloc.giveSize();
    IntArray grloc(rsize), gcloc(csize);
    for ( i = 1; i <= rsize; i++ ) {
        grloc.at(i) = rloc.at(i) - 1;
    }

    for ( i = 1; i <= csize; i++ ) {
        gcloc.at(i) = cloc.at(i) - 1;
    }

    // To allow the insertion of values using MatSetValues in column major order
    MatSetOption(mtrx, MAT_ROW_ORIENTED, PETSC_FALSE);
    MatSetValues(this->mtrx, rsize, grloc.givePointer(),
                 csize, gcloc.givePointer(), mat.givePointer(), ADD_VALUES);

#ifdef __PARALLEL_MODE
}
#endif
    // increment version
    this->version++;

    return 1;
}

int
PetscSparseMtrx :: assembleBegin()
{
    return MatAssemblyBegin(this->mtrx, MAT_FINAL_ASSEMBLY);
}

int
PetscSparseMtrx :: assembleEnd()
{
    return MatAssemblyEnd(this->mtrx, MAT_FINAL_ASSEMBLY);
}




SparseMtrx *
PetscSparseMtrx :: zero()
{
    // test if receiver is aslready assembled
    PetscTruth assembled;
    MatAssembled(this->mtrx, & assembled);
    if ( assembled ) {
        MatZeroEntries(this->mtrx);
    }

    return this;
}

double &
PetscSparseMtrx :: at(int i, int j)
{
    static double a;
    OOFEM_ERROR("PetscSparseMtrx::at(i,j) - unsupported");
    return a;
}

double
PetscSparseMtrx :: at(int i, int j) const
{
    OOFEM_ERROR("PetscSparseMtrx::at(i,j) - unsupported");
    return 0;
}

void
PetscSparseMtrx :: toFloatMatrix(FloatMatrix &answer) const
{
    OOFEM_ERROR("PetscSparseMtrx::toFloatMatrix() - unsupported");
}

void
PetscSparseMtrx :: printStatistics() const
{
    PetscViewerSetFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_INFO);
    MatView(this->mtrx, PETSC_VIEWER_STDOUT_SELF);
}

void
PetscSparseMtrx :: printYourself() const
{
    PetscViewerSetFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_DENSE);
    MatView(this->mtrx, PETSC_VIEWER_STDOUT_SELF);
}

FloatArray
PetscSparseMtrx :: trans_mult(const FloatArray &x) const
{
    OOFEM_ERROR("PetscSparseMtrx::trans_mult() - unsupported");
    return x; // to supress compiler warning
}

} // end namespace oofem
#endif //ifdef __PETSC_MODULE





