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
 #include "petscvec.h"
#endif

namespace oofem {
SparseMtrx *
PetscSparseMtrx :: GiveCopy() const
{
  
  PetscSparseMtrx* answer = new PetscSparseMtrx(nRows, nColumns);
  if (answer) {
    MatDuplicate(this->mtrx, MAT_COPY_VALUES, &(answer->mtrx));
    answer->symmFlag = this->symmFlag;
    answer->mType    = this->mType;
    answer->leqs     = this->leqs;
    answer->geqs     = this->geqs;
    answer->di       = this->di;
    answer->ut       = this->ut;
    answer->emodel   = this->emodel;
    answer->kspInit  = true; // force ksp to be initialized
    answer->newValues= this->newValues;

    return answer;

  } else {
    OOFEM_FATAL("PetscSparseMtrx :: GiveCopy allocation failed");
    return NULL;
  }
}

void
PetscSparseMtrx :: times(const FloatArray &x, FloatArray &answer) const
{
    if ( this->giveNumberOfColumns() != x.giveSize() ) {
        OOFEM_ERROR("Dimension mismatch");
    }

#ifdef __PARALLEL_MODE
    OOFEM_ERROR("PetscSparseMtrx :: times - Not implemented");
#else
    Vec globX, globY;
    VecCreateSeqWithArray(PETSC_COMM_SELF, x.giveSize(), x.givePointer(), & globX);
    VecCreate(PETSC_COMM_SELF, & globY);
    VecSetType(globY, VECSEQ);
    VecSetSizes(globY, PETSC_DECIDE, this->nRows);

    MatMult(this->mtrx, globX, globY);
    double *ptr;
    VecGetArray(globY, & ptr);
    answer.resize(this->nRows);
    for ( int i = 0; i < this->nRows; i++ ) {
        answer(i) = ptr [ i ];
    }

    VecRestoreArray(globY, & ptr);
    VecDestroy(globX);
    VecDestroy(globY);
#endif
}

void
PetscSparseMtrx :: timesT(const FloatArray &x, FloatArray &answer) const
{
    if ( this->giveNumberOfRows() != x.giveSize() ) {
        OOFEM_ERROR("Dimension mismatch");
    }

#ifdef __PARALLEL_MODE
    OOFEM_ERROR("PetscSparseMtrx :: times - Not implemented");
#else
    Vec globX, globY;
    VecCreateSeqWithArray(PETSC_COMM_SELF, x.giveSize(), x.givePointer(), & globX);
    VecCreate(PETSC_COMM_SELF, & globY);
    VecSetType(globY, VECSEQ);
    VecSetSizes(globY, PETSC_DECIDE, this->nColumns);

    MatMultTranspose(this->mtrx, globX, globY);
    double *ptr;
    VecGetArray(globY, & ptr);
    answer.resize(this->nColumns);
    for ( int i = 0; i < this->nColumns; i++ ) {
        answer(i) = ptr [ i ];
    }

    VecRestoreArray(globY, & ptr);
    VecDestroy(globX);
    VecDestroy(globY);
#endif
}


void
PetscSparseMtrx :: times(const FloatMatrix &B, FloatMatrix &answer) const
{
    if ( this->giveNumberOfColumns() != B.giveNumberOfRows() ) {
        OOFEM_ERROR("Dimension mismatch");
    }

#ifdef __PARALLEL_MODE
    OOFEM_ERROR("PetscSparseMtrx :: times - Not implemented");
#else
    // I'm opting to work with a set of vectors, as i think it might be faster and more robust. / Mikael

    int nr = this->giveNumberOfRows();
    int nc = B.giveNumberOfColumns();
    answer.resize(nr, nc);
    double *aptr = answer.givePointer();

    /*
     *      // Approach using several vectors. Not sure if it is optimal, but it includes petsc calls which i suspect are inefficient. / Mikael
     *      // UNTESTED!
     *  Vec globX, globY;
     *  VecCreate(PETSC_COMM_SELF, &globY);
     *  VecSetType(globY, VECSEQ);
     *  VecSetSizes(globY, PETSC_DECIDE, nr);
     *  int nrB = B.giveNumberOfRows();
     *  for (int k = 0; k < nc; k++) {
     *      double colVals[nrB];
     *      for (int i = 0; i < nrB; i++) colVals[i] = B(i,k); // B.copyColumn(Bk,k);
     *      VecCreateSeqWithArray(PETSC_COMM_SELF, nrB, colVals, &globX);
     *      MatMult(this->mtrx, globX, globY );
     *              double *ptr;
     *              VecGetArray(globY, &ptr);
     *              for (int i = 0; i < nr; i++) *aptr++ = ptr[i]; // answer.setColumn(Ak,k);
     *              VecRestoreArray(globY, &ptr);
     *              VecDestroy(globX);
     *  }
     *  VecDestroy(globY);
     */

    Mat globB, globC;
    MatCreateSeqDense(PETSC_COMM_SELF, B.giveNumberOfRows(), B.giveNumberOfColumns(), B.givePointer(), & globB);
    MatMatMult(this->mtrx, globB, MAT_INITIAL_MATRIX, PETSC_DEFAULT, & globC);
    const double *vals;
    for ( int r = 0; r < nr; r++ ) {
        MatGetRow(globC, r, NULL, NULL, & vals);
        for ( int i = 0, i2 = r; i < nc; i++, i2 += nr ) {
            aptr [ i2 ] = vals [ i ];
        }
    }

    MatDestroy(globB);
    MatDestroy(globC);
#endif
}

void
PetscSparseMtrx :: timesT(const FloatMatrix &B, FloatMatrix &answer) const
{
    if ( this->giveNumberOfRows() != B.giveNumberOfRows() ) {
        OOFEM_ERROR("Dimension mismatch");
    }

#ifdef __PARALLEL_MODE
    OOFEM_ERROR("PetscSparseMtrx :: times - Not implemented");
#else
    int nr = this->giveNumberOfColumns();
    int nc = B.giveNumberOfColumns();
    answer.resize(nr, nc);
    double *aptr = answer.givePointer();

    // For some reason SEQAIJ and SEQDENSE are incompatible with each other for MatMatMultTranspose (MatMatMult is OK). I don't know why.
    // Either way, this is not to bad, except for an extra transposition.

    Mat globB, globC;
    FloatMatrix BT;
    BT.beTranspositionOf(B);
    MatCreateSeqDense(PETSC_COMM_SELF, BT.giveNumberOfRows(), BT.giveNumberOfColumns(), BT.givePointer(), & globB);
    MatMatMult(globB, this->mtrx, MAT_INITIAL_MATRIX, PETSC_DEFAULT, & globC);
    const double *vals;
    for ( int r = 0; r < nc; r++ ) {
        MatGetRow(globC, r, NULL, NULL, & vals);
        for ( int i = 0; i < nr; i++ ) {
            * aptr++ = vals [ i ];
        }
    }

    MatDestroy(globB);
    MatDestroy(globC);
#endif
}

void
PetscSparseMtrx :: times(double x)
{
    MatScale(this->mtrx, x);
}

// NOTE! I haven't looked at the parallel code yet (lack of time right now, and i want to see it work first). / Mikael
int
PetscSparseMtrx :: buildInternalStructure(EngngModel *eModel, int di, EquationID ut, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
#ifdef __PARALLEL_MODE
    OOFEM_ERROR("PetscSparseMtrx :: buildInternalStructure - Not implemented");
#else
    Domain *domain = eModel->giveDomain(di);
    int nelem;

    if ( mtrx ) {
        MatDestroy(mtrx);
    }

    this->ut = ut;
    // This should be based on the numberingscheme. Also, geqs seems redundant.

    /*
     * // This could simplify things.
     * if (r_s.isPrescribed())
     *  nRows = geqs = leqs = eModel->giveNumberOfPrescribedEquations(ut);
     * else
     *  nRows = geqs = leqs = eModel->giveNumberOfEquations(ut);
     *
     * if (c_s.isPrescribed())
     *  nColumns = eModel->giveNumberOfPrescribedEquations(ut);
     * else
     *  nColumns = eModel->giveNumberOfEquations(ut);
     */
    int totalEquations = eModel->giveNumberOfEquations(ut) + eModel->giveNumberOfPrescribedEquations(ut);

    //determine nonzero structure of matrix
    int i, ii, j, jj, n;
    Element *elem;
    IntArray r_loc, c_loc;
    std :: vector< std :: set< int > >rows(totalEquations);
    nRows = nColumns = 0;

    nelem = domain->giveNumberOfElements();
    for ( n = 1; n <= nelem; n++ ) {
        elem = domain->giveElement(n);
        elem->giveLocationArray(r_loc, ut, r_s);
        elem->giveLocationArray(c_loc, ut, c_s);
        for ( i = 1; i <= r_loc.giveSize(); i++ ) {
            if ( ( ii = r_loc.at(i) ) ) {
                if ( ii > nRows ) {
                    nRows = ii;
                }
                for ( j = 1; j <= c_loc.giveSize(); j++ ) {
                    if ( ( jj = c_loc.at(j) ) ) {
                        rows [ ii - 1 ].insert(jj - 1);
                    }
                }
            }
        }
        // Column counting should be like this. The matrix should in all reasonable situations get
        // the full corresponding size, regardless if elements are assembled there.
        // Mixed matrices like this are not likely used for solving, so a zero row/column is not unreasonable.
        for ( j = 1; j <= c_loc.giveSize(); j++ ) {
            if ( (jj = c_loc.at(j)) > nColumns ) {
                nColumns = jj;
            }
        }
    }

    geqs = leqs = nRows;

    IntArray d_nnz(leqs);
    for ( i = 0; i < leqs; i++ ) {
        d_nnz(i) = rows [ i ].size();
    }

    // create PETSc mat
    /*
     * MatCreateSeqAIJ(PETSC_COMM_SELF,leqs,leqs,0,d_nnz.givePointer(),&mtrx);
     */
    MatCreate(PETSC_COMM_WORLD, & mtrx);
    // Trying to do this with leqs and geqs doesn't make any sense. They can only mean lines if i understand this correctly (local vs global). / Mikael.
    // Besides, geqs is meaningless, as it will always be the number of rows. Why not lRows and lColumns instead?
    MatSetSizes(mtrx, nRows, nColumns, nRows, nColumns);
    MatSetType(mtrx, MATSEQAIJ);
    MatSetFromOptions(mtrx);
    MatSeqAIJSetPreallocation( mtrx, 0, d_nnz.givePointer() );

#endif
    this->newValues = true;
    return true;
}

int
PetscSparseMtrx :: buildInternalStructure(EngngModel *eModel, int di, EquationID ut, const UnknownNumberingScheme &s)
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
        std :: vector< std :: set< int > >rows(leqs);  // ??

        IntArray d_nnz(leqs), lloc;

        //fprintf (stderr,"[%d] n2l map: ",rank);
        //for (n=1; n<=n2l.giveN2Lmap()->giveSize(); n++) fprintf (stderr, "%d ", n2l.giveN2Lmap()->at(n));

        nelem = domain->giveNumberOfElements();
        for ( n = 1; n <= nelem; n++ ) {
            //fprintf (stderr, "(elem %d) ", n);
            elem = domain->giveElement(n);
            elem->giveLocationArray(loc, ut, s);
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
    std :: vector< std :: set< int > >rows(leqs);  // ??
    IntArray d_nnz(leqs);

    nelem = domain->giveNumberOfElements();
    for ( n = 1; n <= nelem; n++ ) {
        elem = domain->giveElement(n);
        elem->giveLocationArray(loc, ut, s);
        for ( i = 1; i <= loc.giveSize(); i++ ) {
            if ( ( ii = loc.at(i) ) ) {
                for ( j = 1; j <= loc.giveSize(); j++ ) {
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

    //MatCreateSeqAIJ(PETSC_COMM_SELF,leqs,leqs,0,d_nnz.givePointer(),&mtrx);
    MatCreate(PETSC_COMM_WORLD, & mtrx);
    MatSetSizes(mtrx, leqs, leqs, geqs, geqs);
    MatSetType(mtrx, MATSEQAIJ);
    MatSetFromOptions(mtrx);
    MatSeqAIJSetPreallocation( mtrx, 0, d_nnz.givePointer() );

#ifdef __PARALLEL_MODE
}
#endif

    nRows = nColumns = geqs;
    this->newValues = true;
    return true;
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
    this->newValues = true;
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
    this->newValues = true;
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
    this->newValues = true;
    return MatAssemblyEnd(this->mtrx, MAT_FINAL_ASSEMBLY);
}


void
PetscSparseMtrx :: zero()
{
    // test if receiver is already assembled
    PetscTruth assembled;
    MatAssembled(this->mtrx, & assembled);
    if ( assembled ) {
        MatZeroEntries(this->mtrx);
    }
    this->newValues = true;
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
} // end namespace oofem
#endif //ifdef __PETSC_MODULE
