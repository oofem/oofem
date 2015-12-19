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

#include "petscsparsemtrx.h"
#include "engngm.h"
#include "activebc.h"
#include "element.h"
#include "dofmanager.h"
#include "sparsemtrxtype.h"
#include "classfactory.h"

#include <vector>
#include <petscksp.h>
#include <petscvec.h>
#include <petscmat.h>

namespace oofem {
REGISTER_SparseMtrx(PetscSparseMtrx, SMT_PetscMtrx);


PetscSparseMtrx :: PetscSparseMtrx(int n, int m) : SparseMtrx(n, m),
    mtrx(NULL), symmFlag(false), leqs(0), geqs(0), blocksize(1), di(0), kspInit(false), newValues(true), localIS(NULL), globalIS(NULL) { }


PetscSparseMtrx :: PetscSparseMtrx() : SparseMtrx(),
    mtrx(NULL), symmFlag(false), leqs(0), geqs(0), blocksize(1), di(0), kspInit(false), newValues(true), localIS(NULL), globalIS(NULL){ }


PetscSparseMtrx :: ~PetscSparseMtrx()
{
    MatDestroy(& this->mtrx);
    if ( this->kspInit ) {
        KSPDestroy(& this->ksp);
    }
    if ( localIS ) {
        ISDestroy(& localIS);
        ISDestroy(& globalIS);
    }
}


SparseMtrx *
PetscSparseMtrx :: GiveCopy() const
{
    PetscSparseMtrx *answer = new PetscSparseMtrx(nRows, nColumns);
    MatDuplicate( this->mtrx, MAT_COPY_VALUES, & ( answer->mtrx ) );
    answer->symmFlag = this->symmFlag;
    answer->mType    = this->mType;
    answer->leqs     = this->leqs;
    answer->geqs     = this->geqs;
    answer->di       = this->di;
    answer->emodel   = this->emodel;
    answer->kspInit  = false;
    answer->newValues = this->newValues;

    return answer;
}

void
PetscSparseMtrx :: times(const FloatArray &x, FloatArray &answer) const
{
    if ( emodel->isParallel() ) {
        answer.resize(x.giveSize());

        Vec globX;
        Vec globY;

        /*
         * scatter and gather x to global representation
         */
        this->createVecGlobal(& globX);
        this->scatterL2G(x, globX);

        VecDuplicate(globX, & globY);

        MatMult(this->mtrx, globX, globY);

        this->scatterG2L(globY, answer);

        VecDestroy(& globX);
        VecDestroy(& globY);
    } else {
        if ( this->giveNumberOfColumns() != x.giveSize() ) {
            OOFEM_ERROR("Dimension mismatch");
        }

        Vec globX, globY;
        VecCreateSeqWithArray(PETSC_COMM_SELF, 1, x.giveSize(), x.givePointer(), & globX);
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
        VecDestroy(& globX);
        VecDestroy(& globY);
    }
}

void
PetscSparseMtrx :: timesT(const FloatArray &x, FloatArray &answer) const
{
    if ( this->giveNumberOfRows() != x.giveSize() ) {
        OOFEM_ERROR("Dimension mismatch");
    }

    if ( emodel->isParallel() ) {
        OOFEM_ERROR("Not implemented");
    }

    Vec globX, globY;
    VecCreateSeqWithArray(PETSC_COMM_SELF, 1, x.giveSize(), x.givePointer(), & globX);
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
    VecDestroy(& globX);
    VecDestroy(& globY);
}


void
PetscSparseMtrx :: times(const FloatMatrix &B, FloatMatrix &answer) const
{
    if ( this->giveNumberOfColumns() != B.giveNumberOfRows() ) {
        OOFEM_ERROR("Dimension mismatch");
    }

    if ( emodel->isParallel() ) {
        OOFEM_ERROR("Not implemented");
    }

    // I'm opting to work with a set of vectors, as i think it might be faster and more robust. / Mikael

    int nr = this->giveNumberOfRows();
    int nc = B.giveNumberOfColumns();
    answer.resize(nr, nc);

#if 1
    // Using several vectors are necessary because PETSc is missing most of it operations on symmetric matrices.
    // If PETSc ever adds the missing "MatMatMult_***_***" implementations, we can use the second approach below
    Vec globX, globY;
    VecCreate(PETSC_COMM_SELF, & globY);
    VecSetType(globY, VECSEQ);
    VecSetSizes(globY, PETSC_DECIDE, nr);
    int nrB = B.giveNumberOfRows();
    FloatArray colVals(nrB);
    double *resultsPtr;
    for ( int k = 0; k < nc; ++k ) {
        B.copyColumn(colVals, k + 1);
        VecCreateSeqWithArray(PETSC_COMM_SELF, 1, nrB, colVals.givePointer(), & globX);
        MatMult(this->mtrx, globX, globY);
        VecGetArray(globY, & resultsPtr);
        for ( int i = 0; i < nr; ++i ) {
            answer(i, k) = resultsPtr [ i ];
        }
        VecRestoreArray(globY, & resultsPtr);
        VecDestroy(& globX);
    }
    VecDestroy(& globY);
#else
    double *aptr = answer.givePointer();
    Mat globB, globC;
    MatCreateSeqDense(PETSC_COMM_SELF, B.giveNumberOfRows(), B.giveNumberOfColumns(), B.givePointer(), & globB);
    MatMatMult(this->mtrx, globB, MAT_INITIAL_MATRIX, PETSC_DEFAULT, & globC);
    const double *vals;
    for ( int r = 0; r < nr; r++ ) {
        MatGetRow(globC, r, NULL, NULL, & vals);
        for ( int i = 0, i2 = r; i < nc; i++, i2 += nr ) {
            aptr [ i2 ] = vals [ i ];
        }
        MatRestoreRow(globC, r, NULL, NULL, & vals);
    }

    MatDestroy(& globB);
    MatDestroy(& globC);
#endif
}

void
PetscSparseMtrx :: timesT(const FloatMatrix &B, FloatMatrix &answer) const
{
    if ( this->giveNumberOfRows() != B.giveNumberOfRows() ) {
        OOFEM_ERROR("Dimension mismatch");
    }

    if ( emodel->isParallel() ) {
        OOFEM_ERROR("Not implemented");
    }

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
        MatRestoreRow(globC, r, NULL, NULL, & vals);
    }

    MatDestroy(& globB);
    MatDestroy(& globC);
}

void
PetscSparseMtrx :: times(double x)
{
    MatScale(this->mtrx, x);
}


void
PetscSparseMtrx :: add(double x, SparseMtrx &m)
{
    PetscSparseMtrx *M = dynamic_cast< PetscSparseMtrx * >( &m );
    MatAXPY(*this->giveMtrx(), x, *M->giveMtrx(), SAME_NONZERO_PATTERN);
}


void
PetscSparseMtrx :: addDiagonal(double x, FloatArray &m)
{
    OOFEM_WARNING("Calling function that has not been tested (remove this message when it is verified)");
    Vec globM;
    if ( emodel->isParallel() ) {
        this->createVecGlobal(& globM);
        this->scatterL2G(m, globM);
    } else {
        VecCreateSeqWithArray(PETSC_COMM_SELF, 1, m.giveSize(), m.givePointer(), & globM);
    }
    MatDiagonalSet(this->mtrx, globM, ADD_VALUES);
    VecDestroy(& globM);
}

///@todo I haven't looked at the parallel code yet (lack of time right now, and i want to see it work first). / Mikael
int
PetscSparseMtrx :: buildInternalStructure(EngngModel *eModel, int di, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    IntArray loc;
    Domain *domain = eModel->giveDomain(di);
    blocksize = 1;

    if ( mtrx ) {
        MatDestroy(& mtrx);
    }

    if ( this->kspInit ) {
        KSPDestroy(& ksp);
        this->kspInit  = false; // force ksp to be initialized
    }

    this->emodel = eModel;
    this->di = di;

    if ( emodel->isParallel() ) {
        OOFEM_ERROR("Not implemented");
    }

    nRows = eModel->giveNumberOfDomainEquations(di, r_s);
    nColumns = eModel->giveNumberOfDomainEquations(di, c_s);

    geqs = leqs = nRows;

    int total_nnz;
    IntArray d_nnz(leqs), d_nnz_sym(leqs);

    {
        //determine nonzero structure of matrix
        IntArray r_loc, c_loc;
        std :: vector< IntArray > rows_upper(nRows), rows_lower(nRows);

        for ( auto &elem : domain->giveElements() ) {
            elem->giveLocationArray(r_loc, r_s);
            elem->giveLocationArray(c_loc, c_s);
            for ( int ii : r_loc ) {
                if ( ii > 0 ) {
                    for ( int jj : c_loc ) {
                        if ( jj > 0 ) {
                            if ( jj >= ii ) {
                                rows_upper [ ii - 1 ].insertSortedOnce(jj - 1, c_loc.giveSize() / 2);
                            } else {
                                rows_lower [ ii - 1 ].insertSortedOnce(jj - 1, c_loc.giveSize() / 2);
                            }
                        }
                    }
                }
            }
        }
        // Structure from active boundary conditions.
        std :: vector< IntArray >r_locs, c_locs;
        for ( auto &gbc : domain->giveBcs() ) {
            ActiveBoundaryCondition *activebc = dynamic_cast< ActiveBoundaryCondition * >( gbc.get() );
            if ( activebc ) {
                activebc->giveLocationArrays(r_locs, c_locs, TangentStiffnessMatrix, r_s, c_s);
                for ( std :: size_t k = 0; k < r_locs.size(); k++ ) {
                    IntArray &krloc = r_locs [ k ];
                    IntArray &kcloc = c_locs [ k ];
                    for ( int ii : krloc ) {
                        if ( ii > 0 ) {
                            for ( int jj : kcloc ) {
                                if ( jj > 0 ) {
                                    if ( jj >= ii ) {
                                        rows_upper [ ii - 1 ].insertSortedOnce(jj - 1, kcloc.giveSize() / 2);
                                    } else {
                                        rows_lower [ ii - 1 ].insertSortedOnce(jj - 1, kcloc.giveSize() / 2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        total_nnz = 0;
        for ( int i = 0; i < leqs; i++ ) {
            d_nnz(i) = rows_upper [ i ].giveSize() + rows_lower [ i ].giveSize();
        }
    }

    // create PETSc mat
    MatCreate(PETSC_COMM_SELF, & mtrx);
    MatSetSizes(mtrx, nRows, nColumns, nRows, nColumns);
    MatSetFromOptions(mtrx);

    if ( total_nnz > nRows * nColumns / 10 ) { // More than 10% nnz, then we just force the dense matrix.
        MatSetType(mtrx, MATDENSE);
    } else {
        MatSetType(mtrx, MATSEQAIJ);
    }

    //The incompatible preallocations are ignored automatically.
    MatSetUp(mtrx);
    MatSeqAIJSetPreallocation( mtrx, 0, d_nnz.givePointer() );

    MatSetOption(mtrx, MAT_ROW_ORIENTED, PETSC_FALSE); // To allow the insertion of values using MatSetValues in column major order
    MatSetOption(mtrx, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

    this->newValues = true;
    return true;
}

int
PetscSparseMtrx :: buildInternalStructure(EngngModel *eModel, int di, const UnknownNumberingScheme &s)
{
    IntArray loc;
    Domain *domain = eModel->giveDomain(di);
    blocksize = 1;

    // Delete old stuff first;
    if ( mtrx ) {
        MatDestroy(& mtrx);
    }

    if ( this->kspInit ) {
        KSPDestroy(& ksp);
        this->kspInit  = false; // force ksp to be initialized
    }

    if ( localIS ) {
        ISDestroy(& localIS);
        ISDestroy(& globalIS);
        localIS = NULL;
        globalIS = NULL;
    }

    this->emodel = eModel;
    this->di = di;


#ifdef __PARALLEL_MODE
    if ( eModel->isParallel() ) {
        Natural2GlobalOrdering *n2g;
        Natural2LocalOrdering *n2l;
        ParallelContext *context = eModel->giveParallelContext(di);
        n2g = context->giveN2Gmap();
        n2l = context->giveN2Lmap();

        n2l->init(eModel, di, s);
        n2g->init(eModel, di, s);

 #ifdef __VERBOSE_PARALLEL
        VERBOSEPARALLEL_PRINT( "PetscSparseMtrx:: buildInternalStructure", "", eModel->giveRank() );
 #endif

        leqs = n2g->giveNumberOfLocalEqs();
        geqs = n2g->giveNumberOfGlobalEqs();

        //printf("%d, leqs = %d, geqs = %d\n", this->emodel->giveRank(), leqs, geqs);

 #ifdef __VERBOSE_PARALLEL
        OOFEM_LOG_INFO( "[%d]PetscSparseMtrx:: buildInternalStructure: l_eqs = %d, g_eqs = %d, n_eqs = %d\n", eModel->giveRank(), leqs, geqs, eModel->giveNumberOfDomainEquations(di, s) );
 #endif

        IntArray d_nnz(leqs), o_nnz(leqs), d_nnz_sym(leqs), o_nnz_sym(leqs);

        {
            // determine nonzero structure of a "local (maintained)" part of matrix, and the off-diagonal part
            // allocation map
            ///@todo Split this into upper and lower triangular part like for the sequential matrix (we can then use SBAIJ which is a huge performance boost)
            std :: vector< IntArray > d_rows_upper(leqs), d_rows_lower(leqs);   // diagonal sub-matrix allocation
            std :: vector< IntArray > o_rows_upper(leqs), o_rows_lower(leqs);   // off-diagonal allocation

            IntArray lloc, gloc;

            //fprintf (stderr,"[%d] n2l map: ",rank);
            //for (n=1; n<=n2l.giveN2Lmap()->giveSize(); n++) fprintf (stderr, "%d ", n2l.giveN2Lmap()->at(n));

            for ( auto &elem : domain->giveElements() ) {
                //fprintf (stderr, "(elem %d) ", n);
                elem->giveLocationArray(loc, s);
                n2l->map2New(lloc, loc, 0); // translate natural->local numbering (remark, 1-based indexing)
                n2g->map2New(gloc, loc, 0); // translate natural->global numbering (remark, 0-based indexing)
                // See the petsc manual for details on how this allocation is constructed.
                int ii, jj;
                for ( int i = 1; i <= lloc.giveSize(); i++ ) {
                    if ( ( ii = lloc.at(i) ) ) {
                        for ( int j = 1; j <= lloc.giveSize(); j++ ) {
                            if ( ( jj = gloc.at(j) ) >= 0 ) { // if negative, then it is prescribed
                                if ( lloc.at(j) ) { // if true, then its the local part (the diagonal block matrix)
                                    if ( jj >= ( ii - 1 ) ) { // NOTE: ii (the rows) is in 1-based indexing, jj is in 0-base (!)
                                        d_rows_upper [ ii - 1 ].insertSortedOnce(jj, loc.giveSize() / 2);
                                    } else {
                                        d_rows_lower [ ii - 1 ].insertSortedOnce(jj, loc.giveSize() / 2);
                                    }
                                } else { // Otherwise it must be off-diagonal
                                    if ( jj >= ( ii - 1 ) ) {
                                        o_rows_upper [ ii - 1 ].insertSortedOnce(jj, loc.giveSize() / 2);
                                    } else {
                                        o_rows_lower [ ii - 1 ].insertSortedOnce(jj, loc.giveSize() / 2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            ///@todo Take into account active boundary conditions here.

            // Diagonal must always be allocated; this code ensures that for every local line, it adds the global column number
            IntArray *n2gmap = n2g->giveN2Gmap();
            IntArray *n2lmap = n2l->giveN2Lmap();
            for ( int n = 1; n <= n2lmap->giveSize(); ++n ) {
                if ( n2lmap->at(n) ) {
                    d_rows_upper [ n2lmap->at(n) - 1 ].insertSortedOnce( n2gmap->at(n) );
                }
            }

            for ( int i = 0; i < leqs; i++ ) {
                d_nnz(i) = d_rows_upper [ i ].giveSize() + d_rows_lower [ i ].giveSize();
                o_nnz(i) = o_rows_upper [ i ].giveSize() + o_rows_lower [ i ].giveSize();

                d_nnz_sym(i) = d_rows_upper [ i ].giveSize();
                o_nnz_sym(i) = o_rows_upper [ i ].giveSize();
            }
        }

        //fprintf (stderr,"\n[%d]PetscSparseMtrx: Profile ...",rank);
        //for (i=0; i<leqs; i++) fprintf(stderr, "%d ", d_nnz(i));
        //fprintf (stderr,"\n[%d]PetscSparseMtrx: Creating MPIAIJ Matrix ...\n",rank);

        // create PETSc mat
        MatCreate(this->emodel->giveParallelComm(), & mtrx);
        MatSetSizes(mtrx, leqs, leqs, geqs, geqs);
        MatSetType(mtrx, MATMPIAIJ);
        MatSetFromOptions(mtrx);
        MatMPIAIJSetPreallocation( mtrx, 0, d_nnz.givePointer(), 0, o_nnz.givePointer() );
        //MatMPIBAIJSetPreallocation( mtrx, 1, 0, d_nnz.givePointer(), 0, o_nnz.givePointer() );
        //MatMPISBAIJSetPreallocation( mtrx, 1, 0, d_nnz_sym.givePointer(), 0, o_nnz_sym.givePointer() );
        //MatXAIJSetPreallocation( mtrx, 1, d_nnz.givePointer(), o_nnz.givePointer(), d_nnz_sym.givePointer(), o_nnz_sym.givePointer());

        // Creates scatter context for PETSc.
        ISCreateGeneral(this->emodel->giveParallelComm(), context->giveNumberOfNaturalEqs(), n2g->giveN2Gmap()->givePointer(), PETSC_USE_POINTER, & globalIS);
        ISCreateStride(this->emodel->giveParallelComm(), context->giveNumberOfNaturalEqs(), 0, 1, & localIS);

 #ifdef __VERBOSE_PARALLEL
        VERBOSEPARALLEL_PRINT("PetscSparseMtrx:: buildInternalStructure", "done", eModel->giveRank());
 #endif
    } else
#endif
    {
        leqs = geqs = eModel->giveNumberOfDomainEquations(di, s);
        IntArray d_nnz, d_nnz_sym, indices, rowstart;

        MatType type = MATSEQBAIJ;

        MatCreate(PETSC_COMM_SELF, & mtrx);
        MatSetSizes(mtrx, leqs, leqs, geqs, geqs);
        MatSetType(mtrx, type);
        MatSetFromOptions(mtrx);
        MatGetType( mtrx, &type );

        if ( strcmp(type, MATDENSE) != 0 ) {
            // Dry assembly (computes the nonzero structure) (and the blocks)
            std :: vector< IntArray > rows_upper(leqs), blocks;

            for ( auto &elem : domain->giveElements() ) {
                elem->giveLocationArray(loc, s);
                for ( int ii : loc ) {
                    if ( ii > 0 ) {
                        for ( int jj : loc ) {
                            if ( jj >= ii ) {
                                rows_upper [ ii - 1 ].insertSortedOnce( jj - 1, loc.giveSize() );
                            }
                        }
                    }
                }
            }

            // Structure from active boundary conditions.
            std :: vector< IntArray > r_locs, c_locs;
            for ( auto &gbc : domain->giveBcs() ) {
                ActiveBoundaryCondition *activebc = dynamic_cast< ActiveBoundaryCondition * >( gbc.get() );
                if ( activebc ) {
                    activebc->giveLocationArrays(r_locs, c_locs, TangentStiffnessMatrix, s, s);
                    for ( std :: size_t k = 0; k < r_locs.size(); k++ ) {
                        for ( int ii : r_locs [ k ] ) {
                            if ( ii > 0 ) {
                                for ( int jj : c_locs [ k ] ) {
                                    if ( jj >= ii ) {
                                        rows_upper [ ii - 1 ].insertSortedOnce( jj - 1, loc.giveSize() );
                                    }
                                }
                            }
                        }
                    }
                }
            }


            int nnz = 0;
            for ( auto &row_upper : rows_upper ) {
                nnz += row_upper.giveSize();
            }

            if ( strcmp(type, MATSEQAIJ) != 0 ) {
                // Compute optimal block size (test various sizes and compute the "wasted" zeros).
                int maxblock = domain->giveNumberOfDofManagers() > 0 ? domain->giveDofManager(1)->giveNumberOfDofs() : 0;
                for ( int bsize = maxblock; bsize > 1; --bsize ) {
                    int nblocks = ceil(rows_upper.size() / (double)bsize);
                    blocks.clear();
                    blocks.resize(nblocks);
                    for ( int i = 0; i < leqs; i++ ) {
                        int blockrow = i / bsize;
                        for ( int j : rows_upper [ i ] ) {
                            int blockcol = j / bsize;
                            blocks [ blockrow ].insertSortedOnce( blockcol );
                        }
                    }
                    // See how much space (and operations) we risk wasting with this block size:
                    int bnnz = 0;
                    for ( auto &block : blocks ) {
                        bnnz += block.giveSize() * bsize * bsize;
                    }
                    OOFEM_LOG_DEBUG("Block size %d resulted in nnz = %d -> %d using %d^2 blocks\n", bsize, nnz, bnnz, nblocks);
                    if ( bnnz <=  nnz * 1.2 ) {
                        blocksize = bsize;
                        break;
                    }
                }
            }

            if ( blocksize == 1 ) {
                blocks = rows_upper;
            }

            int nblocks = ceil(rows_upper.size() / (double)blocksize);
            d_nnz_sym.resize(nblocks);
            d_nnz.resize(nblocks);
            for ( int i = 0; i < nblocks; i++ ) {
                d_nnz_sym[i] = blocks [ i ].giveSize();
                // We can optimize to use only upper half by using the fact that the problem has symmetric nonzero-structure.
                d_nnz[i] += d_nnz_sym[i];
                for ( int jj : blocks [ i ] ) {
                    if (  jj != i ) {
                        d_nnz[ jj ]++;
                    }
                }
            }
        }

        // Based on the user selected type, determine the block size;
        if ( strcmp(type, MATSEQAIJ) == 0 ) {
            MatSeqAIJSetPreallocation( mtrx, 0, d_nnz.givePointer() );
        } else if ( strcmp(type, MATSEQBAIJ) == 0 ) {
            MatSeqBAIJSetPreallocation( mtrx, blocksize, 0, d_nnz.givePointer() );
        } else if ( strcmp(type, MATSEQSBAIJ) == 0 )  {
            MatSeqSBAIJSetPreallocation( mtrx, blocksize, 0, d_nnz_sym.givePointer() );
        }
    }

    MatSetOption(mtrx, MAT_ROW_ORIENTED, PETSC_FALSE); // To allow the insertion of values using MatSetValues in column major order
    MatSetOption(mtrx, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    MatSetOption(mtrx, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE);

    nRows = nColumns = geqs;
    this->newValues = true;
    return true;
}

int
PetscSparseMtrx :: assemble(const IntArray &loc, const FloatMatrix &mat)
{
    int ndofe = mat.giveNumberOfRows();
    IntArray gloc(ndofe);

#ifdef __PARALLEL_MODE
    if ( emodel->isParallel() ) {
        // translate local code numbers to global ones
        emodel->giveParallelContext(this->di)->giveN2Gmap()->map2New(gloc, loc, 0);

        //fprintf (stderr, "[?] gloc=");
        //for (int i=1; i<=ndofe; i++) fprintf (stderr, "%d ", gloc.at(i));
    } else
#endif
    {
        for ( int i = 1; i <= ndofe; i++ ) {
            gloc.at(i) = loc.at(i) - 1;
        }
    }

    this->version++;
    this->newValues = true;

    // Block'ify the matrix in order to use MatSetValuesBlocked (doing this work really is faster)
    if ( this->blocksize > 1 && ndofe % this->blocksize == 0 ) {
        int nblocke = ndofe / this->blocksize;
        IntArray bloc(nblocke);
        for ( int b = 0; b < nblocke; ++b ) {
            int i = b * blocksize;
            // Have to check alignment inside the block:
            if ( gloc[i] != -1 && gloc[i] % this->blocksize != 0 ) {
                MatSetValues(this->mtrx, ndofe, gloc.givePointer(), ndofe, gloc.givePointer(), mat.givePointer(), ADD_VALUES);
                return 1;
            }
            // Have to verify that this is indeed a series of increments in blocks:
            for ( int k = 1; k < this->blocksize; ++k ) {
                bool ok_seq = gloc[i] + k == gloc[i+k];
                bool ok_bc = gloc[i] == -1 && gloc[i+k] == -1;
                if ( !(ok_bc || ok_seq) ) {
                    // Then we can't use block assembling for this element.
                    MatSetValues(this->mtrx, ndofe, gloc.givePointer(), ndofe, gloc.givePointer(), mat.givePointer(), ADD_VALUES);
                    return 1;
                }
            }
            bloc[b] = gloc[i] == -1 ? -1 : gloc[i] / blocksize;
        }
        MatSetValuesBlocked(this->mtrx, nblocke, bloc.givePointer(), nblocke, bloc.givePointer(), mat.givePointer(), ADD_VALUES);
    } else {
        MatSetValues(this->mtrx, ndofe, gloc.givePointer(), ndofe, gloc.givePointer(), mat.givePointer(), ADD_VALUES);
    }

    //MatView(this->mtrx,PETSC_VIEWER_STDOUT_SELF);
    return 1;
}

int
PetscSparseMtrx :: assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat)
{
#ifdef __PARALLEL_MODE
    if ( emodel->isParallel() ) {
        // translate eq numbers
        IntArray grloc( rloc.giveSize() ), gcloc( cloc.giveSize() );
        emodel->giveParallelContext(this->di)->giveN2Gmap()->map2New(grloc, rloc, 0);
        emodel->giveParallelContext(this->di)->giveN2Gmap()->map2New(gcloc, cloc, 0);

        MatSetValues(this->mtrx, grloc.giveSize(), grloc.givePointer(),
                     gcloc.giveSize(), gcloc.givePointer(), mat.givePointer(), ADD_VALUES);
    } else {
#endif
    int rsize = rloc.giveSize(), csize = cloc.giveSize();
    IntArray grloc(rsize), gcloc(csize);
    for ( int i = 1; i <= rsize; i++ ) {
        grloc.at(i) = rloc.at(i) - 1;
    }

    for ( int i = 1; i <= csize; i++ ) {
        gcloc.at(i) = cloc.at(i) - 1;
    }

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
    PetscBool assembled;
    MatAssembled(this->mtrx, & assembled);
    if ( assembled ) {
        MatZeroEntries(this->mtrx);
    }
    this->newValues = true;
}

double
PetscSparseMtrx :: computeNorm() const
{
    double norm;
    MatNorm(this->mtrx, NORM_1, & norm);
    return norm;
}

double &
PetscSparseMtrx :: at(int i, int j)
{
    static double a;
    OOFEM_ERROR("unsupported");
    return a;
}

double
PetscSparseMtrx :: at(int i, int j) const
{
    OOFEM_ERROR("unsupported");
    return 0;
    //double value;
    //int row = i-1, col = j-1;
    //MatGetValues(this->mtrx, 1, &row, 1, &col, &value);
    //return value;
}

void
PetscSparseMtrx :: toFloatMatrix(FloatMatrix &answer) const
{
    OOFEM_ERROR("unsupported");
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

void
PetscSparseMtrx :: printMatlab() const
{
    PetscViewerSetFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_MATLAB);
    MatView(this->mtrx, PETSC_VIEWER_STDOUT_SELF);
}

void
PetscSparseMtrx :: writeToFile(const char *fname) const
{
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, fname, & viewer);
    MatView(this->mtrx, viewer);
    PetscViewerDestroy(& viewer);
}


void
PetscSparseMtrx :: createVecGlobal(Vec *answer) const
{
#ifdef __PARALLEL_MODE
    if ( emodel->isParallel() ) {
        VecCreate(this->emodel->giveParallelComm(), answer);
        VecSetSizes(* answer, this->leqs, this->geqs);
        VecSetFromOptions(* answer);
    } else {
#endif
    VecCreateSeq(PETSC_COMM_SELF, this->giveNumberOfRows(), answer);
#ifdef __PARALLEL_MODE
}
#endif
}


int
PetscSparseMtrx :: scatterG2L(Vec src, FloatArray &dest) const
{
    PetscScalar *ptr;

#ifdef __PARALLEL_MODE
    if ( emodel->isParallel() ) {
        ParallelContext *context = emodel->giveParallelContext(di);
        int neqs = context->giveNumberOfNaturalEqs();
        Vec locVec;
        VecCreateSeq(PETSC_COMM_SELF, neqs, & locVec);

        VecScatter n2gvecscat;
        VecScatterCreate(locVec, localIS, src, globalIS, & n2gvecscat);
        VecScatterBegin(n2gvecscat, src, locVec, INSERT_VALUES, SCATTER_REVERSE);
        VecScatterEnd(n2gvecscat, src, locVec, INSERT_VALUES, SCATTER_REVERSE);
        VecScatterDestroy(& n2gvecscat);

        dest.resize(neqs);
        VecGetArray(locVec, & ptr);
        for ( int i = 0; i < neqs; i++ ) {
            dest.at(i + 1) = ptr [ i ];
        }

        VecRestoreArray(locVec, & ptr);
        VecDestroy(& locVec);
    } else {
#endif
    int neqs = this->giveNumberOfRows();
    dest.resize(neqs);
    VecGetArray(src, & ptr);
    for ( int i = 0; i < neqs; i++ ) {
        dest.at(i + 1) = ptr [ i ];
    }
    VecRestoreArray(src, & ptr);
#ifdef __PARALLEL_MODE
}
#endif
    return 1;
}


int
PetscSparseMtrx :: scatterL2G(const FloatArray &src, Vec dest) const
{
    const PetscScalar *ptr = src.givePointer();

#ifdef __PARALLEL_MODE
    if ( emodel->isParallel() ) {
        ParallelContext *context = this->emodel->giveParallelContext(di);
        int size = src.giveSize();

        Natural2LocalOrdering *n2l = context->giveN2Lmap();
        Natural2GlobalOrdering *n2g = context->giveN2Gmap();
        for ( int i = 0; i < size; i++ ) {
            if ( n2l->giveNewEq(i + 1) ) {
                int eqg = n2g->giveNewEq(i + 1);
                VecSetValue(dest, eqg, ptr [ i ], INSERT_VALUES);
            }
        }

        VecAssemblyBegin(dest);
        VecAssemblyEnd(dest);
    } else {
#endif

    int size = src.giveSize();
    for ( int i = 0; i < size; i++ ) {
        //VecSetValues(dest, 1, & i, ptr + i, mode);
        VecSetValue(dest, i, ptr [ i ], INSERT_VALUES);
    }

    VecAssemblyBegin(dest);
    VecAssemblyEnd(dest);
#ifdef __PARALLEL_MODE
}
#endif
    return 1;
}


SparseMtrxType PetscSparseMtrx :: giveType() const
{
    return SMT_PetscMtrx;
}

bool PetscSparseMtrx :: isAsymmetric() const
{
    return !symmFlag;
}

Mat *PetscSparseMtrx :: giveMtrx()
{
    return & this->mtrx;
}

bool PetscSparseMtrx :: giveSymmetryFlag() const
{
    return symmFlag;
}

int PetscSparseMtrx :: setOption(MatOption op, PetscBool flag)
{
    return MatSetOption(this->mtrx, op, flag);
}

int PetscSparseMtrx :: giveLeqs()
{
    return leqs;
}

int PetscSparseMtrx :: giveDomainIndex() const
{
    return di;
}
} // end namespace oofem
