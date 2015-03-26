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

#include "spoolessolver.h"
#include "spoolessparsemtrx.h"
#include "floatarray.h"
#include "verbose.h"
#include "timer.h"
#include "classfactory.h"

namespace oofem {
REGISTER_SparseLinSolver(SpoolesSolver, ST_Spooles);

SpoolesSolver :: SpoolesSolver(Domain *d, EngngModel *m) : SparseLinearSystemNM(d, m)
{
    Lhs = NULL;
    msglvl = 0;
    msgFile = NULL;
    msgFileCloseFlag = 0;

    frontmtx = NULL;
    oldToNewIV = newToOldIV = NULL;
    frontETree = NULL;
    adjIVL = symbfacIVL = NULL;
    mtxmanager = NULL;
    graph = NULL;
}


SpoolesSolver :: ~SpoolesSolver()
{
    if ( msgFileCloseFlag ) {
        fclose(msgFile);
    }

    if ( frontmtx ) {
        FrontMtx_free(frontmtx);
    }

    if ( newToOldIV ) {
        IV_free(newToOldIV);
    }

    if ( oldToNewIV ) {
        IV_free(oldToNewIV);
    }

    if ( frontETree ) {
        ETree_free(frontETree);
    }

    if ( symbfacIVL ) {
        IVL_free(symbfacIVL);
    }

    if ( mtxmanager ) {
        SubMtxManager_free(mtxmanager);
    }

    if ( graph ) {
        Graph_free(graph);
    }
}

IRResultType
SpoolesSolver :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    int val;
    std :: string msgFileName;

    val = -3;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_SpoolesSolver_msglvl);
    msglvl = val;
    IR_GIVE_OPTIONAL_FIELD(ir, msgFileName, _IFT_SpoolesSolver_msgfile);
    if ( !msgFileName.empty() ) {
        msgFile = fopen(msgFileName.c_str(), "w");
        msgFileCloseFlag = 1;
    } else {
        msgFile = stdout;
        msgFileCloseFlag = 0;
    }

    return IRRT_OK;
}


NM_Status
SpoolesSolver :: solve(SparseMtrx &A, FloatArray &b, FloatArray &x)
{
    int errorValue, mtxType, symmetryflag;
    int seed = 30145, pivotingflag = 0;
    int *oldToNew, *newToOld;
    double droptol = 0.0, tau = 1.e300;
    double cpus [ 10 ];
    int stats [ 20 ];

    ChvManager *chvmanager;
    Chv *rootchv;
    InpMtx *mtxA;
    DenseMtx *mtxY, *mtxX;

    x.resize(b.giveSize();

    Timer timer;
    timer.startTimer();

    SpoolesSparseMtrx *As = dynamic_cast< SpoolesSparseMtrx * >(&A);
    if ( !As ) {
        OOFEM_ERROR("SpoolesSparseMtrx Expected");
    }

    mtxA = As->giveInpMtrx();
    mtxType = As->giveValueType();
    symmetryflag = As->giveSymmetryFlag();

    int i;
    int neqns = A->giveNumberOfRows();
    int nrhs = 1;
    /* convert right-hand side to DenseMtx */
    mtxY = DenseMtx_new();
    DenseMtx_init(mtxY, mtxType, 0, 0, neqns, nrhs, 1, neqns);
    DenseMtx_zero(mtxY);
    for ( i = 0; i < neqns; i++ ) {
        DenseMtx_setRealEntry( mtxY, i, 0, b->at(i + 1) );
    }

    if ( ( Lhs != &A ) || ( this->lhsVersion != A.giveVersion() ) ) {
        //
        // lhs has been changed -> new factorization
        //
        ///@todo These factorizations should be kept in the matrix itself, rather than the solver.
        Lhs = &A;
        this->lhsVersion = A.giveVersion();

        if ( frontmtx ) {
            FrontMtx_free(frontmtx);
        }

        if ( newToOldIV ) {
            IV_free(newToOldIV);
        }

        if ( oldToNewIV ) {
            IV_free(oldToNewIV);
        }

        if ( frontETree ) {
            ETree_free(frontETree);
        }

        if ( symbfacIVL ) {
            IVL_free(symbfacIVL);
        }

        if ( mtxmanager ) {
            SubMtxManager_free(mtxmanager);
        }

        if ( graph ) {
            Graph_free(graph);
        }

        /*
         * -------------------------------------------------
         * STEP 3 : find a low-fill ordering
         * (1) create the Graph object
         * (2) order the graph using multiple minimum degree
         * -------------------------------------------------
         */
        int nedges;
        graph = Graph_new();
        adjIVL = InpMtx_fullAdjacency(mtxA);
        nedges = IVL_tsize(adjIVL);
        Graph_init2(graph, 0, neqns, 0, nedges, neqns, nedges, adjIVL,
                    NULL, NULL);
        if ( msglvl > 2 ) {
            fprintf(msgFile, "\n\n graph of the input matrix");
            Graph_writeForHumanEye(graph, msgFile);
            fflush(msgFile);
        }

        frontETree = orderViaMMD(graph, seed, msglvl, msgFile);
        if ( msglvl > 0 ) {
            fprintf(msgFile, "\n\n front tree from ordering");
            ETree_writeForHumanEye(frontETree, msgFile);
            fflush(msgFile);
        }

        /*
         * ----------------------------------------------------
         * STEP 4: get the permutation, permute the front tree,
         * permute the matrix and right hand side, and
         * get the symbolic factorization
         * ----------------------------------------------------
         */
        oldToNewIV = ETree_oldToNewVtxPerm(frontETree);
        oldToNew   = IV_entries(oldToNewIV);
        newToOldIV = ETree_newToOldVtxPerm(frontETree);
        newToOld   = IV_entries(newToOldIV);
        ETree_permuteVertices(frontETree, oldToNewIV);
        InpMtx_permute(mtxA, oldToNew, oldToNew);
        if (  symmetryflag == SPOOLES_SYMMETRIC ||
              symmetryflag == SPOOLES_HERMITIAN ) {
            InpMtx_mapToUpperTriangle(mtxA);
        }

        InpMtx_changeCoordType(mtxA, INPMTX_BY_CHEVRONS);
        InpMtx_changeStorageMode(mtxA, INPMTX_BY_VECTORS);
        symbfacIVL = SymbFac_initFromInpMtx(frontETree, mtxA);
        if ( msglvl > 2 ) {
            fprintf(msgFile, "\n\n old-to-new permutation vector");
            IV_writeForHumanEye(oldToNewIV, msgFile);
            fprintf(msgFile, "\n\n new-to-old permutation vector");
            IV_writeForHumanEye(newToOldIV, msgFile);
            fprintf(msgFile, "\n\n front tree after permutation");
            ETree_writeForHumanEye(frontETree, msgFile);
            fprintf(msgFile, "\n\n input matrix after permutation");
            InpMtx_writeForHumanEye(mtxA, msgFile);
            fprintf(msgFile, "\n\n symbolic factorization");
            IVL_writeForHumanEye(symbfacIVL, msgFile);
            fflush(msgFile);
        }

        Tree_writeToFile(frontETree->tree, ( char * ) "haggar.treef");
        /*--------------------------------------------------------------------*/
        /*
         * ------------------------------------------
         * STEP 5: initialize the front matrix object
         * ------------------------------------------
         */
        frontmtx   = FrontMtx_new();
        mtxmanager = SubMtxManager_new();
        SubMtxManager_init(mtxmanager, NO_LOCK, 0);
        FrontMtx_init(frontmtx, frontETree, symbfacIVL, mtxType, symmetryflag,
                      FRONTMTX_DENSE_FRONTS, pivotingflag, NO_LOCK, 0, NULL,
                      mtxmanager, msglvl, msgFile);
        /*--------------------------------------------------------------------*/
        /*
         * -----------------------------------------
         * STEP 6: compute the numeric factorization
         * -----------------------------------------
         */
        chvmanager = ChvManager_new();
        ChvManager_init(chvmanager, NO_LOCK, 1);
        DVfill(10, cpus, 0.0);
        IVfill(20, stats, 0);
        rootchv = FrontMtx_factorInpMtx(frontmtx, mtxA, tau, droptol,
                                        chvmanager, & errorValue, cpus, stats, msglvl, msgFile);
        ChvManager_free(chvmanager);
        if ( msglvl > 0 ) {
            fprintf(msgFile, "\n\n factor matrix");
            FrontMtx_writeForHumanEye(frontmtx, msgFile);
            fflush(msgFile);
        }

        if ( rootchv != NULL ) {
            fprintf(msgFile, "\n\n matrix found to be singular\n");
            exit(-1);
        }

        if ( errorValue >= 0 ) {
            fprintf(msgFile, "\n\n error encountered at front %d", errorValue);
            exit(-1);
        }

        /*--------------------------------------------------------------------*/
        /*
         * --------------------------------------
         * STEP 7: post-process the factorization
         * --------------------------------------
         */
        FrontMtx_postProcess(frontmtx, msglvl, msgFile);
        if ( msglvl > 2 ) {
            fprintf(msgFile, "\n\n factor matrix after post-processing");
            FrontMtx_writeForHumanEye(frontmtx, msgFile);
            fflush(msgFile);
        }

        /*--------------------------------------------------------------------*/
    }

    /*
     * ----------------------------------------------------
     * STEP 4: permute the right hand side
     * ----------------------------------------------------
     */
    DenseMtx_permuteRows(mtxY, oldToNewIV);
    if ( msglvl > 2 ) {
        fprintf(msgFile, "\n\n right hand side matrix after permutation");
        DenseMtx_writeForHumanEye(mtxY, msgFile);
    }

    /*
     * -------------------------------
     * STEP 8: solve the linear system
     * -------------------------------
     */
    mtxX = DenseMtx_new();
    DenseMtx_init(mtxX, mtxType, 0, 0, neqns, nrhs, 1, neqns);
    DenseMtx_zero(mtxX);
    FrontMtx_solve(frontmtx, mtxX, mtxY, mtxmanager,
                   cpus, msglvl, msgFile);
    if ( msglvl > 2 ) {
        fprintf(msgFile, "\n\n solution matrix in new ordering");
        DenseMtx_writeForHumanEye(mtxX, msgFile);
        fflush(msgFile);
    }

    /*--------------------------------------------------------------------*/
    /*
     * -------------------------------------------------------
     * STEP 9: permute the solution into the original ordering
     * -------------------------------------------------------
     */
    DenseMtx_permuteRows(mtxX, newToOldIV);
    if ( msglvl > 0 ) {
        fprintf(msgFile, "\n\n solution matrix in original ordering");
        DenseMtx_writeForHumanEye(mtxX, msgFile);
        fflush(msgFile);
    }

    // DenseMtx_writeForMatlab(mtxX, "x", msgFile) ;
    /*--------------------------------------------------------------------*/
    /* fetch data to oofem vectors */
    double *xptr = x.givePointer();
    for ( i = 0; i < neqns; i++ ) {
        DenseMtx_realEntry(mtxX, i, 0, xptr + i);
        // printf ("x(%d) = %e\n", i+1, *(xptr+i));
    }

    // DenseMtx_copyRowIntoVector(mtxX, 0, x->givePointer());

    timer.stopTimer();
    OOFEM_LOG_DEBUG( "SpoolesSolver info: user time consumed by solution: %.2fs\n", timer.getUtime() );

    /*
     * -----------
     * free memory
     * -----------
     */
    DenseMtx_free(mtxX);
    DenseMtx_free(mtxY);
    /*--------------------------------------------------------------------*/
    return ( 1 );
}
} // end namespace oofem
