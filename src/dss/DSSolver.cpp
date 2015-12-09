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

#include "DSSolver.h"
#include "SparseGridMtxLU.h"
#include "SparseGridMtxLL.h"
#include "SkyLineMtxLDL.h"
#include "MathTracer.h"

DSS_NAMESPASE_BEGIN

DSSolver :: DSSolver(MathTracer *pMT)
{
    SetMT(pMT);
    matrix = NULL;
    matrixPD = NULL;
    blockSize = 3;
    n_blocks = 0;
    neq = 0;
    decompid = 0;
    m_eState = ISolver :: None;
    tmpR = NULL;
    dom_order = NULL;
    mcn = NULL;
    MatrixType = eDSSparseMatrix;

    orig_matrix = NULL;

    fixed = NULL;
    lncn = NULL;

    SolverType = eDSSFactorizationLDLT;
    run_code = 0;

    OrderingType = Ordering :: ApproxMinimumDegree;

    eMT->Writeln("");
    eMT->Writeln("___FemCAD--DirectSparse_Solver______");
    eMT->Writeln("   Richard Vondracek  (c) 2001-2003 ");
}

DSSolver :: ~DSSolver()
{
    Dispose();
}

long DSSolver :: Initialize(unsigned char run_code, eDSSolverType solverType, eDSMatrixType matrixType)
{
    Dispose();
    eMT->Writeln("");
    eMT->Writeln("Sparse Direct initialized");
    this->SolverType = solverType;
    this->MatrixType = matrixType;
    this->run_code = run_code;

    switch ( MatrixType ) {
    case eDSSparseMatrix:
    {
        switch ( SolverType ) {
        case eDSSFactorizationLLTIncomplete:
            OrderingType = Ordering :: ApproxMinimumDegreeIncomplete;
            //OrderingType = Ordering::ApproxMinimumDegree;
            break;
        case eDSSFactorizationLDLTIncomplete:
            OrderingType = Ordering :: ApproxMinimumDegreeIncomplete;
            //OrderingType = Ordering::ApproxMinimumDegree;
            break;
        case eDSSFactorizationLDLT:
        case eDSSFactorizationLLT:
        case eDSSFactorizationLU:
            OrderingType = Ordering :: ApproxMinimumDegree;
            //OrderingType = Ordering::None;
            break;
        case eDSSFastCG:
            OrderingType = Ordering :: None;
            break;
        default:
            OrderingType = Ordering :: None;
            break;
        }

        break;
    }
    case eDSSkylineMatrix:
    {
        switch ( SolverType ) {
        case eDSSFactorizationLLTIncomplete:
            OrderingType = Ordering :: ReverseCuthillMcKee;
            break;
        default:
            OrderingType = Ordering :: None;
            break;
        }
    }
    break;
    }

    m_eState = ISolver :: Initialized;
    return 1;
}

void DSSolver :: SetMT(MathTracer *pMT)
{
    if ( pMT ) {
        eMT = pMT;
    } else {
        eMT = & MT;
    }
}

void DSSolver :: Dispose()
{
    if ( matrix ) {
        delete matrix;
        matrix = NULL;
    }

    if ( matrixPD ) {
        delete matrixPD;
        matrixPD = NULL;
    }

    if ( orig_matrix ) {
        delete orig_matrix;
        orig_matrix = NULL;
    }

    if ( tmpR ) {
        delete [] tmpR;
        tmpR = NULL;
    }

    if ( dom_order ) {
        delete dom_order;
        dom_order = NULL;
    }

    if ( mcn ) {
        delete mcn;
        mcn = NULL;
    }

    if ( fixed ) {
        delete fixed;
        fixed = NULL;
    }

    if ( lncn ) {
        delete lncn;
        lncn = NULL;
    }
}

bool DSSolver :: SetOrderingType(Ordering :: Type otype)
{
    switch ( otype ) {
    case Ordering :: None:
    case Ordering :: MinimumDegree:
    case Ordering :: ApproxMinimumDegree:
    case Ordering :: ReverseCuthillMcKee:
    case Ordering :: CuthillMcKee:
    case Ordering :: NestedGraphBisection:
#ifdef _LINK_METIS_
    case Ordering :: MetisND:
#endif
#ifdef _LINK_COLAMD_
    case Ordering :: ColAMD:
#endif
        OrderingType = otype;
        return true;

    default:
        eMT->Writeln("Selected ordering type is not implemented! Default value is used.");
        break;
    }

    return false;
}


bool DSSolver :: LoadMatrix(unsigned long neq, unsigned char block_size, double *a, unsigned long *ci, unsigned long *adr)
{
    SparseMatrixF loc_sm(neq, a, ci, adr);
    return LoadMatrix(& loc_sm, block_size);
}

bool DSSolver :: LoadMatrix(SparseMatrixF *smt, unsigned char block_size)
{
    SetMatrixPattern(smt, block_size);
    sm.CreateLocalCopy();
    m_eState = ISolver :: Allocated;
    return 1;
}

bool DSSolver :: SetMatrixPattern(SparseMatrixF *smt, unsigned char block_size)
{
    Dispose();
    sm = * smt;
    blockSize = block_size;
    neq = ( long ) sm.neq;

    eMT->Writeln("Loading sparse matrix");
    sprintf(str, " number of equations : %ld ", sm.neq);
    eMT->Writeln(str);
    sprintf( str, " number of nonzeros  : %ld nonzeros", sm.Nonzeros() );
    eMT->Writeln(str);

    if ( SolverType == eDSSDiagonalScaling ) {
        matrixD.Init(neq);
        sm.ReadDiagonal( matrixD.DataPtr() );
    }

    return 1;
}

bool DSSolver :: decomp()
{
    return decompid;
}
void DSSolver :: changedecomp()
{
    if ( decompid == 0 ) {
        decompid = 1;
    } else {
        decompid = 0;
    }
}

bool DSSolver :: IsFactorized()
{
    return m_eState == ISolver :: Factorized;
}

bool DSSolver :: Factorize()
{
    return StartSolver();
}

SparseGridMtx *DSSolver :: CreateNewSparseGridMtx(IntArrayList *fixed)
{
    long no_nonzeros = sm.Nonzeros();
    eMT->CS();

    SparseConectivityMtxII *block_conmtx = new SparseConectivityMtxII(sm, this->mcn, blockSize);
    block_conmtx->eMT = eMT;

    long no_init_blocks = ( block_conmtx->Nonzeros() ) / 2 + block_conmtx->N();
    n_blocks = block_conmtx->N();

    Ordering *order = block_conmtx->GetPermutationAndPattern(OrderingType, fixed);
    delete block_conmtx;
    eMT->Writeln( eMT->MC_() );
    block_conmtx = NULL;

    SparseGridMtx *matrix = NULL;

    eMT->Writeln("Allocating block sparse matrix");
    switch ( SolverType ) {
    case eDSSFactorizationLDLTIncomplete:
        matrix = new SparseGridMtxLDL(sm, blockSize, order, mcn, eMT, true);
        break;
    case eDSSFactorizationLLTIncomplete:
        matrix = new SparseGridMtxLL(sm, blockSize, order, mcn, eMT, true);
        break;
    case eDSSFactorizationLDLT:
        matrix = new SparseGridMtxLDL(sm, blockSize, order, mcn, eMT, true);
        break;
    case eDSSFactorizationLLT:
        matrix = new SparseGridMtxLL(sm, blockSize, order, mcn, eMT, true);
        break;
    case eDSSFactorizationLU:
        matrix = new SparseGridMtxLU(sm, blockSize, order, mcn, eMT, true);
        break;
    case eDSSFastCG:
        matrix = new SparseGridMtxLDL(sm, blockSize, order, mcn, eMT, true);
        break;
    default:
        eMT->Writeln("Unknown solver type.");
        break;
    }

    if ( matrix == NULL ) {
        return NULL;
    }

    order = NULL;

    matrix->WriteStatistics(no_init_blocks, no_nonzeros);

    if ( fixed != NULL ) {
        sprintf(str, " Schur complement    : %ldx%ld\n", fixed->Count, fixed->Count);
        eMT->Write(str);
    }

    if ( tmpR ) {
        delete [] tmpR;
        tmpR = NULL;
    }

    tmpR = new double [ n_blocks * blockSize ];
    Array :: Clear(tmpR, 0, n_blocks * blockSize);
    return matrix;
}

bool DSSolver :: PreFactorize()
{
    if ( sm.neq == 0 || neq == 0 ) {
        eMT->Write("Can't factorize empty matrix!");
        return false;
    }

    switch ( MatrixType ) {
    case eDSSparseMatrix:
    {
        this->matrix = CreateNewSparseGridMtx();
        if ( this->matrix == NULL ) {
            return false;
        }
    }
    break;

    case eDSSkylineMatrix:
        switch ( SolverType ) {
        case eDSSFactorizationLDLT:
        //this->matrix = new SkyLineMtxLDL(sm,NULL,eMT);
        //break;
        case eDSSFactorizationLDLTIncomplete:
        case eDSSFactorizationLLTIncomplete:
        case eDSSFactorizationLLT:
        case eDSSFactorizationLU:
        case eDSSFastCG:
            eMT->Writeln("Skyline matrix doesn't support this solver type.");
            break;
        default:
            eMT->Writeln("Unknown solver type.");
            break;
        }

        break;
    default:
        eMT->Writeln("Unknown matrix storage type.");
        break;
    }

    if ( this->matrix == NULL ) {
        return false;
    }

    return true;
}

void DSSolver :: LoadZeros()
{
    if ( matrix == NULL ) {
        eMT->Writeln("First allocate the matrix.");
        return;
    }

    if ( matrix ) {
        matrix->LoadZeros();
    }

    m_eState = ISolver :: Allocated;
}

double &DSSolver :: ElementAt(int i, int j)
{
    return matrix->ElementAt(i, j);
}

bool DSSolver :: LoadNumbers(SparseMatrixF *sm)
{
    if ( matrix == NULL ) {
        eMT->Writeln("First load matrix.");
        return 0;
    }

    if ( this->sm.a == NULL ) {
        this->sm = * sm;
    }

    eMT->CS();
    eMT->Write("Loading matrix data");
    matrix->LoadMatrixNumbers(* sm);
    eMT->Writeln( eMT->MC_() );

    m_eState = ISolver :: Allocated;
    return 1;
}

void DSSolver :: WriteFactorizationInfo()
{
    switch ( SolverType ) {
    case eDSSFactorizationLDLTIncomplete:
        eMT->Write("Incomplete LDL factorization : ");
        break;
    case eDSSFactorizationLLTIncomplete:
        eMT->Write("Incomplete LL factorization : ");
        break;
    case eDSSFactorizationLDLT:
        eMT->Write("Numerical LDL factorization : ");
        break;
    case eDSSFactorizationLLT:
        eMT->Write("Numerical LLT factorization : ");
        break;
    case eDSSFactorizationLU:
        eMT->Write("Numerical LU factorization  : ");
        break;
    case eDSSFastCG:
        break;
    default:
        eMT->Writeln("Unknown solver type.");
        break;
    }
}

bool DSSolver :: ReFactorize()
{
    if ( matrix == NULL && matrixPD == NULL ) {
        eMT->Writeln("The matrix has to be loaded prior to the factorization.");
        return false;
    }

    eMT->CS();

    DenseMatrixArithmetics :: zero_pivots = 0;

    WriteFactorizationInfo();
    switch ( SolverType ) {
    case eDSSFactorizationLDLTIncomplete:
    case eDSSFactorizationLLTIncomplete:
    case eDSSFactorizationLDLT:
    case eDSSFactorizationLLT:
    case eDSSFactorizationLU:
        if ( matrixPD != NULL ) {
            matrixPD->SchurComplementFactorization();
        } else {
            matrix->Factorize();
        }

        break;
    case eDSSFastCG:
        break;
    default:
        eMT->Writeln("Unknown solver type.");
        break;
    }

    eMT->Writeln( eMT->MC_() );
    if ( DenseMatrixArithmetics :: zero_pivots > 0 ) {
        sprintf(str, "Warning: %ld zero pivot", DenseMatrixArithmetics :: zero_pivots);
        eMT->Write(str);
        if ( DenseMatrixArithmetics :: zero_pivots == 1 ) {
            eMT->Writeln("has been found during the factorization !!!");
        } else {
            eMT->Writeln("s have been found during the factorization !!!");
        }
    }

    m_eState = ISolver :: Factorized;
    return true;
}


bool DSSolver :: Solve(double *r, double *f)
{
    if ( m_eState == ISolver :: Allocated ) {
        if ( !Factorize() ) {
            return false;
        }
    }

    if ( matrix == NULL ) {
        if ( matrixD.N() ) {
            matrixD.DiagonalSolve(f, r);
        } else
        if ( r != f ) {
            Array :: Copy(f, r, neq);
        }

        return false;
    }

    if ( mcn ) {
        long i;
        long *perm = mcn->perm->Items;
        for ( i = 0; i < neq; i++ ) {
            tmpR [ perm [ i ] ] = f [ i ];
        }

        matrix->Solve(tmpR, tmpR);

        for ( i = 0; i < neq; i++ ) {
            r [ i ] = tmpR [ perm [ i ] ];
        }
    } else   {
        if ( neq % blockSize ) { // the vectors are not alligned
            Array :: Copy(f, tmpR, neq);
            matrix->Solve(tmpR, tmpR);
            Array :: Copy(tmpR, r, neq);
        } else   {
            matrix->Solve(f, r);
        }
    }

    return true;
}


long DSSolver :: Close()
{
    Dispose();
    return 0;
}

void DSSolver :: StartSolverWriteInfo()
{
    switch ( SolverType ) {
    case eDSSFactorizationLDLT:
        eMT->Writeln("Sparse matrix LDL factorization");
        break;
    case eDSSFactorizationLLT:
        eMT->Writeln("Sparse matrix LL factorization");
        break;
    case eDSSFactorizationLU:
        eMT->Writeln("Sparse matrix LU factorization");
        break;
    case eDSSFactorizationLDLTIncomplete:
        eMT->Writeln("Incomplete matrix LDL factorization");
        break;
    case eDSSFactorizationLLTIncomplete:
        eMT->Writeln("Incomplete matrix LL factorization");
        break;
    default:
        ;
    }

    sprintf( str, "Starting on          : %s", eMT->NowString() );
    eMT->Write(str);
}

void DSSolver :: EndSolverWriteInfo()
{
    eMT->Write("Ending on            : ");
    eMT->Write( eMT->NowString() );
    sprintf( str, "Computational effort        : %ld block multiplications", matrix->No_Multiplications() );
    eMT->Writeln(str);
}

bool DSSolver :: StartSolver()
{
    if ( sm.neq == 0 ) {
        eMT->Write("Can't factorize empty matrix.");
        return false;
    }

    //long usedmem =  0;    //GC.GetTotalMemory(true);
    clock_t solution_start = clock();

    StartSolverWriteInfo();
    if ( !PreFactorize() ) {
        return false;
    }

    if ( !ReFactorize() ) {
        return false;
    }

    clock_t solution_end = clock();
    EndSolverWriteInfo();

    double solution_duration = ( double ) ( solution_end - solution_start ) / CLOCKS_PER_SEC;
    sprintf(str, "Whole solution (end-start)  : %0.3f s", solution_duration);
    eMT->Writeln(str);

    //usedmem =  0;    //GC.GetTotalMemory(true);
    //sprintf(str,"%ld of used memory",usedmem);
    //eMT->Writeln(str);
    return true;
}

void DSSolver :: SetSM(SparseMatrixF *sm)
{
    this->sm = * sm;
}

double DSSolver :: GetFactorizationError()
{
    if ( neq == 0 ) {
        return 0.0;
    }

    if ( sm.a == NULL ) {
        return -1;
    }

    LargeVector r(neq);
    LargeVector b(neq);
    LargeVector x(neq);
    for ( int i = 0; i < neq; i++ ) {
        b.DataPtr() [ i ] = 1.0;
    }

    Solve( x.DataPtr(), b.DataPtr() );

    //r = b - Ax
    MulMatrixByVector( x.DataPtr(), r.DataPtr() );
    r.LinComb(-1, r, b);

    //delta_new = r r;
    double err = LargeVector :: InnerProduct(r.DataPtr(), r.DataPtr(), neq);
    return sqrt(err) / neq;
}


// counts number of block which combine both fixed and nonfixed DOFs.
// Each of these must be doubled to create only nonfixed and only fixed block
void DSSolver :: ExpandMCN(IntArrayList &mcn)
{
    long nb = mcn.Count / blockSize;
    for ( long b = 0; b < nb; b++ ) {
        bool fixed = false;
        bool nonfixed = false;
        for ( long i = 0; i < blockSize; i++ ) {
            long v = mcn [ b * blockSize + i ];
            if ( v < 0 ) {
                fixed = true;
            }

            if ( v >= 0 ) {
                nonfixed = true;
            }
        }

        if ( fixed && nonfixed ) {
            for ( long i = 0; i < blockSize; i++ ) {
                long cfix = -1;
                long v = mcn [ b * blockSize + i ];
                if ( v < 0 ) {
                    cfix = v;
                    mcn [ b * blockSize + i ] = -1;
                }

                mcn.Add(cfix);
            }
        }
    }
}

// Loads the "MatrixCodeNumbers" vector[n_blocks*block_size]
// each entry in the vector[i] can be :
// vector[i] >=  0  normal unknown                      (means the corresponding row in SparseMatrixF)
// vector[i] == -1  this entry is not used      (no corresponding row in SparseMatrixF)
bool DSSolver :: LoadMCN(unsigned long n_blocks, unsigned char block_size, long *mcn)
{
    this->blockSize = block_size;
    IntArrayList *mcn_order = new IntArrayList(n_blocks * block_size);
    mcn_order->Alloc();
    Array :: Copy(mcn, mcn_order->Items, n_blocks * block_size);
    return LoadMCN_int(mcn_order);
}

bool DSSolver :: LoadMCN_int(IntArrayList *mcn_order)
{
    ExpandMCN(* mcn_order);
    this->n_blocks = mcn_order->Count / blockSize;

    long i = 0;
    IntArrayList *mcn_perm = new IntArrayList(mcn_order->Count);
    mcn_perm->Alloc();

    long no_noncondensed_DOFs = 0;

    for ( i = 0; i < mcn_perm->Count; i++ ) {
        mcn_perm->Items [ i ] = -1;
    }

    for ( i = 0; i < mcn_perm->Count; i++ ) {
        long oi = mcn_order->Items [ i ];
        if ( oi >= 0 ) {
            mcn_perm->Items [ oi ] = i;
        } else if ( oi == -1 ) {} else {
            mcn_perm->Items [ -( oi + 2 ) ] = i;
            no_noncondensed_DOFs++;
        }
    }

    if ( this->mcn ) {
        delete this->mcn;
    }

    this->mcn = new Ordering(mcn_perm, mcn_order);

    return CreateFixedArray(no_noncondensed_DOFs);
}

bool DSSolver :: LoadMCN(IntArrayList &mcn)
{
    IntArrayList *mcn_new = new IntArrayList( ( ( mcn.Count - 1 ) / blockSize + 1 ) * blockSize );
    mcn_new->Alloc();
    mcn_new->Initialize(-1);
    Array :: Copy(mcn.Items, 0, mcn_new->Items, 0, mcn.Count);
    return LoadMCN_int(mcn_new);
}

bool DSSolver :: CreateFixedArray(long no_noncondensed_DOFs)
{
    fixed = new IntArrayList();
    lncn = new IntArrayList(no_noncondensed_DOFs);
    //lncn->Count = no_noncondensed_DOFs;
    lncn->Alloc();
    long lastFixed = -1;

    long i = 0;
    for ( i = 0; i < mcn->perm->Count; i++ ) {
        long oi = mcn->order->Items [ i ];
        if ( oi < -1 ) {
            fixed->AddIfNotLast(i / blockSize);
            //lncn->Add(-(oi+2));
            // Now I know that code numbers of noncondensed
            // nodes form a serie [(n-no_noncondensed)..(n)]
            long cn = -( oi + 2 );
            if ( cn < lastFixed ) {
                eMT->Writeln("Error: Array of fixed nodes is not monotonously growing !!");
                m_eState = ErrorInMCN;
            } else
            if ( cn >= neq ) {
                eMT->Writeln("Error: Array of fixed nodes contains index higher than (neq) !!");
                m_eState = ErrorInMCN;
            } else if ( cn < ( neq - no_noncondensed_DOFs ) ) {
                eMT->Writeln("Error: Array of fixed nodes contains index lower than (neq-no_fixed_DOFs) !!");
                m_eState = ErrorInMCN;
            } else {
                lncn->Items [ cn - ( neq - no_noncondensed_DOFs ) ] = lastFixed = cn;
            }
        }
    }

    return m_eState != ErrorInMCN;
}

/*
 * void MakeNewDOForder(long *lncn,long nbn)
 * {
 *      // The task of this method is to chose appropriate division of dofs into
 * // blocks. The mcm is the original division of dofs into nodes.
 * // The result should respect
 *
 * long b,i=0;
 *
 * IntArrayList* DOForder = new IntArrayList(n_blocks*blockSize);
 *
 * for (b=0; b<n_blocks;b++)
 * {
 *              bool normal = false;
 *              bool dirichlet = false;
 *              for (i=0; i<blockSize; i++)
 *              {
 *                      long cn = mcn->order->Items[b*blockSize+i];
 *                      if (cn==-4) dirichlet = true;
 *                      if (cn>=0) normal = true;
 *              }
 *
 *              // the whole block is eliminated
 *
 *              // the block is normal
 *
 *              // the whole block is dirichlet condition
 *
 *              // the block is divided into two new blocks
 *
 * }
 *
 * IntArrayList dofs_in_blocks(n_blocks);
 * dofs_in_blocks.Alloc();
 * dofs_in_blocks.Initialize();
 *
 * // compute the amount of dofs in each block
 * for (dof=0; dof<mcn->order->Count; dof++)
 *  if (mcn->order->Items[dof]>=0)
 *    dofs_in_blocks.Items[dof/blockSize]++;
 *
 * IntArrayList free_dofs_in_blocks(dofs_in_blocks);
 *
 * IntArrayList fixed_dofs(mcn->order->Count);
 * fixed_dofs.Alloc();
 * fixed_dofs.Initialize();
 *
 * IntArrayList fixed_blocks;
 * IntArrayList doubled_blocks;
 *
 * for (i=0; i<nbn; i++)
 * {
 *  long fixed_dof = lncn[i];
 *  long b_no = fixed_dof / blockSize;
 *
 *  if (--free_dofs_in_blocks[b_no] == 0) {
 *    fixed_blocks.Add(b_no);
 *  } else {
 *    if (dofs_in_blocks.Items[b_no]==1+free_dofs_in_blocks[b_no])
 *      doubled_blocks.Add(b_no);
 *  }
 * }
 *
 * for (i=0; i<doubled_blocks.Count; i++)
 * {
 *              long pdb = doubled_blocks.Items[i];
 *              if (free_dofs_in_blocks.Items[pdb]>0)
 * }
 * //this->mcn->
 *
 * //Ordering* new_node_oed
 * }
 */

bool DSSolver :: FactorizeSchur()
{
    if ( sm.neq == 0 || mcn == 0 ) {
        eMT->Writeln("Can't factorize empty matrix!");
        return false;
    }

    if ( mcn == 0 ) {
        eMT->Writeln("The MCN vector was not set!");
        return false;
    }

    StoreFixedLastPermutation_dom_order();

    SparseGridMtx *matrix = CreateNewSparseGridMtx(fixed);
    matrixPD = new SparseGridMtxPD(matrix, fixed->Count);

    ReFactorize();
    return true;
}

// Sorts the domain internal DOFs first and the fixed DOFs at the end
// Mostly gives dom_order=indentity as JK sorts all the fixed DOFs at the end already
void DSSolver :: StoreFixedLastPermutation_dom_order()
{
    // store "fixed-last" permutation
    if ( dom_order ) {
        delete dom_order;
    }

    dom_order = new IntArrayList(sm.neq);
    dom_order->Alloc();

    int lastdo = 0;

    long *dmap = new long [ neq ];
    Array :: Clear(dmap, 0, neq);

    if ( lncn->Count ) {
        long *pI = lncn->Items;
        long idx = lncn->Count - 1;
        for ( long *p = pI + idx; p >= pI; p--, idx-- ) {
            dmap [ * p ] = ~( idx + ( neq - lncn->Count ) );
        }
    }

    for ( long d = 0; d < neq; d++ ) {
        if ( dmap [ d ] == 0 ) {
            dmap [ d ] = ~lastdo;
            dom_order->Items [ lastdo++ ] = d;
        }
    }

    Array :: Copy(lncn->Items, dom_order->Items + neq - lncn->Count, lncn->Count);
    delete dmap;
    dmap = NULL;
}

void DSSolver :: condense(double *a, double *lhs, double *rhs, long tc)
{
    if ( tc == 1 ) {
        if ( m_eState != ISolver :: Factorized || matrixPD == NULL ) {
            if ( !FactorizeSchur() ) {
                return;
            }
        }

        matrixPD->WriteCondensedMatrixA22(a, mcn, lncn);

        long i = 0;
        long *perm = mcn->perm->Items;
        long *dor = dom_order->Items;
        for ( i = 0; i < neq; i++ ) {
            tmpR [ perm [ dor [ i ] ] ] = rhs [ i ];
        }

        matrixPD->Sub_A21_A11inv(tmpR);

        for ( i = 0; i < neq; i++ ) {
            rhs [ i ] = tmpR [ perm [ dor [ i ] ] ];
        }
    } else if ( tc == 2 )       {
        long i;
        long *perm = mcn->perm->Items;
        long *dor = dom_order->Items;
        for ( i = 0; i < neq - lncn->Count; i++ ) {
            tmpR [ perm [ dor [ i ] ] ] = rhs [ i ];
        }

        for ( i = neq - lncn->Count; i < neq; i++ ) {
            tmpR [ perm [ dor [ i ] ] ] = lhs [ i ];
        }

        matrixPD->Sub_A11inv_A12(tmpR);

        for ( i = 0; i < neq; i++ ) {
            lhs [ i ] = tmpR [ perm [ dor [ i ] ] ];
        }
    } else if ( tc == 3 )       {
        long i;
        long *perm = mcn->perm->Items;
        long *dor = dom_order->Items;
        for ( i = 0; i < neq - lncn->Count; i++ ) {
            tmpR [ perm [ dor [ i ] ] ] = rhs [ i ];
        }

        matrixPD->SolveA11(tmpR);

        for ( i = 0; i < neq - lncn->Count; i++ ) {
            lhs [ i ] = tmpR [ perm [ dor [ i ] ] ];
        }
    }
}

void DSSolver :: GetA12block(double *pA12)
{
    sm.GetA12block(pA12, lncn->Count);
}


// c = A * b
void DSSolver :: MulMatrixByVector(double *b, double *c)
{
    if ( SolverType == eDSSFastCG ) {
        LargeVectorAttach vb(neq, b);
        LargeVectorAttach vx(neq, c);
        matrix->MultiplyByVector(vb, vx);
    }

    if ( SolverType == eDSSFactorizationLU ) {
        sm.MulNonsymMatrixByVector(b, c);
    } else {
        sm.MulMatrixByVector(b, c);
    }

    //sm.mxv_scr(b,c);
}

//multiply by a scalar
void DSSolver :: times(double x) {
    matrix->times(x);
}

// positive number n is success after n iterations
// negative number n is failed after -n iterations
int DSSolver :: PreCG(double *b_, double *x_, double epsilon, int max_iter)
{
    double epsilon2 = epsilon * epsilon;

    LargeVectorAttach b(neq, b_);
    LargeVectorAttach x(neq, x_);

    //r = b - Ax
    LargeVector r(neq);
    MulMatrixByVector( x.DataPtr(), r.DataPtr() );
    //  for (int i=0; i<neq; i++)
    //    r[i] = b[i] - r[i];
    r.LinComb(-1, r, b);

    //d = M r
    LargeVector d(neq);
    Solve( d.DataPtr(), r.DataPtr() );

    //delta_new = r d;
    double delta_new = LargeVector :: InnerProduct(r.DataPtr(), d.DataPtr(), neq);

    if ( delta_new == 0 ) {
        eMT->Writeln("Zero right hand side vector.");
        return 0;
    }

    double delta_0 = delta_new;

    eMT->Writeln("Starting PCG ");
    eMT->Write(" epsilon : ");
    eMT->Write(epsilon);
    eMT->Writeln();
    eMT->Write(" residual : ");
    eMT->Writeln();

    LargeVector q(neq);
    LargeVector s(neq);
    double alfa = 0;
    int iter;
    for ( iter = 1; iter < max_iter; iter++ ) {
        // q = A.d
        MulMatrixByVector( d.DataPtr(), q.DataPtr() );

        // a1 = q.d
        double a1 = LargeVector :: InnerProduct(d.DataPtr(), q.DataPtr(), neq);
        alfa = delta_new / a1;

        //x = x+alfa*d;
        x.AddMult(alfa, d);
        //LargeVector::AddMul(x,alfa,d.DataPtr(),neq);

        if ( iter % 50 == 0 ) {
            // r = b - Ax
            MulMatrixByVector( x.DataPtr(), r.DataPtr() );
            r.LinComb(-1, r, b);
            //for (int i=0; i<neq; i++)
            //  r[i] = b[i] - r[i];
        } else   {
            // r = r - alfa*q
            r.AddMult(-alfa, q);
        }

        // s = M r
        Solve( s.DataPtr(), r.DataPtr() );

        double delta_old = delta_new;

        delta_new = LargeVector :: InnerProduct(r.DataPtr(), s.DataPtr(), neq);

        double beta = delta_new / delta_old;

        // d = s + beta*d
        //for (int i=0; i<neq; i++)
        //  d[i] = s[i] + beta*d[i];
        d.LinComb(beta, d, s);

        eMT->Write("  ");
        eMT->Write( sqrt(delta_new / delta_0) );
        eMT->Writeln();

        //if (x_i) //je dost presne then konec
        if ( fabs(delta_new) < fabs(epsilon2 * delta_0) ) {
            eMT->Write("Success after ");
            eMT->Write(iter);
            eMT->Writeln(" iterations");
            return iter;
        }
    }

    eMT->Write("Failed after ");
    eMT->Write(iter);
    eMT->Writeln(" iterations");
    return -iter;
}

// positive number n is success after n iterations
// negative number n is failed after -n iterations
int DSSolver :: CG(double *b_, double *x_, double epsilon, int max_iter)
{
    double epsilon2 = epsilon * epsilon;

    LargeVectorAttach b(neq, b_);
    LargeVectorAttach x(neq, x_);

    //r = b - Ax
    LargeVector r(neq);
    MulMatrixByVector( x.DataPtr(), r.DataPtr() );

    //for (int i=0; i<neq; i++)
    //  r[i] = b[i] - r[i];
    r.LinComb(-1, r, b);

    //d = r
    LargeVector d(r);

    //delta_new = r r;
    double delta_new = LargeVector :: InnerProduct(r.DataPtr(), r.DataPtr(), neq);

    if ( delta_new == 0 ) {
        eMT->Writeln("Zero right hand side vector.");
        //return 0;
    }

    double delta_0 = delta_new;

    eMT->Writeln("Starting PCG ");
    eMT->Write(" epsilon : ");
    eMT->Write(epsilon);
    eMT->Writeln();
    eMT->Write(" residual : ");
    eMT->Writeln();

    LargeVector q(neq);
    LargeVector s(neq);
    double alfa = 0;
    int iter;
    for ( iter = 1; iter < max_iter; iter++ ) {
        // q = A.d
        MulMatrixByVector( d.DataPtr(), q.DataPtr() );

        // a1 = q.d
        double a1 = LargeVector :: InnerProduct(d.DataPtr(), q.DataPtr(), neq);
        alfa = delta_new / a1;

        //x = x+alfa*d;
        x.AddMult(alfa, d);
        //LargeVector::AddMul(x.DataPtr(),alfa,d.DataPtr(),neq);

        if ( iter % 50 == 0 ) {
            // r = b - Ax
            MulMatrixByVector( x.DataPtr(), r.DataPtr() );
            r.LinComb(-1, r, b);
            //for (int i=0; i<neq; i++)
            //  r[i] = b[i] - r[i];
        } else   {
            // r = r - alfa*q
            r.AddMult(-alfa, q);
        }

        double delta_old = delta_new;

        delta_new = LargeVector :: InnerProduct(r.DataPtr(), r.DataPtr(), neq);

        double beta = delta_new / delta_old;
        double conv = delta_new / delta_0;

        eMT->Write("  ");
        eMT->Write( sqrt(conv) );
        //    eMT->Write(sqrt(delta_new));
        eMT->Writeln();

        // d = r + beta*d
        //for (int i=0; i<neq; i++)
        //  d[i] = r[i] + beta*d[i];
        d.LinComb(beta, d, r);

        //if (x_i) //je dost presne then konec
        if ( fabs(conv) < fabs(epsilon2) ) {
            return iter;
        }
    }

    return -iter;
}




DSS_NAMESPASE_END
