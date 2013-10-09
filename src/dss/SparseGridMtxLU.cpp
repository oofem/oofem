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

#include "SparseGridMtxLU.h"

DSS_NAMESPASE_BEGIN

// Allocates new space according to bskl and reads old matrix with respect
// to permutation blockP
SparseGridMtxLU :: SparseGridMtxLU(SparseMatrixF &sm, long block_size, Ordering *block_order, MathTracer *eMT, bool load_data) :
    SparseGridMtx(sm, block_size, block_order, eMT)
{
    IConectMatrix *bskl = block_order->cm;

    this->AlocateMemoryByPattern(bskl);
    ComputeBlocks();

    if ( sm.a != NULL && load_data ) {
        LoadMatrixNumbers(sm);
    }
}

// Allocates new space according to bskl and reads old matrix with respect
// to permutation blockP
SparseGridMtxLU :: SparseGridMtxLU(SparseMatrixF &sm, long block_size, Ordering *block_order, Ordering *node_order, MathTracer *eMT, bool load_data) :
    SparseGridMtx(sm, block_size, block_order, node_order, eMT)
{
    IConectMatrix *bskl = block_order->cm;

    this->AlocateMemoryByPattern(bskl);
    ComputeBlocks();

    if ( sm.a != NULL && load_data ) {
        LoadMatrixNumbers(sm);
    }
}

double &SparseGridMtxLU :: ElementAt(int i, int j)
{
    long *blockP = this->block_order->perm->Items;
    long *nodeP = this->node_order ? node_order->perm->Items : NULL;

    long aux_bi_idx = 0;        // index of the last found block in the column
    long aux_bj_idx = 0;

    long nj = ( nodeP == NULL ) ? j : nodeP [ j ];
    long sj = nj % block_size;
    long bj = nj / block_size;
    long nbj = ( blockP == NULL ) ? bj : blockP [ bj ];

    SparseGridColumn *Column_nbj = Columns [ nbj ];
    SparseGridColumn *Column_nbi = NULL;

    double *column_nbj_bi_data = NULL;
    double *row_nbi_bj_data = NULL;

    long ni = ( nodeP == NULL ) ? i : nodeP [ i ];
    long si = ni % block_size;
    long bi = ni / block_size;
    long nbi = ( blockP == NULL ) ? bi : blockP [ bi ];

    Column_nbi = Columns [ nbi ];
    aux_bj_idx = aux_bi_idx = -1;
    if ( nbi < nbj ) {
        aux_bi_idx = Column_nbj->FindExistingBlockIndex(nbi);
        column_nbj_bi_data = Columns_data + Column_nbj->column_start_idx + aux_bi_idx * block_storage;
    } else if ( nbi > nbj ) {
        aux_bj_idx = Column_nbi->FindExistingBlockIndex(nbj);
        row_nbi_bj_data = Rows_data + Column_nbi->column_start_idx + aux_bj_idx * block_storage;
    }

    if ( aux_bi_idx < 0 && aux_bj_idx < 0 ) {   //diagonal index
        return Diagonal_data [ nbi * block_storage + si + sj * block_size ];
    } else if ( aux_bi_idx >= 0 ) { // write into column bj
        return column_nbj_bi_data [ block_size * sj + si ];
    } else {                                    // write into row bi
        return row_nbi_bj_data [ block_size * si + sj ];
    }
}

void SparseGridMtxLU :: LoadZeros()
{
    memset( this->Columns_data, 0, columns_data_length * sizeof( double ) );
    memset( this->Rows_data, 0, columns_data_length * sizeof( double ) );
    memset( this->Diagonal_data, 0, n_blocks * block_storage * sizeof( double ) );

    for ( long bj = 0; bj < n_blocks; bj++ ) {
        for ( long j = 0; j < block_size; j++ ) {
            Diagonal_data [ bj * block_storage + block_size * j + j ] = 1.0;
        }
    }

    long *blockP = this->block_order->perm->Items;
    long *nodeP = this->node_order ? node_order->perm->Items : NULL;

    long neq = n_blocks * block_size - noDummyDOFs;
    for ( long i = 0; i < neq; i++ ) {
        long ni = ( nodeP == NULL ) ? i : nodeP [ i ];
        long si = ni % block_size;
        long bi = ni / block_size;
        long nbi = ( blockP == NULL ) ? bi : blockP [ bi ];
        Diagonal_data [ nbi * block_storage + si + si * block_size ] = 0.0;
    }
}

void SparseGridMtxLU :: LoadMatrixNumbers(SparseMatrixF &sm)
{
    LoadZeros();
    long *blockP = this->block_order->perm->Items;
    long *nodeP = this->node_order ? node_order->perm->Items : NULL;

    long aux_bi_idx = 0;        // index of the last found block in the column
    long aux_bj_idx = 0;

    for ( unsigned long j = 0; j < sm.neq; j++ ) {
        long nj = ( nodeP == NULL ) ? j : nodeP [ j ];

        long sj = nj % block_size;
        long bj = nj / block_size;
        long nbj = ( blockP == NULL ) ? bj : blockP [ bj ];

        SparseGridColumn *Column_nbj = Columns [ nbj ];
        SparseGridColumn *Column_nbi = NULL;

        double *column_nbj_bi_data = NULL;
        double *row_nbi_bj_data = NULL;

        long old_nbi = -1;

        for ( unsigned long ad = sm.Adr(j); ad < sm.Adr(j + 1); ad++ ) {
            long i = ( long ) sm.Ci(ad);
            long ni = ( nodeP == NULL ) ? i : nodeP [ i ];
            double val = sm.a [ ad ];

            long si = ni % block_size;
            long bi = ni / block_size;
            long nbi = ( blockP == NULL ) ? bi : blockP [ bi ];

            if ( old_nbi != nbi ) {
                Column_nbi = Columns [ nbi ];
                aux_bj_idx = aux_bi_idx = -1;
                if ( nbi < nbj ) {
                    aux_bi_idx = Column_nbj->FindExistingBlockIndex(nbi);
                    column_nbj_bi_data = Columns_data + Column_nbj->column_start_idx + aux_bi_idx * block_storage;
                } else if ( nbi > nbj ) {
                    aux_bj_idx = Column_nbi->FindExistingBlockIndex(nbj);
                    row_nbi_bj_data = Rows_data + Column_nbi->column_start_idx + aux_bj_idx * block_storage;
                }

                old_nbi = nbi;
            }

            if ( aux_bi_idx < 0 && aux_bj_idx < 0 ) {           //diagonal index
                Diagonal_data [ nbi * block_storage + si + sj * block_size ] = val;
            } else if ( aux_bi_idx >= 0 ) {     // write into column bj
                column_nbj_bi_data [ block_size * sj + si ] = val;
            } else {                                            // write into row bi
                row_nbi_bj_data [ block_size * si + sj ] = val;
            }

            nonzeros++;
        }
    }

    //long bl = n_blocks-1;
    //long nbl = (blockP==NULL)?bl:blockP[bl];

    // Dummy DOFs
    //for (long i=block_size-1; i>=block_size-noDummyDOFs; i--)
    //Columns_data[nbl*block_storage + block_size*i + i] = 1.0;
    //this->SetValue(nbl,nbl,i,i,1.0,aux_bi_idx,aux_bj_idx);

    /*long dzeros = 0;
     *      //Check diagonal for zeros
     *      for (long bj=0; bj<n_blocks; bj++)
     *              for (long j=0; j<block_size; j++)
     *                      if (Columns_data[bj*block_storage + block_size*j + j] == 0.0)
     *                      {
     *                              Columns_data[bj*block_storage + block_size*j + j] = 1.0;
     *                              dzeros++;
     *                      }*/
    //if (dzeros>0)
    //{
    //char str[128];
    //sprintf(str,"Matrix had %ld zeros on diagonal!",dzeros);
    //eMT->Writeln(str);
    //}
}

SparseGridMtxLU :: ~SparseGridMtxLU()
{
    if ( Columns_data ) {
        delete Columns_data;
        Columns_data = NULL;
    }

    if ( Rows_data ) {
        delete Rows_data;
        Rows_data = NULL;
    }

    if ( Diagonal_data ) {
        delete Diagonal_data;
        Diagonal_data = NULL;
    }
}

void SparseGridMtxLU :: AlocateMemoryByPattern(IConectMatrix *bskl)
{
    clock_t ts = clock();
    this->n_blocks = bskl->N();
    columns_data_length = 0;
    SparseGridColumn *column;

    for ( long bj = 0; bj < n_blocks; bj++ ) {
        IntArrayList *indexes = bskl->DetachIndexesAboveDiagonalInColumn(bj);
        Columns [ bj ] = column = new SparseGridColumn(indexes);
        column->column_start_idx = columns_data_length;
        columns_data_length += column->Entries * block_storage;
    }

    if ( columns_data_length == 0 ) {   // This array seems to need nonzero size, needed for the fixed statement
        columns_data_length = 1;
    }

    long memn = ( ( long ) columns_data_length ) * sizeof( double ) / 1024;
    char str [ 256 ];
    //Write(" aloc "+memn.ToString("### ### ##0")+ "kB..");
    sprintf(str, " sparse matrix size  : %ld kB..", memn);
    Write(str);
    this->Columns_data = new double [ columns_data_length ];
    this->Rows_data = new double [ columns_data_length ];
    this->Diagonal_data = new double [ n_blocks * block_storage ];

    clock_t end = clock();
    double duration = ( double ) ( end - ts ) / CLOCKS_PER_SEC;
    sprintf(str, "%.3f s", duration);
    Writeln(str);
    //Write(duration.TotalSeconds.ToString("g4")+" s");
}

// y = Ax
void SparseGridMtxLU :: MultiplyByVector(const LargeVectorAttach &x, LargeVectorAttach &y)
{
    /*y.Initialize();
     *
     *      double* py = y.DataPtr();
     *      double* px = x.DataPtr();
     *      for (long bj=0; bj<n_blocks; bj++)
     *      {
     *              SparseGridColumn& columnJ = *Columns[bj];
     *              //Columns[bj].MultBlockByVectorSym(BlockArith,x,ref y,block_size,bj);
     *              long cnt = columnJ.Entries;
     *              for (long idx = 0; idx<cnt; idx++)
     *              {
     *
     *                      long bi = columnJ.IndexesUfa->Items[idx];
     *                      //DenseMatrix block = this->Blocksfa[idx];
     *                      long block = columnJ.column_start_idx + idx*block_storage;
     *                      BlockArith->AddMultBlockByVectorSym(Columns_data+block,px,py,block_size*bi,bj*block_size);
     *              }
     *              //BlockArith.MultDiagonalBlockByVector(DiagonalBlocks[bj],x, ref y,block_size,bj*block_size);
     *              BlockArith->MultDiagonalBlockByVector(Columns_data+bj*block_storage,px+bj*block_size, py+bj*block_size);
     *      }*/
    // x;y;
} //MultiplyByVector

void SparseGridMtxLU :: Factorize()
{
    long bi, bj, min_bi_J = 0;

    // This is a pattern of blocks in J-th column
    long *p_blockJ_pattern = new long [ n_blocks + 1 ];
    memset( p_blockJ_pattern, 0, ( n_blocks + 1 ) * sizeof( long ) );

    long newdotcnount = ( long ) ceil( ( double ) n_blocks / 24.0 );
    //long newdotcnount = n / 24;
    BlockArith->zero_pivots = 0;

    double *cd = this->Columns_data;
    double *rd = this->Rows_data;
    double *dd = this->Diagonal_data;
    {
        //double* idd = dd;
        long Djj = 0;
        eMT->act_block = 0;
        for ( bj = 0; bj < n_blocks; bj++ ) {
            SparseGridColumn &columnJ = * Columns [ bj ];
            long noJentries = columnJ.Entries;
            if ( noJentries > 0 ) {
                long *columnJentries = columnJ.IndexesUfa->Items;

                double *pBkj = cd + columnJ.column_start_idx;
                double *pBij = pBkj;
                double *pAkj = rd + columnJ.column_start_idx;
                double *pAij = pAkj;


                //columnJ.DrawColumnPattern(p_blockJ_pattern,ref min_bi_J,bj);
                min_bi_J = bj;
                for ( long i = noJentries - 1; i >= 0; i-- ) {
                    p_blockJ_pattern [ min_bi_J = columnJentries [ i ] ] = ~( i * block_storage );
                }

                // eliminate above diagonal
                for ( long idx_J = 0; idx_J < noJentries; idx_J++ ) {
                    bi = columnJentries [ idx_J ];

                    SparseGridColumn &columnI = * Columns [ bi ];
                    long noIentries = columnI.Entries;

                    if ( noIentries > 0 ) {
                        double *pBki = cd + columnI.column_start_idx + ( noIentries - 1 ) * block_storage;
                        double *pAki = rd + columnI.column_start_idx + ( noIentries - 1 ) * block_storage;

                        long *columnIentries = columnI.IndexesUfa->Items;
                        for ( long *columnIentry = columnIentries + noIentries - 1; columnIentry >= columnIentries; pAki -= block_storage, pBki -= block_storage ) {
                            long idx_K = p_blockJ_pattern [ * columnIentry-- ];
                            if ( idx_K == 0 ) {
                                continue;
                            }

                            BlockArith->SubATBproduct(pBij, pAki, pBkj + ~idx_K);
                            BlockArith->SubATBproduct(pAij, pBki, pAkj + ~idx_K);
                        }
                    }

                    BlockArith->ULT_BlockSolve(dd + bi * block_storage, pAij);
                    pAij += block_storage;
                    pBij += block_storage;
                }

                // compute the diagonal and divide by it
                //DenseMatrix Djj = DiagonalBlocks[bj];// Diagonal
                for ( long idx = noJentries - 1; idx >= 0; idx-- ) {
                    bi = columnJentries [ idx ];
                    //Clear pattern
                    p_blockJ_pattern [ bi ] = 0;

                    int ij = columnJ.column_start_idx + idx * block_storage;
                    BlockArith->SubATBproduct(dd + Djj, rd + ij, cd + ij);
                }
            }

            // Factorize diagonal block
            BlockArith->LU_Decomposition(dd + Djj);

            if ( ( bj % newdotcnount ) == 0 ) {
                Write(".");
            }

            Djj += block_storage;
            eMT->act_block += block_size;
            if ( eMT->break_flag ) {
                break;
            }
        }
    }
    delete [] p_blockJ_pattern;
    ComputeBlocks();
}

// This functions computes LDL' decomposition of first (nblocks-fixed_bn) columns
// The resulting submatrix will contain the schur complement A22 - A21 * A11^(-1) * A12
void SparseGridMtxLU :: SchurComplementFactorization(int fixed_blocks)
{
    long bi, bj, min_bi_J = 0;

    long blocks_to_factor = n_blocks - fixed_blocks;

    // This is a pattern of blocks in J-th column
    long *p_blockJ_pattern = new long [ n_blocks + 1 ];
    memset( p_blockJ_pattern, 0, ( n_blocks + 1 ) * sizeof( long ) );

    //long newdotcnount = (long)ceil((double)n_blocks / 24.0);
    //long newdotcnount = n / 24;
    BlockArith->zero_pivots = 0;

    double *cd = this->Columns_data;
    double *rd = this->Rows_data;
    double *dd = this->Diagonal_data;
    {
        //double* idd = dd;
        long Djj = 0;
        eMT->act_block = 0;
        for ( bj = 0; bj < n_blocks; bj++ ) {
            SparseGridColumn &columnJ = * Columns [ bj ];
            long noJentries = columnJ.Entries;
            if ( noJentries > 0 ) {
                long *columnJentries = columnJ.IndexesUfa->Items;

                double *pBkj = cd + columnJ.column_start_idx;
                double *pBij = pBkj;
                double *pAkj = rd + columnJ.column_start_idx;
                double *pAij = pAkj;

                //columnJ.DrawColumnPattern(p_blockJ_pattern,ref min_bi_J,bj);
                min_bi_J = bj;

                if ( bj < blocks_to_factor ) {                  // This is normal scalar product
                    for ( long i = noJentries - 1; i >= 0; i-- ) {
                        p_blockJ_pattern [ min_bi_J = columnJentries [ i ] ] = ~( i * block_storage );
                    }
                } else {                                                      // Here we would read from the A22 matrix, therefore the scalar product is shorter (bi<blocks_to_factor)
                    for ( long i = noJentries - 1; i >= 0; i-- ) {
                        if ( columnJentries [ i ] < blocks_to_factor ) {
                            p_blockJ_pattern [ min_bi_J = columnJentries [ i ] ] = ~( i * block_storage );
                        }
                    }
                }

                // eliminate above diagonal
                for ( long idx_J = 0; idx_J < noJentries; idx_J++ ) {
                    bi = columnJentries [ idx_J ];

                    SparseGridColumn &columnI = * Columns [ bi ];
                    long noIentries = columnI.Entries;

                    if ( noIentries > 0 ) {
                        double *pBki = cd + columnI.column_start_idx + ( noIentries - 1 ) * block_storage;
                        double *pAki = rd + columnI.column_start_idx + ( noIentries - 1 ) * block_storage;

                        long *columnIentries = columnI.IndexesUfa->Items;
                        for ( long *columnIentry = columnIentries + noIentries - 1; columnIentry >= columnIentries; pAki -= block_storage, pBki -= block_storage ) {
                            long idx_K = p_blockJ_pattern [ * columnIentry-- ];
                            if ( idx_K == 0 ) {
                                continue;
                            }

                            BlockArith->SubATBproduct(pBij, pAki, pBkj + ~idx_K);
                            BlockArith->SubATBproduct(pAij, pBki, pAkj + ~idx_K);
                        }
                    }

                    if ( bi < blocks_to_factor ) {
                        BlockArith->ULT_BlockSolve(dd + bi * block_storage, pAij);
                    }

                    pAij += block_storage;
                    pBij += block_storage;
                }

                // compute the diagonal and divide by it
                //DenseMatrix Djj = DiagonalBlocks[bj];// Diagonal
                for ( long idx = noJentries - 1; idx >= 0; idx-- ) {
                    bi = columnJentries [ idx ];
                    //Clear pattern
                    p_blockJ_pattern [ bi ] = 0;

                    if ( bi >= blocks_to_factor ) {
                        break;
                    }

                    int ij = columnJ.column_start_idx + idx * block_storage;
                    BlockArith->SubATBproduct(dd + Djj, rd + ij, cd + ij);
                }
            }

            // Factorize diagonal block
            if ( bj < blocks_to_factor ) {
                BlockArith->LU_Decomposition(dd + Djj);
            }

            //if ((bj % newdotcnount) == 0)
            //Write(".");

            Djj += block_storage;
            eMT->act_block += block_size;
            if ( eMT->break_flag ) {
                break;
            }
        }
    }
    delete [] p_blockJ_pattern;
    ComputeBlocks();
}

void SparseGridMtxLU :: Solve(double *b, double *x)
{
    if ( x != b ) {
        Array :: Copy(b, x, n);
    }

    SolveLU(x);
}

void SparseGridMtxLU :: SolveLV(const LargeVector &b, LargeVector &x)
{
    if ( & x != & b ) {
        x.Initialize(b, 0, n);
    }

    SolveLU( x.DataPtr() );
}

// x = L^(-1) * x
void SparseGridMtxLU :: ForwardSubstL(double *x, long fixed_blocks)
{
    if ( this->N() == 0 ) {
        return;
    }

    long blocks_to_factor = n_blocks - fixed_blocks;

    int bi;
    long *ord = this->block_order->order->Items;

    // forward substitution L z = f  --> z
    for ( bi = 1; bi < blocks_to_factor; bi++ ) {
        SparseGridColumn &rowI = * Columns [ bi ];
        int no = rowI.Entries;
        if ( no > 0 ) {
            double *dst = x + block_size * ord [ bi ];
            double *Aij = Rows_data + rowI.column_start_idx;

            long *idxs = rowI.IndexesUfa->Items;
            //r[i] -= Lji * r[j]
            for ( int idx = 0; idx < no; idx++, Aij += block_storage ) {
                BlockArith->SubMultTBlockByVector(Aij, x + block_size * ord [ * ( idxs++ ) ], dst);
            }
        }
    }
}

void SparseGridMtxLU :: BackSubstU(double *x, long fixed_blocks)
{
    int bi;
    long *ord = this->block_order->order->Items;
    if ( this->N() == 0 ) {
        return;
    }

    double *Dii = Diagonal_data + block_storage * ( n_blocks - fixed_blocks - 1 );
    // back substitution U r = z'
    for ( bi = n_blocks - fixed_blocks - 1; bi >= 0; bi--, Dii -= block_storage ) {
        BlockArith->LU_Solve(Dii, x + block_size * ord [ bi ]);

        SparseGridColumn &columnI = * Columns [ bi ];
        int no = columnI.Entries;
        if ( no > 0 ) {
            double *src = x + block_size * ord [ bi ];
            double *Aij = Columns_data + columnI.column_start_idx;

            long *idxs = columnI.IndexesUfa->Items;
            //r[j] -= Uij * r[i]
            for ( int idx = 0; idx < no; idx++, Aij += block_storage ) {
                BlockArith->SubMultBlockByVector(Aij, src, x + block_size * ord [ * ( idxs++ ) ]);
            }
        }
    }
}

// x = A^(-1) * b
void SparseGridMtxLU :: SolveLU(double *x, long fixed_blocks)
{
    ForwardSubstL(x, fixed_blocks);
    BackSubstU(x, fixed_blocks);
}

/// <summary>  y -= U12 * x  </summary>
void SparseGridMtxLU :: SubMultU12(double *px, double *py, long fixed_blocks)
{
    long schur_shift = n_blocks - fixed_blocks;
    if ( fixed_blocks == 0 ) {
        return;
    }

    long *ord = block_order->order->Items;
    double *cd = this->Columns_data;
    for ( long bj = schur_shift; bj < n_blocks; bj++ ) {
        SparseGridColumn &columnJ = * Columns [ bj ];

        long no = columnJ.Entries;
        if ( no == 0 ) {
            continue;
        }

        double *src = px + block_size * ord [ bj ];
        double *Aij = cd +  columnJ.column_start_idx;

        long *idxs = columnJ.IndexesUfa->Items;
        //y[i] -= Lji^T * x[j]
        for ( long idx = 0; idx < no && * idxs < schur_shift; idx++, Aij += block_storage ) {
            BlockArith->SubMultBlockByVector(Aij, src, py + block_size * ord [ * ( idxs++ ) ]);
        }
    }
}

/// <summary>  y -= L12 * x  </summary>
void SparseGridMtxLU :: SubMultL21(double *px, double *py, long fixed_blocks)
{
    long schur_shift = n_blocks - fixed_blocks;
    if ( fixed_blocks == 0 ) {
        return;
    }

    long *ord = block_order->order->Items;
    for ( long bj = schur_shift; bj < n_blocks; bj++ ) {
        SparseGridColumn &columnJ = * Columns [ bj ];

        long no = columnJ.Entries;
        if ( no == 0 ) {
            continue;
        }

        double *dst = py + block_size * ord [ bj ];
        double *Aij = Rows_data  +  columnJ.column_start_idx;

        long *idxs = columnJ.IndexesUfa->Items;
        //y[i] -= Lji * x[j]
        for ( long idx = 0; idx < no && * idxs < schur_shift; idx++, Aij += block_storage ) {
            BlockArith->SubMultTBlockByVector(Aij, px + block_size * ord [ * ( idxs++ ) ], dst);
        }
    }
}


void SparseGridMtxLU :: SolveA11(double *x, long fixed_blocks)
{
    SolveLU(x, fixed_blocks);
}

//   f2 - A21*A11inv*f1 =
// = f2 - A21*U11inv*L11inv*f1 =
// = f2 - L21*L11inv*f1
void SparseGridMtxLU :: Sub_A21_A11inv(double *x, long fixed_blocks)
{
    ForwardSubstL(x, fixed_blocks);
    SubMultL21(x, x, fixed_blocks);
}

//   A11inv(f1-A12*r1) =
// = U11inv(L11inv*f1 - L11inv*A12*r2) =
// = U11inv(L11inv*f1 - U12*r2)
void SparseGridMtxLU :: Sub_A11inv_A12(double *x, long fixed_blocks)
{
    SubMultU12(x, x, fixed_blocks);
    BackSubstU(x, fixed_blocks);
}

//resulting matrix a is stored by columns
void SparseGridMtxLU :: WriteCondensedMatrixA22(double *a, Ordering *mcn, IntArrayList *lncn)
{
    long blockSize = BlockSize();
    long i, j;
    long nbn = lncn->Count;
    long aux_bi_idx = 0, aux_bj_idx = 0;
    for ( j = 0; j < lncn->Count; j++ ) {
        long nj = mcn->perm->Items [ lncn->Items [ j ] ];
        long nbj = nj / blockSize;
        long nsj = nj % blockSize;
        long obj = block_order->perm->Items [ nbj ];

        for ( i = 0; i < lncn->Count; i++ ) {
            long ni = mcn->perm->Items [ lncn->Items [ i ] ];
            long nbi = ni / blockSize;
            long nsi = ni % blockSize;
            long obi = block_order->perm->Items [ nbi ];

            double val = GetValue(obi, obj, nsi, nsj, aux_bi_idx, aux_bj_idx);
            a [ ( j ) * nbn + ( i ) ] = val;
        }
    }
}

DSS_NAMESPASE_END
