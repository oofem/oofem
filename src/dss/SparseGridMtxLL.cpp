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

#include "SparseGridMtxLL.h"

DSS_NAMESPASE_BEGIN

// Allocates new space according to bskl and reads old matrix with respect
// to permutation blockP
SparseGridMtxLL :: SparseGridMtxLL(SparseMatrixF &sm, long block_size, Ordering *block_order, MathTracer *eMT, bool load_data) :
    SparseGridMtx(sm, block_size, block_order, eMT)
{
    this->BlockArith->prefered_decomposition = eLL_decomposition;
    IConectMatrix *bskl = block_order->cm;

    this->AlocateMemoryByPattern(bskl);
    ComputeBlocks();

    if ( sm.a != NULL && load_data ) {
        LoadMatrixNumbers(sm);
    }
}

// Allocates new space according to bskl and reads old matrix with respect
// to permutation blockP
SparseGridMtxLL :: SparseGridMtxLL(SparseMatrixF &sm, long block_size, Ordering *block_order, Ordering *node_order, MathTracer *eMT, bool load_data) :
    SparseGridMtx(sm, block_size, block_order, node_order, eMT)
{
    this->BlockArith->prefered_decomposition = eLL_decomposition;
    IConectMatrix *bskl = block_order->cm;

    this->AlocateMemoryByPattern(bskl);
    ComputeBlocks();

    if ( sm.a != NULL && load_data ) {
        LoadMatrixNumbers(sm);
    }
}

double &SparseGridMtxLL :: ElementAt(int i, int j)
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

    double *column_nbj_bi_data = NULL;
    double *column_nbi_bj_data = NULL;

    long ni = ( nodeP == NULL ) ? i : nodeP [ i ];
    long si = ni % block_size;
    long bi = ni / block_size;
    long nbi = ( blockP == NULL ) ? bi : blockP [ bi ];

    SparseGridColumn *Column_nbi = Columns [ nbi ];

    aux_bj_idx = aux_bi_idx = -1;
    if ( nbi < nbj ) {
        aux_bi_idx = Column_nbj->FindExistingBlockIndex(nbi);
        column_nbj_bi_data = Columns_data + Column_nbj->column_start_idx + aux_bi_idx * block_storage;
    } else if ( nbi > nbj ) {
        aux_bj_idx = Column_nbi->FindExistingBlockIndex(nbj);
        column_nbi_bj_data = Columns_data + Column_nbi->column_start_idx + aux_bj_idx * block_storage;
    }

    if ( aux_bi_idx < 0 && aux_bj_idx < 0 ) {   //diagonal index
        if ( si < sj ) {
            return Columns_data [ nbi * block_storage + si + sj * block_size ];
        } else {
            return Columns_data [ nbi * block_storage + sj + si * block_size ];
        }
    } else if ( aux_bi_idx >= 0 ) { // write into column bj
        return column_nbj_bi_data [ block_size * sj + si ];
    } else {                                    // write into column bi
        return column_nbi_bj_data [ block_size * si + sj ];
    }
}

void SparseGridMtxLL :: LoadZeros()
{
    memset( this->Columns_data, 0, columns_data_length * sizeof( double ) );
    for ( long bj = 0; bj < n_blocks; bj++ ) {
        for ( long j = 0; j < block_size; j++ ) {
            Columns_data [ bj * block_storage + block_size * j + j ] = 1.0;
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
        Columns_data [ nbi * block_storage + si + si * block_size ] = 0.0;
    }
}

void SparseGridMtxLL :: LoadMatrixNumbers(SparseMatrixF &sm)
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
        double *column_nbi_bj_data = NULL;

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
                    column_nbi_bj_data = Columns_data + Column_nbi->column_start_idx + aux_bj_idx * block_storage;
                }

                old_nbi = nbi;
            }

            if ( aux_bi_idx < 0 && aux_bj_idx < 0 ) {           //diagonal index
                if ( si < sj ) {
                    Columns_data [ nbi * block_storage + si + sj * block_size ] = val;
                } else {
                    Columns_data [ nbi * block_storage + sj + si * block_size ] = val;
                }
            } else if ( aux_bi_idx >= 0 ) { // write into column bj
                column_nbj_bi_data [ block_size * sj + si ] = val;
            } else {                                            // write into column bi
                column_nbi_bj_data [ block_size * si + sj ] = val;
            }

            //this->SetValue(nbi,nbj,si,sj,(double)ColumnData[idx],ref aux_bi_idx,ref aux_bj_idx);

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

SparseGridMtxLL :: ~SparseGridMtxLL()
{
    if ( Columns_data ) {
        delete Columns_data;
        Columns_data = NULL;
    }
}

void SparseGridMtxLL :: AlocateMemoryByPattern(IConectMatrix *bskl)
{
    clock_t ts = clock();
    this->n_blocks = bskl->N();
    columns_data_length = n_blocks * block_storage;
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

    clock_t end = clock();
    double duration = ( double ) ( end - ts ) / CLOCKS_PER_SEC;
    sprintf(str, "%.3f s", duration);
    Writeln(str);
    //Write(duration.TotalSeconds.ToString("g4")+" s");
}

// y = Ax
void SparseGridMtxLL :: MultiplyByVector(const LargeVectorAttach &x, LargeVectorAttach &y)
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

void SparseGridMtxLL :: Factorize()
{
    //double* Atmp = new double[block_storage];
    long bi, bj, min_bi_J = 0;

    // This is a pattern of blocks in J-th column
    long *p_blockJ_pattern = new long [ n_blocks + 1 ];
    memset( p_blockJ_pattern, 0, ( n_blocks + 1 ) * sizeof( long ) );

    long newdotcnount = ( long ) ceil( ( double ) n_blocks / 24.0 );
    //long newdotcnount = n / 24;
    BlockArith->zero_pivots = 0;

    double *cd = this->Columns_data;
    //double* atmp = Atmp;
    {
        double *dd = cd;
        //double* idd = cd;
        long Djj = 0;
        eMT->act_block = 0;
        for ( bj = 0; bj < n_blocks; bj++ ) {
            SparseGridColumn &columnJ = * Columns [ bj ];
            long noJentries = columnJ.Entries;
            if ( noJentries > 0 ) {
                long *columnJentries = columnJ.IndexesUfa->Items;
                double *pAkj = cd + columnJ.column_start_idx;
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
                        double *pAki = cd + columnI.column_start_idx + ( noIentries - 1 ) * block_storage;
                        long *columnIentries = columnI.IndexesUfa->Items;
                        for ( long *columnIentry = columnIentries + noIentries - 1; columnIentry >= columnIentries; pAki -= block_storage ) {
                            long idx_K = p_blockJ_pattern [ * columnIentry-- ];
                            if ( idx_K == 0 ) {
                                continue;
                            }

                            BlockArith->SubATBproduct(pAij, pAki, pAkj + ~idx_K);
                        }
                    }

                    BlockArith->L_BlockSolve(dd + bi * block_storage, pAij);
                    pAij += block_storage;
                }

                // compute the diagonal and divide by it
                //DenseMatrix Djj = DiagonalBlocks[bj];// Diagonal
                for ( long idx = noJentries - 1; idx >= 0; idx-- ) {
                    bi = columnJentries [ idx ];
                    //Clear pattern
                    p_blockJ_pattern [ bi ] = 0;

                    int ij = columnJ.column_start_idx + idx * block_storage;
                    BlockArith->SubATBproduct(dd + Djj, cd + ij, cd + ij);
                }
            }

            // Factorize diagonal block
            BlockArith->LL_Decomposition(dd + Djj);

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
    //delete [] Atmp;
    ComputeBlocks();
}

void SparseGridMtxLL :: Factorize_Incomplete()
{
    //double* Atmp = new double[block_storage];
    long bi, bj, min_bi_J = 0;

    // This is a pattern of blocks in J-th column
    long *p_blockJ_pattern = new long [ n_blocks + 1 ];
    memset( p_blockJ_pattern, 0, ( n_blocks + 1 ) * sizeof( long ) );

    long newdotcnount = ( long ) ceil( ( double ) n_blocks / 24.0 );
    //long newdotcnount = n / 24;
    BlockArith->zero_pivots = 0;

    double *cd = this->Columns_data;
    //double* atmp = Atmp;
    {
        double *dd = cd;
        //double* idd = cd;
        long Djj = 0;
        eMT->act_block = 0;
        for ( bj = 0; bj < n_blocks; bj++ ) {
            SparseGridColumn &columnJ = * Columns [ bj ];
            long noJentries = columnJ.Entries;
            if ( noJentries > 0 ) {
                long *columnJentries = columnJ.IndexesUfa->Items;
                double *pAkj = cd + columnJ.column_start_idx;
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
                        double *pAki = cd + columnI.column_start_idx + ( noIentries - 1 ) * block_storage;
                        long *columnIentries = columnI.IndexesUfa->Items;
                        for ( long *columnIentry = columnIentries + noIentries - 1; columnIentry >= columnIentries; pAki -= block_storage ) {
                            long idx_K = p_blockJ_pattern [ * columnIentry-- ];
                            if ( idx_K == 0 ) {
                                continue;
                            }

                            BlockArith->SubATBproduct(pAij, pAki, pAkj + ~idx_K);
                        }
                    }

                    BlockArith->L_BlockSolve(dd + bi * block_storage, pAij);
                    pAij += block_storage;
                }

                // compute the diagonal and divide by it
                //DenseMatrix Djj = DiagonalBlocks[bj];// Diagonal
                for ( long idx = noJentries - 1; idx >= 0; idx-- ) {
                    bi = columnJentries [ idx ];
                    //Clear pattern
                    p_blockJ_pattern [ bi ] = 0;

                    int ij = columnJ.column_start_idx + idx * block_storage;
                    BlockArith->SubATBproduct(dd + Djj, cd + ij, cd + ij);
                }


                // Maximal number of entries in this column, which can be left
                //long maxJentries = (long)(columnJ.Entries*1.0);
                long maxJentries = columnJ.Entries - 1;
                if ( maxJentries < 2 ) {
                    maxJentries = 2;
                }

                if ( columnJ.Entries > maxJentries ) {
                    Array :: Copy(columnJ.IndexesUfa->Items + ( columnJ.Entries - maxJentries ),
                                  columnJ.IndexesUfa->Items,
                                  maxJentries);
                    columnJ.Entries = maxJentries;
                }
            }

            // Factorize diagonal block
            BlockArith->LL_Decomposition(dd + Djj);

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
    //delete [] Atmp;
    ComputeBlocks();
}

// This functions computes LDL' decomposition of first (nblocks-fixed_bn) columns
// The resulting submatrix will contain the schur complement A22 - A21 * A11^(-1) * A12
void SparseGridMtxLL :: SchurComplementFactorization(int fixed_blocks)
{
    //double* Atmp = new double[block_storage];
    long bi, bj, min_bi_J = 0;

    long blocks_to_factor = n_blocks - fixed_blocks;

    // This is a pattern of blocks in J-th column
    long *p_blockJ_pattern = new long [ n_blocks + 1 ];
    memset( p_blockJ_pattern, 0, ( n_blocks + 1 ) * sizeof( long ) );

    //long newdotcnount = (long)ceil((double)n_blocks / 24.0);
    //long newdotcnount = n / 24;
    BlockArith->zero_pivots = 0;

    double *cd = this->Columns_data;
    //double* atmp = Atmp;
    {
        double *dd = cd;
        //double* idd = cd;
        long Djj = 0;
        eMT->act_block = 0;
        for ( bj = 0; bj < n_blocks; bj++ ) {
            SparseGridColumn &columnJ = * Columns [ bj ];
            long noJentries = columnJ.Entries;
            if ( noJentries > 0 ) {
                long *columnJentries = columnJ.IndexesUfa->Items;
                double *pAkj = cd + columnJ.column_start_idx;
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
                        double *pAki = cd + columnI.column_start_idx + ( noIentries - 1 ) * block_storage;
                        long *columnIentries = columnI.IndexesUfa->Items;
                        for ( long *columnIentry = columnIentries + noIentries - 1; columnIentry >= columnIentries; pAki -= block_storage ) {
                            long idx_K = p_blockJ_pattern [ * columnIentry-- ];
                            if ( idx_K == 0 ) {
                                continue;
                            }

                            BlockArith->SubATBproduct(pAij, pAki, pAkj + ~idx_K);
                        }
                    }

                    BlockArith->L_BlockSolve(dd + bi * block_storage, pAij);
                    pAij += block_storage;
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
                    BlockArith->SubATBproduct(dd + Djj, cd + ij, cd + ij);
                }
            }

            // Factorize diagonal block
            if ( bj < blocks_to_factor ) {
                BlockArith->LL_Decomposition(dd + Djj);
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
    //delete [] Atmp;
    ComputeBlocks();
}

void SparseGridMtxLL :: Solve(double *b, double *x)
{
    if ( x != b ) {
        Array :: Copy(b, x, n);
    }

    SolveLL(x);
}

void SparseGridMtxLL :: SolveLV(const LargeVector &b, LargeVector &x)
{
    if ( & x != & b ) {
        x.Initialize(b, 0, n);
    }

    SolveLL( x.DataPtr() );
}

// x = L^(-1) * x
void SparseGridMtxLL :: ForwardSubstL(double *x, long fixed_blocks)
{
    if ( this->N() == 0 ) {
        return;
    }

    long blocks_to_factor = n_blocks - fixed_blocks;

    int bi;
    long *ord = this->block_order->order->Items;

    // forward substitution L z = f  --> z
    for ( bi = 0; bi < blocks_to_factor; bi++ ) {
        SparseGridColumn &rowI = * Columns [ bi ];
        int no = rowI.Entries;
        if ( no > 0 ) {
            double *dst = x + block_size * ord [ bi ];
            double *Aij = Columns_data + rowI.column_start_idx;

            long *idxs = rowI.IndexesUfa->Items;
            //r[i] -= Lji * r[j]
            for ( int idx = 0; idx < no; idx++, Aij += block_storage ) {
                BlockArith->SubMultTBlockByVector(Aij, x + block_size * ord [ * ( idxs++ ) ], dst);
            }
        }

        // Diagonal solve
        BlockArith->SubstSolveL(Columns_data + block_storage * bi, x + block_size * ord [ bi ]);
    }
}

void SparseGridMtxLL :: BackSubstLT(double *x, long fixed_blocks)
{
    if ( this->N() == 0 ) {
        return;
    }

    long blocks_to_factor = n_blocks - fixed_blocks;

    int bi;
    long *ord = this->block_order->order->Items;
    //double* Dii = dd + block_storage*(n_blocks-1);
    // back substitution U r = z'
    for ( bi = blocks_to_factor - 1; bi >= 0; bi-- ) {
        BlockArith->SubstSolveLT(Columns_data + block_storage * bi, x + block_size * ord [ bi ]);

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
void SparseGridMtxLL :: SolveLL(double *x, long fixed_blocks)
{
    ForwardSubstL(x, fixed_blocks);
    BackSubstLT(x, fixed_blocks);
}

void SparseGridMtxLL :: SolveA11(double *x, long fixed_blocks)
{
    if ( x == 0 || fixed_blocks == 0 ) {
        x = 0;
    }
}

void SparseGridMtxLL :: Sub_A21_A11inv(double *x, long fixed_blocks)
{
    if ( x == 0 || fixed_blocks == 0 ) {
        x = 0;
    }
}

void SparseGridMtxLL :: Sub_A11inv_A12(double *x, long fixed_blocks)
{
    if ( x == 0 || fixed_blocks == 0 ) {
        x = 0;
    }
}

void SparseGridMtxLL :: WriteCondensedMatrixA22(double *a, Ordering *mcn, IntArrayList *lncn)
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

        for ( i = 0; i <= j; i++ ) {
            long ni = mcn->perm->Items [ lncn->Items [ i ] ];
            long nbi = ni / blockSize;
            long nsi = ni % blockSize;
            long obi = block_order->perm->Items [ nbi ];

            double val = GetValue(obi, obj, nsi, nsj, aux_bi_idx, aux_bj_idx);
            a [ ( j ) + ( i ) * nbn ] = val;
            a [ ( i ) + ( j ) * nbn ] = val;
        }
    }
}


DSS_NAMESPASE_END
