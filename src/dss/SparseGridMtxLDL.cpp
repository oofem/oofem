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

#include "SparseGridMtxLDL.h"

DSS_NAMESPASE_BEGIN

// Allocates new space according to bskl and reads old matrix with respect
// to permutation blockP
SparseGridMtxLDL :: SparseGridMtxLDL(SparseMatrixF &sm, long block_size, Ordering *block_order, MathTracer *eMT, bool load_data) :
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
SparseGridMtxLDL :: SparseGridMtxLDL(SparseMatrixF &sm, long block_size, Ordering *block_order, Ordering *node_order, MathTracer *eMT, bool load_data) :
    SparseGridMtx(sm, block_size, block_order, node_order, eMT)
{
    IConectMatrix *bskl = block_order->cm;
    this->AlocateMemoryByPattern(bskl);
    ComputeBlocks();

    if ( sm.a != NULL && load_data ) {
        LoadMatrixNumbers(sm);
    }
}

void SparseGridMtxLDL :: LoadZeros()
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

void SparseGridMtxLDL :: LoadMatrixNumbers(SparseMatrixF &sm)
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

SparseGridMtxLDL :: ~SparseGridMtxLDL()
{
    if ( Columns_data ) {
        delete [] Columns_data;
        Columns_data = NULL;
    }
}


void SparseGridMtxLDL :: AlocateMemoryByPattern(IConectMatrix *bskl)
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
    if ( Columns_data == NULL ) {
        Writeln("Out of memory!");
        return;
    }

    clock_t end = clock();
    double duration = ( double ) ( end - ts ) / CLOCKS_PER_SEC;
    sprintf(str, "%.3f s", duration);
    Writeln(str);
    //Write(duration.TotalSeconds.ToString("g4")+" s");
}

double &SparseGridMtxLDL :: ElementAt(int i, int j)
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

/*
 * // Procedures using sealed data
 * void SparseGridMtxLDL::SetValue(long bi,long bj,long si, long sj,double val,long& aux_bi_idx,long& aux_bj_idx)
 * {
 *      if (bi==bj)
 *              this->Columns_data[bi*block_storage+si+sj*block_size] = val;
 *      else
 *              if (bj>bi)
 *                      Columns[bj]->SetValue(block_size,bi,si,sj,val,Columns_data,aux_bj_idx);
 *              else
 *                      Columns[bi]->SetValue(block_size,bj,sj,si,val,Columns_data,aux_bi_idx);
 * }
 *
 * // Procedures using sealed data
 * double SparseGridMtxLDL::GetValue(long bi,long bj,long si, long sj,long& aux_bi_idx,long& aux_bj_idx)
 * {
 *      if (bi==bj)
 *              return this->Columns_data[bi*block_storage+si+sj*block_size];
 *      else
 *              if (bj>bi)
 *                      return Columns[bj]->GetValue(block_size,bi,si,sj,Columns_data,aux_bj_idx);
 *              else
 *                      return Columns[bi]->GetValue(block_size,bj,sj,si,Columns_data,aux_bi_idx);
 * }
 *
 * void SparseGridMtxLDL::AddValue(long bi,long bj,long si, long sj,double val)
 * {
 *      if (bi==bj)
 *              this->Columns_data[bi*block_storage+si+sj*block_size] += val;
 *      else
 *              if (bj>bi)
 *                      Columns[bj]->AddValue(block_size,bi,si,sj,val,Columns_data);
 *              else
 *                      Columns[bi]->AddValue(block_size,bj,sj,si,val,Columns_data);
 * } */

// y = Ax
void SparseGridMtxLDL :: MultiplyByVector(const LargeVectorAttach &x, LargeVectorAttach &y)
{
    y.Initialize();

    double *py = y.DataPtr();
    double *px = x.DataPtr();
    for ( long bj = 0; bj < n_blocks; bj++ ) {
        SparseGridColumn &columnJ = *Columns [ bj ];
        //Columns[bj].MultBlockByVectorSym(BlockArith,x,ref y,block_size,bj);
        long cnt = columnJ.Entries;
        for ( long idx = 0; idx < cnt; idx++ ) {
            long bi = columnJ.IndexesUfa->Items [ idx ];
            //DenseMatrix block = this->Blocksfa[idx];
            long block = columnJ.column_start_idx + idx * block_storage;
            BlockArith->AddMultBlockByVectorSym(Columns_data + block, px, py, block_size * bi, bj * block_size);
        }

        //BlockArith.MultDiagonalBlockByVector(DiagonalBlocks[bj],x, ref y,block_size,bj*block_size);
        BlockArith->MultDiagonalBlockByVector(Columns_data + bj * block_storage, px + bj * block_size, py + bj * block_size);
    }
} //MultiplyByVector

void SparseGridMtxLDL :: times(double x) {
//     long neq = n_blocks * block_size - noDummyDOFs;
//     for ( long i = 0; i < neq; i++ ) {
//         for ( long j = 0; j < neq; j++ ) {
//             printf("%ld %ld %lf\n", i,j, ElementAt(i, j));
//         }
//     }

//     printf("n_blocks %ld block_size %ld columns_data_length %ld \n", n_blocks, block_size, columns_data_length );

    for ( long i = 0; i < columns_data_length; i++ ) {
        Columns_data[i] *= x;
    }

}

void SparseGridMtxLDL :: Factorize()
{
    double *Atmp = new double [ block_storage ];
    long bi, bj, min_bi_J = 0;

    // This is a pattern of blocks in J-th column
    long *p_blockJ_pattern = new long [ n_blocks + 1 ];
    memset( p_blockJ_pattern, 0, ( n_blocks + 1 ) * sizeof( long ) );

    long newdotcnount = ( long ) ceil( ( double ) n_blocks / 24.0 );
    //long newdotcnount = n / 24;
    BlockArith->zero_pivots = 0;

    no_multiplications = 0;

    double *cd = this->Columns_data;
    double *atmp = Atmp;
    {
        double *dd = cd;
        double *idd = cd;
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
                for ( long idx_J = 1; idx_J < noJentries; idx_J++ ) {
                    pAij += block_storage;
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
                            //no_multiplications++;
                        }
                    }
                }

                // compute the diagonal and divide by it
                //DenseMatrix Djj = DiagonalBlocks[bj];// Diagonal
                for ( long idx = noJentries - 1; idx >= 0; idx-- ) {
                    bi = columnJentries [ idx ];
                    //Clear pattern
                    p_blockJ_pattern [ bi ] = 0;

                    //DenseMatrix Aij = columnJ.Blocksfa[idx];
                    long Aij = columnJ.column_start_idx + idx * block_storage;

                    //Aij.CopyTo(ref Atmp,block_size);
                    Array :: Copy(this->Columns_data, Aij, Atmp, 0, block_storage);

                    //L12 = D1(-1) * A12
                    BlockArith->SubstSolveBlock(idd + bi * block_storage, cd + Aij);

                    // Atmp = D1 * L12
                    // D2 = A22 - L12(T) * D1 * L12
                    // D2 = A22 - L12(T) * Atmp
                    BlockArith->SubATBproduct(dd + Djj, atmp, cd + Aij);
                }
            }

            // Factorize diagonal block
            BlockArith->FactorizeBlock(dd + Djj);

            if ( ( bj % newdotcnount ) == 0 ) {
           //     Write("."); //JB
            }

            Djj += block_storage;
            eMT->act_block += block_size;
            if ( eMT->break_flag ) {
                break;
            }
        }
    }
    delete [] p_blockJ_pattern;
    delete [] Atmp;
    ComputeBlocks();
}

void SparseGridMtxLDL :: Factorize_Incomplete()
{
    double *Atmp = new double [ block_storage ];
    long bi, bj, min_bi_J = 0;

    // This is a pattern of blocks in J-th column
    long *p_blockJ_pattern = new long [ n_blocks + 1 ];
    memset( p_blockJ_pattern, 0, ( n_blocks + 1 ) * sizeof( long ) );

    long newdotcnount = ( long ) ceil( ( double ) n_blocks / 24.0 );
    //long newdotcnount = n / 24;
    BlockArith->zero_pivots = 0;

    double *cd = this->Columns_data;
    double *atmp = Atmp;
    {
        double *dd = cd;
        double *idd = cd;
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
                for ( long idx_J = 1; idx_J < noJentries; idx_J++ ) {
                    pAij += block_storage;
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
                }

                // compute the diagonal and divide by it
                //DenseMatrix Djj = DiagonalBlocks[bj];// Diagonal
                for ( long idx = noJentries - 1; idx >= 0; idx-- ) {
                    bi = columnJentries [ idx ];
                    //Clear pattern
                    p_blockJ_pattern [ bi ] = 0;

                    //DenseMatrix Aij = columnJ.Blocksfa[idx];
                    long Aij = columnJ.column_start_idx + idx * block_storage;

                    //Aij.CopyTo(ref Atmp,block_size);
                    Array :: Copy(this->Columns_data, Aij, Atmp, 0, block_storage);

                    //L12 = D1(-1) * A12
                    BlockArith->SubstSolveBlock(idd + bi * block_storage, cd + Aij);

                    // Atmp = D1 * L12
                    // D2 = A22 - L12(T) * D1 * L12
                    // D2 = A22 - L12(T) * Atmp
                    BlockArith->SubATBproduct(dd + Djj, atmp, cd + Aij);
                }

                // Maximal number of entries in this column, which can be left
                //long maxJentries = (long)(columnJ.Entries*1.0);
                long maxJentries = columnJ.Entries - 1;
                if ( maxJentries < 2 ) {
                    maxJentries = 2;
                }

                if ( columnJ.Entries > maxJentries ) {
                    columnJ.Entries = maxJentries;
                }

                // Do the Incomplete magic. Scan this column and remove some entries.
                /*for (long idx=noJentries-1; idx>=0; idx--)
                 * {
                 *      bi = columnJentries[idx];
                 *      //Clear pattern
                 *      p_blockJ_pattern[bi] = 0;
                 *
                 *      //DenseMatrix Aij = columnJ.Blocksfa[idx];
                 *      long Aij = columnJ.column_start_idx + idx*block_storage;
                 *
                 *      //Aij.CopyTo(ref Atmp,block_size);
                 *      Array::Copy(this->Columns_data,Aij,Atmp,0,block_storage);
                 *
                 *      //L12 = D1(-1) * A12
                 *      BlockArith->SubstSolveBlock(idd+bi*block_storage,cd+Aij);
                 *
                 *      // Atmp = D1 * L12
                 *      // D2 = A22 - L12(T) * D1 * L12
                 *      // D2 = A22 - L12(T) * Atmp
                 *      BlockArith->SubATBproduct(dd+Djj,atmp,cd+Aij);
                 * }*/
            }

            // Factorize diagonal block
            BlockArith->FactorizeBlock(dd + Djj);

            if ( ( bj % newdotcnount ) == 0 ) {
                // Write("."); // JB
            }

            Djj += block_storage;
            eMT->act_block += block_size;
            if ( eMT->break_flag ) {
                break;
            }
        }
    }
    delete [] p_blockJ_pattern;
    delete [] Atmp;
    ComputeBlocks();
}


// This functions computes LDL' decomposition of first (nblocks-fixed_bn) columns
// The resulting submatrix will contain the schur complement A22 - A21 * A11^(-1) * A12
void SparseGridMtxLDL :: SchurComplementFactorization(int fixed_blocks)
{
    double *Atmp = new double [ block_storage ];
    long bi, bj, min_bi_J = 0;

    long blocks_to_factor = n_blocks - fixed_blocks;

    // This is a pattern of blocks in J-th column
    long *p_blockJ_pattern = new long [ n_blocks + 1 ];
    memset( p_blockJ_pattern, 0, ( n_blocks + 1 ) * sizeof( long ) );

    //long newdotcnount = (long)ceil((double)n_blocks / 24.0);
    //long newdotcnount = n / 24;
    BlockArith->zero_pivots = 0;

    double *cd = this->Columns_data;
    double *atmp = Atmp;
    {
        double *dd = cd;
        double *idd = cd;
        long Djj = 0;
        for ( bj = 0; bj < n_blocks; bj++ ) {
            SparseGridColumn &columnJ = * Columns [ bj ];
            long noJentries = columnJ.Entries;


            long *columnJentries = columnJ.IndexesUfa->Items;
            double *pAkj = cd + columnJ.column_start_idx;
            double *pAij = pAkj;

            //columnJ.DrawColumnPattern(p_blockJ_pattern,ref min_bi_J,bj);
            min_bi_J = bj;

            if ( bj < blocks_to_factor ) {              // This is normal scalar product
                for ( long i = noJentries - 1; i >= 0; i-- ) {
                    p_blockJ_pattern [ min_bi_J = columnJentries [ i ] ] = ~( i * block_storage );
                }
            } else {                                                  // Here we would read from the A22 matrix, therefore the scalar product is shorter (bi<blocks_to_factor)
                for ( long i = noJentries - 1; i >= 0; i-- ) {
                    if ( columnJentries [ i ] < blocks_to_factor ) {
                        p_blockJ_pattern [ min_bi_J = columnJentries [ i ] ] = ~( i * block_storage );
                    }
                }
            }

            // eliminate above diagonal
            for ( long idx_J = 1; idx_J < noJentries; idx_J++ ) {
                pAij += block_storage;
                bi = columnJentries [ idx_J ];

                SparseGridColumn &columnI = * Columns [ bi ];
                long noIentries = columnI.Entries;

                if ( noIentries > 0 ) { // This is normal scalar product
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
            }

            // compute the diagonal and divide by it
            //DenseMatrix Djj = DiagonalBlocks[bj];// Diagonal
            for ( long idx = 0; idx < noJentries; idx++ ) {
                bi = columnJentries [ idx ];
                //Clear pattern
                p_blockJ_pattern [ bi ] = 0;

                if ( bi >= blocks_to_factor ) {
                    break;
                }

                //DenseMatrix Aij = columnJ.Blocksfa[idx];
                long Aij = columnJ.column_start_idx + idx * block_storage;

                //Aij.CopyTo(ref Atmp,block_size);
                Array :: Copy(this->Columns_data, Aij, Atmp, 0, block_storage);

                //L12 = D1(-1) * A12
                BlockArith->SubstSolveBlock(idd + bi * block_storage, cd + Aij);

                // Atmp = D1 * L12
                // D2 = A22 - L12(T) * D1 * L12
                // D2 = A22 - L12(T) * Atmp
                BlockArith->SubATBproduct(dd + Djj, atmp, cd + Aij);
            }

            if ( bj < blocks_to_factor ) {
                BlockArith->FactorizeBlock(dd + Djj);
            }

            //if ((bj % newdotcnount) == 0)
            //Write(".");

            Djj += block_storage;
        }
    }

    delete [] p_blockJ_pattern;
    delete [] Atmp;
    ComputeBlocks();
}



void SparseGridMtxLDL :: Solve(double *b, double *x)
{
    if ( x != b ) {
        Array :: Copy(b, x, n);
    }

    SolveLDL(x);
}

void SparseGridMtxLDL :: SolveLV(const LargeVector &f, LargeVector &r)
{
    SolveLDL_node_perm(f, r);
}

void SparseGridMtxLDL :: SolveLDL_node_perm(const LargeVector &b, LargeVector &x)
{
    if ( node_order != NULL ) {
        b.GetPermuted_1Vector(tmp_vector_BS_nodes, node_order->perm->Items);
        SolveLDL_block_perm(* tmp_vector_BS_nodes, * tmp_vector_BS_nodes);
        tmp_vector_BS_nodes->GetPermutedVector(& x, node_order->perm->Items);
    } else {
        SolveLDL_block_perm(b, x);
    }
}

void SparseGridMtxLDL :: SolveLDL_block_perm(const LargeVector &b, LargeVector &x)
{
    Solve( x.DataPtr(), b.DataPtr() );
}

// x = A^(-1) * b
void SparseGridMtxLDL :: SolveLDL(double *x, long fixed_blocks)
{
    ForwardSubstL(x, fixed_blocks);
    SolveD(x, fixed_blocks);
    BackSubstLT(x, fixed_blocks);
}

/// <summary>  y -= L12^T * x  </summary>
void SparseGridMtxLDL :: SubMultL12T(double *px, double *py, long fixed_blocks)
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

        double *dst = py + block_size * ord [ bj ];
        double *Aij = cd +  columnJ.column_start_idx;

        long *idxs = columnJ.IndexesUfa->Items;
        //y[i] -= Lji^T * x[j]
        for ( long idx = 0; idx < no && * idxs < schur_shift; idx++, Aij += block_storage ) {
            BlockArith->SubMultTBlockByVector(Aij, px + block_size * ord [ * ( idxs++ ) ], dst);
        }
    }
}

/// <summary>  y -= L12 * x  </summary>
void SparseGridMtxLDL :: SubMultL12(double *px, double *py, long fixed_blocks)
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
        //y[i] -= Lji * x[j]
        for ( long idx = 0; idx < no && * idxs < schur_shift; idx++, Aij += block_storage ) {
            BlockArith->SubMultBlockByVector(Aij, src, py + block_size * ord [ * ( idxs++ ) ]);
        }
    }
}

// x = L^(-1) * b
void SparseGridMtxLDL :: ForwardSubstL(double *x, long fixed_blocks)
{
    if ( this->N() == 0 ) {
        return;
    }

    long blocks_to_factor = n_blocks - fixed_blocks;
    long *ord = block_order->order->Items;
    // forward substitution L z = f  --> z
    for ( long bi = 1; bi < blocks_to_factor; bi++ ) {
        SparseGridColumn &rowI = * Columns [ bi ];
        long no = rowI.Entries;
        if ( no == 0 ) {
            continue;
        }

        double *dst = x + block_size * ord [ bi ];
        double *Aij = Columns_data + rowI.column_start_idx;

        long *idxs = rowI.IndexesUfa->Items;
        //r[i] -= Lji^T * r[j]
        for ( long idx = 0; idx < no; idx++, Aij += block_storage ) {
            BlockArith->SubMultTBlockByVector(Aij, x + block_size * ord [ * ( idxs++ ) ], dst);
        }
    }
}

void SparseGridMtxLDL :: SolveD(double *x, long fixed_blocks)
{
    if ( this->N() == 0 ) {
        return;
    }

    long *ord = block_order->order->Items;
    // Diagonal solution D z' = z
    double *Dii = Columns_data;
    long blocks_to_factor = n_blocks - fixed_blocks;
    for ( long bi = 0; bi < blocks_to_factor; bi++, Dii += block_storage ) {
        BlockArith->SubstSolve(Dii, x + block_size * ord [ bi ]);
    }
}

void SparseGridMtxLDL :: BackSubstLT(double *x, long fixed_blocks)
{
    if ( this->N() == 0 ) {
        return;
    }

    long *ord = block_order->order->Items;
    double *cd = this->Columns_data;
    // back substitution L^T r = z'
    for ( long bi = n_blocks - fixed_blocks - 1; bi >= 0; bi-- ) {
        SparseGridColumn &columnI = * Columns [ bi ];
        long no = columnI.Entries;
        if ( no == 0 ) {
            continue;
        }

        double *src = x + block_size * ord [ bi ];
        double *Aij = cd + columnI.column_start_idx;
        long *idxs = columnI.IndexesUfa->Items;
        //r[j] -= Lij * r[i]
        for ( long idx = 0; idx < no; idx++, Aij += block_storage ) {
            BlockArith->SubMultBlockByVector(Aij, src, x + block_size * ord [ * ( idxs++ ) ]);
        }
    }
}

void SparseGridMtxLDL :: SolveA11(double *x, long fixed_blocks)
{
    SolveLDL(x, fixed_blocks);
}

void SparseGridMtxLDL :: Sub_A21_A11inv(double *x, long fixed_blocks)
{
    ForwardSubstL(x, fixed_blocks);
    SubMultL12T(x, x, fixed_blocks);
    SolveD(x, fixed_blocks);
}

void SparseGridMtxLDL :: Sub_A11inv_A12(double *x, long fixed_blocks)
{
    SubMultL12(x, x, fixed_blocks);
    BackSubstLT(x, fixed_blocks);
}

void SparseGridMtxLDL :: WriteCondensedMatrixA22(double *a, Ordering *mcn, IntArrayList *lncn)
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
