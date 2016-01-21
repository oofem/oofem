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

#include "SparseGridMtx.h"

DSS_NAMESPASE_BEGIN

// Allocates new space according to bskl and reads old matrix with respect
// to permutation blockP
SparseGridMtx :: SparseGridMtx(SparseMatrixF &sm, long block_size, Ordering *block_order, MathTracer *eMT)
{
    ///tmp_vector_BS = NULL;
    //tmp_vector_BS_nodes = NULL;

    //long* blockO = order->order->Items;
    long *blockP = block_order->perm->Items;
    nonzeros = 0;

    this->eMT = eMT;
    this->block_size = block_size;
    this->n_blocks = ( long ) sm.neq / block_size + ( ( ( sm.neq % block_size ) == 0 ) ? 0 : 1 );
    this->n = ( blockP == NULL ) ? ( long ) sm.neq : n_blocks * block_size;
    this->noDummyDOFs = n_blocks * block_size - ( long ) sm.neq;
    this->block_storage = block_size * block_size;

    this->BlockArith = DenseMatrixArithmetics :: NewArithmetics(block_size);
    this->BlockArith->eMT = eMT;
    this->block_order = block_order;
    this->node_order = NULL;
    this->Columns = new SparseGridColumn * [ n_blocks ];
    this->no_multiplications = 0;
}

// Allocates new space according to bskl and reads old matrix with respect
// to permutation blockP
SparseGridMtx :: SparseGridMtx(SparseMatrixF &sm, long block_size, Ordering *block_order, Ordering *node_order, MathTracer *eMT)
{
    //tmp_vector_BS = NULL;
    //tmp_vector_BS_nodes = NULL;

    //long* blockO = order->order->Items;
    long *blockP = block_order->perm->Items;
    nonzeros = 0;

    this->eMT = eMT;
    this->block_size = block_size;
    if ( node_order ) {
        this->n_blocks = ( long ) ( node_order->order->Count ) / block_size + ( ( ( sm.neq % block_size ) == 0 ) ? 0 : 1 );
    } else {
        this->n_blocks = ( long ) sm.neq / block_size + ( ( ( sm.neq % block_size ) == 0 ) ? 0 : 1 );
    }

    this->n = ( blockP == NULL ) ? ( long ) sm.neq : n_blocks * block_size;
    this->noDummyDOFs = n_blocks * block_size - ( long ) sm.neq;
    this->block_storage = block_size * block_size;

    this->BlockArith = DenseMatrixArithmetics :: NewArithmetics(block_size);
    this->BlockArith->eMT = eMT;
    this->block_order = block_order;
    this->node_order = node_order;
    this->Columns = new SparseGridColumn * [ n_blocks ];
    this->no_multiplications = 0;
}

SparseGridMtx :: ~SparseGridMtx()
{
    if ( BlockArith ) {
        delete BlockArith;
        BlockArith = NULL;
    }

    if ( Columns ) {
        for ( long i = 0; i < n_blocks; i++ ) {
            if ( Columns [ i ] ) {
                delete Columns [ i ];
                Columns [ i ] = NULL;
            }
        }

        delete [] Columns;
        Columns = NULL;
    }

    // This array is not owned by the matrix
    //if (node_order){delete node_order;node_order = NULL;}
    if ( block_order ) {
        delete block_order;
        block_order = NULL;
    }

    //if (tmp_vector_BS){delete tmp_vector_BS;tmp_vector_BS = NULL;}
    //if (tmp_vector_BS_nodes){delete tmp_vector_BS_nodes;tmp_vector_BS_nodes = NULL;}
}

void SparseGridMtx :: ComputeBlocks()
{
    this->blocks = n_blocks;
    for ( long bi = 0; bi < n_blocks; bi++ ) {
        blocks += Columns [ bi ]->Entries;
    }
}

double SparseGridMtx :: GetWaste()
{
    return 1.0 - ( double ) nonzeros / ( block_storage * blocks );
}

void SparseGridMtx :: WriteStatistics(long no_init_blocks, long no_nonzeros)
{
    char str [ 512 ];
    long blockSize = BlockSize();
    double Waste = 100.0 * ( 1.0 - ( double ) no_nonzeros / ( no_init_blocks * blockSize * blockSize ) );

    sprintf(str, " blocks size         : %ldx%ld", blockSize, blockSize);
    eMT->Writeln(str);
    sprintf( str, " number of blocks    : %ld", Blocks() );
    eMT->Writeln(str);
    sprintf( str, " number of nonzeros  : %ld", Nonzeros() );
    eMT->Writeln(str);
    sprintf(str, " allocation waste    : %.2f", Waste);
    eMT->Write(str);
    eMT->Writeln("%% ");
}

void SparseGridMtx :: MultiplyByVector(const LargeVectorAttach &x, LargeVectorAttach &y)
{
    printf("SparseGridMtx::MultiplyByVector not implemented, file %s, line %d\n", __FILE__, __LINE__);
    exit(1);
}

void SparseGridMtx :: times(double x){
    printf("SparseGridMtx::times not implemented, file %s, line %d\n", __FILE__, __LINE__);
    exit(1);
}

DSS_NAMESPASE_END
