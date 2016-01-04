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

// Class SymCompCol

// inspired by
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*             ********   ***                                 SparseLib++    */
/*          *******  **  ***       ***      ***               v. 1.5c        */
/*           *****      ***     ******** ********                            */
/*            *****    ***     ******** ********              R. Pozo        */
/*       **  *******  ***   **   ***      ***                 K. Remington   */
/*        ********   ********                                 A. Lumsdaine   */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                                                                           */
/*                                                                           */
/*                     SparseLib++ : Sparse Matrix Library                   */
/*                                                                           */
/*               National Institute of Standards and Technology              */
/*                        University of Notre Dame                           */
/*              Authors: R. Pozo, K. Remington, A. Lumsdaine                 */
/*                                                                           */
/*                                 NOTICE                                    */
/*                                                                           */
/* Permission to use, copy, modify, and distribute this software and         */
/* its documentation for any purpose and without fee is hereby granted       */
/* provided that the above notice appear in all copies and supporting        */
/* documentation.                                                            */
/*                                                                           */
/* Neither the Institutions (National Institute of Standards and Technology, */
/* University of Notre Dame) nor the Authors make any representations about  */
/* the suitability of this software for any purpose.  This software is       */
/* provided ``as is'' without expressed or implied warranty.                 */
/*                                                                           */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*          Compressed column symmetric sparse matrix (0-based)          */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "symcompcol.h"
#include "floatarray.h"
#include "engngm.h"
#include "domain.h"
#include "element.h"
#include "sparsemtrxtype.h"
#include "activebc.h"
#include "classfactory.h"

#include <set>

namespace oofem {
REGISTER_SparseMtrx(SymCompCol, SMT_SymCompCol);

SymCompCol :: SymCompCol(void) : CompCol()
{ }


SymCompCol :: SymCompCol(int n) : CompCol(n)
{ }


/*****************************/
/*  Copy constructor         */
/*****************************/

SymCompCol :: SymCompCol(const SymCompCol &S) : CompCol(S)
{ }


SparseMtrx *SymCompCol :: GiveCopy() const
{
    SymCompCol *result = new SymCompCol(*this);
    return result;
}


#define MAP(i, j) map [ ( j ) * neq - ( j ) * ( ( j ) + 1 ) / 2 + ( i ) ]

int SymCompCol :: buildInternalStructure(EngngModel *eModel, int di, const UnknownNumberingScheme &s)
{
    IntArray loc;
    Domain *domain = eModel->giveDomain(di);
    int neq = eModel->giveNumberOfDomainEquations(di, s);
    int indx;
    // allocation map
    std :: vector< std :: set< int > > columns(neq);

    this->nz_ = 0;

    for ( auto &elem : domain->giveElements() ) {
        elem->giveLocationArray(loc, s);

        for ( int ii : loc ) {
            if ( ii > 0 ) {
                for ( int jj : loc ) {
                    if ( jj > 0 && ii >= jj ) {
                        columns [ jj - 1 ].insert(ii - 1);
                    }
                }
            }
        }
    }

    // loop over active boundary conditions
    std :: vector< IntArray >r_locs;
    std :: vector< IntArray >c_locs;

    for ( auto &gbc : domain->giveBcs() ) {
        ActiveBoundaryCondition *bc = dynamic_cast< ActiveBoundaryCondition * >( gbc.get() );
        if ( bc != NULL ) {
            bc->giveLocationArrays(r_locs, c_locs, UnknownCharType, s, s);
            for ( std :: size_t k = 0; k < r_locs.size(); k++ ) {
                IntArray &krloc = r_locs [ k ];
                IntArray &kcloc = c_locs [ k ];
                for ( int ii : krloc ) {
                    if ( ii > 0 ) {
                        for ( int jj : kcloc ) {
                            if ( jj > 0 && ii >= jj ) {
                                columns [ jj - 1 ].insert(ii - 1);
                            }
                        }
                    }
                }
            }
        }
    }



    for ( auto &val : columns ) {
        this->nz_ += val.size();
    }

    rowind_.resize(nz_);
    colptr_.resize(neq + 1);
    indx = 0;

    for ( int j = 0; j < neq; j++ ) { // column loop
        colptr_(j) = indx;
        for ( int row: columns [ j ] ) { // row loop
            rowind_(indx++) = row;
        }
    }

    colptr_(neq) = indx;

    // allocate value array
    val_.resize(nz_);
    val_.zero();

    OOFEM_LOG_INFO("SymCompCol info: neq is %d, nwk is %d\n", neq, nz_);

    dim_ [ 0 ] = dim_ [ 1 ] = nColumns = nRows = neq;

    // increment version
    this->version++;

    return true;
}



void SymCompCol :: times(const FloatArray &x, FloatArray &answer) const
{
    int M = dim_ [ 0 ];
    int N = dim_ [ 1 ];

#if DEBUG
    if ( x.giveSize() != N ) {
        OOFEM_ERROR("incompatible dimensions");
    }
#endif

    answer.resize(M);
    answer.zero();

    for ( int j = 0; j < N; j++ ) {
        double rhs = x(j);
        double sum = 0.0;
        for ( int t = colptr_(j) + 1; t < colptr_(j + 1); t++ ) {
            answer( rowind_(t) ) += val_(t) * rhs; // column loop
            sum += val_(t) * x( rowind_(t) ); // row loop
        }

        answer(j) += sum;
        answer(j) += val_( colptr_(j) ) * rhs; // diagonal
    }
}

void SymCompCol :: times(double x)
{
    val_.times(x);

    // increment version
    this->version++;
}



int SymCompCol :: assemble(const IntArray &loc, const FloatMatrix &mat)
{
    int dim = mat.giveNumberOfRows();

#  ifdef DEBUG
    if ( dim != loc.giveSize() ) {
        OOFEM_ERROR("dimension of 'k' and 'loc' mismatch");
    }
#  endif

    for ( int j = 0; j < dim; j++ ) {
        int jj = loc[j];
        if ( jj ) {
            int cstart = colptr_[jj - 1];
            int t = cstart;
            int last_ii = this->nRows + 1; // Ensures that t is set correctly the first time.
            for ( int i = 0; i < dim; i++ ) {
                int ii = loc[i];
                if ( ii >= jj ) { // assemble only lower triangular part
                    // Some heuristics that speed up most cases ( benifits are large for incremental sub-blocks, e.g. locs = [123, 124, 125, 245, 246, 247] ):
                    if ( ii < last_ii )
                        t = cstart;
                    else if ( ii > last_ii )
                        t++;
                    for ( ; rowind_[t] < ii - 1; t++ ) {
#  ifdef DEBUG
                        if ( t >= colptr_[jj] )
                            OOFEM_ERROR("Couldn't find row %d in the sparse structure", ii);
#  endif
                    }
                    val_[t] += mat(i, j);
                    last_ii = ii;
                }
            }
        }
    }

    // increment version
    this->version++;

    return 1;
}

int SymCompCol :: assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat)
{
    int dim1, dim2;

    dim1 = mat.giveNumberOfRows();
    dim2 = mat.giveNumberOfColumns();

    for ( int j = 0; j < dim2; j++ ) {
        int jj = cloc[j];
        if ( jj ) {
            int cstart = colptr_[jj - 1];
            int t = cstart;
            int last_ii = this->nRows + 1; // Ensures that t is set correctly the first time.
            for ( int i = 0; i < dim1; i++ ) {
                int ii = rloc[i];
                if ( ii >= jj ) { // assemble only lower triangular part
                    // Some heuristics that speed up most cases ( benifits are large for incremental sub-blocks, e.g. locs = [123, 124, 125, 245, 246, 247] ):
                    if ( ii < last_ii )
                        t = cstart;
                    else if ( ii > last_ii )
                        t++;
                    for ( ; rowind_[t] < ii - 1; t++ ) {
#  ifdef DEBUG
                        if ( t >= colptr_[jj] )
                            OOFEM_ERROR("Couldn't find row %d in the sparse structure", ii);
#  endif
                    }
                    val_[t] += mat(i, j);
                    last_ii = ii;
                }
            }
        }
    }

    // increment version
    this->version++;

    return 1;
}

void SymCompCol :: zero()
{
    val_.zero();

    // increment version
    this->version++;
}

/*********************/
/*   Array access    */
/*********************/

double &SymCompCol :: at(int i, int j)
{
    int ii = i, jj = j;
    if ( ii < jj ) {
        ii = j;
        jj = i;
    }

    // increment version
    this->version++;

    for ( int t = colptr_(jj - 1); t < colptr_(jj); t++ ) {
        if ( rowind_(t) == ( ii - 1 ) ) {
            return val_(t);
        }
    }

    OOFEM_ERROR("Array accessing exception -- out of bounds");
    return val_(0); // return to suppress compiler warning message
}


double SymCompCol :: at(int i, int j) const
{
    int ii = i, jj = j;
    if ( ii < jj ) {
        ii = j;
        jj = i;
    }

    for ( int t = colptr_(jj - 1); t < colptr_(jj); t++ ) {
        if ( rowind_(t) == ( ii - 1 ) ) {
            return val_(t);
        }
    }

    if ( i <= dim_ [ 0 ] && j <= dim_ [ 1 ] ) {
        return 0.0;
    } else {
        OOFEM_ERROR("Array accessing exception -- index out of bounds (%d,%d)", i, j);
        return ( 0 ); // return to suppress compiler warning message
    }
}

double SymCompCol :: operator() (int i, int j)  const
{
    int ii = i, jj = j;
    if ( ii < jj ) {
        ii = j;
        jj = i;
    }

    for ( int t = colptr_(jj); t < colptr_(jj + 1); t++ ) {
        if ( rowind_(t) == ii ) {
            return val_(t);
        }
    }

    if ( i < dim_ [ 0 ] && j < dim_ [ 1 ] ) {
        return 0.0;
    } else {
        OOFEM_ERROR("Array accessing exception, index out of bounds (%d,%d)", i, j);
        return ( 0 ); // return to suppress compiler warning message
    }
}

double &SymCompCol :: operator() (int i, int j)
{
    int ii = i, jj = j;
    if ( ii < jj ) {
        ii = j;
        jj = i;
    }

    // increment version
    this->version++;

    for ( int t = colptr_(jj); t < colptr_(jj + 1); t++ ) {
        if ( rowind_(t) == ii ) {
            return val_(t);
        }
    }

    OOFEM_ERROR("Array element (%d,%d) not in sparse structure -- cannot assign", i, j);
    return val_(0); // return to suppress compiler warning message
}
} // end namespace oofem
