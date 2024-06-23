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
/*          Compressed column sparse matrix (0-based)                    */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "compcol.h"
#include "floatarray.h"
#include "engngm.h"
#include "domain.h"
#include "element.h"
#include "sparsemtrxtype.h"
#include "activebc.h"
#include "classfactory.h"

#include <set>

namespace oofem {
REGISTER_SparseMtrx(CompCol, SMT_CompCol);


CompCol :: CompCol(int n) : SparseMtrx(n, n),
    val(0),
    rowind(0),
    colptr(n),
    base(0),
    nz(0)
{}


CompCol :: CompCol(const CompCol &S) : SparseMtrx(S.nRows, S.nColumns),
    val(S.val),
    rowind(S.rowind),
    colptr(S.colptr),
    base(S.base),
    nz(S.nz)
{}


CompCol &CompCol :: operator = ( const CompCol & C )
{
    nRows   = C.nRows;
    nColumns = C.nColumns;

    base   = C.base;
    nz     = C.nz;
    val    = C.val;
    rowind = C.rowind;
    colptr = C.colptr;
    this->version = C.version;

    return * this;
}


std::unique_ptr<SparseMtrx> CompCol :: clone() const
{
    return std::make_unique<CompCol>(*this);
}


void CompCol :: times(const FloatArray &x, FloatArray &answer) const
{
    if ( x.giveSize() != this->giveNumberOfColumns() ) {
        OOFEM_ERROR("incompatible dimensions");
    }

    answer.resize(this->giveNumberOfRows());
    answer.zero();

    for ( int j = 0; j < this->giveNumberOfColumns(); j++ ) {
        double rhs = x[j];
        for ( int t = colptr[j]; t < colptr[j + 1]; t++ ) {
            answer[ rowind[t] ] += val[t] * rhs;
        }
    }
}


void CompCol :: timesT(const FloatArray &x, FloatArray &answer) const
{
    if ( x.giveSize() != this->giveNumberOfRows() ) {
        OOFEM_ERROR("Error in CompCol -- incompatible dimensions");
    }

    answer.resize(this->giveNumberOfColumns());
    answer.zero();

    for ( int i = 0; i < this->giveNumberOfColumns(); i++ ) {
        double r = 0.0;
        for ( int t = colptr[i]; t < colptr[i + 1]; t++ ) {
            r += val[t] * x[ rowind[t] ];
        }

        answer[i] = r;
    }
}


void CompCol :: times(double x)
{
    val.times(x);

    this->version++;
}


int CompCol :: buildInternalStructure(EngngModel *eModel, int di, const UnknownNumberingScheme &s)
{
    IntArray loc;
    Domain *domain = eModel->giveDomain(di);
    int neq = eModel->giveNumberOfDomainEquations(di, s);
    // allocation map
    std :: vector< std :: set< int > > columns(neq);

    this->nz = 0;

    for ( auto &elem : domain->giveElements() ) {
        elem->giveLocationArray(loc, s);

        for ( int ii : loc ) {
            if ( ii > 0 ) {
                for ( int jj : loc ) {
                    if ( jj > 0 ) {
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
                    if ( ii ) {
                        for ( int jj : kcloc ) {
                            if ( jj ) {
                                columns [ jj - 1 ].insert(ii - 1);
                            }
                        }
                    }
                }
            }
        }
    }

    for ( int i = 0; i < neq; i++ ) {
        this->nz += columns [ i ].size();
    }

    rowind.resize(nz);
    colptr.resize(neq + 1);
    int indx = 0;

    for ( int j = 0; j < neq; j++ ) { // column loop
        colptr[j] = indx;
        for ( int row: columns [ j ] ) { // row loop
            rowind[indx++] = row;
        }
    }

    colptr[neq] = indx;

    // allocate value array
    val.resize(nz);
    val.zero();

    OOFEM_LOG_DEBUG("CompCol info: neq is %d, nwk is %d\n", neq, nz);

    nColumns = nRows = neq;

    this->version++;

    return true;
}


int CompCol :: assemble(const IntArray &loc, const FloatMatrix &mat)
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
            int cstart = colptr[jj - 1];
            int t = cstart;
            int last_ii = this->nRows + 1; // Ensures that t is set correctly the first time.
            for ( int i = 0; i < dim; i++ ) {
                int ii = loc[i];
                if ( ii ) {
                    // Some heuristics that speed up most cases ( benifits are large for incremental sub-blocks, e.g. locs = [123, 124, 125, 245, 246, 247] ):
                    if ( ii < last_ii )
                        t = cstart;
                    else if ( ii > last_ii )
                        t++;
                    for ( ; rowind[t] < ii - 1; t++ ) {
#  ifdef DEBUG
                        if ( t >= colptr[jj] )
                            OOFEM_ERROR("Couldn't find row %d in the sparse structure", ii);
#  endif
                    }
                    val[t] += mat(i, j);
                    last_ii = ii;
                }
            }
        }
    }

    // increment version
    this->version++;

    return 1;
}

int CompCol :: assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat)
{
    int dim1, dim2;

    dim1 = mat.giveNumberOfRows();
    dim2 = mat.giveNumberOfColumns();

    for ( int j = 0; j < dim2; j++ ) {
        int jj = cloc[j];
        if ( jj ) {
            int cstart = colptr[jj - 1];
            int t = cstart;
            int last_ii = this->nRows + 1; // Ensures that t is set correctly the first time.
            for ( int i = 0; i < dim1; i++ ) {
                int ii = rloc[i];
                if ( ii ) {
                    // Some heuristics that speed up most cases ( benifits are large for incremental sub-blocks, e.g. locs = [123, 124, 125, 245, 246, 247] ):
                    if ( ii < last_ii )
                        t = cstart;
                    else if ( ii > last_ii )
                        t++;
                    for ( ; rowind[t] < ii - 1; t++ ) {
#  ifdef DEBUG
                        if ( t >= colptr[jj] )
                            OOFEM_ERROR("Couldn't find row %d in the sparse structure", ii);
#  endif
                    }
                    val[t] += mat(i, j);
                    last_ii = ii;
                }
            }
        }
    }

    this->version++;

    return 1;
}

void CompCol :: zero()
{
    val.zero();

    this->version++;
}


void CompCol :: toFloatMatrix(FloatMatrix &answer) const
{ }

void CompCol :: printYourself() const
{ }


double &CompCol :: at(int i, int j)
{
    this->version++;

    for ( int t = colptr[j - 1]; t < colptr[j]; t++ ) {
        if ( rowind[t] == i - 1 ) {
            return val[t];
        }
    }

    OOFEM_ERROR("Array accessing exception -- (%d,%d) out of bounds", i, j);
}


double CompCol :: at(int i, int j) const
{
    for ( int t = colptr[j - 1]; t < colptr[j]; t++ ) {
        if ( rowind[t] == i - 1 ) {
            return val[t];
        }
    }

    if ( i <= this->giveNumberOfRows() && j <= this->giveNumberOfColumns() ) {
        return 0.0;
    } else {
        OOFEM_ERROR("Array accessing exception -- (%d,%d) out of bounds", i, j);
    }
    // return 0.; // return to suppress compiler warning message
}

double CompCol :: operator() (int i, int j)  const
{
    for ( int t = colptr[j]; t < colptr[j + 1]; t++ ) {
        if ( rowind[t] == i ) {
            return val[t];
        }
    }

    if ( i < this->giveNumberOfRows() && j < this->giveNumberOfColumns() ) {
        return 0.0;
    } else {
        OOFEM_ERROR("Array accessing exception -- (%d,%d) out of bounds", i, j);
    }
    // return 0.; // return to suppress compiler warning message
}

double &CompCol :: operator() (int i, int j)
{
    this->version++;

    for ( int t = colptr[j]; t < colptr[j + 1]; t++ ) {
        if ( rowind[t] == i ) {
            return val[t];
        }
    }

    OOFEM_ERROR("Array element (%d,%d) not in sparse structure -- cannot assign", i, j);
}

} // end namespace oofem
