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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

// inspired by SL++

#include "dyncompcol.h"
#include "floatarray.h"
#include "engngm.h"
#include "domain.h"
#include "mathfem.h"
#include "element.h"
#include "sparsemtrxtype.h"
#include "activebc.h"
#include "classfactory.h"
#ifdef __MPM_MODULE
#include "../mpm/integral.h"
#endif

namespace oofem {
REGISTER_SparseMtrx(DynCompCol, SMT_DynCompCol);


DynCompCol :: DynCompCol(int n) : SparseMtrx(n, n),
    base(0)
{}


DynCompCol :: DynCompCol(const DynCompCol &S) : SparseMtrx(S.nRows, S.nColumns),
    columns(S.columns),
    rowind(S.rowind),
    base(S.base)
{
    this->version = S.version;
}


DynCompCol &DynCompCol :: operator = ( const DynCompCol & C )
{
    columns = C.columns;
    rowind = C.rowind;
    base = C.base;

    nRows   = C.nRows;
    nColumns = C.nColumns;
    version = C.version;

    return * this;
}

std::unique_ptr<SparseMtrx> DynCompCol :: clone() const
{
    return std::make_unique<DynCompCol>(*this);
}

void DynCompCol :: times(const FloatArray &x, FloatArray &answer) const
{
    if ( x.giveSize() != nColumns ) {
        OOFEM_ERROR("incompatible dimensions");
    }

    answer.resize(nRows);
    answer.zero();

    for ( int j = 0; j < nColumns; j++ ) {
        double rhs = x[j];
        for ( int t = 1; t <= columns[ j ].giveSize(); t++ ) {
            answer[ rowind[ j ].at(t) ] += columns[ j ].at(t) * rhs;
        }
    }
}

void DynCompCol :: times(double x)
{
    for ( auto &column : columns ) {
        column.times(x);
    }

    this->version++;
}

int DynCompCol :: buildInternalStructure(EngngModel *eModel, int di, const UnknownNumberingScheme &s)
{
    int neq = eModel->giveNumberOfDomainEquations(di, s);

    IntArray loc;
    Domain *domain = eModel->giveDomain(di);

    nColumns = nRows = neq;

    rowind.clear();
    columns.clear();
    this->growTo(neq);

    for ( auto &elem : domain->giveElements() ) {
        elem->giveLocationArray(loc, s);

        for ( int ii : loc ) {
            if ( ii > 0 ) {
                for ( int jj : loc ) {
                    if ( jj > 0 ) {
                        this->insertRowInColumn(ii - 1, jj - 1);
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
        if ( bc ) {
            bc->giveLocationArrays(r_locs, c_locs, UnknownCharType, s, s);
            for ( std :: size_t k = 0; k < r_locs.size(); k++ ) {
                IntArray &krloc = r_locs [ k ];
                IntArray &kcloc = c_locs [ k ];
                for ( int ii : krloc ) {
                    if ( ii > 0 ) {
                        for ( int jj : kcloc ) {
                            if ( jj > 0 ) {
                                this->insertRowInColumn(jj - 1, ii - 1);
                            }
                        }
                    }
                }
            }
        }
    }

#ifdef __MPM_MODULE
    IntArray locr, locc;
    // loop over integrals 
    for (auto &in: eModel->giveIntegralList()) {
        // loop over integral domain
        for (auto &elem: in->set->giveElementList()) {
            // get code numbers for integral.term on element
            in->getElementTermCodeNumbers (locr, locc, domain->giveElement(elem), *in->term, s) ;
            for ( int ii : locr ) {
                if ( ii > 0 ) {
                    for ( int jj : locc ) {
                        if ( jj > 0 ) {
                            this->insertRowInColumn(jj - 1, ii - 1);
                        }
                    }
                }
            }
        }
    }
#endif

    int nz_ = 0;
    for ( int j = 0; j < neq; j++ ) {
        nz_ += this->rowind[ j ].giveSize();
    }

    OOFEM_LOG_DEBUG("DynCompCol info: neq is %d, nelem is %d\n", neq, nz_);

    this->version++;

    return true;
}


int DynCompCol :: assemble(const IntArray &loc, const FloatMatrix &mat)
{
    int dim = mat.giveNumberOfRows();

#  ifdef DEBUG
    if ( dim != loc.giveSize() ) {
        OOFEM_ERROR("dimension of 'k' and 'loc' mismatch");
    }
#  endif

    for ( int j = 1; j <= dim; j++ ) {
        int jj = loc.at(j);
        if ( jj ) {
            for ( int i = 1; i <= dim; i++ ) {
                int ii = loc.at(i);
                if ( ii ) {
                    this->at(ii, jj) += mat.at(i, j);
                }
            }
        }
    }

    this->version++;

    return 1;
}

int DynCompCol :: assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat)
{
#if 0
     // slow and safe
     this->checkSizeTowards(rloc, cloc);
     int dim1 = mat.giveNumberOfRows();
     int dim2 = mat.giveNumberOfColumns();
     for (int i = 1 ; i <= dim1; i++) {
        int ii = rloc.at(i);
        if (ii)
            for (j=1 ; j<= dim2; j++) {
                int jj = cloc.at(j);
                if (jj) this->at(ii,jj) += mat.at(i,j);
            }
     }
     return 1;
#endif

    // optimized low-end implementation
    IntArray rowsToAdd( rloc.giveSize() );
    int rsize = rloc.giveSize();
    int csize = cloc.giveSize();

#if 0
     // adjust the size of receiver
     int maxid = 0;
     for (int i = 0; i < rsize; i++) maxid = max(maxid, rloc[i]);
     for (int i = 0; i < csize; i++) maxid = max(maxid, cloc[i]);
     this->growTo(maxid);
#endif

    for ( int i = 0; i < csize; i++ ) {
        int ii = cloc[i];
        if ( ii ) {
            int ii1 = ii - 1;
            for ( int j = 0; j < rsize; j++ ) {
                int jj = rloc[j];
                if ( jj ) {
                    int jj1 = jj - 1;

                    int rowindx = this->insertRowInColumn(ii1, jj1);
                    columns[ ii1 ].at(rowindx) += mat(j, i);
                }
            }
        }
    }

    this->version++;

    return 1;
}


void DynCompCol :: zero()
{
    for ( auto &column : columns ) {
        column.zero();
    }

    this->version++;
}


void DynCompCol :: printStatistics() const
{
    int nz_ = 0;
    for ( auto &row : rowind ) {
        nz_ += row.giveSize();
    }

    OOFEM_LOG_DEBUG("DynCompCol info: neq is %d, nelem is %d\n", nColumns, nz_);
}


double &DynCompCol :: at(int i, int j)
{
    this->version++;

    int rowIndx = this->giveRowIndx(j - 1, i - 1);
    if ( rowIndx ) {
        return columns[ j - 1 ].at(rowIndx);
    }

    OOFEM_ERROR("Array accessing exception -- (%d,%d) out of bounds", i, j);
    return columns[ 0 ].at(1); // return to suppress compiler warning message
}


double DynCompCol :: at(int i, int j) const
{
    int rowIndx = this->giveRowIndx(j - 1, i - 1);
    if ( rowIndx ) {
        return columns[ j - 1 ].at(rowIndx);
    }

    if ( i <= nRows && j <= nColumns ) {
        return 0.0;
    } else {
        OOFEM_ERROR("Array accessing exception -- (%d,%d) out of bounds", i, j);
        return columns[ 0 ].at(1); // return to suppress compiler warning message
    }
}

double DynCompCol :: operator() (int i, int j)  const
{
    int rowIndx;
    if ( ( rowIndx = this->giveRowIndx(j, i) ) ) {
        return columns[ j ].at(rowIndx);
    }

    if ( i < nRows && j < nColumns ) {
        return 0.0;
    } else {
        OOFEM_ERROR("Array accessing exception -- (%d,%d) out of bounds", i, j);
        return columns[ 0 ].at(1); // return to suppress compiler warning message
    }
}

double &DynCompCol :: operator() (int i, int j)
{
    this->version++;

    int rowIndx = this->giveRowIndx(j, i);
    if ( rowIndx ) {
        return columns[ j ].at(rowIndx);
    }

    OOFEM_ERROR("Array element (%d,%d) not in sparse structure -- cannot assign", i, j);
    return columns[ 0 ].at(1); // return to suppress compiler warning message
}


void DynCompCol :: timesT(const FloatArray &x, FloatArray &answer) const
{
    if ( x.giveSize() != nRows ) {
        OOFEM_ERROR("Error in CompCol -- incompatible dimensions");
    }
    answer.resize(nColumns);
    answer.zero();

    for ( int i = 0; i < nColumns; i++ ) {
        double r = 0.0;
        for ( int t = 1; t <= columns[ i ].giveSize(); t++ ) {
            r += columns[ i ].at(t) * x[ rowind[ i ].at(t) ];
        }

        answer[i] = r;
    }
}


void DynCompCol :: checkSizeTowards(IntArray &loc)
{
    int maxid = 0;
    int size = loc.giveSize();
    // adjust the size of receiver
    for ( int i = 0; i < size; i++ ) {
        maxid = max( maxid, loc[i] );
    }

    this->growTo(maxid);

    for ( int i = 0; i < size; i++ ) {
        int ii = loc[i];
        if ( ii ) {
            for ( int j = 0; j < size; j++ ) {
                int jj = loc[j];
                if ( jj ) {
                    this->insertRowInColumn(ii - 1, jj - 1);
                }
            }
        }
    }
}


void DynCompCol :: checkSizeTowards(const IntArray &rloc, const IntArray &cloc)
{
    int maxid = 0;
    int rsize = rloc.giveSize();
    int csize = cloc.giveSize();

    // adjust the size of receiver
    for ( int i = 0; i < rsize; i++ ) {
        maxid = max( maxid, rloc[i] );
    }

    for ( int i = 0; i < csize; i++ ) {
        maxid = max( maxid, cloc[i] );
    }

    this->growTo(maxid);

    IntArray rowsToAdd( rloc.giveSize() );

    for ( int i = 0; i < csize; i++ ) {
        int ii = cloc[i];
        if ( ii ) {
            for ( int j = 0; j < rsize; j++ ) {
                int jj = rloc[j];
                if ( jj ) {
                    insertRowInColumn(ii - 1, jj - 1);
                }
            }
        }
    }
}



void DynCompCol :: growTo(int ns)
{
    if ( ns > nColumns ) {
        columns.resize(ns);
        rowind.resize(ns);
        nColumns = nRows = ns;
    }
}


int DynCompCol :: giveRowIndx(int col, int row) const
{
    // fast row indx search, based on assumption, that row indices are sorted
    int left = 1, right = this->rowind[ col ].giveSize();
    int middle = ( left + right ) / 2;
    int middleVal;

    if ( right == 0 ) {
        return 0;
    }

    if ( this->rowind[ col ].at(right) == row ) {
        return right;
    }

    while ( !( ( ( middleVal = this->rowind[ col ].at(middle) ) == row ) || ( middle == left ) ) ) {
        if ( row > middleVal ) {
            left = middle;
        } else {
            right = middle;
        }

        middle = ( left + right ) / 2;
    }

    if ( middleVal == row ) {
        return middle;
    }

    return 0;
}


int
DynCompCol :: insertRowInColumn(int col, int row)
{
    // insert row into column, preserving order of row indexes.
    int oldsize = this->rowind[ col ].giveSize();
    int left = 1, right = oldsize;
    int middle = ( left + right ) / 2;
    int middleVal;

    if ( oldsize == 0 ) {
        rowind[ col ].resizeWithValues(1, DynCompCol_CHUNK);
        columns[ col ].resizeWithValues(1, DynCompCol_CHUNK);
        columns[ col ].at(1) = 0.0;
        rowind[ col ].at(1) = row;
        return 1;
    }

    if ( this->rowind[ col ].at(right) == row ) {
        return right;
    }

    while ( !( ( ( middleVal = this->rowind[ col ].at(middle) ) == row ) || ( middle == left ) ) ) {
        if ( row > middleVal ) {
            left = middle;
        } else {
            right = middle;
        }

        middle = ( left + right ) / 2;
    }

    if ( middleVal == row ) {
        return middle;
    }

    // we have to insert new row entry
    if ( row > this->rowind[ col ].at(oldsize) ) {
        right = oldsize + 1;
    } else if ( row < this->rowind[ col ].at(1) ) {
        right = 1;
    }

    // insert row at middle+1 position
    rowind[ col ].resizeWithValues(oldsize + 1, DynCompCol_CHUNK);
    columns[ col ].resizeWithValues(oldsize + 1, DynCompCol_CHUNK);

    for ( int i = oldsize; i >= right; i-- ) {
        rowind[ col ].at(i + 1) = rowind[ col ].at(i);
        columns[ col ].at(i + 1) = columns[ col ].at(i);
    }

    columns[ col ].at(right) = 0.0;
    rowind[ col ].at(right) = row;
    return right;
}

} // end namespace oofem
