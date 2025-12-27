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

#include "dyncomprow.h"
#include "floatarray.h"
#include "engngm.h"
#include "domain.h"
#include "mathfem.h"
#include "verbose.h"
#include "element.h"
#include "sparsemtrxtype.h"
#include "activebc.h"
#include "classfactory.h"
#ifdef __MPM_MODULE
#include "../mpm/integral.h"
#endif

#ifdef TIME_REPORT
 #include "timer.h"
#endif

namespace oofem {
REGISTER_SparseMtrx(DynCompRow, SMT_DynCompRow);


DynCompRow :: DynCompRow(int n) : SparseMtrx(n, n),
    base(0)
{}


DynCompRow :: DynCompRow(const DynCompRow &S) : SparseMtrx(S.nRows, S.nColumns),
    rows(S.rows),
    colind(S.colind),
    diag(S.diag),
    base(S.base)
{
    this->version = S.version;
}


DynCompRow &DynCompRow :: operator = ( const DynCompRow & C )
{
    rows = C.rows;
    colind = C.colind;
    diag = C.diag;
    base = C.base;

    nRows   = C.nRows;
    nColumns = C.nColumns;
    version = C.version;
    return * this;
}

std::unique_ptr<SparseMtrx> DynCompRow :: clone() const
{
    return std::make_unique<DynCompRow>(*this);
}


void DynCompRow :: times(const FloatArray &x, FloatArray &answer) const
{
    if ( x.giveSize() != nColumns ) {
        OOFEM_ERROR("Error in CompRow -- incompatible dimensions");
    }

    answer.resize(nRows);
    answer.zero();

    for ( int j = 0; j < nRows; j++ ) {
        double r = 0.0;
        for ( int t = 1; t <= rows [ j ].giveSize(); t++ ) {
            r += rows [ j ].at(t) * x[ colind [ j ].at(t) ];
        }

        answer(j) = r;
    }
}

void DynCompRow :: times(double x)
{
    for ( auto &row : rows ) {
        row.times(x);
    }

    this->version++;
}

int DynCompRow :: buildInternalStructure(EngngModel *eModel, int di, const UnknownNumberingScheme &s)
{
    int neq = eModel->giveNumberOfDomainEquations(di, s);

    IntArray loc;
    Domain *domain = eModel->giveDomain(di);

#ifdef TIME_REPORT
    Timer timer;
    timer.startTimer();
#endif

    this->colind.clear();
    this->rows.clear();
    this->growTo(neq);

    for ( auto &elem : domain->giveElements() ) {
        elem->giveLocationArray(loc, s);
        for ( int i = 1; i <= loc.giveSize(); i++ ) { // row indx
            int ii;
            if ( ( ii = loc.at(i) ) ) {
                for ( int j = 1; j <= loc.giveSize(); j++ ) { // col indx
                    int jj;
                    if ( ( jj = loc.at(j) ) ) {
                        this->insertColInRow(ii - 1, jj - 1);
                    }
                }
            }
        }
    }


    std :: vector< IntArray >r_locs;
    std :: vector< IntArray >c_locs;

    for ( auto &gbc : domain->giveBcs() ) {
        ActiveBoundaryCondition *bc = dynamic_cast< ActiveBoundaryCondition * >( gbc.get() );
        if ( bc ) {
            bc->giveLocationArrays(r_locs, c_locs, UnknownCharType, s, s);
            for ( std :: size_t k = 0; k < r_locs.size(); k++ ) {
                IntArray &krloc = r_locs [ k ];
                IntArray &kcloc = c_locs [ k ];
                for ( int i = 1; i <= krloc.giveSize(); i++ ) {
                    int ii;
                    if ( ( ii = krloc.at(i) ) ) {
                        for ( int j = 1; j <= kcloc.giveSize(); j++ ) {
                            int jj;
                            if ( ( jj = kcloc.at(j) ) ) {
                                this->insertColInRow(ii - 1, jj - 1);
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
            for ( int i = 1; i <= locr.giveSize(); i++ ) {
                int ii;
                if ( ( ii = locr.at(i) ) ) {
                    for ( int j = 1; j <= locc.giveSize(); j++ ) {
                        int jj;
                        if ( ( jj = locc.at(j) ) ) {
                            this->insertColInRow(ii - 1, jj - 1);
                        }
                    }
                }
            }
        }
    }
#endif


    int nz_ = 0;
    for ( auto &col : this->colind ) {
        nz_ += col.giveSize();
    }

    OOFEM_LOG_DEBUG("DynCompRow info: neq is %d, nelem is %d\n", neq, nz_);

    this->version++;
#ifdef TIME_REPORT
    timer.stopTimer();
    OOFEM_LOG_DEBUG( "DynCompRow::buildInternalStructure: user time consumed: %.2fs\n", timer.getUtime() );
#endif

    return true;
}


int DynCompRow :: assemble(const IntArray &loc, const FloatMatrix &mat)
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

int DynCompRow :: assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat)
{
    // optimized low-end implementation
    IntArray colsToAdd( rloc.giveSize() );
    int rsize = rloc.giveSize();
    int csize = cloc.giveSize();

    for ( int i = 0; i < rsize; i++ ) {
        int ii;
        if ( ( ii = rloc[i] ) ) {
            int ii1 = ii - 1;
            for ( int j = 0; j < csize; j++ ) {
                int jj;
                if ( ( jj = cloc[j] ) ) {
                    int jj1 = jj - 1;

                    int colindx = this->insertColInRow(ii1, jj1);
                    rows [ ii1 ].at(colindx) += mat(i, j);
                }
            }
        }
    }

    this->version++;
    return 1;
}

void DynCompRow :: zero()
{
    for ( auto &row: rows ) {
        row.zero();
    }

    this->version++;
}


void DynCompRow :: printStatistics() const
{
    int nz = 0;
    for ( auto &row : rows ) {
        nz += row.giveSize();
    }

    printf("\nDynCompRow info: neq is %d, nelem is %d\n", nRows, nz);
}


double &DynCompRow :: at(int i, int j)
{
    int colIndx;

    this->version++;
    if ( ( colIndx = this->giveColIndx(i - 1, j - 1) ) ) {
        return rows [ i - 1 ].at(colIndx);
    }

    OOFEM_ERROR("Array accessing exception -- (%d,%d) out of bounds", i, j);
    return rows [ 0 ].at(1); // return to suppress compiler warning message
}


double DynCompRow :: at(int i, int j) const
{
    int colIndx;
    if ( ( colIndx = this->giveColIndx(i - 1, j - 1) ) ) {
        return rows [ i - 1 ].at(colIndx);
    }

    if ( i <= nRows && j <= nColumns ) {
        return 0.0;
    } else {
        OOFEM_ERROR("Array accessing exception -- (%d,%d) out of bounds", i, j);
        return rows [ 0 ].at(1); // return to suppress compiler warning message
    }
}

double DynCompRow :: operator() (int i, int j)  const
{
    int colIndx;
    if ( ( colIndx = this->giveColIndx(i, j) ) ) {
        return rows [ i ].at(colIndx);
    }

    if ( i < nRows && j < nColumns ) {
        return 0.0;
    } else {
        OOFEM_ERROR("Array accessing exception -- (%d,%d) out of bounds", i, j);
        return rows [ 0 ].at(1); // return to suppress compiler warning message
    }
}

double &DynCompRow :: operator() (int i, int j)
{
    int colIndx;

    // increment version
    this->version++;

    if ( ( colIndx = this->giveColIndx(i, j) ) ) {
        return rows [ i ].at(colIndx);
    }

    OOFEM_ERROR("Array element (%d,%d) not in sparse structure -- cannot assign", i, j);
    return rows [ 0 ].at(1); // return to suppress compiler warning message
}

void DynCompRow :: timesT(const FloatArray &x, FloatArray &answer) const
{
    //      Check for compatible dimensions:
    if ( x.giveSize() != nRows ) {
        OOFEM_ERROR("Error in CompRow -- incompatible dimensions");
    }

    answer.resize(nColumns);
    answer.zero();

    for ( int i = 0; i < nColumns; i++ ) {
        double r = x[i];
        for ( int t = 1; t <= rows [ i ].giveSize(); t++ ) {
            answer[ colind [ i ].at(t) ] += rows [ i ].at(t) * r;
        }
    }
}



void DynCompRow :: checkSizeTowards(IntArray &loc)
{
    int maxid = 0;
    int size = loc.giveSize();
    // adjust the size of receiver
    for ( int i = 0; i < size; i++ ) {
        maxid = max( maxid, loc[i] );
    }

    this->growTo(maxid);

    for ( int i = 0; i < size; i++ ) {
        int ii;
        if ( ( ii = loc[i] ) ) {
            for ( int j = 0; j < size; j++ ) {
                int jj;
                if ( ( jj = loc[j] ) ) {
                    this->insertColInRow(ii - 1, jj - 1);
                }
            }
        }
    }
}


void DynCompRow :: checkSizeTowards(const IntArray &rloc, const IntArray &cloc)
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

    for ( int i = 0; i < rsize; i++ ) {
        int ii;
        if ( ( ii = rloc[i] ) ) {
            for ( int j = 0; j < csize; j++ ) {
                int jj;
                if ( ( jj = cloc[j] ) ) {
                    this->insertColInRow(ii - 1, jj - 1);
                }
            }
        }
    }
}


void DynCompRow :: growTo(int ns)
{
    if ( ns > nRows ) {
        this->rows.resize(ns);
        this->colind.resize(ns);
        nColumns = nRows = ns;
    }
}


int DynCompRow :: giveColIndx(int row, int col) const
{
    // fast col indx search, based on assumption, that col indices are sorted
    int left = 1, right = this->colind [ row ].giveSize();
    int middle = ( left + right ) / 2;
    int middleVal;

    if ( right == 0 ) {
        return 0;
    }

    if ( this->colind [ row ].at(right) == col ) {
        return right;
    }

    while ( !( ( ( middleVal = this->colind [ row ].at(middle) ) == col ) || ( middle == left ) ) ) {
        if ( col > middleVal ) {
            left = middle;
        } else {
            right = middle;
        }

        middle = ( left + right ) / 2;
    }

    if ( middleVal == col ) {
        return middle;
    }

    return 0;
}


int
DynCompRow :: insertColInRow(int row, int col)
{
    // insert col entry into row, preserving order of col indexes.
    int oldsize = this->colind [ row ].giveSize();
    int left = 1, right = oldsize;
    int middle = ( left + right ) / 2;
    int middleVal;

    if ( oldsize == 0 ) {
        colind [ row ].resizeWithValues(1, DynCompRow_CHUNK);
        rows [ row ].resizeWithValues(1, DynCompRow_CHUNK);
        rows [ row ].at(1) = 0.0;
        colind [ row ].at(1) = col;
        return 1;
    }

    if ( this->colind [ row ].at(right) == col ) {
        return right;
    }

    while ( !( ( ( middleVal = this->colind [ row ].at(middle) ) == col ) || ( middle == left ) ) ) {
        if ( col > middleVal ) {
            left = middle;
        } else {
            right = middle;
        }

        middle = ( left + right ) / 2;
    }

    if ( middleVal == col ) {
        return middle;
    }

    // we have to insert new row entry
    if ( col > this->colind [ row ].at(oldsize) ) {
        right = oldsize + 1;
    } else if ( col < this->colind [ row ].at(1) ) {
        right = 1;
    }

    // insert col at middle+1 position
    colind [ row ].resizeWithValues(oldsize + 1, DynCompRow_CHUNK);
    rows [ row ].resizeWithValues(oldsize + 1, DynCompRow_CHUNK);

    for ( int i = oldsize; i >= right; i-- ) {
        colind [ row ].at(i + 1) = colind [ row ].at(i);
        rows [ row ].at(i + 1) = rows [ row ].at(i);
    }

    rows [ row ].at(right) = 0.0;
    colind [ row ].at(right) = col;
    return right;
}


//#define ILU_0
//#define ILU_DROP_TOL 1.e-8
#define ILU_ROW_CHUNK 10
//#define ILU_PART_FILL 5

void
DynCompRow :: ILUPYourself(int part_fill, double drop_tol)
{
    IntArray irw(nColumns), iw;
    FloatArray w;
    diag.resize(nRows);

#ifdef TIME_REPORT
    Timer timer;
    timer.startTimer();
#endif

    for ( int i = 0; i < nRows; i++ ) { // row loop
        if ( ( diag[i] = giveColIndx(i, i) ) == 0 ) { // giveColIndx returns 1-based indexing
            OOFEM_ERROR("zero value on diagonal");
        }
    }

    /* FACTOR MATRIX */

    for ( int i = 1; i < nRows; i++ ) { // loop  over rows
        double inorm = 0.0;
        for ( int ii = 1; ii <= rows [ i ].giveSize(); ii++ ) {
            double val = rows [ i ].at(ii);
            inorm += val * val;
        }

        inorm = sqrt(inorm);

        w.resizeWithValues(rows [ i ].giveSize(), ILU_ROW_CHUNK);
        iw.resizeWithValues(rows [ i ].giveSize(), ILU_ROW_CHUNK);
        for ( int kk = 1; kk <= rows [ i ].giveSize(); kk++ ) {
            irw[ colind [ i ].at(kk) ] = kk;
            iw[kk - 1] = colind [ i ].at(kk);
            w[kk - 1] = rows [ i ].at(kk);
        }

        //for (k=0; k < (diag(i)-1); k++) { // loop 1,...,i-1 for (i,k) \in NZ(A)
        int k = 0;
        while ( iw.at(k + 1) < i ) {
            // initialize k-th row indexes
            int krow = iw.at(k + 1);
            //multiplier = (rows[i].at(k+1) /= rows[krow].at(diag(krow)));
            double multiplier = ( w.at(k + 1) /= rows [ krow ].at( diag[krow] ) );

#ifndef ILU_0
            // first dropping rule for aik
            if ( fabs(multiplier) >= drop_tol * inorm )
#endif
            { // first drop rule
                for ( int j = 0; j < colind [ krow ].giveSize(); j++ ) {
                    int jcol = colind [ krow ].at(j + 1);
                    if ( jcol > krow ) {
                        if ( irw[jcol] ) {
                            //rows[i].at(irw(jcol)) -= multiplier*rows[krow].at(j+1);
                            w.at( irw[jcol] ) -= multiplier * rows [ krow ].at(j + 1);
                        } else {
#ifndef ILU_0
                            // insert new entry
                            int newsize = w.giveSize() + 1;
                            w.resizeWithValues(newsize, ILU_ROW_CHUNK);
                            iw.resizeWithValues(newsize, ILU_ROW_CHUNK);

                            iw.at(newsize) = jcol;
                            w.at(newsize) = -multiplier * rows [ krow ].at(j + 1);
                            irw[jcol] = newsize;
#endif

                            /*
                             * ipos = insertColInRow (i,jcol) ;
                             * for (kk=ipos+1;  kk<= rows[i].giveSize(); kk++)
                             * irw(colind[i].at(kk))++;
                             *
                             * ipos = insertColInRow (i,jcol) ;
                             * rows[i].at(ipos) = -multiplier*rows[krow].at(j+1);
                             * irw(jcol) = ipos;
                             * if (jcol < i) diag(i)++;
                             */
                        }
                    }
                }
            }

            // scan iw to find closest index to krow
            int ck = nColumns + 1;
            for ( int kk = 0; kk < iw.giveSize(); kk++ ) {
                if ( ( ( iw[kk] - krow ) > 0 ) && ( ( iw[kk] - krow ) < ( ck - krow ) ) ) {
                    ck = iw[kk];
                }
            }

            k = irw[ck] - 1;
        }

#ifndef ILU_0

        int end = iw.giveSize();
        int curr = 1;
        // second drop rule
        while ( curr <= end ) {
            if ( ( fabs( w.at(curr) ) < drop_tol * inorm ) && ( iw.at(curr) != i ) ) {
                // remove entry
                w.at(curr) = w.at(end);
                irw[ iw.at(curr) ] = 0;
                iw.at(curr) = iw.at(end);
                if ( curr != end ) {
                    irw[ iw.at(curr) ] = curr;
                }

                end--;
            } else {
                curr++;
            }
        }

        // cutt off
        w.resize(end);
        iw.resize(end);

        int count = end;

        // select only the p-largest w values
        this->qsortRow(iw, irw, w, 0, iw.giveSize() - 1);
        //
        int lsizeLimit = diag[i] - 1;
        int usizeLimit = rows [ i ].giveSize() - lsizeLimit;

        lsizeLimit += part_fill;
        usizeLimit += part_fill;

        int lnums = 0;
        int unums = 0;
        count = 0;
        for ( int kk = 1; kk <= iw.giveSize(); kk++ ) {
            if ( iw.at(kk) < i ) { // lpart
                if ( ++lnums > lsizeLimit ) {
                    irw[ iw.at(kk) ] = 0;
                } else {
                    count++;
                }
            } else if ( iw.at(kk) > i ) { // upart
                if ( ++unums > usizeLimit ) {
                    irw[ iw.at(kk) ] = 0;
                } else {
                    count++;
                }
            } else { // diagonal is always kept
                count++;
            }
        }

#else
        int count = iw.giveSize();
#endif
        rows [ i ].resize(count);
        colind [ i ].resize(count);

        int icount = 1;
        int indx, idist, previndx = -1;
        int kkend = iw.giveSize();

        for ( int kk = 1; kk <= count; kk++ ) {
            idist = nColumns + 2;
            indx = 0;

            for ( int kki = 1; kki <= kkend; kki++ ) {
                if ( ( irw[ iw.at(kki) ] != 0 ) && ( iw.at(kki) > previndx ) && ( ( iw.at(kki) - previndx ) < idist ) ) {
                    idist = iw.at(kki) - previndx;
                    indx = kki;
                }
            }

            if ( indx == 0 ) {
                OOFEM_ERROR("internal error");
            }

            previndx = iw.at(indx);
            rows [ i ].at(icount) = w.at(indx);
            colind [ i ].at(icount) = iw.at(indx);
            if ( colind [ i ].at(icount) == i ) {
                diag[i] = icount;
            }

            icount++;

            // exclude the indx entry from search by moving it to the end of list
            irw[ iw.at(indx) ] = 0;
            iw.at(indx) = iw.at(kkend);
            w.at(indx)  = w.at(kkend);
            if ( irw[ iw.at(indx) ] != 0 ) {
                irw[ iw.at(indx) ] = indx;
            }

            kkend--;

            // exclude the indx entry from search by moving it to the end of list
#if 0
            std::swap(irw[iw.at(indx)], irw[iw.at(kkend)]);
            std::swap(iw.at(indx), iw.at(kkend));
            std::swap(w.at(indx), w.at(kkend));
            kkend--;
#endif
        }


#if 0
        int icount = 1;
        for (kk=1;  kk<= nColumns; kk++) {
            if ( irw.at(kk) > 0 ) {
                rows[i].at(icount) = w.at(abs(irw[kk-1]));
                colind[i].at(icount) = iw.at(abs(irw[kk-1]));
                if (colind[i].at(icount) == i) diag[i] = icount;
                icount++;
            }
        }
#endif
        if ( ( icount - count ) != 1 ) {
            OOFEM_ERROR("%d - row errorr (%d,%d)\n", i, icount, count);
        }

        //Refresh all iw enries to zero
        for ( int kk = 1; kk <= iw.giveSize(); kk++ ) {
            irw[ iw.at(kk) ] = 0;
        }

        //irw.zero();
    }

#ifdef TIME_REPORT
    timer.stopTimer();
    OOFEM_LOG_DEBUG( "\nILUT(%d,%e): user time consumed by factorization: %.2fs\n", part_fill, drop_tol, timer.getUtime() );
#endif

    this->version++;
}



void
DynCompRow :: ILUPsolve(const FloatArray &x, FloatArray &y) const
{
    int M = x.giveSize();
    FloatArray work(M);

    // solve Lw=x
    for ( int i = 0; i < M; i++ ) {
        double r = x[i];
        for ( int t = 0; t < ( diag[i] - 1 ); t++ ) {
            r -= rows [ i ].at(t + 1) * work( colind [ i ].at(t + 1) );
        }

        work[i] = r;
    }

    y.resize(M);
    y.zero();
    // solve Uy=w
    for ( int i = M - 1; i >= 0; i-- ) {
        double r = work[i];
        for ( int t = diag[i]; t < rows [ i ].giveSize(); t++ ) {
            r -= rows [ i ].at(t + 1) * y( colind [ i ].at(t + 1) );
        }

        y[i] = r / rows [ i ].at( diag[i] );
    }
}


void
DynCompRow :: ILUPtrans_solve(const FloatArray &x, FloatArray &y) const
{
    int M = x.giveSize();
    FloatArray work(M);

    y.resize(M);

    // solve for U^Tw = x
    for ( int i = 0; i < M; i++ ) {
        work[i] = ( x[i] + work[i] ) / rows [ i ].at( diag[i] );
        for ( int t = diag[i]; t < rows [ i ].giveSize(); t++ ) {
            work( colind [ i ].at(t + 1) ) -= rows [ i ].at(t + 1) * work[i];
        }
    }

    // solve for L^T y = w
    for ( int i = M - 1; i >= 0; i-- ) {
        y[i] = work[i];
        for ( int t = 1; t < diag[i]; t++ ) {
            work( colind [ i ].at(t) ) -= rows [ i ].at(t) * y[i];
        }
    }
}


void
DynCompRow :: qsortRow(IntArray &ind, IntArray &ir, FloatArray &val, int l, int r)
{
    if ( r <= l ) {
        return;
    }

    int i = qsortRowPartition(ind, ir, val, l, r);
    qsortRow(ind, ir, val, l, i - 1);
    qsortRow(ind, ir, val, i + 1, r);
}


int
DynCompRow :: qsortRowPartition(IntArray &ind, IntArray &ir, FloatArray &val, int l, int r)
{
    int i = l - 1, j = r;
    double v = fabs( val[r] );

    for ( ; ; ) {
        while ( ( fabs( val(++i) ) >  v ) ) {
            ;
        }

        while ( ( v > fabs( val(--j) ) ) ) {
            if ( j == l ) {
                break;
            }
        }

        if ( i >= j ) {
            break;
        }

        std::swap(ir[ ind[i] ], ir[ ind[j] ]);
        std::swap(ind[i], ind[j]);
        std::swap(val[i], val[j]);
    }

    std::swap(ir[ ind[i] ], ir[ ind[r] ]);
    std::swap(ind[i], ind[r]);
    std::swap(val[i], val[r]);

    return i;
}
} // end namespace oofem
