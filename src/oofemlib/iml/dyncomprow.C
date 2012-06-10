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
 *               Copyright (C) 1993 - 2012   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */
// Class DynCompRow

// inspired by SL++

#include "dyncomprow.h"
#include "flotarry.h"
#include "engngm.h"
#include "domain.h"
#include "mathfem.h"
#include "verbose.h"
#include "element.h"

#ifdef TIME_REPORT
 #include "clock.h"
#endif

namespace oofem {
DynCompRow :: DynCompRow(void) : SparseMtrx(), base_(0)
{
    rows_ = NULL;
    colind_  = NULL;
}


DynCompRow :: DynCompRow(int n) : SparseMtrx(n, n), base_(0)
{
    rows_ = NULL;
    colind_  = NULL;

    nRows = nColumns = n;
}


/*****************************/
/*  Copy constructor         */
/*****************************/

DynCompRow :: DynCompRow(const DynCompRow &S) : SparseMtrx(S.nRows, S.nColumns), base_(S.base_)
{
    int i;
    if ( S.rows_ ) {
        this->rows_ = new FloatArray * [ S.nRows ];
        for ( i = 0; i < S.nRows; i++ ) {
            this->rows_ [ i ] = new FloatArray(*S.rows_ [ i ]);
        }
    } else {
        this->rows_ = NULL;
    }

    if ( S.colind_ ) {
        this->colind_ = new IntArray * [ S.nRows ];
        for ( i = 0; i < S.nRows; i++ ) {
            this->colind_ [ i ] = new IntArray(*S.colind_ [ i ]);
        }
    } else {
        this->colind_ = NULL;
    }

    this->nRows = S.nRows;
    this->nColumns = S.nColumns;
    this->version = S.version;
}


// Destructor
DynCompRow :: ~DynCompRow()
{
    int i;

    if ( this->rows_ ) {
        for ( i = 0; i < nRows; i++ ) {
            delete this->rows_ [ i ];
        }

        delete this->rows_;
    }

    if ( this->colind_ ) {
        for ( i = 0; i < nRows; i++ ) {
            delete this->colind_ [ i ];
        }

        delete this->colind_;
    }
}



/***************************/
/* Assignment operator...  */
/***************************/

DynCompRow &DynCompRow :: operator = ( const DynCompRow & C )
{
    base_   = C.base_;

    int i;

    if ( this->rows_ ) {
        for ( i = 0; i < nRows; i++ ) {
            delete this->rows_ [ i ];
        }

        delete this->rows_;
    }

    if ( C.rows_ ) {
        this->rows_ = new FloatArray * [ C.nRows ];
        for ( i = 0; i < C.nRows; i++ ) {
            this->rows_ [ i ] = new FloatArray(*C.rows_ [ i ]);
        }
    } else {
        this->rows_ = NULL;
    }


    if ( this->colind_ ) {
        for ( i = 0; i < nRows; i++ ) {
            delete this->colind_ [ i ];
        }

        delete this->colind_;
    }

    if ( C.colind_ ) {
        this->colind_ = new IntArray * [ C.nRows ];
        for ( i = 0; i < C.nRows; i++ ) {
            this->colind_ [ i ] = new IntArray(*C.colind_ [ i ]);
        }
    } else {
        this->colind_ = NULL;
    }

    nRows   = C.nRows;
    nColumns = C.nColumns;
    version = C.version;
    return * this;
}

SparseMtrx *DynCompRow :: GiveCopy() const
{
    DynCompRow *result = new DynCompRow(*this);
    return result;
}


void DynCompRow :: times(const FloatArray &x, FloatArray &answer) const
{
    //      Check for compatible dimensions:
    if ( x.giveSize() != nColumns ) {
        OOFEM_ERROR("DynCompRow::times: Error in CompRow -- incompatible dimensions");
    }

    answer.resize(nRows);
    answer.zero();

    int j, t;
    double r;

    for ( j = 0; j < nRows; j++ ) {
        r = 0.0;
        for ( t = 1; t <= rows_ [ j ]->giveSize(); t++ ) {
            r += rows_ [ j ]->at(t) * x( colind_ [ j ]->at(t) );
        }

        answer(j) = r;
    }
}

void DynCompRow :: times(double x)
{
    for ( int j = 0; j < nRows; j++ ) {
        rows_ [ j ]->times(x);
    }

    // increment version
    this->version++;
}

int DynCompRow :: buildInternalStructure(EngngModel *eModel, int di, EquationID ut, const UnknownNumberingScheme &s)
{
    /*
     * int neq = eModel -> giveNumberOfDomainEquations (di);
     *
     * IntArray  loc;
     * Domain* domain = eModel->giveDomain(di);
     * int nelem = domain -> giveNumberOfElements() ;
     * int i,ii,j,jj,n, indx;
     * Element* elem;
     * IntArray rowItems(neq);
     * // allocation map
     * char* map = new char[neq*neq]; // row by row storage
     * if (map == NULL) {
     * printf ("CompRow::buildInternalStructure - map creation failed");
     * exit (1);
     * }
     *
     * for (i=0; i<neq*neq; i++)
     * map[i]=0;
     *
     *
     * for (n=1 ; n<=nelem ; n++) {
     * elem = domain -> giveElement(n);
     * elem -> giveLocationArray (loc) ;
     *
     * for (i=1 ; i <= loc.giveSize() ; i++) { // row indx
     * if ((ii = loc.at(i))) {
     *  for (j=1; j <= loc.giveSize() ; j++) { // col indx
     *   if ((jj=loc.at(j)))
     *    if (map[(ii-1)*neq+jj-1] == 0) {
     *     map[(ii-1)*neq+jj-1] = 1;
     *     rowItems.at(ii) ++;
     *    }
     *  }
     * }
     * }
     * }
     *
     * if (colind_) {
     * for (i=0; i< nRows; i++) delete this->colind_[i];
     * delete this->colind_;
     * }
     * colind_ = (IntArray**) new (IntArray*)[neq];
     * for (j=0; j<neq; j++) colind_[j] = new IntArray(rowItems(j));
     *
     * indx = 1;
     * for (i=0; i<neq; i++) { // row loop
     * indx = 1;
     * for (j=0; j<neq; j++) { // column loop
     * if (map[i*neq+j]) {
     *  colind_[i]->at(indx) = j;
     *  indx++;
     * }
     * }
     * }
     *
     * // delete map
     * delete (map);
     *
     * // allocate value array
     * if (rows_) {
     * for (i=0; i< nRows; i++) delete this->rows_[i];
     * delete this->rows_;
     * }
     * rows_= (FloatArray**) new (FloatArray*)[neq];
     * int nz_ = 0;
     * for (j=0; j<neq; j++) {
     * rows_[j] = new FloatArray (rowItems(j));
     * nz_ += rowItems(j);
     * }
     *
     * printf ("\nDynCompRow info: neq is %d, nelem is %d\n",neq,nz_);
     * nColumns = nRows = neq;
     *
     * // increment version
     * this->version++;
     * return true;
     */

    int neq = eModel->giveNumberOfDomainEquations(di, ut);

    IntArray loc;
    Domain *domain = eModel->giveDomain(di);
    int nelem = domain->giveNumberOfElements();
    int i, ii, j, jj, n;
    Element *elem;

#ifdef TIME_REPORT
    //clock_t tstart = clock();
    oofem_timeval tstart;
    getUtime(tstart);
#endif

    nColumns = nRows = neq;

    if ( colind_ ) {
        for ( i = 0; i < nRows; i++ ) {
            delete this->colind_ [ i ];
        }

        delete this->colind_;
    }

    colind_ = ( IntArray ** ) new IntArray * [ neq ];
    for ( j = 0; j < neq; j++ ) {
        colind_ [ j ] = new IntArray();
    }

    // allocate value array
    if ( rows_ ) {
        for ( i = 0; i < nRows; i++ ) {
            delete this->rows_ [ i ];
        }

        delete this->rows_;
    }

    rows_ = ( FloatArray ** ) new FloatArray * [ neq ];
    for ( j = 0; j < neq; j++ ) {
        rows_ [ j ] = new FloatArray();
    }

    for ( n = 1; n <= nelem; n++ ) {
        elem = domain->giveElement(n);
        elem->giveLocationArray(loc, ut, s);

        for ( i = 1; i <= loc.giveSize(); i++ ) { // row indx
            if ( ( ii = loc.at(i) ) ) {
                for ( j = 1; j <= loc.giveSize(); j++ ) { // col indx
                    if ( ( jj = loc.at(j) ) ) {
                        this->insertColInRow(ii - 1, jj - 1);
                    }
                }
            }
        }
    }

    int nz_ = 0;
    for ( j = 0; j < neq; j++ ) {
        nz_ += this->colind_ [ j ]->giveSize();
    }

    OOFEM_LOG_DEBUG("DynCompRow info: neq is %d, nelem is %d\n", neq, nz_);

    // increment version
    this->version++;
#ifdef TIME_REPORT
    oofem_timeval tfin;
    getRelativeUtime(tfin, tstart);
    OOFEM_LOG_DEBUG( "DynCompRow::buildInternalStructure: user time consumed: %.2fs\n",
                    ( double ) ( tfin.tv_sec + tfin.tv_usec / ( double ) OOFEM_USEC_LIM ) );
#endif

    return true;
}


int DynCompRow :: assemble(const IntArray &loc, const FloatMatrix &mat)
{
    int i, j, ii, jj, dim;

#  ifdef DEBUG
    dim = mat.giveNumberOfRows();
    if ( dim != loc.giveSize() ) {
        OOFEM_ERROR("DynCompRow::assemble : dimension of 'k' and 'loc' mismatch");
    }

#  endif

    dim = mat.giveNumberOfRows();

    for ( j = 1; j <= dim; j++ ) {
        jj = loc.at(j);
        if ( jj ) {
            for ( i = 1; i <= dim; i++ ) {
                ii = loc.at(i);
                if ( ii ) {
                    this->at(ii, jj) += mat.at(i, j);
                }
            }
        }
    }

    // increment version
    this->version++;
    return 1;
}

int DynCompRow :: assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat)
{
    /// optimized low-end implementation
    IntArray colsToAdd( rloc.giveSize() );
    int i, ii, ii1, j, jj, jj1, colindx;
    int rsize = rloc.giveSize();
    int csize = cloc.giveSize();

    for ( i = 0; i < rsize; i++ ) {
        if ( ( ii = rloc(i) ) ) {
            ii1 = ii - 1;
            for ( j = 0; j < csize; j++ ) {
                if ( ( jj = cloc(j) ) ) {
                    jj1 = jj - 1;

                    colindx = this->insertColInRow(ii1, jj1);
                    //if (rowind_[j-1]->at(t) == (i-1)) return columns_[j-1]->at(t);
                    rows_ [ ii1 ]->at(colindx) += mat(i, j);
                }
            }
        }
    }

    // increment version
    this->version++;
    return 1;
}

void DynCompRow :: zero()
{
    for ( int j = 0; j < nRows; j++ ) {
        rows_ [ j ]->zero();
    }

    // increment version
    this->version++;
}


void DynCompRow :: printStatistics() const
{
    int nz_ = 0;
    for ( int j = 0; j < nRows; j++ ) {
        nz_ += rows_ [ j ]->giveSize();
    }

    printf("\nDynCompRow info: neq is %d, nelem is %d\n", nRows, nz_);
}


/*********************/
/*   Array access    */
/*********************/

double &DynCompRow :: at(int i, int j)
{
    int colIndx;

    // increment version
    this->version++;
    if ( ( colIndx = this->giveColIndx(i - 1, j - 1) ) ) {
        return rows_ [ i - 1 ]->at(colIndx);
    }

    OOFEM_ERROR3("DynCompRow::operator at(): Array accessing exception -- (%d,%d) out of bounds", i, j);
    return rows_ [ 0 ]->at(1); // return to suppress compiler warning message
}


double DynCompRow :: at(int i, int j) const
{
    int colIndx;
    if ( ( colIndx = this->giveColIndx(i - 1, j - 1) ) ) {
        return rows_ [ i - 1 ]->at(colIndx);
    }

    if ( i <= nRows && j <= nColumns ) {
        return 0.0;
    } else {
        OOFEM_ERROR3("DynCompRow::operator at(): Array accessing exception -- (%d,%d) out of bounds", i, j);
        return rows_ [ 0 ]->at(1); // return to suppress compiler warning message
    }
}

double DynCompRow :: operator() (int i, int j)  const
{
    int colIndx;
    if ( ( colIndx = this->giveColIndx(i, j) ) ) {
        return rows_ [ i ]->at(colIndx);
    }

    if ( i < nRows && j < nColumns ) {
        return 0.0;
    } else {
        OOFEM_ERROR3("DynCompRow::operator(): Array accessing exception -- (%d,%d) out of bounds", i, j);
        return rows_ [ 0 ]->at(1); // return to suppress compiler warning message
    }
}

double &DynCompRow :: operator() (int i, int j)
{
    int colIndx;

    // increment version
    this->version++;

    if ( ( colIndx = this->giveColIndx(i, j) ) ) {
        return rows_ [ i ]->at(colIndx);
    }

    OOFEM_ERROR3("DynCompRow::operator(): Array element (%d,%d) not in sparse structure -- cannot assign", i, j);
    return rows_ [ 0 ]->at(1); // return to suppress compiler warning message
}

void DynCompRow :: timesT(const FloatArray &x, FloatArray &answer) const
{
    //      Check for compatible dimensions:
    if ( x.giveSize() != nRows ) {
        OOFEM_ERROR("DynCompRow::trans_mult: Error in CompRow -- incompatible dimensions");
    }

    answer.resize(nColumns);
    answer.zero();

    int i, t;
    double r;

    for ( i = 0; i < nColumns; i++ ) {
        r = x(i);
        for ( t = 1; t <= rows_ [ i ]->giveSize(); t++ ) {
            answer( colind_ [ i ]->at(t) ) += rows_ [ i ]->at(t) * r;
        }
    }
}



void DynCompRow :: checkSizeTowards(IntArray &loc)
{
    int i, maxid = 0;
    int size = loc.giveSize();
    // adjust the size of receiver
    for ( i = 0; i < size; i++ ) {
        maxid = max( maxid, loc(i) );
    }

    this->growTo(maxid);

    int ii, j, jj;

    for ( i = 0; i < size; i++ ) {
        if ( ( ii = loc(i) ) ) {
            for ( j = 0; j < size; j++ ) {
                if ( ( jj = loc(j) ) ) {
                    this->insertColInRow(ii - 1, jj - 1);
                }
            }
        }
    }
}


void DynCompRow :: checkSizeTowards(const IntArray &rloc, const IntArray &cloc)
{
    int i, maxid = 0;
    int rsize = rloc.giveSize();
    int csize = cloc.giveSize();

    // adjust the size of receiver
    for ( i = 0; i < rsize; i++ ) {
        maxid = max( maxid, rloc(i) );
    }

    for ( i = 0; i < csize; i++ ) {
        maxid = max( maxid, cloc(i) );
    }

    this->growTo(maxid);

    int ii, j, jj;

    for ( i = 0; i < rsize; i++ ) {
        if ( ( ii = rloc(i) ) ) {
            for ( j = 0; j < csize; j++ ) {
                if ( ( jj = cloc(j) ) ) {
                    this->insertColInRow(ii - 1, jj - 1);
                }
            }
        }
    }
}



void DynCompRow :: growTo(int ns)
{
    if ( ns > nRows ) {
        FloatArray **newrows_ = new FloatArray * [ ns ];
        IntArray **newcolind_ = new IntArray * [ ns ];

        // copy existing columns
        for ( int i = 0; i < nRows; i++ ) {
            newrows_ [ i ]   = rows_ [ i ];
            newcolind_ [ i ] = colind_ [ i ];
        }

        delete rows_;
        delete colind_;

        rows_ = newrows_;
        colind_  = newcolind_;

        nColumns = nRows = ns;
    }
}


int DynCompRow :: giveColIndx(int row, int col) const
{
    // fast col indx search, based on assumption, that col indices are sorted
    int left = 1, right = this->colind_ [ row ]->giveSize();
    int middle = ( left + right ) / 2;
    int middleVal;

    if ( right == 0 ) {
        return 0;
    }

    if ( this->colind_ [ row ]->at(right) == col ) {
        return right;
    }

    while ( !( ( ( middleVal = this->colind_ [ row ]->at(middle) ) == col ) || ( middle == left ) ) ) {
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
    int i, oldsize = this->colind_ [ row ]->giveSize();
    int left = 1, right = oldsize;
    int middle = ( left + right ) / 2;
    int middleVal;

    if ( oldsize == 0 ) {
        colind_ [ row ]->resize(1, DynCompRow_CHUNK);
        rows_ [ row ]->resize(1, DynCompRow_CHUNK);
        rows_ [ row ]->at(1) = 0.0;
        colind_ [ row ]->at(1) = col;
        return 1;
    }

    if ( this->colind_ [ row ]->at(right) == col ) {
        return right;
    }

    while ( !( ( ( middleVal = this->colind_ [ row ]->at(middle) ) == col ) || ( middle == left ) ) ) {
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
    if ( col > this->colind_ [ row ]->at(oldsize) ) {
        right = oldsize + 1;
    } else if ( col < this->colind_ [ row ]->at(1) ) {
        right = 1;
    }

    // insert col at middle+1 position
    colind_ [ row ]->resize(oldsize + 1, DynCompRow_CHUNK);
    rows_ [ row ]->resize(oldsize + 1, DynCompRow_CHUNK);

    for ( i = oldsize; i >= right; i-- ) {
        colind_ [ row ]->at(i + 1) = colind_ [ row ]->at(i);
        rows_ [ row ]->at(i + 1) = rows_ [ row ]->at(i);
    }

    rows_ [ row ]->at(right) = 0.0;
    colind_ [ row ]->at(right) = col;
    return right;
}





/* reference ILU(0) version
 * void
 * DynCompRow::ILUPYourself ()
 * {
 * int i, j, k, kk, krow, colk;
 * double multiplier;
 *
 * IntArray iw (nColumns);
 * diag_rowptr_.resize(nRows);
 *
 *#ifndef DynCompRow_USE_STL_SETS
 * for (i = 0; i < nRows; i++) // row loop
 *  if ((diag_rowptr_(i) = giveColIndx (i, i)) == 0) { // giveColIndx returns 1-based indexing
 * printf ("DynCompRow:: Zero diagonal member\n");
 * exit(1);
 * }
 *
 *#else
 * std::map<int, double>::iterator pos;
 * for (i = 0; i < nRows; i++) {// row loop
 * pos  = this->rows[i]->find(i);
 * if (pos != this->rows[i-1]->end())
 * diag_rowptr_(i) = *pos->first;
 * else {
 * printf ("DynCompRow:: Zero diagonal member\n");
 * exit(1);
 * }
 * }
 *#endif
 *
 * // FACTOR MATRIX //
 *#ifndef DynCompRow_USE_STL_SETS
 *
 * for (i=1; i< nRows; i++) { // loop  over rows
 * for (k=0; k < (diag_rowptr_(i)-1); k++) { // loop 1,...,i-1 for (i,k) \in NZ(A)
 * // initialize k-th row indexes
 * krow = colind_[i]->at(k+1);
 * for (kk=1;  kk<= rows_[krow]->giveSize(); kk++)
 *  iw(colind_[krow]->at(kk)) = kk;
 * multiplier = (rows_[i]->at(k+1) /= rows_[krow]->at(diag_rowptr_(krow)));
 * for (j=k+1; j < rows_[i]->giveSize(); j++) { // loop over k+1,...,n for (i,j) \in NZ(A)
 *  colk = iw (colind_[i]->at(j+1));   // column position of a(i,j) at row k
 *  if (colk) rows_[i]->at(j+1) -= multiplier * rows_[krow]->at(colk); // aij=aij-aik*akj
 * }
 * // Refresh all iw enries to zero
 * for (kk=1;  kk<= rows_[krow]->giveSize(); kk++)
 *  iw(colind_[krow]->at(kk)) = 0;
 * }
 * }
 *
 *#else
 * NOT IMPLEMENTED NOW
 *#endif
 * }
 */


//#define ILU_0
//#define ILU_DROP_TOL 1.e-8
#define ILU_ROW_CHUNK 10
//#define ILU_PART_FILL 5

void
DynCompRow :: ILUPYourself(int part_fill, double drop_tol)
{
    int i, ii, j, jcol, k, kk, krow, ck;
    int end, curr;
    double multiplier, inorm, val;

    IntArray irw(nColumns), iw;
    FloatArray w;
    diag_rowptr_.resize(nRows);

#ifdef TIME_REPORT
    //clock_t tstart = clock();
    oofem_timeval tstart;
    getUtime(tstart);
#endif

    for ( i = 0; i < nRows; i++ ) { // row loop
        if ( ( diag_rowptr_(i) = giveColIndx(i, i) ) == 0 ) { // giveColIndx returns 1-based indexing
            OOFEM_ERROR("DynCompRow::ILUPYourself : zero value on diagonal");
        }
    }

    /* FACTOR MATRIX */

    for ( i = 1; i < nRows; i++ ) { // loop  over rows
        inorm = 0.0;
        for ( ii = 1; ii <= rows_ [ i ]->giveSize(); ii++ ) {
            val = rows_ [ i ]->at(ii);
            inorm += val * val;
        }

        inorm = sqrt(inorm);

        w.resize(rows_ [ i ]->giveSize(), ILU_ROW_CHUNK);
        iw.resize(rows_ [ i ]->giveSize(), ILU_ROW_CHUNK);
        for ( kk = 1; kk <= rows_ [ i ]->giveSize(); kk++ ) {
            irw( colind_ [ i ]->at(kk) ) = kk;
            iw(kk - 1) = colind_ [ i ]->at(kk);
            w(kk - 1) = rows_ [ i ]->at(kk);
        }

        //for (k=0; k < (diag_rowptr_(i)-1); k++) { // loop 1,...,i-1 for (i,k) \in NZ(A)
        k = 0;
        while ( iw.at(k + 1) < i ) {
            // initialize k-th row indexes
            krow = iw.at(k + 1);
            //multiplier = (rows_[i]->at(k+1) /= rows_[krow]->at(diag_rowptr_(krow)));
            multiplier = ( w.at(k + 1) /= rows_ [ krow ]->at( diag_rowptr_(krow) ) );


#ifndef ILU_0
            // first dropping rule for aik
            if ( fabs(multiplier) >= drop_tol * inorm )
#endif
            { // first drop rule
                for ( j = 0; j < colind_ [ krow ]->giveSize(); j++ ) {
                    jcol = colind_ [ krow ]->at(j + 1);
                    if ( jcol > krow ) {
                        if ( irw(jcol) ) {
                            //rows_[i]->at(irw(jcol)) -= multiplier*rows_[krow]->at(j+1);
                            w.at( irw(jcol) ) -= multiplier * rows_ [ krow ]->at(j + 1);
                        } else {
#ifndef ILU_0
                            // insert new entry
                            int newsize = w.giveSize() + 1;
                            w.resize(newsize, ILU_ROW_CHUNK);
                            iw.resize(newsize, ILU_ROW_CHUNK);

                            iw.at(newsize) = jcol;
                            w.at(newsize) = -multiplier * rows_ [ krow ]->at(j + 1);
                            irw(jcol) = newsize;
#endif

                            /*
                             * ipos = insertColInRow (i,jcol) ;
                             * for (kk=ipos+1;  kk<= rows_[i]->giveSize(); kk++)
                             * irw(colind_[i]->at(kk))++;
                             *
                             * ipos = insertColInRow (i,jcol) ;
                             * rows_[i]->at(ipos) = -multiplier*rows_[krow]->at(j+1);
                             * irw(jcol) = ipos;
                             * if (jcol < i) diag_rowptr_(i)++;
                             */
                        }
                    }
                }
            }

            // scan iw to find closest index to krow
            ck = nColumns + 1;
            for ( kk = 0; kk < iw.giveSize(); kk++ ) {
                if ( ( ( iw(kk) - krow ) > 0 ) && ( ( iw(kk) - krow ) < ( ck - krow ) ) ) {
                    ck = iw(kk);
                }
            }

            k = irw(ck) - 1;
        }

#ifndef ILU_0

        end = iw.giveSize();
        curr = 1;
        // second drop rule
        while ( curr <= end ) {
            if ( ( fabs( w.at(curr) ) < drop_tol * inorm ) && ( iw.at(curr) != i ) ) {
                // remove entry
                w.at(curr) = w.at(end);
                irw( iw.at(curr) ) = 0;
                iw.at(curr) = iw.at(end);
                if ( curr != end ) {
                    irw( iw.at(curr) ) = curr;
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
        int lsizeLimit = diag_rowptr_(i) - 1;
        int usizeLimit = rows_ [ i ]->giveSize() - lsizeLimit;

        lsizeLimit += part_fill;
        usizeLimit += part_fill;

        int lnums = 0;
        int unums = 0;
        count = 0;
        for ( kk = 1; kk <= iw.giveSize(); kk++ ) {
            if ( iw.at(kk) < i ) { // lpart
                if ( ++lnums > lsizeLimit ) {
                    irw( iw.at(kk) ) = 0;
                } else {
                    count++;
                }
            } else if ( iw.at(kk) > i ) { // upart
                if ( ++unums > usizeLimit ) {
                    irw( iw.at(kk) ) = 0;
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
        rows_ [ i ]->resize(count);
        colind_ [ i ]->resize(count);

        int icount = 1;
        int kki, indx, idist, previndx = -1;
        int kkend = iw.giveSize();

        for ( kk = 1; kk <= count; kk++ ) {
            idist = nColumns + 2;
            indx = 0;

            for ( kki = 1; kki <= kkend; kki++ ) {
                if ( ( irw( iw.at(kki) ) != 0 ) && ( iw.at(kki) > previndx ) && ( ( iw.at(kki) - previndx ) < idist ) ) {
                    idist = ( iw.at(kki) - previndx );
                    indx = kki;
                }
            }

            if ( indx == 0 ) {
                OOFEM_ERROR("DynCompRow::ILUPYourself : internal error");
            }

            previndx = iw.at(indx);
            rows_ [ i ]->at(icount) = w.at(indx);
            colind_ [ i ]->at(icount) = iw.at(indx);
            if ( colind_ [ i ]->at(icount) == i ) {
                diag_rowptr_(i) = icount;
            }

            icount++;


            // exclude the indx entry from search by moving it to the end of list
            irw( iw.at(indx) ) = 0;
            iw.at(indx) = iw.at(kkend);
            w.at(indx)  = w.at(kkend);
            if ( irw( iw.at(indx) ) != 0 ) {
                irw( iw.at(indx) ) = indx;
            }

            kkend--;

            // exclude the indx entry from search by moving it to the end of list
            //swap = irw(iw.at(indx)); irw(iw.at(indx)) = irw(iw.at(kkend)); irw(iw.at(kkend)) = swap;
            //swap = iw.at(indx); iw.at(indx) = iw.at(kkend); iw.at(kkend) = swap;
            //dswap= w.at(indx); w.at(indx) = w.at(kkend); w.at(kkend) = dswap;

            //kkend--;
        }


        /*
         * int icount = 1;
         * for (kk=1;  kk<= nColumns; kk++) {
         * if ( irw.at(kk) > 0 ) {
         * rows_[i]->at(icount) = w.at(abs(irw(kk-1)));
         * colind_[i]->at(icount) = iw.at(abs(irw(kk-1)));
         * if (colind_[i]->at(icount) == i) diag_rowptr_(i) = icount;
         * icount++;
         * }
         * }
         */
        if ( ( icount - count ) != 1 ) {
            OOFEM_ERROR4("DynCompRow::ILUPYourself : %d - row errorr (%d,%d)\n", i, icount, count);
        }

        //Refresh all iw enries to zero
        for ( kk = 1; kk <= iw.giveSize(); kk++ ) {
            irw( iw.at(kk) ) = 0;
        }

        //irw.zero();
    }

#ifdef TIME_REPORT
    oofem_timeval ut;
    getRelativeUtime(ut, tstart);
    OOFEM_LOG_DEBUG( "\nILUT(%d,%e): user time consumed by factorization: %.2fs\n", part_fill, drop_tol, ( double ) ( ut.tv_sec + ut.tv_usec / ( double ) OOFEM_USEC_LIM ) );
#endif

    // increment version
    this->version++;
}



void
DynCompRow :: ILUPsolve(const FloatArray &x, FloatArray &y) const
{
    int M = x.giveSize();
    double r;
    FloatArray work(M);

    int i, t;

    y.resize(M);

    work.zero();
    // solve Lw=x
    for ( i = 0; i < M; i++ ) {
        r = x(i);
        for ( t = 0; t < ( diag_rowptr_(i) - 1 ); t++ ) {
            r -= rows_ [ i ]->at(t + 1) * work( colind_ [ i ]->at(t + 1) );
        }

        work(i) = r;
    }

    y.zero();
    // solve Uy=w
    for ( i = M - 1; i >= 0; i-- ) {
        r = work(i);
        for ( t = diag_rowptr_(i); t < rows_ [ i ]->giveSize(); t++ ) {
            r -= rows_ [ i ]->at(t + 1) * y( colind_ [ i ]->at(t + 1) );
        }

        y(i) = r / rows_ [ i ]->at( diag_rowptr_(i) );
    }

    //return y;
}


void
DynCompRow :: ILUPtrans_solve(const FloatArray &x, FloatArray &y) const
{
    int M = x.giveSize();
    FloatArray work(M);

    int i, t;

    y.resize(M);
    //work.zero();
    // solve for U^Tw = x
    for ( i = 0; i < M; i++ ) {
        work(i) = ( x(i) + work(i) ) / rows_ [ i ]->at( diag_rowptr_(i) );
        for ( t = diag_rowptr_(i); t < rows_ [ i ]->giveSize(); t++ ) {
            work( colind_ [ i ]->at(t + 1) ) -= rows_ [ i ]->at(t + 1) * work(i);
        }
    }

    // solve for L^T y = w
    for ( i = M - 1; i >= 0; i-- ) {
        y(i) = work(i);
        for ( t = 1; t < diag_rowptr_(i); t++ ) {
            work( colind_ [ i ]->at(t) ) -= rows_ [ i ]->at(t) * y(i);
        }
    }

    //return y;
}


/*
 * void
 * DynCompRow::qsortRow (IntArray& ind, IntArray& ir, FloatArray& val, int l, int r)
 * {
 * if (r<=l) return;
 * int i = qsortRowPartition (ind, ir, val, l, r);
 * qsortRow (ind, ir, val, l, i-1);
 * qsortRow (ind, ir, val, i+1, r);
 * }
 *
 *
 * int
 * DynCompRow::qsortRowPartition (IntArray& ind, IntArray& ir, FloatArray& val, int l, int r)
 * {
 * int i=l-1, j=r;
 * double v = fabs(val(r));
 * int swap;
 * double dswap;
 *
 * for (;;) {
 * while (( fabs(val(++i)) >  v ));
 * while (( v > fabs(val(--j)))) if (j==l) break;
 * if (i >= j) break;
 * swap = ir(ind(i)); ir(ind(i)) = ir(ind(j)); ir(ind(j)) = swap;
 * swap = ind(i); ind(i) = ind(j); ind(j) = swap;
 * dswap= val(i); val(i) = val(j); val(j) = dswap;
 * }
 * swap = ir(ind(i)); ir(ind(i)) = ir(ind(r)); ir(ind(r)) = swap;
 * swap = ind(i); ind(i) = ind(r); ind(r) = swap;
 * dswap= val(i); val(i) = val(r); val(r) = dswap;
 *
 * return i;
 * }
 */

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
    double v = fabs( val(r) );
    int swap;
    double dswap;

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

        swap = ir( ind(i) );
        ir( ind(i) ) = ir( ind(j) );
        ir( ind(j) ) = swap;
        swap = ind(i);
        ind(i) = ind(j);
        ind(j) = swap;
        dswap = val(i);
        val(i) = val(j);
        val(j) = dswap;
    }

    swap = ir( ind(i) );
    ir( ind(i) ) = ir( ind(r) );
    ir( ind(r) ) = swap;
    swap = ind(i);
    ind(i) = ind(r);
    ind(r) = swap;
    dswap = val(i);
    val(i) = val(r);
    val(r) = dswap;

    return i;
}
} // end namespace oofem
