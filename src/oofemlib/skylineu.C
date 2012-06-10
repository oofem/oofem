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

#include "skylineu.h"
#include "rowcol.h"
#include "flotmtrx.h"
#include "intarray.h"
#include "domain.h"
#include "engngm.h"
#include "element.h"
#include "mathfem.h"
#include "verbose.h"
#include "error.h"

#ifdef TIME_REPORT
 #include "clock.h"
#endif

namespace oofem {
SkylineUnsym :: SkylineUnsym(int n) : SparseMtrx(n, n)
    // Constructor. Creates an empty skyline unsymmetric.

{
    size         = n;
    rowColumns   = NULL;
    isFactorized = false;
}

SkylineUnsym :: SkylineUnsym() : SparseMtrx()
{
    // Constructor. Creates a skyline of size 0.
    // nRows = nColumns = 0;  // set by SparseMtrx constructor
    size         = 0;
    rowColumns   = NULL;
    isFactorized = false;
}

SkylineUnsym :: ~SkylineUnsym()
// Destructor.

{
    RowColumn **p;
    int i;

    if ( size ) {
        i = size;
        p = rowColumns;
        while ( i-- ) {
            delete *p++;
        }

        delete rowColumns;
    }
}

void
SkylineUnsym :: toFloatMatrix(FloatMatrix &answer) const
{
    int i, j, start;

    answer.resize(size, size);
    answer.zero();
    for ( j = 1; j <= size; j++ ) {
        start = this->giveRowColumn(j)->giveStart();
        for ( i = start; i <= j; i++ ) {
            answer.at(i, j) = this->at(i, j);
            answer.at(j, i) = this->at(j, i);
        }
    }
}


double &
SkylineUnsym :: at(int i, int j)
// Returns the (i,j) coefficient of the receiver. Not very efficient.

{
    // increment version
    this->version++;

    if ( i < j ) {
        return this->giveRowColumn(j)->atU(i);
    } else if ( i > j ) {
        return this->giveRowColumn(i)->atL(j);
    } else {
        return this->giveRowColumn(i)->atDiag();
    }
}

double
SkylineUnsym :: at(int i, int j) const
// Returns the (i,j) coefficient of the receiver. Not very efficient.

{
    if ( i < j ) {
        return this->giveRowColumn(j)->atU(i);
    } else if ( i > j ) {
        return this->giveRowColumn(i)->atL(j);
    } else {
        return this->giveRowColumn(i)->atDiag();
    }
}


int
SkylineUnsym :: assemble(const IntArray &loc, const FloatMatrix &mat)
// Assembles the elemental matrix 'k' to the receiver, using 'loc' as a
// location array. The values in k corresponding to a zero coefficient
// in loc are not assembled.
// Warning : k is not supposed to be an instance of DiagonalMatrix.

{
    int i, j, ii, jj, dim;
    //   RowColumn* rowColumnJJ ; // 13 November 1996 - not needed anymore

#  ifdef DEBUG
    dim = mat.giveNumberOfRows();
    if ( dim != loc.giveSize() ) {
        OOFEM_ERROR("SkylineUnsym :: assemble : dimension of 'k' and 'loc' mismatch");
    }

    this->checkSizeTowards(loc);
#  endif

    // added to make it work for nonlocal model !!!!!!!!!!!!!!!!!!!!!!!!!!
    // checkSizeTowards(loc) ;


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


int
SkylineUnsym :: assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat)
// Assembles the elemental matrix 'k' to the receiver, using 'loc1' as a
// location array for rows and 'loc2' as a location array for columns.
// The values in 'k' corresponding to a zero coefficient
// in loc are not assembled.
// Warning : k is not supposed to be an instance of DiagonalMatrix.

{
    int i, j, ii, jj, dim1, dim2;

    this->checkSizeTowards(rloc, cloc);

    dim1 = mat.giveNumberOfRows();
    dim2 = mat.giveNumberOfColumns();
    for ( i = 1; i <= dim1; i++ ) {
        ii = rloc.at(i);
        if ( ii ) {
            for ( j = 1; j <= dim2; j++ ) {
                jj = cloc.at(j);
                if ( jj ) {
                    this->at(ii, jj) += mat.at(i, j);
                }
            }
        }
    }

    // increment version
    this->version++;

    return 1;
}

void
SkylineUnsym :: checkSizeTowards(const IntArray &rloc, const IntArray &cloc)
// Increases the number of columns of the receiver if 'rloc, cloc' points to
// out-of-range columns.
{
    int maxCol, i, j, ii, jj, dim1, dim2;

    maxCol = 0;                         // the largest coefficient in 'loc'
    dim1    = rloc.giveSize();
    dim2    = cloc.giveSize();

    for ( i = 1; i <= dim1; i++ ) {
        maxCol = max( maxCol, rloc.at(i) );
    }

    for ( j = 1; j <= dim2; j++ ) {
        maxCol = max( maxCol, cloc.at(j) );
    }

    if ( maxCol > size ) {              // enlarge the matrix
        growTo(maxCol);
    }

    for ( i = 1; i <= dim1; i++ ) {
        ii = rloc.at(i);
        if ( ii ) {
            giveRowColumn(ii)->checkSizeTowards(cloc);
        }
    }

    for ( j = 1; j <= dim2; j++ ) {
        jj = cloc.at(j);
        if ( jj ) {
            giveRowColumn(jj)->checkSizeTowards(rloc);
        }
    }
}



int
SkylineUnsym :: buildInternalStructure(EngngModel *eModel, int di, EquationID ut, const UnknownNumberingScheme &s)
{
    // Instanciates the profile of the receiver and initializes all coefficients to zero.
    // Warning : case diagonal (lumped) matrix to expected.

    IntArray loc;
    int n, i, first, ii, nonlocal;
    Domain *domain = eModel->giveDomain(di);
    int neq = eModel->giveNumberOfDomainEquations(di, ut);
    int nelem = domain->giveNumberOfElements();
    Element *elem;
    //nonlocal = aDomain -> giveAlgorithm() -> useNonlocalStiffness();
    nonlocal = eModel->useNonlocalStiffnessOption();

    // clear receiver if exist

    if ( size ) {
        RowColumn **_p;
        int _i;
        _i = size;
        _p = rowColumns;
        while ( _i-- ) {
            delete *_p++;
        }

        delete rowColumns;
        rowColumns = NULL;
        size = 0;
    }

    this->growTo(neq); // from now on, size = MaxIndex

    // Set up the array with indices of first nonzero elements in each row
    IntArray firstIndex(neq);
    for ( i = 1; i <= neq; i++ ) {
        firstIndex.at(i) = i;
    }

    for ( n = 1; n <= nelem; n++ ) {
        elem = domain->giveElement(n);
        //elem -> giveLocationArray (loc) ;

        if ( !nonlocal ) {
            elem->giveLocationArray(loc, ut, s);
        }
        //else ((StructuralElement*)elem) -> giveNonlocalLocationArray(loc) ;
        else {
            elem->giveLocationArray(loc, ut, s);
        }

        //    Find 'first', the smallest positive number in LocArray
        first = neq;
        for ( i = 1; i <= loc.giveSize(); i++ ) {
            ii = loc.at(i);
            if ( ii && ii < first ) {
                first = ii;
            }
        }

        //    Make sure that the FirstIndex is not larger than 'first'
        for ( i = 1; i <= loc.giveSize(); i++ ) {
            ii = loc.at(i);
            if ( ii && ( first < firstIndex.at(ii) ) ) {
                firstIndex.at(ii) = first;
            }
        }
    }

    // Enlarge the rowcolumns
    for ( i = 1; i <= neq; i++ ) {
        this->giveRowColumn(i)->growTo( firstIndex.at(i) );
    }

    // Print out basic information.
    this->printStatistics();

    nRows = nColumns = neq;
    // increment version
    this->version++;

    return true;
}


int SkylineUnsym :: setInternalStructure(IntArray *adr1)
{
    // allocates and built structure according to given
    // array of maximal column heights
    //
    IntArray mht, firstIndex;
    int i;
    int n = adr1->giveSize();
    this->growTo(n - 1);
    size = n - 1; // check

    // build first index array
    mht.resize(n - 1);
    mht.zero();
    firstIndex.resize(n - 1);
    firstIndex.zero();
    for ( i = 1; i <= n - 1; i++ ) {
        mht.at(i) = adr1->at(i + 1) - adr1->at(i);
        firstIndex.at(i) = i - mht.at(i) + 1;
    }

    // Enlarge the rowcolumns
    for ( i = 1; i <= n - 1; i++ ) {
        this->giveRowColumn(i)->growTo( firstIndex.at(i) );
    }

    // Print out basic information.
    this->printStatistics();

    nRows = nColumns = n - 1;
    // increment version
    this->version++;

    return true;
}



void
SkylineUnsym :: checkSizeTowards(const IntArray &loc)
// Increases the number of columns of the receiver if 'loc' points to
// out-of-range columns.
{
    int i, ii, dim, maxCol;

    maxCol = 0;                         // the largest coefficient in 'loc'
    dim    = loc.giveSize();
    for ( i = 1; i <= dim; i++ ) {
        maxCol = max( maxCol, loc.at(i) );
    }

    if ( maxCol > size ) {              // enlarge the matrix
        growTo(maxCol);
    }

    for ( i = 1; i <= dim; i++ ) {
        ii = loc.at(i);
        if ( ii ) {
            giveRowColumn(ii)->checkSizeTowards(loc);
        }
    }
}

/*
 * int
 * SkylineUnsym :: computeNumberNegativeEigenvalue()
 * // Compute the number of negative eigenvalue, equivalent to the number of
 * // negative diagonal terms.
 *
 * {
 * int j, count;
 *
 * if (! isFactorized) factorized();
 *
 * count = 0;
 * for (j=1 ; j<=size ; j++) {
 * if (at(j,j) <= 0.)
 * count = count + 1;
 * }
 * return count;
 * }
 */

SparseMtrx *
SkylineUnsym :: factorized()
// Returns the receiver in L.D.U factorized form. From Golub & van Loan,
// 1rst edition, pp 83-84.

{
    RowColumn *rowColumnK, *rowColumnI;
    FloatArray r, w;
    double diag;
    int k, i, p, start, startK, startI;
#ifdef TIME_REPORT
    //clock_t tstart = clock();
    oofem_timeval tstart;
    getUtime(tstart);
#endif

    if ( isFactorized ) {
        return this;
    }

    if ( !size ) {
        OOFEM_WARNING("SkylineUnsym::factorized : null-sized matrix factorized");
        isFactorized = 1;
        return this;
    }

    for ( k = 1; k <= size; k++ ) {
        rowColumnK = this->giveRowColumn(k);
        startK   = rowColumnK->giveStart();

        // compute vectors r and w
        r.resize(k - 1);
        r.zero();
        w.resize(k - 1);
        w.zero();
        for ( p = startK; p < k; p++ ) {
            diag     = this->giveRowColumn(p)->atDiag();
            r.at(p) = diag * rowColumnK->atU(p);
            w.at(p) = diag * rowColumnK->atL(p);
        }

        // compute diagonal coefficient of rowColumn k
        rowColumnK->atDiag() -= rowColumnK->dot(r, 'R', startK, k - 1);
        diag = rowColumnK->atDiag();

        // test pivot not too small
        if ( fabs(diag) < SkylineUnsym_TINY_PIVOT ) {
            rowColumnK->atDiag() = diag = SkylineUnsym_TINY_PIVOT;
            OOFEM_LOG_DEBUG("SkylineUnsym :: factorized: zero pivot %d artificially set to a small value", k);
        }

        // compute off-diagonal coefficients of rowColumns i>k
        for ( i = k + 1; i <= size; i++ ) {
            rowColumnI = giveRowColumn(i);
            startI   = rowColumnI->giveStart();
            if ( startI <= k ) {
                start = max(startI, startK);
                rowColumnI->atL(k) -= rowColumnI->dot(r, 'R', start, k - 1);
                rowColumnI->atL(k) /= diag;
                rowColumnI->atU(k) -= rowColumnI->dot(w, 'C', start, k - 1);
                rowColumnI->atU(k) /= diag;
            }
        }
    }

    isFactorized = true;

#ifdef TIME_REPORT
    //printf ("\nSkylineU info: user time consumed by factorization: %.2gs", (clock()-tstart)/(double)CLOCKS_PER_SEC);
    oofem_timeval ut;
    getRelativeUtime(ut, tstart);
    OOFEM_LOG_DEBUG( "SkylineU info: user time consumed by factorization: %.2fs\n", ( double ) ( ut.tv_sec + ut.tv_usec / ( double ) OOFEM_USEC_LIM ) );
#endif
    // increment version
    //this->version++;

    return this;
}



FloatArray *
SkylineUnsym :: backSubstitutionWith(FloatArray &y) const
// Returns the solution x of the system U.x = y , where U is the upper
// half of the receiver. Note : x overwrites y.

{
    int i, k, start;
    double yK, diag;
    RowColumn *rowColumnK;


    // forwardReductionWith
    // if (! isFactorized) this -> factorized();

    if ( !size ) {
        return & y;                               // null size system
    }

    if ( y.giveSize() != size ) {
        OOFEM_ERROR("SkylineUnsym::backSubstitutionWith : size mismatch");
    }

    for ( k = 1; k <= size; k++ ) {
        rowColumnK  = this->giveRowColumn(k);
        start     = rowColumnK->giveStart();
        y.at(k) -= rowColumnK->dot(y, 'R', start, k - 1);
    }

    // diagonalScaling
    /*  uprava diagonalniho prvku  */
    for ( k = 1; k <= size; k++ ) {
        diag = this->giveRowColumn(k)->atDiag();
#     ifdef DEBUG
        if ( fabs(diag) < SkylineUnsym_TINY_PIVOT ) {
            OOFEM_ERROR2("SkylineUnsym::backSubstitutionWith : diagScaling: pivot %d is small", k);
        }

#     endif
        y.at(k) /= diag;
    }


    for ( k = size; k > 0; k-- ) {
        rowColumnK = this->giveRowColumn(k);
        yK       = y.at(k);
        start    = rowColumnK->giveStart();
        for ( i = start; i < k; i++ ) {
            y.at(i) -= rowColumnK->atU(i) * yK;
        }
    }

    return & y;
}

SparseMtrx *
SkylineUnsym :: GiveCopy() const
{
    int i;
    RowColumn **newRowColumns, *rowColumni;
    newRowColumns = new RowColumn * [ size ];
    SkylineUnsym *answer;

    for ( i = 1; i <= size; i++ ) {
        rowColumni  = this->giveRowColumn(i);
        if ( rowColumni ) {
            newRowColumns [ i - 1 ] = rowColumni->GiveCopy();
        } else {
            newRowColumns [ i - 1 ] = NULL;
        }
    }

    answer = new SkylineUnsym(newRowColumns, this->size, this->isFactorized);
    return answer;
}

void
SkylineUnsym :: times(const FloatArray &x, FloatArray &answer) const
{
    int starti, i, j;
    RowColumn *rowColumni;
    //
    // first check sizes
    //
    if ( this->size != x.giveSize() ) {
        OOFEM_ERROR("SkylineUnsym::times : size mismatch");
    }

    answer.resize(this->size);
    answer.zero();

    for ( i = 1; i <= size; i++ ) {
        rowColumni  = this->giveRowColumn(i);
        starti     = rowColumni->giveStart();
        answer.at(i) += rowColumni->dot(x, 'R', starti, i - 1);
        answer.at(i) += rowColumni->atDiag() * x.at(i);

        for ( j = starti; j <= i - 1; j++ ) {
            answer.at(j) += rowColumni->atU(j) * x.at(i);
        }
    }
}

void SkylineUnsym :: times(double x)
{
    int starti, i, j;
    RowColumn *rowColumni;

    for ( i = 1; i <= size; i++ ) {
        rowColumni  = this->giveRowColumn(i);
        starti     = rowColumni->giveStart();
        for ( j = starti; j <= i - 1; j++ ) {
            rowColumni->atU(j) *= x;
            rowColumni->atL(j) *= x;
        }

        rowColumni->atDiag() *= x;
    }

    // increment version
    this->version++;
}


RowColumn *
SkylineUnsym :: giveRowColumn(int j) const
// Returns the j-th rowColumn of the receiver. Creates it if necessary.

{
    if ( size < j ) {
        // this -> growTo(j) ;
        OOFEM_ERROR("SkylineUnsym::giveRowColumn : size mismatch");
    }

    if ( !rowColumns [ j - 1 ] ) {
        rowColumns [ j - 1 ] = new RowColumn(j, j);
    }

    return rowColumns [ j - 1 ];
}


void
SkylineUnsym :: growTo(int n)
// Enlarges the receiver to n columns.
{
    RowColumn **newRowColumns, **p1, **p2;
    int i;

    if ( n == size ) {
        return;
    }

#  ifdef DEBUG
    else if ( n <= size ) {
        OOFEM_ERROR3("SkylineUnsym::growTo : cannot grow from %d to %d", size, n);
    }
#  endif

    newRowColumns = new RowColumn * [ n ];
    p1          = rowColumns;
    p2          = newRowColumns;

    i = size;
    while ( i-- ) {
        * p2++ = * p1++;
    }

    i = n - size;
    while ( i-- ) {
        * p2++ = NULL;
    }

    if ( rowColumns ) {
        delete rowColumns;
    }

    rowColumns = newRowColumns;
    size     = n;
}

void
SkylineUnsym :: printStatistics() const
{
    int nelem = 0;
    for ( int i = 1; i <= size; i++ ) {
        nelem += this->giveRowColumn(i)->giveSize();
    }

    OOFEM_LOG_INFO("Skylineu info: neq is %d, nwk is %d\n", size, nelem);
}


void
SkylineUnsym :: writeToFile(const char* fname) const
{
    FILE *file = fopen(fname,"w");
    FloatMatrix copy;
    this->toFloatMatrix(copy);
    for ( int i = 1; i <= nRows; ++i ) {
        for ( int j = 1; j <= nColumns; ++j ) {
            fprintf(file, "%.16e ", copy.at(i, j) );
        }
        fprintf(file, "\n");
    }
    fclose(file);
}


void
SkylineUnsym :: printYourself() const
// Prints the receiver.
{
    FloatMatrix copy;

    this->toFloatMatrix(copy);
    copy.printYourself();
}

void
SkylineUnsym :: zero()
{
    // Returns the receiver with all coefficients set to zero.
    int j;

    for ( j = 1; j <= size; j++ ) {
        this->giveRowColumn(j)->zero();
    }

    isFactorized = false;

    // increment version
    this->version++;
}

/*
 * SkylineSym* SkylineUnsym :: giveSymmetricPart()
 * {
 * int i, j, rowcolsize, colsize;
 * SkylineSym* answer = new SkylineSym();
 * answer -> growTo(size);
 * for (i=1; i<=size; i++){
 * rowcolsize = rowColumns[i-1] -> giveSize();
 * colsize = (rowcolsize+1)/2;
 * answer -> giveColumn(i) -> growTo(colsize);
 * answer->at(i,i) = at(i,i);
 * j=i+1-colsize;
 * while(j<i){
 * answer->at(j,i) = (at(i,j)+at(j,i))/2.;
 * j++;
 * }
 * }
 * return answer;
 * }
 *
 *
 * int  SkylineUnsym :: computeNumberNegativePivotsOfSymPart()
 * {
 * SkylineSym* SymPart = giveSymmetricPart();
 * int answer = SymPart -> computeNumberNegativeEigenvalue();
 * delete SymPart;
 * return answer;
 * }
 *
 */

SkylineUnsym :: SkylineUnsym(RowColumn **newRowCol, int newSize, int isFact) : SparseMtrx(newSize, newSize)
{
    size         = newSize;
    rowColumns   = newRowCol;
    isFactorized = isFact;
}

void SkylineUnsym :: timesT(const FloatArray &x, FloatArray &answer) const
{
    int starti, i, j;
    RowColumn *rowColumni;

    // first check sizes
    if ( this->size != x.giveSize() ) {
        OOFEM_ERROR("SkylineUnsym::trans_mult : size mismatch");
    }

    answer.resize(this->size);
    answer.zero();

    for ( i = 1; i <= size; i++ ) {
        rowColumni  = this->giveRowColumn(i);
        starti     = rowColumni->giveStart();
        answer.at(i) += rowColumni->dot(x, 'C', starti, i - 1);
        answer.at(i) += rowColumni->atDiag() * x.at(i);

        for ( j = starti; j <= i - 1; j++ ) {
            answer.at(j) += rowColumni->atL(j) * x.at(i);
        }
    }
}
} // end namespace oofem
