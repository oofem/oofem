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

// Class CompCol

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
#include "flotarry.h"
#include "engngm.h"
#include "domain.h"
#include "element.h"

#include <set>

namespace oofem {
CompCol :: CompCol(void) : SparseMtrx(), val_(0), rowind_(0), colptr_(0), base_(0), nz_(0)
{
    dim_ [ 0 ] = 0;
    dim_ [ 1 ] = 0;
}


CompCol :: CompCol(int n) : SparseMtrx(n, n), val_(0), rowind_(0), colptr_(n), base_(0), nz_(0)
{
    dim_ [ 0 ] = n;
    dim_ [ 1 ] = n;
}


/*****************************/
/*  Copy constructor         */
/*****************************/

CompCol :: CompCol(const CompCol &S) : SparseMtrx(S.nRows, S.nColumns),
    val_(S.val_), rowind_(S.rowind_), colptr_(S.colptr_),
    base_(S.base_), nz_(S.nz_)
{
    dim_ [ 0 ] = S.dim_ [ 0 ];
    dim_ [ 1 ] = S.dim_ [ 1 ];
    nRows = S.nRows;
    nColumns = S.nColumns;
}



/***************************/
/* Assignment operator...  */
/***************************/

CompCol &CompCol :: operator = ( const CompCol & C )
{
    dim_ [ 0 ] = C.dim_ [ 0 ];
    dim_ [ 1 ] = C.dim_ [ 1 ];
    nRows   = C.nRows;
    nColumns = C.nColumns;

    base_   = C.base_;
    nz_     = C.nz_;
    val_    = C.val_;
    rowind_ = C.rowind_;
    colptr_ = C.colptr_;
    this->version = C.version;

    return * this;
}



SparseMtrx *CompCol :: GiveCopy() const
{
    CompCol *result = new CompCol(*this);
    return result;
}


void CompCol :: times(const FloatArray &x, FloatArray &answer) const
{
    int M = dim_ [ 0 ];
    int N = dim_ [ 1 ];

    //      Check for compatible dimensions:
    if ( x.giveSize() != N ) {
        OOFEM_ERROR("CompCol::times: Error in CompCol -- incompatible dimensions");
    }

    answer.resize(M);
    answer.zero();

    int j, t;
    double rhs;

    for ( j = 0; j < N; j++ ) {
        rhs = x(j);
        for ( t = colptr_(j); t < colptr_(j + 1); t++ ) {
            answer( rowind_(t) ) += val_(t) * rhs;
        }
    }
}

void CompCol :: times(double x)
{
    for ( int t = 0; t < nz_; t++ ) {
        val_(t) *= x;
    }

    // increment version
    this->version++;
}

int CompCol :: buildInternalStructure(EngngModel *eModel, int di, EquationID ut, const UnknownNumberingScheme &s)
{
    /*
     * IntArray  loc;
     * Domain* domain = eModel->giveDomain(di);
     * int neq = eModel -> giveNumberOfDomainEquations (di);
     * int nelem = domain -> giveNumberOfElements() ;
     * int i,ii,j,jj,n, indx;
     * Element* elem;
     * // allocation map
     * char* map = new char[neq*neq];
     * if (map == NULL) {
     * printf ("CompCol::buildInternalStructure - map creation failed");
     * exit (1);
     * }
     *
     * for (i=0; i<neq*neq; i++)
     * map[i]=0;
     *
     *
     * this->nz_ = 0;
     *
     * for (n=1 ; n<=nelem ; n++) {
     * elem = domain -> giveElement(n);
     * elem -> giveLocationArray (loc) ;
     *
     * for (i=1 ; i <= loc.giveSize() ; i++) {
     * if ((ii = loc.at(i))) {
     *  for (j=1; j <= loc.giveSize() ; j++) {
     *   if ((jj=loc.at(j))) {
     *    if (map[(ii-1)*neq+jj-1] == 0) {
     *     map[(ii-1)*neq+jj-1] = 1;
     *     this->nz_ ++;
     *    }
     *   }
     *  }
     * }
     * }
     * }
     *
     * rowind_.resize (nz_);
     * colptr_.resize (neq+1);
     * indx = 0;
     * for (j=0; j<neq; j++) { // column loop
     * colptr_(j) = indx;
     * for (i=0; i<neq; i++) { // row loop
     *  if (map[i*neq+j]) {
     *   rowind_(indx) = i;
     *   indx++;
     *  }
     * }
     * }
     * colptr_(neq) = indx;
     *
     * // delete map
     * delete (map);
     *
     * // allocate value array
     * val_.resize(nz_);
     * val_.zero();
     *
     * printf ("\nCompCol info: neq is %d, nwk is %d\n",neq,nz_);
     *
     * dim_[0] = dim_[1] = nColumns = nRows = neq;
     *
     * // increment version
     * this->version++;
     *
     * return true;
     */
    IntArray loc;
    Domain *domain = eModel->giveDomain(di);
    int neq = eModel->giveNumberOfDomainEquations(di, ut);
    int nelem = domain->giveNumberOfElements();
    int i, ii, j, jj, n, indx;
    Element *elem;
    // allocation map
    std :: vector< std :: set< int > >columns(neq);
    /*
     * std::set<int> **columns = new std::set<int>*[neq];
     * for (j=0; j<neq; j++) {
     * columns[j] = new std::set<int>;
     * }
     */

    this->nz_ = 0;

    for ( n = 1; n <= nelem; n++ ) {
        elem = domain->giveElement(n);
        elem->giveLocationArray(loc, ut, s);

        for ( i = 1; i <= loc.giveSize(); i++ ) {
            if ( ( ii = loc.at(i) ) ) {
                for ( j = 1; j <= loc.giveSize(); j++ ) {
                    if ( ( jj = loc.at(j) ) ) {
                        columns [ jj - 1 ].insert(ii - 1);
                    }
                }
            }
        }
    }

    for ( i = 0; i < neq; i++ ) {
        this->nz_ += columns [ i ].size();
    }

    rowind_.resize(nz_);
    colptr_.resize(neq + 1);
    indx = 0;

    std :: set< int > :: iterator pos;
    for ( j = 0; j < neq; j++ ) { // column loop
        colptr_(j) = indx;
        for ( pos = columns [ j ].begin(); pos != columns [ j ].end(); ++pos ) { // row loop
            rowind_(indx++) = * pos;
        }
    }

    colptr_(neq) = indx;

    /*
     * // delete map
     * for (i=0; i< neq; i++) {columns[i]->clear(); delete columns[i];}
     * delete columns;
     */

    // allocate value array
    val_.resize(nz_);
    val_.zero();

    OOFEM_LOG_DEBUG("CompCol info: neq is %d, nwk is %d\n", neq, nz_);

    dim_ [ 0 ] = dim_ [ 1 ] = nColumns = nRows = neq;

    // increment version
    this->version++;

    return true;
}


int CompCol :: assemble(const IntArray &loc, const FloatMatrix &mat)
{
    int i, j, ii, jj, dim;
    //   RowColumn* rowColumnJJ ; // 13 November 1996 - not needed anymore

#  ifdef DEBUG
    dim = mat.giveNumberOfRows();
    if ( dim != loc.giveSize() ) {
        OOFEM_ERROR("CompCol::assemble : dimension of 'k' and 'loc' mismatch");
    }

    //this -> checkSizeTowards(loc) ;
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

int CompCol :: assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat)
{
    int i, j, ii, jj, dim1, dim2;

    // this->checkSizeTowards(rloc, cloc);

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

void CompCol :: zero()
{
    for ( int t = 0; t < nz_; t++ ) {
        val_(t) = 0.0;
    }

    // increment version
    this->version++;
}


void CompCol :: toFloatMatrix(FloatMatrix &answer) const
{ }

void CompCol :: printYourself() const
{ }

/*********************/
/*   Array access    */
/*********************/

double &CompCol :: at(int i, int j)
{
    // increment version
    this->version++;

    for ( int t = colptr_(j - 1); t < colptr_(j); t++ ) {
        if ( rowind_(t) == i - 1 ) {
            return val_(t);
        }
    }

    OOFEM_ERROR3("CompCol::operator(): Array accessing exception -- (%d,%d) out of bounds", i, j);
    return val_(0); // return to suppress compiler warning message
}


double CompCol :: at(int i, int j) const
{
    for ( int t = colptr_(j - 1); t < colptr_(j); t++ ) {
        if ( rowind_(t) == i - 1 ) {
            return val_(t);
        }
    }

    if ( i <= dim_ [ 0 ] && j <= dim_ [ 1 ] ) {
        return 0.0;
    } else {
        OOFEM_ERROR3("CompCol::operator(): Array accessing exception -- (%d,%d) out of bounds", i, j);
        return ( 0 ); // return to suppress compiler warning message
    }
}

double CompCol :: operator() (int i, int j)  const
{
    for ( int t = colptr_(j); t < colptr_(j + 1); t++ ) {
        if ( rowind_(t) == i ) {
            return val_(t);
        }
    }

    if ( i < dim_ [ 0 ] && j < dim_ [ 1 ] ) {
        return 0.0;
    } else {
        OOFEM_ERROR3("CompCol::operator(): Array accessing exception -- (%d,%d) out of bounds", i, j);
        return ( 0 ); // return to suppress compiler warning message
    }
}

double &CompCol :: operator() (int i, int j)
{
    // increment version
    this->version++;

    for ( int t = colptr_(j); t < colptr_(j + 1); t++ ) {
        if ( rowind_(t) == i ) {
            return val_(t);
        }
    }

    OOFEM_ERROR3("CompCol::operator(): Array element (%d,%d) not in sparse structure -- cannot assign", i, j);
    return val_(0); // return to suppress compiler warning message
}

void CompCol :: timesT(const FloatArray &x, FloatArray &answer) const
{
    int M = dim_ [ 0 ];
    int N = dim_ [ 1 ];

    //      Check for compatible dimensions:
    if ( x.giveSize() != M ) {
        OOFEM_ERROR("CompCol::timesT: Error in CompCol -- incompatible dimensions");
    }

    answer.resize(N);
    answer.zero();
    int i, t;
    double r;

    for ( i = 0; i < N; i++ ) {
        r = 0.0;
        for ( t = colptr_(i); t < colptr_(i + 1); t++ ) {
            r += val_(t) * x( rowind_(t) );
        }

        answer(i) = r;
    }
}
} // end namespace oofem
