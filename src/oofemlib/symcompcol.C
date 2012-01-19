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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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
#include "flotarry.h"
#include "engngm.h"
#include "domain.h"
#include "element.h"
#ifndef __MAKEDEPEND
 #include <set>
#endif

namespace oofem {
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

int SymCompCol :: buildInternalStructure(EngngModel *eModel, int di, EquationID ut, const UnknownNumberingScheme &s)
{
    /*
     * IntArray  loc;
     * Domain* domain = eModel->giveDomain(di);
     * int neq = eModel -> giveNumberOfDomainEquations (di);
     * int nelem = domain -> giveNumberOfElements() ;
     * int i,ii,j,jj,n, indx;
     * Element* elem;
     * // allocation map
     * int nmap = neq*neq-neq*(neq-1)/2;
     * char* map = new char[nmap];
     * if (map == NULL) {
     * printf ("CompCol::buildInternalStructure - map creation failed");
     * exit (1);
     * }
     *
     * for (i=0; i<nmap; i++)
     * map[i]=0;
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
     *   if ((jj=loc.at(j)) && (ii>=jj)) {
     *    if (MAP(ii-1,jj-1) == 0) {
     *     MAP(ii-1,jj-1) = 1;
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
     * for (i=j; i<neq; i++) { // row loop
     *  if (MAP(i,j)) {
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
     * printf ("\nSymCompCol info: neq is %d, nwk is %d\n",neq,nz_);
     *
     * dim_[0] = dim_[1] = nColumns = nRows = neq;
     *
     * // increment version
     * this->version++;
     *
     * return TRUE;
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
                    if ( ( jj = loc.at(j) ) && ( ii >= jj ) ) {
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

    OOFEM_LOG_INFO("SymCompCol info: neq is %d, nwk is %d\n", neq, nz_);

    dim_ [ 0 ] = dim_ [ 1 ] = nColumns = nRows = neq;

    // increment version
    this->version++;

    return TRUE;
}



void SymCompCol :: times(const FloatArray &x, FloatArray &answer) const
{
    int M = dim_ [ 0 ];
    int N = dim_ [ 1 ];

    //      Check for compatible dimensions:
    if ( x.giveSize() != N ) {
        OOFEM_ERROR("SymCompCol::times: Error in CompCol -- incompatible dimensions");
    }

    answer.resize(M);
    answer.zero();

    int j, t;
    double rhs, sum;

    for ( j = 0; j < N; j++ ) {
        rhs = x(j);
        sum = 0.0;
        for ( t = colptr_(j) + 1; t < colptr_(j + 1); t++ ) {
            answer( rowind_(t) ) += val_(t) * rhs; // column loop
            sum += val_(t) * x( rowind_(t) ); // row loop
        }

        answer(j) += sum;
        answer(j) += val_( colptr_(j) ) * rhs; // diagonal
    }

    return;
}

void SymCompCol :: times(double x)
{
    for ( int t = 0; t < nz_; t++ ) {
        val_(t) *= x;
    }

    // increment version
    this->version++;
}



int SymCompCol :: assemble(const IntArray &loc, const FloatMatrix &mat)
{
    int i, j, ii, jj, dim;
    //   RowColumn* rowColumnJJ ; // 13 November 1996 - not needed anymore

#  ifdef DEBUG
    dim = mat.giveNumberOfRows();
    if ( dim != loc.giveSize() ) {
        OOFEM_ERROR("SymCompCol::assemble error : dimension of 'k' and 'loc' mismatch");
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
                if ( ii && ( ii >= jj ) ) { // assemble only lower triangular part
                    this->at(ii, jj) += mat.at(i, j);
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
    int i, j, ii, jj, dim1, dim2;

    // this->checkSizeTowards(rloc, cloc);

    dim1 = mat.giveNumberOfRows();
    dim2 = mat.giveNumberOfColumns();
    for ( i = 1; i <= dim1; i++ ) {
        ii = rloc.at(i);
        if ( ii ) {
            for ( j = 1; j <= dim2; j++ ) {
                jj = cloc.at(j);
                if ( jj && ( ii >= jj ) ) {
                    this->at(ii, jj) += mat.at(i, j);
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
    for ( int t = 0; t < nz_; t++ ) {
        val_(t) = 0.0;
    }

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

    OOFEM_ERROR("SymCompCol::operator(): Array accessing exception -- out of bounds");
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
        OOFEM_ERROR3("SymCompCol::operator(): Array accessing exception -- index out of bounds (%d,%d)", i, j);
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
        OOFEM_ERROR3("SymCompCol::operator(): Array accessing exception, index out of bounds (%d,%d)", i, j);
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

    OOFEM_ERROR3("SymCompCol::operator(): Array element (%d,%d) not in sparse structure -- cannot assign", i, j);
    return val_(0); // return to suppress compiler warning message
}
} // end namespace oofem
