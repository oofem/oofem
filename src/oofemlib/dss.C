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

#ifdef __DSS_MODULE

 #include "dss.h"
 #include "error.h"
 #include "flotarry.h"
 #include "engngm.h"
 #include "domain.h"
 #include "element.h"
 #include "dofmanager.h"
 #include "DSSolver.h"
 #ifndef __MAKEDEPEND
  #include <set>
 #endif

namespace oofem {
DSSMatrix :: DSSMatrix(dssType _t) : SparseMtrx()
{
    eDSSolverType _st = eDSSFactorizationLDLT;
    _dss = new DSSolver();
    _type = _t;
    if ( _t == sym_LDL ) {
        _st = eDSSFactorizationLDLT;
    } else if ( _t == sym_LL ) {
        _st = eDSSFactorizationLLT;
    } else if ( _t == unsym_LU ) {
        _st = eDSSFactorizationLU;
    } else {
        OOFEM_ERROR("DSSMatrix::DSSMatrix() -- unknown dssType");
    }

    _dss->Initialize(0, _st);
    _sm = NULL;
    isFactorized = false;
}


DSSMatrix :: DSSMatrix(dssType _t, int n) : SparseMtrx(n, n)
{
    eDSSolverType _st = eDSSFactorizationLDLT;
    _dss = new DSSolver();
    _type = _t;
    if ( _t == sym_LDL ) {
        _st = eDSSFactorizationLDLT;
    } else if ( _t == sym_LL ) {
        _st = eDSSFactorizationLLT;
    } else if ( _t == unsym_LU ) {
        _st = eDSSFactorizationLU;
    } else {
        OOFEM_ERROR("DSSMatrix::DSSMatrix() -- unknown dssType");
    }

    _dss->Initialize(0, _st);
    _sm = NULL;
    isFactorized = false;
}

DSSMatrix :: ~DSSMatrix()
{
    delete _dss;
    delete _sm;
}

/*****************************/
/*  Copy constructor         */
/*****************************/

DSSMatrix :: DSSMatrix(const DSSMatrix &S) : SparseMtrx(S.nRows, S.nColumns)
{
    OOFEM_ERROR("DSSMatrix::DSSMatrix(const DSSMatrix &S) -- not implemented");
}

SparseMtrx *DSSMatrix :: GiveCopy() const
{
    OOFEM_ERROR("DSSMatrix::GiveCopy -- not implemented");
    return NULL;
}

void DSSMatrix :: times(const FloatArray &x, FloatArray &answer) const
{
    OOFEM_ERROR("DSSMatrix::times -- not implemented");
}

void DSSMatrix :: times(double x)
{
    OOFEM_ERROR("DSSMatrix::times -- not implemented");
}

int DSSMatrix :: buildInternalStructure(EngngModel *eModel, int di, EquationID ut, const UnknownNumberingScheme &s)
{
    IntArray loc;
    Domain *domain = eModel->giveDomain(di);
    int neq = eModel->giveNumberOfDomainEquations(di, ut);
    int nelem = domain->giveNumberOfElements();
    int i, ii, j, jj, n;
    unsigned long indx;
    Element *elem;
    // allocation map
    std :: vector< std :: set< int > >columns(neq);

    unsigned long nz_ = 0;

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
        nz_ += columns [ i ].size();
    }

    unsigned long *rowind_ = new unsigned long [ nz_ ];
    unsigned long *colptr_ = new unsigned long [ neq + 1 ];
    if ( ( rowind_ == NULL ) || ( colptr_ == NULL ) ) {
        OOFEM_ERROR("DSSMatrix::buildInternalStructure: free store exhausted, exiting");
    }

    indx = 0;

    std :: set< int > :: iterator pos;
    for ( j = 0; j < neq; j++ ) { // column loop
        colptr_ [ j ] = indx;
        for ( pos = columns [ j ].begin(); pos != columns [ j ].end(); ++pos ) { // row loop
            rowind_ [ indx++ ] = * pos;
        }
    }

    colptr_ [ neq ] = indx;

    if ( _sm ) {
        delete _sm;
    }

    if ( ( _sm = new SparseMatrixF(neq, NULL, rowind_, colptr_, 0, 0, true) ) == NULL ) {
        OOFEM_ERROR("DSSMatrix::buildInternalStructure: free store exhausted, exiting");
    }

    int bsize = eModel->giveDomain(1)->giveDefaultNodeDofIDArry().giveSize();
    /*
     *  Assemble block to equation mapping information
     */

    bool _succ = true;
    int _ndofs, _neq, ndofmans = domain->giveNumberOfDofManagers();
    long *mcn = new long [ ndofmans * bsize ];
    long _c = 0;
    DofManager *dman;

    if ( mcn == NULL ) {
        OOFEM_ERROR("DSSMatrix::buildInternalStructure: free store exhausted, exiting");
    }

    for ( n = 1; n <= ndofmans; n++ ) {
        dman = domain->giveDofManager(n);
        _ndofs = dman->giveNumberOfDofs();
        if ( _ndofs > bsize ) {
            _succ = false;
            break;
        }

        for ( i = 1; i <= _ndofs; i++ ) {
            if ( dman->giveDof(i)->isPrimaryDof() ) {
                _neq = dman->giveDof(i)->giveEquationNumber(s);
                if ( _neq > 0 ) {
                    mcn [ _c++ ] = _neq - 1;
                } else {
                    mcn [ _c++ ] = -1; // no corresponding row in sparse mtrx structure
                }
            }
        }

        for ( i = _ndofs + 1; i <= bsize; i++ ) {
            mcn [ _c++ ] = -1;                         // no corresponding row in sparse mtrx structure
        }
    }

    if ( _succ ) {
        _dss->SetMatrixPattern(_sm, bsize);
        _dss->LoadMCN(ndofmans, bsize, mcn);
    } else {
        OOFEM_LOG_INFO("DSSMatrix: using assumed block structure");
        _dss->SetMatrixPattern(_sm, bsize);
    }

    _dss->PreFactorize();
    // zero matrix, put unity on diagonal with supported dofs
    _dss->LoadZeros();
    delete[] mcn;

    OOFEM_LOG_DEBUG("DSSMatrix info: neq is %d, bsize is %d\n", neq, nz_);

    // increment version
    this->version++;

    return true;
}


int DSSMatrix :: assemble(const IntArray &loc, const FloatMatrix &mat)
{
    int i, j, ii, jj, dim;

 #  ifdef DEBUG
    dim = mat.giveNumberOfRows();
    if ( dim != loc.giveSize() ) {
        OOFEM_ERROR("CompCol::assemble : dimension of 'k' and 'loc' mismatch");
    }

 #  endif

    dim = mat.giveNumberOfRows();

    if ( _type == unsym_LU ) {
        for ( j = 1; j <= dim; j++ ) {
            jj = loc.at(j);
            if ( jj ) {
                for ( i = 1; i <= dim; i++ ) {
                    ii = loc.at(i);
                    if ( ii ) {
                        _dss->ElementAt(ii - 1, jj - 1) += mat.at(i, j);
                    }
                }
            }
        }
    } else { // symmetric pattern
        for ( j = 1; j <= dim; j++ ) {
            jj = loc.at(j);
            if ( jj ) {
                for ( i = 1; i <= dim; i++ ) {
                    ii = loc.at(i);
                    if ( ii ) {
                        if ( jj > ii ) {
                            continue;
                        }

                        _dss->ElementAt(jj - 1, ii - 1) += mat.at(j, i);
                    }
                }
            }
        }
    }

    // increment version
    this->version++;

    return 1;
}

int DSSMatrix :: assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat)
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
                    _dss->ElementAt(ii - 1, jj - 1) += mat.at(i, j);
                }
            }
        }
    }

    // increment version
    this->version++;

    return 1;
}

void DSSMatrix :: zero()
{
    _dss->LoadZeros();

    // increment version
    this->version++;
    isFactorized = false;
}

SparseMtrx *DSSMatrix :: factorized()
{
    if ( isFactorized ) {
        return this;
    }

    _dss->ReFactorize();
    isFactorized = true;
    return this;
}

void DSSMatrix :: solve(FloatArray *b, FloatArray *x)
{
    x->resize( b->giveSize() );
    _dss->Solve( x->givePointer(), b->givePointer() );
}

/*********************/
/*   Array access    */
/*********************/

double &DSSMatrix :: at(int i, int j)
{
    // increment version
    this->version++;
    return _dss->ElementAt(i - 1, j - 1);
}


double DSSMatrix :: at(int i, int j) const
{
    return _dss->ElementAt(i - 1, j - 1);
}

double DSSMatrix :: operator() (int i, int j)  const
{
    return _dss->ElementAt(i, j);
}

double &DSSMatrix :: operator() (int i, int j)
{
    // increment version
    this->version++;
    return _dss->ElementAt(i, j);
}
} // end namespace oofem

#endif // ifdef __DSS_MODULE
