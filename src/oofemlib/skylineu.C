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

#include "skylineu.h"
#include "rowcol.h"
#include "floatmatrix.h"
#include "intarray.h"
#include "domain.h"
#include "engngm.h"
#include "element.h"
#include "mathfem.h"
#include "verbose.h"
#include "error.h"
#include "sparsemtrxtype.h"
#include "activebc.h"
#include "classfactory.h"

#ifdef TIME_REPORT
 #include "timer.h"
#endif

namespace oofem {
REGISTER_SparseMtrx(SkylineUnsym, SMT_SkylineU);

SkylineUnsym :: SkylineUnsym(int n) : SparseMtrx(n, n),
    isFactorized(false)
{
}


SkylineUnsym :: SkylineUnsym(const SkylineUnsym &s) : SparseMtrx(s.nRows, s.nColumns),
    rowColumns(s.rowColumns),
    isFactorized(s.isFactorized)
{
}


std::unique_ptr<SparseMtrx>
SkylineUnsym :: clone() const
{
    return std::make_unique<SkylineUnsym>(*this);
}


void
SkylineUnsym :: toFloatMatrix(FloatMatrix &answer) const
{
    answer.resize(this->giveNumberOfRows(), this->giveNumberOfColumns());
    answer.zero();
    for ( int j = 1; j <= this->giveNumberOfColumns(); j++ ) {
        int start = this->rowColumns[j-1].giveStart();
        for ( int i = start; i <= j; i++ ) {
            answer.at(i, j) = this->at(i, j);
            answer.at(j, i) = this->at(j, i);
        }
    }
}


double &
SkylineUnsym :: at(int i, int j)
{
    this->version++;

    if ( i < j ) {
        return this->rowColumns[j-1].atU(i);
    } else if ( i > j ) {
        return this->rowColumns[i-1].atL(j);
    } else {
        return this->rowColumns[i-1].atDiag();
    }
}


double
SkylineUnsym :: at(int i, int j) const
{
    if ( i < j ) {
        return this->rowColumns[j-1].atU(i);
    } else if ( i > j ) {
        return this->rowColumns[i-1].atL(j);
    } else {
        return this->rowColumns[i-1].atDiag();
    }
}


int
SkylineUnsym :: assemble(const IntArray &loc, const FloatMatrix &mat)
{
    int dim = mat.giveNumberOfRows();

#  ifdef DEBUG
    if ( dim != loc.giveSize() ) {
        OOFEM_ERROR("dimension of 'k' and 'loc' mismatch");
    }

    this->checkSizeTowards(loc);
#  endif

    // added to make it work for nonlocal model !!!!!!!!!!!!!!!!!!!!!!!!!!
    // checkSizeTowards(loc) ;

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


int
SkylineUnsym :: assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat)
{
    this->checkSizeTowards(rloc, cloc);

    int dim1 = rloc.giveSize();
    int dim2 = cloc.giveSize();
    for ( int i = 1; i <= dim1; i++ ) {
        int ii = rloc.at(i);
        if ( ii ) {
            for ( int j = 1; j <= dim2; j++ ) {
                int jj = cloc.at(j);
                if ( jj ) {
                    this->at(ii, jj) += mat.at(i, j);
                }
            }
        }
    }

    this->version++;

    return 1;
}


void
SkylineUnsym :: checkSizeTowards(const IntArray &rloc, const IntArray &cloc)
{
    int maxCol = 0;

    for ( auto r : rloc ) {
        maxCol = max( maxCol, r );
    }

    for ( auto c : cloc ) {
        maxCol = max( maxCol, c );
    }

    if ( maxCol > this->giveNumberOfColumns() ) {
        growTo(maxCol);
    }

    for ( auto ii : rloc ) {
        if ( ii ) {
            rowColumns[ii-1].checkSizeTowards(cloc);
        }
    }

    for ( auto jj : cloc ) {
        if ( jj ) {
            rowColumns[jj-1].checkSizeTowards(rloc);
        }
    }
}


int
SkylineUnsym :: buildInternalStructure(EngngModel *eModel, int di, const UnknownNumberingScheme &s)
{
    // Instanciates the profile of the receiver and initializes all coefficients to zero.
    // Warning : case diagonal (lumped) matrix to expected.
    IntArray loc;
    int first, nonlocal;
    Domain *domain = eModel->giveDomain(di);
    int neq = eModel->giveNumberOfDomainEquations(di, s);
    int nelem = domain->giveNumberOfElements();
    nonlocal = eModel->useNonlocalStiffnessOption();

    // clear receiver if exist
    this->rowColumns.clear();
    this->growTo(neq); // from now on, size = MaxIndex

    // Set up the array with indices of first nonzero elements in each row
    IntArray firstIndex(neq);
    for ( int i = 1; i <= neq; i++ ) {
        firstIndex.at(i) = i;
    }

    for ( int n = 1; n <= nelem; n++ ) {
        auto elem = domain->giveElement(n);
        //elem -> giveLocationArray (loc) ;

        if ( !nonlocal ) {
            elem->giveLocationArray(loc, s);
        }
        //else ((StructuralElement*)elem) -> giveNonlocalLocationArray(loc) ;
        else {
            elem->giveLocationArray(loc, s);
        }

        //    Find 'first', the smallest positive number in LocArray
        first = neq;
        for ( int i = 1; i <= loc.giveSize(); i++ ) {
            int ii = loc.at(i);
            if ( ii && ii < first ) {
                first = ii;
            }
        }

        //    Make sure that the FirstIndex is not larger than 'first'
        for ( int i = 1; i <= loc.giveSize(); i++ ) {
            int ii = loc.at(i);
            if ( ii && ( first < firstIndex.at(ii) ) ) {
                firstIndex.at(ii) = first;
            }
        }
    }

    // loop over active boundary conditions
    int nbc = domain->giveNumberOfBoundaryConditions();
    std :: vector< IntArray >r_locs;
    std :: vector< IntArray >c_locs;

    for ( int i = 1; i <= nbc; ++i ) {
        ActiveBoundaryCondition *bc = dynamic_cast< ActiveBoundaryCondition * >( domain->giveBc(i) );
        if ( bc ) {
            bc->giveLocationArrays(r_locs, c_locs, UnknownCharType, s, s);
            for ( std :: size_t j = 0; j < r_locs.size(); j++ ) {
                IntArray &krloc = r_locs [ j ];
                IntArray &kcloc = c_locs [ j ];
                first = neq;
                for ( int k = 1; k <= kcloc.giveSize(); k++ ) {
                    int kk = kcloc.at(k);
                    if ( kk ) {
                        first = min(first, kk);
                    }
                }
                for ( int k = 1; k <= krloc.giveSize(); k++ ) {
                    int kk = krloc.at(k);
                    if ( kk && ( first < firstIndex.at(kk) ) ) {
                        firstIndex.at(kk) = first;
                    }
                }
            }
        }
    }


    for ( int i = 1; i <= neq; i++ ) {
        this->rowColumns[i-1].growTo( firstIndex.at(i) );
    }

    this->printStatistics();

    nRows = nColumns = neq;
    this->version++;

    return true;
}


int SkylineUnsym :: setInternalStructure(IntArray &adr1)
{
    // allocates and built structure according to given
    // array of maximal column heights
    //
    IntArray mht, firstIndex;
    int n = adr1.giveSize();
    this->growTo(n - 1);

    // build first index array
    mht.resize(n - 1);
    mht.zero();
    firstIndex.resize(n - 1);
    firstIndex.zero();
    for ( int i = 1; i <= n - 1; i++ ) {
        mht.at(i) = adr1.at(i + 1) - adr1.at(i);
        firstIndex.at(i) = i - mht.at(i) + 1;
    }

    // Enlarge the rowcolumns
    for ( int i = 1; i <= n - 1; i++ ) {
        this->rowColumns[i-1].growTo( firstIndex.at(i) );
    }

    this->printStatistics();

    nRows = nColumns = n - 1;

    this->version++;

    return true;
}



void
SkylineUnsym :: checkSizeTowards(const IntArray &loc)
{
    int maxCol = 0; // the largest coefficient in 'loc'
    for ( auto c: loc ) {
        maxCol = max( maxCol, c );
    }

    if ( maxCol > (int)this->rowColumns.size() ) { // enlarge the matrix
        growTo(maxCol);
    }

    for ( auto ii: loc ) {
        if ( ii ) {
            this->rowColumns[ii-1].checkSizeTowards(loc);
        }
    }
}


SparseMtrx *
SkylineUnsym :: factorized()
// Returns the receiver in L.D.U factorized form. From Golub & van Loan,
// 1rst edition, pp 83-84.
{
    FloatArray r, w;
#ifdef TIME_REPORT
    Timer timer;
    timer.startTimer();
#endif

    if ( isFactorized ) {
        return this;
    }

    if ( this->rowColumns.empty() ) {
        OOFEM_WARNING("null-sized matrix factorized");
        isFactorized = 1;
        return this;
    }

    for ( int k = 1; k <= this->giveNumberOfColumns(); k++ ) {
        auto &rowColumnK = this->rowColumns[k-1];
        int startK = rowColumnK.giveStart();

        // compute vectors r and w
        r.resize(k - 1);
        r.zero();
        w.resize(k - 1);
        w.zero();
        for ( int p = startK; p < k; p++ ) {
            double diag = this->rowColumns[p-1].atDiag();
            r.at(p) = diag * rowColumnK.atU(p);
            w.at(p) = diag * rowColumnK.atL(p);
        }

        // compute diagonal coefficient of rowColumn k
        rowColumnK.atDiag() -= rowColumnK.dot(r, 'R', startK, k - 1);
        double diag = rowColumnK.atDiag();

        // test pivot not too small
        if ( fabs(diag) < SkylineUnsym_TINY_PIVOT ) {
            rowColumnK.atDiag() = diag = SkylineUnsym_TINY_PIVOT;
            OOFEM_LOG_DEBUG("SkylineUnsym :: factorized: zero pivot %d artificially set to a small value", k);
        }

        // compute off-diagonal coefficients of rowColumns i>k
        for ( int i = k + 1; i <= this->giveNumberOfColumns(); i++ ) {
            auto &rowColumnI = this->rowColumns[i-1];
            int startI = rowColumnI.giveStart();
            if ( startI <= k ) {
                int start = max(startI, startK);
                rowColumnI.atL(k) -= rowColumnI.dot(r, 'R', start, k - 1);
                rowColumnI.atL(k) /= diag;
                rowColumnI.atU(k) -= rowColumnI.dot(w, 'C', start, k - 1);
                rowColumnI.atU(k) /= diag;
            }
        }
    }

    isFactorized = true;

#ifdef TIME_REPORT
    timer.stopTimer();
    OOFEM_LOG_DEBUG( "SkylineU info: user time consumed by factorization: %.2fs\n", timer.getUtime() );
#endif

    //this->version++;

    return this;
}



FloatArray *
SkylineUnsym :: backSubstitutionWith(FloatArray &y) const
{
    if ( y.giveSize() != this->giveNumberOfColumns() ) {
        OOFEM_ERROR("size mismatch");
    }

    for ( int k = 1; k <= this->giveNumberOfColumns(); k++ ) {
        auto &rowColumnK = this->rowColumns[k-1];
        int start = rowColumnK.giveStart();
        y.at(k) -= rowColumnK.dot(y, 'R', start, k - 1);
    }

    // diagonalScaling
    for ( int k = 1; k <= this->giveNumberOfColumns(); k++ ) {
        double diag = this->rowColumns[k-1].atDiag();
#     ifdef DEBUG
        if ( fabs(diag) < SkylineUnsym_TINY_PIVOT ) {
            OOFEM_ERROR("pivot %d is small", k);
        }

#     endif
        y.at(k) /= diag;
    }


    for ( int k = this->giveNumberOfColumns(); k > 0; k-- ) {
        auto &rowColumnK = this->rowColumns[k-1];
        double yK = y.at(k);
        int start = rowColumnK.giveStart();
        for ( int i = start; i < k; i++ ) {
            y.at(i) -= rowColumnK.atU(i) * yK;
        }
    }

    return & y;
}


void
SkylineUnsym :: times(const FloatArray &x, FloatArray &answer) const
{
    if ( this->giveNumberOfColumns() != x.giveSize() ) {
        OOFEM_ERROR("size mismatch");
    }

    answer.resize(this->giveNumberOfRows());
    answer.zero();

    for ( int i = 1; i <= this->giveNumberOfColumns(); i++ ) {
        auto &rowColumni = this->rowColumns[i-1];
        int starti = rowColumni.giveStart();
        answer.at(i) += rowColumni.dot(x, 'R', starti, i - 1);
        answer.at(i) += rowColumni.atDiag() * x.at(i);

        for ( int j = starti; j <= i - 1; j++ ) {
            answer.at(j) += rowColumni.atU(j) * x.at(i);
        }
    }
}


void SkylineUnsym :: times(double x)
{
    for ( int i = 1; i <= this->giveNumberOfColumns(); i++ ) {
        auto &rowColumni = this->rowColumns[i-1];
        int starti = rowColumni.giveStart();
        for ( int j = starti; j <= i - 1; j++ ) {
            rowColumni.atU(j) *= x;
            rowColumni.atL(j) *= x;
        }

        rowColumni.atDiag() *= x;
    }

    this->version++;
}


void
SkylineUnsym :: growTo(int n)
{
    this->rowColumns.reserve(n);
    for (int i = 1; i <= n; ++i ) {
        this->rowColumns.emplace_back(i);
    }
    this->nRows = n;
    this->nColumns = n;
}

void
SkylineUnsym :: printStatistics() const
{
    int nelem = 0;
    for ( auto &rc : this->rowColumns ) {
        nelem += rc.giveSize();
    }

    OOFEM_LOG_INFO("Skylineu info: neq is %d, nwk is %d\n", this->giveNumberOfRows(), nelem);
}


void
SkylineUnsym :: writeToFile(const char *fname) const
{
    FILE *file = fopen(fname, "w");
    FloatMatrix copy;
    this->toFloatMatrix(copy);
    for ( int i = 1; i <= nRows; ++i ) {
        for ( int j = 1; j <= nColumns; ++j ) {
            fprintf( file, "%.16e ", copy.at(i, j) );
        }
        fprintf(file, "\n");
    }
    fclose(file);
}


void
SkylineUnsym :: printYourself() const
{
    FloatMatrix copy;
    this->toFloatMatrix(copy);
    copy.printYourself();
    for ( auto &rc : this->rowColumns ) {
        rc.printYourself();
    }
}

void
SkylineUnsym :: zero()
{
    for ( auto &rc : this->rowColumns ) {
        rc.zero();
    }

    isFactorized = false;

    this->version++;
}


void SkylineUnsym :: timesT(const FloatArray &x, FloatArray &answer) const
{
    if ( this->giveNumberOfRows() != x.giveSize() ) {
        OOFEM_ERROR("size mismatch");
    }

    answer.resize(this->giveNumberOfColumns());
    answer.zero();

    for ( int i = 1; i <= this->giveNumberOfRows(); i++ ) {
        auto &rowColumni = this->rowColumns[i-1];
        int starti = rowColumni.giveStart();
        answer.at(i) += rowColumni.dot(x, 'C', starti, i - 1);
        answer.at(i) += rowColumni.atDiag() * x.at(i);

        for ( int j = starti; j <= i - 1; j++ ) {
            answer.at(j) += rowColumni.atL(j) * x.at(i);
        }
    }
}
} // end namespace oofem
