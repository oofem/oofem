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

#ifdef __SPOOLES_MODULE

#include "spoolessparsemtrx.h"
#include "engngm.h"
#include "spoolesinterface.h"

namespace oofem {
SparseMtrx *
SpoolesSparseMtrx :: GiveCopy() const
{
    OOFEM_ERROR("SpoolesSparseMtrx::GiveCopy - Not implemented");
    return NULL;
}

void
SpoolesSparseMtrx :: times(const FloatArray &x, FloatArray &answer) const
{
    double alpha = 1.0, beta = 0.0;
    int result;

    answer.resize( this->giveNumberOfColumns() );
    answer.zero();

    if ( sflag == SPOOLES_SYMMETRIC ) {
        result = InpMtx_sym_gmvm( this->mtrx, & beta, 1, answer.givePointer(), & alpha, 1, x.givePointer() );
    } else if ( sflag == SPOOLES_NONSYMMETRIC ) {
        result = InpMtx_nonsym_gmvm( this->mtrx, & beta, 1, answer.givePointer(), & alpha, 1, x.givePointer() );
    } else {
        OOFEM_ERROR("SpoolesSparseMtrx::times - unsupported symmetry flag");
        exit(1);
    }
}

void
SpoolesSparseMtrx :: times(double x)
{
    OOFEM_ERROR("SpoolesSparseMtrx::times(double x) - unsupported");
}

void
SpoolesSparseMtrx :: timesT(const FloatArray &x, FloatArray &answer) const
{
    double alpha = 1.0, beta = 0.0;
    int result;
    answer.resize(this->giveNumberOfRows());
    answer.zero();

    if ( sflag == SPOOLES_SYMMETRIC ) {
        result = InpMtx_sym_gmvm( this->mtrx, & beta, 1, answer.givePointer(), & alpha, 1, x.givePointer() );
    } else if ( sflag == SPOOLES_NONSYMMETRIC ) {
        result = InpMtx_nonsym_gmvm_T( this->mtrx, & beta, 1, answer.givePointer(), & alpha, 1, x.givePointer() );
    } else {
        OOFEM_ERROR("SpoolesSparseMtrx::timesT - unsupported symmetry flag");
    }

    if ( result != 1 ) {
        OOFEM_ERROR2("SpoolesSparseMtrx::timesT - error code from InpMtx_(non)sym_gmvm %d", result);
    }
}

int
SpoolesSparseMtrx :: buildInternalStructure(EngngModel *eModel, int di, EquationID ut, const UnknownNumberingScheme &s)
{
    // Determine number of equations and estimate number of nonzero entries
    int neq = eModel->giveNumberOfDomainEquations(di, ut);
    int nent = neq * 5;


    this->mtrx = InpMtx_new();
    InpMtx_init(this->mtrx, INPMTX_BY_ROWS, type, nent, neq);
    nRows = nColumns = neq;
    return true;
}

int
SpoolesSparseMtrx :: assemble(const IntArray &loc, const FloatMatrix &mat)
{
    int i, j, ac1, ac2, ndofe;

    ndofe = mat.giveNumberOfRows();

    for ( i = 1; i <= ndofe; i++ ) {
        ac1 = loc.at(i);
        if ( ac1 == 0 ) {
            continue;
        }

        for ( j = 1; j <= ndofe; j++ ) {
            ac2 = loc.at(j);
            if ( ac2 == 0 ) {
                continue;
            }

            if ( ac1 > ac2 ) {
                continue;
            }

            InpMtx_inputRealEntry( this->mtrx, ac1 - 1, ac2 - 1, mat.at(i, j) );
        }
    }

    // increment version
    this->version++;
    return 1;
}

int
SpoolesSparseMtrx :: assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat)
{
    int i, j, ii, jj, dim1, dim2;

    dim1 = mat.giveNumberOfRows();
    dim2 = mat.giveNumberOfColumns();
    for ( i = 1; i <= dim1; i++ ) {
        ii = rloc.at(i);
        if ( ii ) {
            for ( j = 1; j <= dim2; j++ ) {
                jj = cloc.at(j);
                if ( jj ) {
                    InpMtx_inputRealEntry( this->mtrx, ii - 1, jj - 1, mat.at(i, j) );
                }
            }
        }
    }

    // increment version
    this->version++;

    return 1;
}

void
SpoolesSparseMtrx :: zero()
{
    InpMtx_clearData(this->mtrx);
}

double &
SpoolesSparseMtrx :: at(int i, int j)
{
    OOFEM_ERROR("SpoolesSparseMtrx::at(i,j) - unsupported");
    abort();
}

double
SpoolesSparseMtrx :: at(int i, int j) const
{
    OOFEM_ERROR("SpoolesSparseMtrx::at(i,j) - unsupported");
    return 0.0;
}

void
SpoolesSparseMtrx :: printStatistics() const
{
    InpMtx_writeStats(this->mtrx, stdout);
}

void
SpoolesSparseMtrx :: printYourself() const
{
    InpMtx_writeForHumanEye(this->mtrx, stdout);
}

} // end namespace oofem
#endif //ifdef __SPOOLES_MODULE
