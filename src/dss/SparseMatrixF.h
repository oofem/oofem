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

#ifndef _SPARSEMATRIXF_H__
#define _SPARSEMATRIXF_H__

#include "DSSAfx.h"

DSS_NAMESPASE_BEGIN

/**
 * @author: Richard Vondracek
 */

struct IConectMatrix;

struct SparseMatrixF
{
public:
    unsigned long neq;
    double *a;

    enum eOrientation {
        eCompressedRows,      //ci - stores column indices [j]
        eCompressedColumns,   //ci - stores row indices [i]
        eSymmetric            // column indices <-> row indices, but only half is stored
    };

private:
    unsigned long *ci;
    unsigned long *adr;
    unsigned long Aoffset;
    unsigned long Coffset;
    bool bJardaConvention;
    bool bLocalCopy;
    bool m_bIsSymmertric;     // true - only diagonal + one half is stored, false - all entries are stored
    eOrientation m_eSparseOrientation;

public:
    SparseMatrixF();
    SparseMatrixF(unsigned long neq, double *a, unsigned long *ci, unsigned long *adr, unsigned long Aofs = 0, unsigned long Cofs = 0, bool JardaConvention = true, bool bIsSymetric = true, eOrientation sparseOri = eCompressedColumns);
    SparseMatrixF(IConectMatrix *pConMtx);
    ~SparseMatrixF();

    void Delete();
    void Detach();
    void CreateLocalCopy();

    long Nonzeros();
    bool IsSymmetric() { return m_bIsSymmertric; }
    eOrientation GetOrientation() { return m_eSparseOrientation; }

    void GetA12block(double *pA12, long c);
    void LoadMatrix(FILE *stream);
    void SaveMatrix(FILE *stream);

    inline unsigned long Adr(int i) const
    {
        if ( bJardaConvention ) {
            return this->adr [ i ] - Aoffset;
        } else {
            if ( i == 0 ) {
                return 0;
            } else {
                return this->adr [ i - 1 ];
            }
        }
    }

    inline unsigned long Ci(int i) const
    {
        return ( this->ci [ i ] ) - Coffset;
    }

    // c = A.b
    void MulMatrixByVector(double *b, double *c);

    // c = A.b
    void MulNonsymMatrixByVector(double *b, double *c);

    // c = A.b
    void MulSymMatrixByVector(double *b, double *c);


    /**
     * function multiplies %matrix by %vector
     *
     * @param b - array containing %vector b
     * @param c - array containing resulting %vector c = A.b
     *
     * JK
     */
    void mxv_scr(double *b, double *c);

    void ReadDiagonal(double *dv);
};

DSS_NAMESPASE_END

#endif //_SPARSEMATRIXF_H__
