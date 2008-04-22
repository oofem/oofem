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
/*
 * Author: Richard Vondracek, <richard.vondracek@seznam.cz>
 */


#ifndef sparsematrixf_h
#define sparsematrixf_h

#include "DSSAfx.h"

DSS_NAMESPASE_BEGIN

struct IConectMatrix;

struct SparseMatrixF
{
public:
    ULONG neq;
    double *a;
private:

    ULONG *ci;
    ULONG *adr;
    ULONG Aoffset;
    ULONG Coffset;
    bool bJardaConvention;
    bool bLocalCopy;


public:
    SparseMatrixF();
    SparseMatrixF(unsigned long neq, double *a, unsigned long *ci, unsigned long *adr, ULONG Aofs = 0, ULONG Cofs = 0, bool JardaConvention = true);
    SparseMatrixF(IConectMatrix *pConMtx);
    ~SparseMatrixF();

    void Delete();
    void Detach();
    void CreateLocalCopy();

    long Nonzeros();

    void GetA12block(double *pA12, long c);
    void LoadMatrix(FILE *stream);
    void SaveMatrix(FILE *stream);

    inline ULONG Adr(int i) const
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

    inline ULONG Ci(int i) const
    {
        return ( this->ci [ i ] ) - Coffset;
    }

    //void MulMatrixByVector(double *b,double *c);

    void MulNonsymMatrixByVector(double *b, double *c);
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

#endif
