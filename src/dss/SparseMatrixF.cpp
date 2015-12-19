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

#include "SparseMatrixF.h"
#include "Array.h"
#include "BigMatrix.h"

DSS_NAMESPASE_BEGIN

SparseMatrixF :: SparseMatrixF()
{
    a = NULL;
    ci = NULL;
    adr = NULL;
    neq = 0;
    Aoffset = 0;
    Coffset = 0;
    bJardaConvention = true;
    bLocalCopy = false;
    m_bIsSymmertric = false;
    m_eSparseOrientation = eCompressedColumns;
}

SparseMatrixF :: SparseMatrixF(unsigned long neq, double *a, unsigned long *ci, unsigned long *adr, unsigned long Aofs, unsigned long Cofs, bool JardaConvention, bool bIsSymetric, eOrientation sparseOri)
{
    this->a = a;
    this->neq = neq;
    this->adr = adr;
    this->ci = ci;
    Aoffset = Aofs;
    Coffset = Cofs;
    bJardaConvention = JardaConvention;
    bLocalCopy = false;
    m_bIsSymmertric = bIsSymetric;
    m_eSparseOrientation = sparseOri;
}

SparseMatrixF :: SparseMatrixF(IConectMatrix *pConMtx)
{
    /*
     * this->neq = pConMtx->N();
     * this->adr = new unsigned long[neq+1];
     *
     * unsigned long i,length = 0;
     * adr[0] = 0;
     * for (i=0; i<neq; i++)
     * {
     *      IntArrayList* plist = pConMtx->GetIndexesAboveDiagonalInColumn(i);
     *      length += plist->Count;
     *      adr[i+1] = length;
     * }
     *
     * this->ci = new unsigned long[length];
     * for (i=0; i<neq; i++)
     * {
     *      IntArrayList* plist = pConMtx->GetIndexesAboveDiagonalInColumn(i);
     *      Array::Copy(plist->Items,ci+adr[i],plist->Count);
     * }*/
}

SparseMatrixF :: ~SparseMatrixF()
{
    if ( bLocalCopy ) {
        Delete();
    }
}

void SparseMatrixF :: Delete()
{
    if ( a ) {
        delete [] a;
        a = NULL;
    }

    if ( ci ) {
        delete [] ci;
        ci = NULL;
    }

    if ( adr ) {
        delete [] adr;
        adr = NULL;
    }

    neq = 0;
    bLocalCopy = false;
}

void SparseMatrixF :: Detach()
{
    a = NULL;
    ci = NULL;
    adr = NULL;
    neq = 0;
    Aoffset = 0;
    Coffset = 0;
}

void SparseMatrixF :: CreateLocalCopy()
{
    if ( neq ) {
        long neql = neq;
        if ( bJardaConvention ) {
            neql++;
        }

        long nnz = Nonzeros();

        unsigned long *old_adr = adr;
        adr = new unsigned long [ neql ];
        Array :: Copy( ( long * ) old_adr, ( long * ) adr, neql );

        unsigned long *old_ci = ci;
        ci = new unsigned long [ nnz ];
        Array :: Copy( ( long * ) old_ci, ( long * ) ci, nnz );

        double *old_a = a;
        a = new double [ nnz ];
        Array :: Copy(old_a, a, nnz);

        bLocalCopy = true;
    }
}

long SparseMatrixF :: Nonzeros()
{
    if ( adr == 0 ) {
        return 0;
    }

    if ( bJardaConvention ) {
        return ( long ) ( adr [ neq ] - adr [ 0 ] );
    } else {
        return ( long ) ( adr [ neq - 1 ] );
    }
}

void SparseMatrixF :: GetA12block(double *pA12, long c)
{
    for ( unsigned long j = neq - c; j < neq; j++ ) {
        for ( unsigned long ad = Adr(j); ad < Adr(j + 1); ad++ ) {
            unsigned long i = Ci(ad);
            if ( i < ( neq - c ) ) {
                pA12 [ i * c + ( j - ( neq - c ) ) ] = a [ ad ];
            }
        }
    }
}

void SparseMatrixF :: MulMatrixByVector(double *b, double *c)
{
    if ( m_bIsSymmertric ) {
        MulSymMatrixByVector(b, c);
    } else {
        MulNonsymMatrixByVector(b, c);
    }
}

void SparseMatrixF :: MulNonsymMatrixByVector(double *b, double *c)
{
    switch ( m_eSparseOrientation ) {
    case eCompressedColumns:
    {
        for ( unsigned long j = 0; j < neq; j++ ) {
            for ( unsigned long ad = Adr(j); ad < Adr(j + 1); ad++ ) {
                double A = a [ ad ];
                unsigned long i = Ci(ad);
                c [ i ] += A * b [ j ];
            }
        }
    }
    break;
    case eCompressedRows:
    {
        for ( unsigned long j = 0; j < neq; j++ ) {
            for ( unsigned long ad = Adr(j); ad < Adr(j + 1); ad++ ) {
                double A = a [ ad ];
                unsigned long i = Ci(ad);
                c [ j ] += A * b [ i ];
            }
        }
    }
    break;
    default:
        fprintf(stderr, "SparseMatrixF::MulNonsymMatrixByVector: unsupported m_eSparseOrientation value\n");
        abort();
    }
}

void SparseMatrixF :: MulSymMatrixByVector(double *b, double *c)
{
    for ( unsigned long j = 0; j < neq; j++ ) {
        for ( unsigned long ad = Adr(j); ad < Adr(j + 1); ad++ ) {
            double A = a [ ad ];
            unsigned long i = Ci(ad);
            c [ i ] += A * b [ j ];
            if ( i != j ) {
                c [ j ] += A * b [ i ];
            }
        }
    }
}

void SparseMatrixF :: mxv_scr(double *b, double *c)
{
    unsigned long i, j, ii, lj, uj;
    double s, d;

    for ( i = 0; i < neq; i++ ) {
        lj = Adr(i);
        uj = Adr(i + 1);
        s = 0.0;
        d = b [ i ];
        for ( j = lj; j < uj; j++ ) {
            ii = Ci(j);
            s += a [ j ] * b [ ii ];
            c [ ii ] += a [ j ] * d;
        }

        c [ i ] = s;
    }
}

void SparseMatrixF :: ReadDiagonal(double *dv)
{
    for ( unsigned long j = 0; j < neq; j++ ) {
        for ( unsigned long ad = Adr(j); ad < Adr(j + 1); ad++ ) {
            unsigned long i = Ci(ad);
            if ( i == j ) {
                dv [ j ] = a [ ad ];
            }
        }
    }
}

void SparseMatrixF :: LoadMatrix(FILE *stream)
{
    neq = 0;
    adr = NULL;
    ci = NULL;
    a = NULL;

    bool low = false, up = false;
    //BinaryReader b = null;
    //try
    {
        //neq = b.ReadUInt32();
        fread(& neq, sizeof( neq ), 1, stream);

        adr = new unsigned long [ neq + 1 ];
        fread(adr, sizeof( unsigned long ), neq + 1, stream);
        //for (long i=0; i<=neq; i++)
        //adr[i] = b.ReadUInt32();

        ci = new unsigned long [ adr [ neq ] ];
        a = new double [ adr [ neq ] ];

        for ( unsigned long ad = 0; ad < adr [ neq ]; ad++ ) {
            //ci[ad] = b.ReadUInt32();
            fread(ci + ad, sizeof( unsigned long ), 1, stream);
            //a[ad] = b.ReadDouble();
            fread(a + ad, sizeof( double ), 1, stream);
            if ( ci [ ad ] < ad ) {
                low = true;
            }

            if ( ci [ ad ] > ad ) {
                up = true;
            }
        }

        if ( adr [ 0 ] == 1 ) {
            Aoffset = 1;
            Coffset = 1;
        }

        bLocalCopy = true;
        m_bIsSymmertric = ( !low || !up );
    }

    //catch(exception e)
    //{
    //throw e;
    //}
}

void SparseMatrixF :: SaveMatrix(FILE *stream)
{
    //BinaryWriter b = null;
    //try
    {
        //b = new BinaryWriter(s);
        //b.Write(neq);
        fwrite(& neq, sizeof( neq ), 1, stream);

        //for (long i=0; i<=neq; i++)
        //b.Write(adr[i]);
        fwrite(adr, sizeof( unsigned long ), neq + 1, stream);

        for ( unsigned long ad = 0; ad < adr [ neq ]; ad++ ) {
            //b.Write(ci[ad]);
            fwrite(ci + ad, sizeof( unsigned long ), 1, stream);
            //b.Write(a[ad]);
            fwrite(a + ad, sizeof( double ), 1, stream);
        }
    }
    //catch(exception ex)
    //{
    //throw ex;
    //}
}

/*
 * void SaveMatrixToFile(string fileName)
 * {
 *      try
 *      {
 *              Stream s = File.Open(fileName,FileMode.Create);
 *              SaveMatrix(s);
 *      }
 *      catch(Exception ex)
 *      {
 *              throw ex;
 *      }
 * }
 *
 * public static void LoadMatrixFromFile(string fileName,out long neq,out long[] adr,out long[] ci,out double[] a)
 * {
 *      try
 *      {
 *              Stream s = File.Open(fileName,FileMode.Open);
 *              LoadMatrix(s,out neq,out adr,out ci,out a);
 *      }
 *      catch(Exception ex)
 *      {
 *              throw ex;
 *      }
 * }
 */

DSS_NAMESPASE_END
