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

#ifndef _SPARSECONMTX_H__
#define _SPARSECONMTX_H__

#include "ColHash.h"
#include "ConList.h"
#include "BigMatrix.h"
#include "SparseMatrixF.h"
#include "Ordering.h"

DSS_NAMESPASE_BEGIN

/**
 * This is a dynamic sparse matrix format
 * There is full set of rows but in each row are stored only nonzeros.
 *
 * @author Richard Vondracek
 */
class SparseConectivityMtxII :
    public TraceableMatrix,
    public IConectMatrix
{
public:
    IntArrayList **ColumnsIndexes;
    long n;
    long nonzeros;

public:
    long N() const { return n; }
    long Nonzeros() const { return nonzeros; }

public:

    SparseConectivityMtxII(const SparseMatrixF &sm, char block_size);
    SparseConectivityMtxII(const SparseMatrixF &sm, Ordering *node_order, char block_size);
    SparseConectivityMtxII(SparseConectivityMtxII &mtx, Ordering *order);
    SparseConectivityMtxII(const SparseConectivityMtxII &mtx);

    virtual ~SparseConectivityMtxII();

    void AddGrow(long i, long j);

    long ColumnLength(long i);

    // Returns compressed row storage, the caller must delete the data
    void GetCmpRows(long * &adr, long * &ci);

    IntArrayList *GetOrder_Cuthill_McKee2();

    Ordering *GenerateMD(IntArrayList *fixed = NULL);
    Ordering *GenerateAMD(IntArrayList *fixed = NULL);
    Ordering *GenerateAMD_AA(IntArrayList *fixed = NULL);
    Ordering *Get_Cuthill_McKee();
    Ordering *Get_Reverse_Cuthill_McKee();
    Ordering *Get_Unity();
    Ordering *Get_RecursiveBiSection();
    Ordering *Get_MetisDiSection();
    Ordering *Get_ColAMD();

    void GenerateFillInPresorted(Ordering *ord);

    Ordering *GetOrdering(Ordering :: Type ord);

    IntArrayList *GetIndexesAboveDiagonalInColumn(long j);

    IntArrayList *DetachIndexesAboveDiagonalInColumn(long j);

    Ordering *GetPermutationAndPattern(Ordering :: Type ord, IntArrayList *fixed = NULL);
};

/// <summary>
/// Computes the minimum degree or the approximate minimum degree ordering
/// using the quotient graph with support for mass elimination and supervariables.
/// </summary>
class MD_Qqraph
{
public:
    SparseConectivityMtxII *mtx;
    long *degrees;
    long *amd_w;
private:
    long MinDegA, MinDegB, MinDegC;
    long vi_Min;

public:
    IntArrayList **A;
    IntArrayList **E;
    IntArrayList **I;           // supervariables
    IntArrayList **L;

    long *pLIdx;
    long *pPIdx;
    long *pfixed;       // is the i-th node fixed or not

    IntLinkArray variables;
    IntArrayList *elements;
    long no_elements;

    bool approximate_degree;
    bool keep_sorted_order;
    bool aggressive_absorbtion;

    long n;
    long n_1;
    ColHash ht;

    MD_Qqraph(SparseConectivityMtxII *mtx);

    ~MD_Qqraph();

    static void Insert(IntArrayList &al, long i);

    // true  - the variable is available for elimination
    // flase - the variable is fixed (FETI border nodes)
    inline bool IsFree(long i)
    {
        return bool( pfixed == NULL || pfixed [ i ] == 0 );
    }

    void ComputeAMD_w_Le_Lp(long i);

    void ClearAMD_w(long i);

    long ApproximateMinumumDegree(long i, IntArrayList &Lp);

    /// <summary> This method computes the exact external degree </summary>
    /// <param name="i">variable</param>
    /// <returns>degree</returns>
    long ExternalDegree(long i);

    /// <summary> Masselimination </summary>
    /// <param name="p">pivot element</param>
    void Eliminate(long p);

    long Hash(long i)
    {
        long hk = A [ i ]->SumElements() + E [ i ]->SumElements();
        return ( hk % n_1 ) + 1;
    }

    bool IsIndistinguishable(long i, long j);

    IntArrayList hash_parents;
    void SupervariablesDetection(IntArrayList &Lp);

    /// <summary>
    /// Returns the permutation according to minimum degree algorithm.
    /// </summary>
    /// <param name="approximate_degree"> true=approximate degree, false=true external degree</param>
    /// <returns>permutation vector</returns>
    IntArrayList *GenerateMD(bool approximate_degree, IntArrayList *fixed = NULL);
}; //class MD_Qqraph


DSS_NAMESPASE_END

#endif //_SPARSECONMTX_H__
