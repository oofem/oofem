/*

                   *****    *****   ******  ******  ***   ***
                 **   **  **   **  **      **      ** *** **
                **   **  **   **  ****    ****    **  *  **
               **   **  **   **  **      **      **     **
              **   **  **   **  **      **      **     **
              *****    *****   **      ******  **     **


               OOFEM : Object Oriented Finite Element Code

                 Copyright (C) 1993 - 2000   Borek Patzak



         Czech Technical University, Faculty of Civil Engineering,
     Department of Structural Mechanics, 166 29 Prague, Czech Republic


    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/
/*
   Author: Richard Vondracek, <richard.vondracek@seznam.cz>
*/


#ifndef _SKYLINEMTX_H__
#define _SKYLINEMTX_H__

#include "SparseMatrixF.h"
#include "Ordering.h"
#include "BigMatrix.h"

DSS_NAMESPASE_BEGIN

/// <summary>
/// Summary description for SkyLineMtx.
/// </summary>
class SkyLineMtx : 
	public TraceableMatrix,
	public ILargeMatrix
{
public:
	SkyLineMtx(SparseMatrixF& sm,Ordering* order,MathTracer* eMT);
	virtual ~SkyLineMtx();
	
public:
	int* column_starts;
	double* columndata;
	double* D;
	Ordering* order;	

protected:
	long n;
	long nonzeros;

	void GrowSkyline(int i, int j,long* column_ns);
	void AllocateMemory(IConectMatrix* spm,int neq);	

public:
	long N() const	{return n;}
	long Nonzeros() const	{return (long)columns_data_length;}

	// This data is used in the Sealed state
	long columns_data_length;

//ILargeMatrix 
	virtual void WriteStatistics(long no_init_blocks,long no_nonzeros);
	virtual long No_Multiplications() {return 0;}

	virtual void Solve(double* b, double* x) = 0;

//ILargeMatrix 
//	virtual void LoadZeros();
//	virtual void LoadMatrixNumbers(SparseMatrixF& sm) PURE;
//	virtual void SolveLV(const LargeVector& b, LargeVector& x) PURE; 
	virtual void Factorize() = 0;
//	virtual void MultiplyByVector(const LargeVectorAttach& x, LargeVectorAttach& y) PURE;

public:
	virtual void SchurComplementFactorization(int fixed_blocks) = 0;
	virtual void SolveA11(double* x,long fixed_blocks) = 0;
	virtual void Sub_A21_A11inv(double* x,long fixed_blocks) = 0;
	virtual void Sub_A11inv_A12(double* x,long fixed_blocks) = 0;
	virtual void WriteCondensedMatrixA22(double* a,Ordering* mcn,IntArrayList* lncn) = 0;

}; //class SkyLineMtx 

DSS_NAMESPASE_END

#endif// _SKYLINEMTX_H__
