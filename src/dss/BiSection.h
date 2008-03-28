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


#ifndef bisection_h
#define bisection_h

#include "SparseConectivityMtx.h"

DSS_NAMESPASE_BEGIN

class CMcKee
{
private:
	long n;

	SparseConectivityMtxII* mtx;
	long* p_node_level;
	long* p_order; //out
	// 0 - unsorted available
	//


public:
	long* nodes;  //in
	long size;
	long domA,domB;

	CMcKee();
	~CMcKee();

	void Init(SparseConectivityMtxII* mtx);

	BOOL IsAvailable(int v);

	void PrepareValid();

	long FindFirstNode();
	void ComputeLevels();
	void DivideByMidLevel();
};


class CBiSection
{
private:
	SparseConectivityMtxII* mtx;
	CMcKee mck;

public:
	CBiSection(SparseConectivityMtxII* mtx);
	void RecurBiSectOrder(IntArrayList* order);

private:
	void RecurBiSect(long* nodes,long size);
	void BiSect(long* nodes,long size, long& domA, long& domB);
};

DSS_NAMESPASE_END

#endif 
