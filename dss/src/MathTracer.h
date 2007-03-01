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

#ifndef _MATH_TRACER_H__
#define _MATH_TRACER_H__

#include "DSSAfx.h"

DSS_NAMESPASE_BEGIN
////////////////// ////////////////////////////////////////////////////////////////////////////////////
// You can modify this class, or inherit a new one . If you want to redirect the output somewhere else.
class MathTracer 
{
private:
	char m_string[128];

public:
	double min_pivot;
	double stabil_pivot;
	int break_flag;
	long act_block;
	long act_row;

	MathTracer();
	virtual ~MathTracer() {}
	virtual void Write(double a);
	virtual void Write(int a);
	virtual void Writeln();
	virtual void Writeln(char* str);
	virtual void Write(char* str);

//	virtual void DrawProgress(double e);

	// true - continue factorization
	// false - break factorization
	virtual bool CallUnstableDialog();

//	virtual void PrintUnstablePivot(long pivot);

	char* NowString();
	void CS(void);
	char* MC_();

	clock_t ClockStart(void);
	char* MeasureClock(clock_t& clock);

protected:
	time_t m_temporary_measure_start;
	clock_t m_clock_start;
};


DSS_NAMESPASE_END

#endif //_MATH_TRACER_H__
