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

#ifndef _MATH_TRACER_H__
#define _MATH_TRACER_H__

#include "DSSAfx.h"

DSS_NAMESPASE_BEGIN

/**
 * @author Richard Vondracek
 */

class MathTracer
{
private:
    char m_string [ 128 ];

public:
    double min_pivot;
    double stabil_pivot;
    int break_flag;
    long act_block;
    long act_row;

    MathTracer();

    virtual void Write(double a);
    virtual void Write(int a);
    virtual void Writeln();
    virtual void Writeln(const char *str);
    virtual void Write(const char *str);

    //virtual void DrawProgress(double e);

    // true - continue factorization
    // false - break factorization
    virtual bool CallUnstableDialog();

    //virtual void PrintUnstablePivot(long pivot);

    char *NowString();
    void CS(void);
    char *MC_();

    clock_t ClockStart(void);
    char *MeasureClock(clock_t &clock);

protected:
    time_t m_temporary_measure_start;
    clock_t m_clock_start;
};


DSS_NAMESPASE_END

#endif //_MATH_TRACER_H__
