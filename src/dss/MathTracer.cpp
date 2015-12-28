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

#include "MathTracer.h"

DSS_NAMESPASE_BEGIN

MathTracer :: MathTracer()
{
    act_block = 0;
    min_pivot = 0;
    stabil_pivot = 1.0;
    break_flag = 0;
}

void MathTracer :: Write(double a)
{
    sprintf(m_string, "%e", a);
    Write(m_string);
}

void MathTracer :: Write(int a)
{
    sprintf(m_string, "%d", a);
    Write(m_string);
}

void MathTracer :: Writeln()
{
    printf("\n");
}

void MathTracer :: Writeln(const char *str)
{
    puts(str);
    //printf("\n");
}

void MathTracer :: Write(const char *str)
{
    puts(str);
}

/*void MathTracer::DrawProgress(double e)
 * {
 *      //e;
 * }*/

// true - continue factorization
// false - break factorization
bool MathTracer :: CallUnstableDialog()
{
    return true;
}

/*void MathTracer::PrintUnstablePivot(long pivot)
 * {
 *      //pivot;
 * }
 */

char *MathTracer :: NowString()
{
    memcpy(m_string, "\0", 1);
    struct tm *today;
    time_t ltime;
    time(& ltime);
    today = localtime(& ltime);
    if ( today ) {
        //return ctime(&ltime);
        //strftime( m_string, 128,"%C", today );
        sprintf( m_string, "%s", ctime(& ltime) );
    }

    return m_string;
}

clock_t MathTracer :: ClockStart(void)
{
    return clock();
}
char *MathTracer :: MeasureClock(clock_t &clock_start)
{
    clock_t end = clock();
    double duration = ( double ) ( end - clock_start ) / CLOCKS_PER_SEC;
    sprintf(m_string, "%0.3f s", duration);
    return m_string;
}

void MathTracer :: CS()
{
    //temporary_measure_start = CTime::GetCurrentTime();
    time(& m_temporary_measure_start);
    m_clock_start = clock();
}

char *MathTracer :: MC_()
{
    clock_t end = clock();
    double duration = ( double ) ( end - m_clock_start ) / CLOCKS_PER_SEC;
    sprintf(m_string, "%0.3f s", duration);
    return m_string;
}

DSS_NAMESPASE_END
