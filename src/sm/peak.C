/* $Header: /home/cvs/bp/oofem/sm/src/peak.C,v 1.3 2003/04/06 14:08:31 bp Exp $ */
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


//   file PEAK.CC

#include "peak.h"
#ifndef __MAKEDEPEND
#include <math.h>
#endif

double PeakFunction :: __at(double time)
// Returns the value of the receiver at time 'time'.
{
    const double precision = 0.000001;

    if ( fabs(t - time) < precision ) {
        return value;
    } else {
        return 0.;
    }
}


/* void  PeakFunction :: getCoefficients ()
 * // Reads the date anf the time increment of the receiver in the data file.
 * {
 * t     = new double[1] ;
 * value = new double[1] ;
 *
 *t     = this -> read("t") ;
 *value = this -> read("f(t)") ;
 * }*/


IRResultType
PeakFunction :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    LoadTimeFunction::initializeFrom(ir);
    IR_GIVE_FIELD(ir, t, IFT_PeakFunction_t, "t"); // Macro
    IR_GIVE_FIELD(ir, value, IFT_PeakFunction_ft, "f(t)"); // Macro

    return IRRT_OK;
}

