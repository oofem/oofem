/* $Header: /home/cvs/bp/oofem/sm/src/heavisideltf.h,v 1.4 2003/04/06 14:08:30 bp Exp $ */
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

//   ******************************************
//   *** CLASS HEAVISIDE LOAD TIME FUNCTION ***
//   ******************************************
#ifndef heavisideltf_h
#define heavisideltf_h

#include "flotarry.h"
#include "loadtime.h"

namespace oofem {
class HeavisideLTF : public LoadTimeFunction
{
    /*
     * This class implements a Heaviside step constant function.
     * DESCRIPTION
     * The function is defined by the origin of step and value.
     * The result is value*H(t-origin),
     * where H(t) is 0 for all t <= 0;
     *             1 for all t > 0.;
     */
private:
    double origin, value;

public:
    HeavisideLTF(int i, Domain *d) : LoadTimeFunction(i, d)
    { origin = value = 0.; }
    ~HeavisideLTF() { }

    //      void    getPoints () ;
    IRResultType initializeFrom(InputRecord *ir);
    int giveInputRecordString(std :: string &str, bool keyword = true);

    classType    giveClassID() const { return HeavisideLTFClass; }
    const char *giveClassName() const { return "HeavisideLTF"; }

    double  __at(double);
};
} // end namespace oofem
#endif // heavisideltf_h
