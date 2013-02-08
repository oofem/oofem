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

#ifndef heavisideltf_h
#define heavisideltf_h

#include "flotarry.h"
#include "loadtime.h"

namespace oofem {
/**
 * This class implements a Heaviside step load time function.
 *
 * The function is defined by the origin of step and value.
 * The result is value*H(t-origin),
 * where
 * @f[
 * H(t) = \begin{cases} 0,& t\leq 0  \\ 1, & t>0 \end{cases}
 * @f]
 */
class HeavisideLTF : public LoadTimeFunction
{
private:
    double origin, value;

public:
    HeavisideLTF(int i, Domain *d) : LoadTimeFunction(i, d)
    { origin = value = 0.; }
    virtual ~HeavisideLTF() { }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);
    virtual classType giveClassID() const { return HeavisideLTFClass; }
    virtual const char *giveClassName() const { return "HeavisideLTF"; }

    virtual double __at(double);
};
} // end namespace oofem
#endif // heavisideltf_h
