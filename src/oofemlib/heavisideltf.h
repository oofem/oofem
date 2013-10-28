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

#ifndef heavisideltf_h
#define heavisideltf_h

#include "floatarray.h"
#include "loadtimefunction.h"

///@name Input fields for HeavisideLTF
//@{
#define _IFT_HeavisideLTF_Name "heavisideltf"
#define _IFT_HeavisideLTF_origin "origin"
#define _IFT_HeavisideLTF_value "value"
//@}

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
class OOFEM_EXPORT HeavisideLTF : public LoadTimeFunction
{
private:
    double origin, value;

public:
    HeavisideLTF(int i, Domain *d) : LoadTimeFunction(i, d)
    { origin = value = 0.; }
    virtual ~HeavisideLTF() { }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);
    virtual classType giveClassID() const { return HeavisideLTFClass; }
    virtual const char *giveClassName() const { return "HeavisideLTF"; }
    virtual const char *giveInputRecordName() const { return _IFT_HeavisideLTF_Name; }

    virtual double __at(double);
};
} // end namespace oofem
#endif // heavisideltf_h
