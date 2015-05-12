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

#include "mixedgradientpressurebc.h"
#include "floatarray.h"
#include "engngm.h"
#include "element.h"
#include "feinterpol.h"
#include "set.h"

namespace oofem {
double MixedGradientPressureBC :: domainSize()
{
    int nsd = this->domain->giveNumberOfSpatialDimensions();
    double domain_size = 0.0;
    // This requires the boundary to be consistent and ordered correctly.
    Set *set = this->giveDomain()->giveSet(this->set);
    const IntArray &boundaries = set->giveBoundaryList();

    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
        Element *e = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
        int boundary = boundaries.at(pos * 2);
        FEInterpolation *fei = e->giveInterpolation();
        domain_size += fei->evalNXIntegral( boundary, FEIElementGeometryWrapper(e) );
    }
    return domain_size / nsd;
}


IRResultType MixedGradientPressureBC :: initializeFrom(InputRecord *ir)
{
    IRResultType result;

    FloatArray devGradient;
    double pressure;

    IR_GIVE_FIELD(ir, devGradient, _IFT_MixedGradientPressure_devGradient);
    IR_GIVE_FIELD(ir, pressure, _IFT_MixedGradientPressure_pressure);

    this->setPrescribedDeviatoricGradientFromVoigt(devGradient);
    this->setPrescribedPressure(pressure);

    return GeneralBoundaryCondition :: initializeFrom(ir);
}
} // end namespace oofem
