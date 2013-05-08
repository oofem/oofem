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

#include "mixedgradientpressurebc.h"
#include "floatarray.h"
#include "engngm.h"
#include "feinterpol.h"
#include "set.h"
#ifdef __FM_MODULE
#include "../fm/line2boundaryelement.h"
#endif

namespace oofem {
double MixedGradientPressureBC :: domainSize()
{
    ///@todo This code is very general, and necessary for all multiscale simulations with pores, so it should be moved into Domain eventually
    int nsd = this->domain->giveNumberOfSpatialDimensions();
    double domain_size = 0.0;
    // This requires the boundary to be consistent and ordered correctly.
#if 1
    Set *set = this->giveDomain()->giveSet(this->set);
    const IntArray &boundaries = set->giveBoundaryList();

    for (int pos = 1; pos <= boundaries.giveSize()/2; ++pos) {
        Element *e = this->giveDomain()->giveElement( boundaries.at(pos*2-1) );
        int boundary = boundaries.at(pos*2);
        FEInterpolation *fei = e->giveInterpolation();
        domain_size += fei->evalNXIntegral( boundary, FEIElementGeometryWrapper(e) );
    }
    return domain_size/nsd;
#else
 #ifdef __FM_MODULE
    for (int i = 1; i <= this->domain->giveNumberOfElements(); ++i) {
        //BoundaryElement *e = dynamic_cast< BoundaryElement* >(d->giveElement(i));
        ///@todo Support more than 2D
        Element *e = dynamic_cast< Line2BoundaryElement* >(this->domain->giveElement(i));
        if (e) {
            domain_size += e->computeNXIntegral();
        }
    }
 #endif
    if  (domain_size == 0.0) { // No boundary elements? Assume full density;
        return this->domain->giveArea(); ///@todo Support more than 2D
    } else {
        return domain_size/nsd;
    }
#endif

}


IRResultType MixedGradientPressureBC :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom";
    IRResultType result;

    GeneralBoundaryCondition :: initializeFrom(ir);

    FloatArray devGradient;
    double pressure;

    IR_GIVE_FIELD(ir, devGradient, _IFT_MixedGradientPressure_devGradient);
    IR_GIVE_FIELD(ir, pressure, _IFT_MixedGradientPressure_pressure);

    this->setPrescribedDeviatoricGradientFromVoigt(devGradient);
    this->setPrescribedPressure(pressure);

    return IRRT_OK;
}

} // end namespace oofem

