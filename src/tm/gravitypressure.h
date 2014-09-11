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

#ifndef gravpress_h
#define gravpress_h

#include "bodyload.h"
#include "bcgeomtype.h"
#include "valuemodetype.h"

/**
 * This class implements a gravity-like pressure for transport models.
 * The inheritted attribute 'componentArray' contains the components of an
 * loading prescribed per unit volume.
 */

///@name Input fields for GravityPressure
//@{
#define _IFT_GravityPressure_Name "gravitypressure"
#define _IFT_GravityPressure_normal "normal"
#define _IFT_GravityPressure_zerolevel "zerolevel"
//@}

namespace oofem {
/**
 * This class implements a gravity-like load.
 * The attribute 'componentArray' contains the components of an acceleration
 * 'a', expected (but not required) to be downwards vertical.
 * Task: returning the body force  rho*a  acting on a given element.
 */
class GravityPressure : public BodyLoad
{
protected:

    double zeroLevel;
    FloatArray normalVector;

public:
    /// Constructor
    GravityPressure(int i, Domain * d) : BodyLoad(i, d) { }

    virtual bcGeomType giveBCGeoType() const { return GravityPressureBGT; }
    virtual void computeValueAt(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, ValueModeType mode);

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual const char *giveClassName() const { return "GravityPressure"; }
    virtual const char *giveInputRecordName() const { return _IFT_GravityPressure_Name; }
};
} // end namespace oofem

#endif // gravpress_h
