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

#ifndef basicelement_h
#define basicelement_h

#include "Elements/structural2delement.h"
#define _IFT_BasicElement_Name "basicelement"

namespace oofem {
class FEI2dTrLin;

/**
 * This class implements a 'basic' triangular three-node plane-stress
 * finite element in the xy-plane. Each node has 2 degrees of freedom.
 * 
 * The current implementation is intended to serve as a simple prototype element in 
 * order to provide new users with a simple example of how a standard element can 
 * be implemented in OOFEM.
 * 
 * This element is a simplified version of the more general element TrPlaneStress2d,
 * which contains many additional (and optional) features. See 'trplanstrss.h' for
 * more info.
 * 
 * The Elements essentially only provides the interpolator, which is used by the 
 * base class PlaneStressElement for constructing the N- and B-matrices.
 *
 * @author Jim Brouzoulis 
 */
class BasicElement : public PlaneStressElement
{
protected:
    static FEI2dTrLin interp;

public:
    /// Constructor
    BasicElement(int n, Domain * d);
    /// Destructor.
    virtual ~BasicElement() { }

    virtual FEInterpolation *giveInterpolation() const;

    // Necessary for reading input files:
    virtual const char *giveInputRecordName() const { return _IFT_BasicElement_Name; }
    // Necessary for debug messages:
    virtual const char *giveClassName() const { return "BasicElement"; }

protected:

    // Support for computing the mass matrix needed for dynamic simulations
    virtual int giveNumberOfIPForMassMtrxIntegration() { return 4; }
};
} // end namespace oofem
#endif // trplanstrss_h
