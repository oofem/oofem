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

#include "sm/Elements/structural2delement.h"
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

    FEInterpolation *giveInterpolation() const override;

    // Necessary for reading input files:
    const char *giveInputRecordName() const override { return _IFT_BasicElement_Name; }
    // Necessary for debug messages:
    const char *giveClassName() const override { return "BasicElement"; }
    Element_Geometry_Type giveGeometryType() const override {return EGT_triangle_1;}


protected:

    // Support for computing the mass matrix needed for dynamic simulations
    int giveNumberOfIPForMassMtrxIntegration() override { return 4; }
};
} // end namespace oofem
#endif // basicelement_h
