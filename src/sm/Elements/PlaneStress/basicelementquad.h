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

#ifndef basicelementquad_h
#define basicelementquad_h

#include "Elements/structural2delement.h"
#define _IFT_BasicElementQuad_Name "basicelementquad"

namespace oofem {
class FEI2dQuadLin;

/**
 * This class implements a 'basic' quadratic four node plane-stress
 * finite element in the xy-plane. Each node has 2 degrees of freedom.
 * 
 * Based on the basicelement
 *
 * /@author Johannes Främby
 */
class BasicElementQuad : public PlaneStressElement
{
    /**
     * All of the following methods need to be implemented by the element 
     * (if not stated otherwise).
     */
    
protected:
    static FEI2dQuadLin interp;           

public:
    /// Constructor
    BasicElementQuad(int n, Domain * d);    
    /// Destructor.
    virtual ~BasicElementQuad() { }         

    virtual FEInterpolation *giveInterpolation() const;    
    virtual const char *giveInputRecordName() const { return _IFT_BasicElementQuad_Name; }
    virtual const char *giveClassName() const { return "BasicElementQuad"; }

protected:
    
    // - Support for computing the mass matrix needed for dynamic simulations
    virtual int giveNumberOfIPForMassMtrxIntegration() { return 4; }    


};
} // end namespace oofem
#endif // basicelementquad_h
