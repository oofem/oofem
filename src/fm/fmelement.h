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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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

#ifndef fmelement_h
#define fmelement_h

#include "element.h"
#include "intarray.h"

namespace oofem {

///@name Declaration of basic boundary codes.
//@{
#define FMElement_PrescribedTractionBC ( 1 << 0 )
#define FMElement_PrescribedUnBC       ( 1 << 1 )
#define FMElement_PrescribedUsBC       ( 1 << 2 )
#define FMElement_PrescribedPressureBC ( 1 << 3 )
//@}

/**
 * This abstract class represent a general base element class for
 * fluid dynamic problems.
 */
class FMElement : public Element
{
protected:
    /// Array of boundary sides.
    IntArray boundarySides;
    /// Boundary sides codes.
    IntArray boundaryCodes;

public:
    FMElement(int n, Domain *aDomain);
    virtual ~FMElement();

    /**
     * Updates the stabilization coefficients used for CBS and SUPG algorithms.
     * @param tStep Active time step.
     */
    virtual void updateStabilizationCoeffs(TimeStep *tStep) { }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual const char *giveClassName() const { return "FMElement"; }
    virtual classType giveClassID() const { return FMElementClass; }
};
} // end namespace oofem
#endif // fmelement_h
