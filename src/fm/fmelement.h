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
public:
    FMElement(int n, Domain * aDomain);
    virtual ~FMElement();

    /**
     * Updates the stabilization coefficients used for CBS and SUPG algorithms.
     * @param tStep Active time step.
     */
    virtual void updateStabilizationCoeffs(TimeStep *tStep) { }

    virtual void computeVectorOfVelocities(ValueModeType mode, TimeStep *tStep, FloatArray &velocities);
    virtual void computeVectorOfPressures(ValueModeType mode, TimeStep *tStep, FloatArray &pressures);
};
} // end namespace oofem
#endif // fmelement_h
