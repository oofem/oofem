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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#ifndef rotatingboundary_h
#define rotatingboundary_h

#include "boundary.h"
#include "flotarry.h"
#include "flotmtrx.h"

namespace oofem {
/**
 * Class implementing rotating boundary surface.
 * This boundary condition is usually attribute of one or more degrees of freedom (DOF).
 * @author Andreas Feymark
 */
class RotatingBoundary : public BoundaryCondition
{
protected:
    /// Rotation matrix.
    FloatMatrix R;

    /// Axis and center of rotation.
    FloatArray axis, center;

    /// Frequency by which the boundary rotates.
    double frequency;

public:
    /**
     * Constructor. Creates boundary condition with given number, belonging to given domain.
     * @param i Boundary condition number.
     * @param d Domain to which new object will belongs.
     */
    RotatingBoundary(int i, Domain *d) : BoundaryCondition(i, d) { }
    /// Destructor.
    virtual ~RotatingBoundary() { }

    virtual double give(Dof *dof, ValueModeType mode, TimeStep *tStep);

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);
    virtual void scale(double s) { }
    const char *giveClassName() const { return "RotatingBoundary"; }
    classType giveClassID() const { return RotatingBoundaryClass; }
    
};
} // end namespace oofem
#endif // rotatingboundary_h
