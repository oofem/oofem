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

#ifndef rotatingboundary_h
#define rotatingboundary_h

#include "boundarycondition.h"
#include "floatarray.h"
#include "floatmatrix.h"

///@name Input fields for RotatingBoundary
//@{
#define _IFT_RotatingBoundary_Name "rotatingboundary"
#define _IFT_RotatingBoundary_axis "axis"
#define _IFT_RotatingBoundary_center "center"
#define _IFT_RotatingBoundary_frequency "frequency" ///< @todo Unused ( But it makes sense that you'd have this, can you check it Andreas? ) / Mikael
//@}

namespace oofem {
/**
 * Class implementing rotating boundary surface.
 * This boundary condition is usually attribute of one or more degrees of freedom (DOF).
 *
 * @author Andreas Feymark
 */
class RotatingBoundary : public BoundaryCondition
{
protected:
    /// Rotation matrix.
    FloatMatrix R;

    /// Axis and center of rotation.
    FloatArray axis, center;

public:
    /**
     * Constructor. Creates boundary condition with given number, belonging to given domain.
     * @param i Boundary condition number.
     * @param d Domain to which new object will belongs.
     */
    RotatingBoundary(int i, Domain * d) : BoundaryCondition(i, d) { }
    /// Destructor.
    virtual ~RotatingBoundary() { }

    double give(Dof *dof, ValueModeType mode, double time) override;

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    void scale(double s) override { }

    const char *giveInputRecordName() const override { return _IFT_RotatingBoundary_Name; }
    const char *giveClassName() const override { return "RotatingBoundary"; }
};
} // end namespace oofem
#endif // rotatingboundary_h
