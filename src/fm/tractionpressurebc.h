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

#ifndef tractionpressurebc_h
#define tractionpressurebc_h

#include "boundarycondition.h"

#define _IFT_TractionPressureBC_Name "prescribedtractionpressurebc"

namespace oofem {
/**
 * Class implementing prescribed pressure bc due to prescribed tractions (Dirichlet boundary condition on DOF).
 * This boundary condition is usually attribute of one or more degrees of freedom (DOF).
 */
class TractionPressureBC : public BoundaryCondition
{
public:
    /**
     * Constructor. Creates boundary condition with given number, belonging to given domain.
     * @param i Boundary condition number.
     * @param d Domain to which new object will belongs.
     */
    TractionPressureBC(int i, Domain * d) : BoundaryCondition(i, d) { }
    /// Destructor.
    virtual ~TractionPressureBC() { }

    virtual double give(Dof *dof, ValueModeType mode, double time);

    virtual void scale(double s) { }
    virtual const char *giveClassName() const { return "TractionPressureBC"; }
    virtual const char *giveInputRecordName() const { return _IFT_TractionPressureBC_Name; }
};
} // end namespace oofem
#endif // tractionpressurebc_h
