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

#ifndef inteactionboundarycondition_h
#define inteactionboundarycondition_h

#include "boundarycondition.h"
#include "bctype.h"
#include "valuemodetype.h"

/**
 * @name Dirichlet boundary condition.
 *
 */
//@{
#define _IFT_InteractionBoundaryCondition_Name "interactionboundarycondition"
#define _IFT_BoundaryCondition_PrescribedValue "prescribedvalue" ///< [rn,optional] Prescribed value of all DOFs
#define _IFT_BoundaryCondition_PrescribedValue_d "d" ///< [rn,optional] Alternative input field
#define _IFT_BoundaryCondition_values "values" ///< [ra,optional] Vector of prescribed values for each respective DOF.
//@}

namespace oofem {
class TimeStep;
class Dof;

/**
 * This class represent a b.c. which is enforced on InteractionPFEMParticles.
 * The attached structural node provides velocity.
 */
class OOFEM_EXPORT InteractionBoundaryCondition : public BoundaryCondition
{
protected:

public:

    InteractionBoundaryCondition(int i, Domain *d) : BoundaryCondition(i, d)
    { }
    /// Destructor
    virtual ~InteractionBoundaryCondition() { }

    virtual double give(Dof *dof, ValueModeType mode, TimeStep *tStep);

    virtual const char *giveClassName() const { return "IntertactionBoundaryCondition"; }
    virtual const char *giveInputRecordName() const { return _IFT_InteractionBoundaryCondition_Name; }
};
} // end namespace oofem
#endif //inteactionboundarycondition_h
