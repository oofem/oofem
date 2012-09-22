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

#ifndef boudary_h
#define boudary_h

#include "generalbc.h"
#include "bctype.h"
#include "valuemodetype.h"
#include "classtype.h"

namespace oofem {
class TimeStep;
class Dof;

/**
 * Class implementing Dirichlet boundary condition on DOF (primary boundary condition).
 * This boundary condition is usually attribute of one or more degrees of freedom (DOF).
 *
 * The type of unknown (physical meaning)
 * is fully determined by corresponding DOF, to which given BC is
 * associated. The previous implementation uses the 'prescribedValueDictionary' to
 * store the unknowns types, but this makes sense, when BC is associated to node,
 * but when associated to BC, the physical meaning of unknown is determined by DOF.
 *
 * Boundary condition can change its value in time using its inherited loadTimeFunction.
 * It can also switch itself on or off depending on nonzero value of 
 * isImposedTimeFunction load time function. Please note, that previous option must be
 * supported by particular engineering model (because equation renumbering is necessary,
 * and for incremental solution schemes DOFs unknown dictionaries must be used). See 
 * particular engineering model documentation for details.
 *
 * The services provided include
 * - Returning a component, i.e., the prescribed value of an unknown
 *   (displacement, velocity, temperature, etc). Determined using
 *   ValueModeType (determines whether
 *   for example total value or its change from previous step is requested) - method give.
 * - Returning a value of isImposedTimeFunction, indicating whether b.c. is imposed or not
 *   in given time (nonzero indicates imposed BC) - method isImposed.
 */
class BoundaryCondition : public GeneralBoundaryCondition
{
protected:
    /// Prescribed value of DOF.
    double prescribedValue;
    
public:
    /**
     * Constructor. Creates boundary condition with given number, belonging to given domain.
     * @param i Boundary condition number.
     * @param d Domain to which new object will belongs.
     */
    BoundaryCondition(int i, Domain *d) : GeneralBoundaryCondition(i, d)
    { }
    /// Destructor
    virtual ~BoundaryCondition() { }

    /**
     * Returns the value of a prescribed unknown, respecting requested mode for given time.
     * Its physical meaning is determined by corresponding DOF.
     * This function should only be used if the BC is imposed.
     * @see isImposed
     * @param dof Determines the dof subjected to receiver BC.
     * @param mode Unknown char type (if total or incremental value is returned).
     * @param tStep Time step to give value for.
     * @return Prescribed value of unknown or zero if not prescribed.
     */
    virtual double give(Dof *dof, ValueModeType mode, TimeStep *tStep);

    /**
     * Set prescribed value at the input record string of receiver
     * @param s prescribed value
     */
    virtual void setPrescribedValue(double s) { prescribedValue = s; }

    // Overloaded methods:
    virtual bcType giveType() const { return DirichletBT; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);
    virtual void scale(double s) { prescribedValue *= s; }
    virtual const char *giveClassName() const { return "BoundaryCondition"; }
    virtual classType giveClassID() const { return BoundaryConditionClass; }
};
} // end namespace oofem
#endif // boudary_h

