/* $Header: /home/cvs/bp/oofem/oofemlib/src/generalbc.h,v 1.4 2003/04/06 14:08:24 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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


//   ****************************************
//   *** CLASS General Boundary Condition ***
//   ****************************************


#ifndef generalbc_h
#define generalbc_h

#include "femcmpnn.h"
#include "domain.h"
#include "flotarry.h"
#include "dictionr.h"
#include "classtype.h"
#include "bcvaltype.h"
#include "bcgeomtype.h"
#include "bctype.h"

/**
 * Abstract base class for all boudary conditions of problem.
 * Boundary condition is an aribute of the domain (it belongs to).
 * General BC is also atrribute of several elements, nodes, Dofs and so on,
 * which are subjected to boundary condition.
 *
 * This base class only declares itself as a base class of all boundary conditions,
 * and declares only very basic services. This base class introduces
 * 'loadTimeFunction' as an atrribute of each boundary condition.
 * 'loadTimeFunction' represent time variation, its value is dependent on time step.
 * The value (or the components) of a boundary condition (load) will be
 * the product of its value by the value of
 * the associated load time function at given time step.
 * The meaning of boundary condition components is dependent on particular boundary condition type,
 * and should be defined in derived classes documentation.
 *
 * This base class introduces also two general services for requesting  boundary condition physical meaning and
 * boundary condition geometrical character (pointwise, acting on element body or edge and so on).
 *
 * Derived classes should represent the base classes for particular boundary condition type (like
 * force load, or boundary condition prescribed directlly on some dof) and should declare
 * the basic common interface.
 */
class GeneralBoundaryCondition : public FEMComponent
{
protected:
    /// Associated load time function
    int loadTimeFunction;
    /// Physical maening of bc value
    bcValType valType;
public:

    /**
     * Constructor. Creates boundary condition with given number, belonging to given domain.
     * @param n boundary condition number
     * @param d domain to which new object will belongs.
     */
    GeneralBoundaryCondition(int, Domain *);    // constructor
    /// Destructor.
    virtual ~GeneralBoundaryCondition()  { }    // destructor

    // computations
    /**
     * Returns associated load time function of receiver.
     */
    LoadTimeFunction *giveLoadTimeFunction();
    /**
     * Returns a newly allocated boundary condition, with type depending on parameter.
     * Creates new object for following classes BoundaryCondition, DeadWeight, InitialCondition and NodalLoad
     * otherwise calls directly CreateUsrDefLoadOfType global function to allocate
     * new instance of boundary condition of given type.
     * @param aClass string with boundary condition name
     * @return newly allocated boundary condition with required type.
     * @see CreateUsrDefLoadOfType function.
     */
    GeneralBoundaryCondition *ofType(char *);

    /**
     * Returns receiver load type. It distinguish particular boundary conditions according to
     * their "physical" meaning (like StructuralTemperatureLoadLT, StructuralLoadLT).
     * Derived classes should always overload, default implementation returns UnknownLT value.
     * See cltypes.h file for details.
     */
    virtual bcValType  giveBCValType() const { return this->valType; }
    /**
     * Returns type of boundary condition. It allows to distinguish bc according its
     * mathematical meaning, ie like Dirichlet, Neunamm, or Newton type */
    virtual bcType     giveType() const { return UnknownBT; }
    /**
     * Returns geometry character of boundary condition. For available values see cltypes.h file.
     * Derived classes should always overload, default implementation returns UnknownLoadGT value.
     */
    virtual bcGeomType giveBCGeoType()   const { return UnknownBGT; }
    /// Initializes receiver acording to object description stored in input record.
    IRResultType initializeFrom(InputRecord *ir);
    /** Setups the input record string of receiver
     * @param str string to be filled by input record
     * @param keyword print record keyword (default true)
     */
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);
    /**
     * Scales the receiver according to given value. Typically used in nondimensional analysis to scale down BCs and ICs.
     */
    virtual void scale(double) { }


    /// Returns classType id of receiver.
    classType    giveClassID() const { return GeneralBoundaryConditionClass; }
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "GeneralBoundaryCondition"; }


protected:
};

#endif // generalbc_h

