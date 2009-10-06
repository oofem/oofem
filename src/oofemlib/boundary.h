/* $Header: /home/cvs/bp/oofem/oofemlib/src/boundary.h,v 1.11 2003/04/06 14:08:23 bp Exp $ */
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

//   ********************************
//   *** CLASS BOUNDARY CONDITION ***
//   ********************************


#ifndef boudary_h
#define boudary_h

#include "generalbc.h"
#include "dictionr.h"

#include "dof.h"
#include "bctype.h"
#include "valuemodetype.h"

#ifndef __MAKEDEPEND
#include <string.h>
#endif
#include "classtype.h"

class TimeStep;
class Dictionary;

/**
 * Class implementing Dirichlet boundary condition on DOF (primary boundary condition).
 * This boundary condition is usually attribute of one or more degrees of freedom (DOF).
 *
 * The type of unknown (physical meaning)
 * is fully determined by corresponding DOF, to which given BC is
 * associated. The previous implementation uses the 'prescribedValueDictionary' to
 * store the unknowns types, but this makes sense, when bc is associated to node,
 * but when associated to bc, the physical meaning of unknown is determined by DOF.
 *
 * Boundary condition can change its value in time using its inheritted loadTimeFunction.
 * It can also switch itself on or off depending on nonzero value of newly intorduced
 * isImposedTimeFunction load time function. Please note, that previous option must be
 * supported by particular engineering model (because equation renumbering is necessary,
 * and for incremental solution schemes DOFs unknown dictionaries must be used). See particular
 * engineering model documantation for details.
 *
 * The services provided include
 * <UL>
 * <LI>
 * Returning a component, i.e., the prescribed value of an unknown
 * (displacement, velocity, temperature, etc). Determined using
 * ValueModeType (determines whether
 * for example total value or its change from previous step is requested) - method give.</LI>
 * <LI>
 * Returning a value of isImposedTimeFunction, indicating whether b.c. is imposed or not
 * in given time (nonzero indicates imposed b.c.) - method isImposed.</LI>
 * </UL>
 */
class BoundaryCondition : public GeneralBoundaryCondition
{
    /*
     * This class implements a kinematic boundary condition. A b.c. is usually
     * attribute of one or more degrees of freedom.
     * DESCRIPTION
     * The inherited attribute 'componentArray is not used.
     * TASKS
     * -returning a component, i.e., the prescribed value of an  unknown
     * (displacement, velocity, temperature, etc).
     * -returning a value of isImposedTimeFunction, indicating whether b.c. is imposed or not
     * in given time (nonzero indicates imposed b.c.).
     * REMARK
     * Like the other Loads, a b.c. possesses a load-time function, which is ty-
     * pically a ConstantFunction.
     */
protected:
    // Dictionary of prescribed values. - not used
    // Dictionary*  prescribedValueDictionary ;

    /// prescribed value
    double prescribedValue;

    /// Load time function indicating if b.c. is imposed or not.
    int isImposedTimeFunction;

public:
    /**
     * Constructor. Creates boundary condition with given number, belonging to given domain.
     * @param n boundary condition number
     * @param d domain to which new object will belongs.
     */
    BoundaryCondition(int i, Domain *d) : GeneralBoundaryCondition(i, d)
    { isImposedTimeFunction = 0; }
    /// Destructor
    ~BoundaryCondition()            { }


    //      double  give (char,TimeStep*) ;
    /**
     * Returns the value of a prescribed unknown, respecting requested mode for given time.
     * Its physical meaning is determined by corresponding DOF.
     * @param dof determines the dof subjected to receiver bc.
     * @param mode unknown char type (if total or incremental value is returned)
     * @return prescribed value of unknown or zero if not prescribed
     */
    virtual double give(Dof *, ValueModeType, TimeStep *);
    /**
     * Returns nonzero if receiver representing b.c. is imposed at given time, otherwise returns zero.
     * @param tStep time step representing time when receiver is tested.
     * @return nonzero if imposed for given time, zero otherwise.
     */
    int     isImposed(TimeStep *);
    /**
     * Returns receiver load type. It distinguish particular boundary conditions according to
     * their "physical" meaning (like StructuralTemperatureLoadLT, StructuralLoadLT).
     * @return returns BoundaryConditionLT value.
     */
    bcType     giveType() const { return DirichletBT; }

    /// Initializes receiver acording to object description stored in input record.
    IRResultType initializeFrom(InputRecord *ir);

    /** Setups the input record string of receiver
     *  @param str string to be filled by input record
     *  @param keyword print record keyword (default true)
     */
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);
    /**
     * Scales the receiver according to given value. Typically used in nondimensional analysis to scale down BCs and ICs.
     */
    virtual void scale(double s) { prescribedValue *= s; }


    /** Set prescribed value at the input record string of receiver
     *  @param s prescribed value
     */
    virtual void setPrescribedValue(double s) { prescribedValue = s; }

    /// Returns class name of the receiver.
    const char *giveClassName() const { return "BoundaryCondition"; }
    /// Returns classType id of receiver.
    classType giveClassID() const { return BoundaryConditionClass; }
};


#endif // boudary_h

