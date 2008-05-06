/* $Header: /home/cvs/bp/oofem/oofemlib/src/matstatus.h,v 1.11 2003/04/06 14:08:25 bp Exp $ */
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


#ifndef matstatus_h
#define matstatus_h

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif
#include "femcmpnn.h"
#include "classtype.h"

/*
 * This class implements a material status information. It is atribute of
 * gaussPoint. This is only an abstract class, for every instance of material class
 * there should be specialized derived class, which handles are history variables.
 *
 * DESCRIPTION
 *
 * A Material status is intended to be a "status" object, for every instance of
 * material class there should exist corresponding instance of base material
 * status class, which provides specialized interface to capture possible
 * state history variables.
 * storing and restoring context is invoked at material(yield criteria)
 * level, so correspondent material(y.c) will create a new corresponding
 * conntext if needed (happen when restoreContext invoked for the first time)
 * Statuses of such objects are stored in GaussPoint.
 *
 * TASKS
 * This is abstract class - only basic functionality is supported like:
 * - storing and restoring status on tape
 * - printingYourself()
 * - updating Yourself after a new equlibrium state has been reached.
 *
 * REMARK
 * Materials statuses are atributes of GaussPoints, they are stored in
 * MatStatus variable of GaussPoint class instance.
 */

class GaussPoint;
class Dictionary;
class Domain;
class NonlocalMaterialStatusExtension;

/**
 * Abstract base class representing  a material status information.
 *
 * To provide oportunity for storing arbitrary material model related history variables
 * in integration points, associated material status class is introduced.
 * Each new material model class should be declared together with its associted status class
 * (derived from MaterialStatus class). This status can be seen as simple container,
 * storing necessary history variables and providing some access and modification methods.
 * Each integration point can contain material status. Material model should create
 * unique copy of its associated status in each integration point.
 * Because integration point is parameter of all mesages to material model
 * class, material model therefore can easily access  all history variables it needs.
 *
 * Generally, two sets of internal history variables are defined inside material status.
 * One set should always refer to previously reached equilibrium state. The second set is used to
 * describe the state during iteration of equilibrium. After the new equilibrium is reached (on structural level) the
 * variables of first set are updated according to variables of second set.
 * On the other hand, if convergence of iteration is not obtained, the variables of second set (so called temp-variables)
 * are initialized according to non-temp variables and the iteration can begin for example with smaller load increment.
 * The temp and non-temp history variables allow simple iteration restart within one step.
 * The restarts to previous steps are supported, but context for these steps must be stored.
 *
 * The general services for staus initialization and update, as well as services for
 * storing and restoring status context are declared. The implementation is left on derived
 * classes.
 *
 * The unique copy of material status class instance corresponding to material model is
 * created and associated with any integration point.
 */
class MaterialStatus : public FEMComponent
{
protected:
    /// Associated intagration point.
    GaussPoint *gp;
public:

    /**
     * Constructor.
     * @param n receiver's number
     * @param d domain to which new status belongs
     * @param g associated integration point
     */
    MaterialStatus(int n, Domain *d, GaussPoint *g);
    /// Destructor.
    ~MaterialStatus() { }
    /// Print receiver's output to given stream.
    void   printOutputAt(FILE *, TimeStep *) { }

    /**
     * Initializes the temporary internal variables, describing the current state according to
     * previously reached equilibrium internal variables.
     */
    virtual void initTempStatus() { }
    /**
     * Update equilibrium history variables according to temp-variables.
     * Invoked, after new equilibrium state has been reached.
     */
    virtual void updateYourself(TimeStep *) { } // update after new equilibrium state reached

    // definition
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "MaterialStatus"; }
    /// Returns classType id of receiver.
    classType                giveClassID() const
    { return MaterialStatusClass; }

    IRResultType initializeFrom(InputRecord *ir) { return IRRT_OK; }
};

#endif // matstatus_h




