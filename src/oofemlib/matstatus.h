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

#ifndef matstatus_h
#define matstatus_h

#include "integrationpointstatus.h"
#include "classtype.h"

namespace oofem {

class GaussPoint;
class Dictionary;
class Domain;
class NonlocalMaterialStatusExtension;

/**
 * Abstract base class representing  a material status information.
 *
 * To provide opportunity for storing arbitrary material model related history variables
 * in integration points, associated material status class is introduced.
 * Each new material model class should be declared together with its associated status class
 * (derived from MaterialStatus class). This status can be seen as simple container,
 * storing necessary history variables and providing some access and modification methods.
 * Each integration point can contain material status. Material model should create
 * unique copy of its associated status in each integration point.
 * Because integration point is parameter of all messages to material model
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
 * The general services for status initialization and update, as well as services for
 * storing and restoring status context are declared. The implementation is left on derived
 * classes.
 *
 * The unique copy of material status class instance corresponding to material model is
 * created and associated with any integration point.
 *
 * Tasks:
 * This is abstract class - only basic functionality is supported like:
 * - storing and restoring status on tape
 * - printingYourself()
 * - updating Yourself after a new equilibrium state has been reached.
 *
 * @note{Materials statuses are attributes of GaussPoints, they are stored in
 * MatStatus variable of GaussPoint class instance.}
 */
class MaterialStatus : public IntegrationPointStatus
{
protected:
public:
    /**
     * Constructor.
     * @param n receiver's number
     * @param d domain to which new status belongs
     * @param g associated integration point
     */
    MaterialStatus(int n, Domain *d, GaussPoint *g) : IntegrationPointStatus (n,d,g) {}
    /// Destructor.
    virtual ~MaterialStatus() { }
    /// Print receiver's output to given stream.
    virtual void printOutputAt(FILE *file, TimeStep *tStep) { }

    /**
     * Initializes the temporary internal variables, describing the current state according to
     * previously reached equilibrium internal variables.
     */
    virtual void initTempStatus() { }
    /**
     * Update equilibrium history variables according to temp-variables.
     * Invoked, after new equilibrium state has been reached.
     */
    virtual void updateYourself(TimeStep *) { }
    /**
     * Returns the value of material model property stored in receiving status.
     * This is typically used when random variation of some material property is considered,
     * in this case the individual values are to be stored in status (they are no longer material constants)
     * Returns true if property is available in status,  false
     */
    virtual bool giveMaterialProperty(int propID, double &value) { return false; }
    /**
     * Allows to set the value of material model property to be stored in receiving status.
     * This is typically used when random variation of some material property is considered,
     * in this case the individual values are to be stored in status (they are no longer material constants)
     */
    virtual void setMaterialProperty(int propID, double value) {}

    /**
     * Allows to set the value of a specific variable, identified by varID.
     * The meaning of varID is defined in each specific implementation
     * of the method depending on the material model.
     * This method can be used to set the initial values of internal
     * variables, stresses, etc., which have been previously determined
     * by another simulation (e.g. of the manufacturing process).
     */
    virtual void setStatusVariable(int varID, double value) {}
    /**
     * Restores consistency of the status, i.e., computes or corrects
     * the values of certain status variables such that the state is admissible.
     * For instance, if the initial values of some internal variables
     * are read from a file, other internal variables are adjusted accordingly.
     */
    virtual void restoreConsistency() {}

    virtual const char *giveClassName() const { return "MaterialStatus"; }
    virtual classType giveClassID() const { return MaterialStatusClass; }
    virtual IRResultType initializeFrom(InputRecord *ir) { return IRRT_OK; }
};
} // end namespace oofem
#endif // matstatus_h
