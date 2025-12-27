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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#ifndef integrationpointstatus_h
#define integrationpointstatus_h

#include "interface.h"
#include "interfacetype.h"
#include "contextioresulttype.h"
#include "contextmode.h"

namespace oofem {
class GaussPoint;
class TimeStep;
class DataStream;

/** Enumeration defining allowed IntegrationPointStatuses Keys.
  *  The point is that there is a need to store multiple statuses in integration points.
  *  Each status is identified by its key. The key is used to access the status. 
  */
 enum IntegrationPointStatusIDType {
    IPSID_Default = 0,
    PIDControllerStatusID = 1,
 };

/**
 * Abstract base class representing  a integration status.
 *
 * To provide opportunity for storing arbitrary data related to each integration point,
 * each integration point can store any data, represented by a class derived from this parent.
 * The history variables associated to a material law are a typical example.
 *
 * Any object that stores its status in integration point is responsible for its creation,
 * initialization, and serialization.
 */
class OOFEM_EXPORT IntegrationPointStatus
{
protected:
    /// Associated integration point.
    GaussPoint *gp;

public:
    /**
     * Constructor.
     * @param g associated integration point
     */
    IntegrationPointStatus(GaussPoint * g) : gp(g) { }
    /// Destructor.
    virtual ~IntegrationPointStatus() = default;
    /// Print receiver's output to given stream.
    virtual void printOutputAt(FILE *file, TimeStep *tStep) const { }
    /**
     * Update equilibrium history variables according to temp-variables.
     * Invoked, after new equilibrium state has been reached.
     */
    virtual void updateYourself(TimeStep *) { }
    /**
     * Allows to set the value of a specific variable, identified by varID.
     * The meaning of varID is defined in each specific implementation
     * of the method depending on the material model.
     * This method can be used to set the initial values of internal
     * variables, stresses, etc., which have been previously determined
     * by another simulation (e.g. of the manufacturing process).
     */
    virtual void setStatusVariable(int varID, double value) { }
    /**
     * Stores receiver state to output stream.
     * @param stream Output stream.
     * @param mode Determines amount of info required in stream (state, definition, ...).
     * @exception throws an ContextIOERR exception if error encountered.
     */
    virtual void saveContext(DataStream &stream, ContextMode mode) { }
    /**
     * Restores the receiver state previously written in stream.
     * @see saveContext
     * @param stream Input stream.
     * @param mode Determines amount of info available in stream (state, definition, ...).
     * @exception throws an ContextIOERR exception if error encountered.
     */
    virtual void restoreContext(DataStream &stream, ContextMode mode) { }

    virtual Interface *giveInterface(InterfaceType t) { return nullptr; }

    virtual const char *giveClassName() const = 0; //{ return "IntegrationPointStatus"; }
};
} // end namespace oofem
#endif // integrationpointstatus_h
