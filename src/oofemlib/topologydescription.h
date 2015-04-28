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


#ifndef topologydescription_h
#define topologydescription_h

#include "error.h"
#include "oofemcfg.h"

namespace oofem {
class SpatialLocalizer;
class TimeStep;
class Domain;
class DataReader;

/**
 * Determines the state of the evolving topology.
 */
enum TopologyState {
    TS_OK,             ///< Indicates that everything is OK with respect to topology.
    TS_NeedsRemeshing, ///< Indicates that the topology has reached a need for remeshing, as the case with merging surfaces.
};

/**
 * Abstract class for topology description.
 * The topology is coupled to a given domain. The main workload for a topology description is to be able to reproduce a finite element mesh from the results.
 *
 * @author Mikael Öhman
 */
class OOFEM_EXPORT TopologyDescription
{
protected:
    /// Domain which topology belongs to.
    Domain *d;

public:
    TopologyDescription(Domain * d) {
        this->d = d;
    }
    virtual ~TopologyDescription() { }
    /**
     * Instanciates itself.
     */
    virtual bool instanciateYourself(DataReader *dr) = 0;

    /**
     * Updates the topology from the FE solution.
     * @param tStep Active time step.
     */
    virtual TopologyState updateYourself(TimeStep *tStep) = 0;

    /**
     * File output of the current state of the topology description.
     * This is not handled by the export modules, since each type of representation can differ.
     * @param tStep Active time step.
     */
    virtual void doOutput(TimeStep *tStep)
    { OOFEM_ERROR( "Not implemented" ); }

    /**
     * Generates the FE components from the bare mesh.
     * Does not map fields or internal variables.
     * @todo Placing it in a new domain is probably preferable.
     */
    virtual void replaceFEMesh() // (Domain *& newDomain)
    { OOFEM_ERROR( "Not implemented" ); }

    /**
     * Changes the connected domain of receiver.
     * @note{Does not delete any existing objects.}
     */
    virtual void setDomain(Domain *newDomain) { d = newDomain; }

    /**
     * Gives the name of the class.
     */
    virtual const char *giveClassName() const = 0;

    /// Returns string for prepending output (used by error reporting macros).
    std :: string errorInfo(const char *func) const {
        return std::string(this->giveClassName()) + "::" + func;
    }
};
} // end namespace oofem
#endif // topologydescription_h
