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


#ifndef topologydescription_h
#define topologydescription_h

#include "domain.h"
#include "classtype.h"
#include "error.h"

namespace oofem {
class SpatialLocalizer;
class TimeStep;

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
class TopologyDescription
{
protected:
    /// Domain which topology belongs to.
    Domain *d;
    /// If the requested property is displacement (false for velocities)
    bool useDisplacements;

public:
    TopologyDescription(Domain *d) {
        this->d = d;
        domainType dt = d->giveDomainType();
        this->useDisplacements = dt == _2dPlaneStressMode || dt == _PlaneStrainMode || dt == _2dPlaneStressRotMode || dt == _3dMode || dt == _3dAxisymmMode;
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
    { OOFEM_ERROR2("%s::doOutput - Not implemented",this->giveClassName()); };

    /**
     * Generates the FE components from the bare mesh.
     * Does not map fields or internal variables.
     * @todo Placing it in a new domain is probably preferable.
     */
    virtual void replaceFEMesh() // (Domain *& newDomain)
    { OOFEM_ERROR2("%s::replaceFEMesh - Not implemented",this->giveClassName()); }

    /**
     * Changes the connected domain of receiver.
     * @note{Does not delete any existing objects.}
     */
    virtual void setDomain(Domain *newDomain) { d = newDomain; }

    /**
     * Gives the name of the class.
     */
    virtual const char *giveClassName() const { return "TopologyDescription"; }
    /**
     * Gives the class ID.
     */
    virtual classType giveClassID() const { return TopologyDescriptionClass; }
};

} // end namespace oofem
#endif // topologydescription_h
