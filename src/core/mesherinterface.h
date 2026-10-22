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

#ifndef mesherinterface_h
#define mesherinterface_h

#include "oofemenv.h"
#include "inputrecord.h"

#include <string>

///@name Input fields for MesherInterface
//@{
#define _IFT_MesherInterface_remeshcmd "remeshcmd"
//@}

namespace oofem {
class Domain;
class TimeStep;

/**
 * The base class representing the interface to mesh generation package.
 * This interface is primarily responsible for two main tasks:
 * - to create input mesher file, containing all information including the mesh density informations
 *   based on informations from remeshing criteria.
 * - possibly to launch the mesher and transform its output to oofem input
 *
 * All mesher specific details are hidden behind this interface. The client (typically adaptive
 * engineering model) only calls createMesh and receives either a new domain (MI_OK)
 * or information, that the mesh density file has been written and the mesh has to be
 * generated externally (MI_NEEDS_EXTERNAL_ACTION).
 *
 * Meshers, that write a mesh density file for an external mesh generator, can run
 * the external mesher and the converter to oofem format themselves, when the remesh command
 * is provided (keyword remeshcmd). The command is expected to create the new domain
 * file (in oofem .din format, i.e. starting with domain record) and may contain following
 * placeholders, that are substituted before execution:
 * - %d mesh density file written by the mesher
 * - %o domain file to be created
 * - %s serial number of the new domain
 * - %n solution step number
 * - %% the percent character
 */
class OOFEM_EXPORT MesherInterface
{
protected:
    Domain *domain;
    /// External remeshing command (template), empty if not provided.
    std :: string remeshCmd;

public:
    enum returnCode { MI_OK, MI_NEEDS_EXTERNAL_ACTION, MI_FAILED };
    /// Constructor
    MesherInterface(Domain * d) {
        domain = d;
    }
    /// Destructor
    virtual ~MesherInterface() { }

    /**
     * Runs the mesh generation, mesh will be written to corresponding domain din file.
     * @param tStep Time step.
     * @param domainNumber New domain number.
     * @param domainSerNum New domain serial number.
     * @param dNew Newly allocated domain, representing new mesh or set to NULL if external generation has to be performed.
     */
    virtual returnCode createMesh(TimeStep *tStep, int domainNumber, int domainSerNum, Domain **dNew) = 0;
    /**
     * Initializes receiver according to object description stored in input record.
     * This function is called immediately after creating object using
     * constructor. Input record can be imagined as data record in component database
     * belonging to receiver. Receiver may use value-name extracting functions
     * to extract particular field from record.
     */
    virtual void initializeFrom(const std::shared_ptr<InputRecord> &ir);
    /// Sets the domain the receiver operates on (used when the domain is replaced by a new one).
    virtual void setDomain(Domain *d) { domain = d; }

protected:
    /**
     * Runs the external remesh command (if provided) and instanciates the new domain from the
     * domain file created by the command.
     * @param densityFile Name of the mesh density file written by the mesher.
     * @param tStep Time step.
     * @param domainNumber New domain number.
     * @param domainSerNum New domain serial number.
     * @param dNew Newly allocated domain or NULL.
     * @return MI_OK if new domain created, MI_NEEDS_EXTERNAL_ACTION if no remesh command given, MI_FAILED otherwise.
     */
    returnCode remeshExternally(const std :: string &densityFile, TimeStep *tStep, int domainNumber, int domainSerNum, Domain **dNew);
};
} // end namespace oofem
#endif // mesherinterface_h
