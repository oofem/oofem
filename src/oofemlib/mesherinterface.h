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

//   ******************************
//   *** CLASS MESHER INTERFACE ***
//   ******************************

#ifndef mesherinterface_h
#define mesherinterface_h

#include "inputrecord.h"

namespace oofem {
class Domain;
class TimeStep;

/**
 * The base class representing the interface to mesh generation package.
 * This interface is primarily responsible for two main tasks:
 * - to create input mesher file, containing all information including the mesh density informations
 *   based on informations from remeshing criteria.
 * - possibly to launch the mesher and transform its output to oofem input
 */
class MesherInterface
{
protected:
    Domain *domain;
public:
    enum returnCode { MI_OK, MI_NEEDS_EXTERNAL_ACTION, MI_FAILED };
    /// Constructor
    MesherInterface(Domain *d) { domain = d; }
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
    virtual IRResultType initializeFrom(InputRecord *ir) { return IRRT_OK; }
};
} // end namespace oofem
#endif // mesherinterface_h
