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

#ifndef outputexportmodule_h_
#define outputexportmodule_h_

#include <vector>

#include "exportmodule.h"

///@name Input fields for OutputExportModule
//@{
#define _IFT_OutputExportModule_Name "output"
#define _IFT_OutputExportModule_nodeSets "dofman_sets"
#define _IFT_OutputExportModule_elementSets "element_sets"
//@}

namespace oofem {
class Domain;
class Element;
class DofManager;

/**
 * Standard output for OOFEM. Most available data is written in plain text.
 * 
 * Adapted from OutputManager.
 *
 * @author Mikael Öhman (and others)
 */
class OOFEM_EXPORT OutputExportModule : public ExportModule
{
protected:
    /// Set which contains nodes which should be exported
    IntArray nodeSets;

    /// Set which contains elements which should be exported
    IntArray elementSets;

    /**
     * Does the dofmanager output.
     * All selected dofmanagers are requested for doing their output using printOutputAt service.
     */
    void doDofManOutput(FILE *file, Domain *domain, TimeStep *tStep);
    /**
     * Does the element output.
     * All selected elements are requested for doing their output using printOutputAt service.
     */
    void doElementOutput(FILE *file, Domain *domain, TimeStep *tStep);

    /**
     * Tests if given dof manager is required to do its output for given time step.
     * @return nonzero if output required.
     */
    int testDofManOutput(DofManager *dman);
    /**
     * Tests if given element is required to do its output for given time step.
     * @return nonzero if output required.
     */
    int testElementOutput(Element *element);

public:
    OutputExportModule(int n, EngngModel * e);
    virtual ~OutputExportModule() {}
    virtual IRResultType initializeFrom(InputRecord *ir);
    FILE *giveOutputStream();

    virtual void doOutput(TimeStep *tStep, bool forcedOutput = false);
    virtual void terminate();

    virtual const char *giveClassName() const { return "OutputExportModule"; }
    virtual const char *giveInputRecordName() const { return _IFT_OutputExportModule_Name; }
};
} // end namespace oofem
#endif // outputexportmodule_h_
