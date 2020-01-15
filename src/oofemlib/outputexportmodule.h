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
#define _IFT_OutputExportModule_nodeSets "node_sets"
#define _IFT_OutputExportModule_elementSets "element_sets"
//@}

namespace oofem {
class Domain;
class Element;
class DofManager;

/**
 * Standard output for OOFEM. Most available data is written in plain text.
 * Implementation simply relies on EngngModel::printOutputAt
 *
 * @author Mikael Öhman
 */
class OOFEM_EXPORT OutputExportModule : public ExportModule
{
protected:
    FILE *outputStream;

    /// Set which contains nodes which should be exported
    IntArray nodeSets;

    /// Set which contains elements which should be exported
    IntArray elementSets;

public:
    OutputExportModule(int n, EngngModel * e);
    virtual ~OutputExportModule() {}

    void initializeFrom(InputRecord &ir) override;
    FILE *giveOutputStream();

    void doOutput(TimeStep *tStep, bool forcedOutput = false) override;
    void terminate() override;

    const char *giveClassName() const override { return "OutputExportModule"; }
    const char *giveInputRecordName() const { return _IFT_OutputExportModule_Name; }
};
} // end namespace oofem
#endif // outputexportmodule_h_
