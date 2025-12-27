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

#ifndef dmexportmodule_h
#define dmexportmodule_h

#include "exportmodule.h"

///@name Input fields for DofManExportModule
//@{
#define _IFT_DofManExportModule_Name "dm"
#define _IFT_DofManExportModule_dmlist "dmlist"
//@}

namespace oofem {
/**
 * Represents DofManager export module.
 * This module writes the coordinates of all dof managers
 * along with the values of displacements
 * for further processing.
 * @author Milan Jirasek
 */
class DofManExportModule : public ExportModule
{
protected:
    IntArray dofManList;

public:
    /// Constructor
    DofManExportModule(int n, EngngModel * e);

    /// Destructor
    virtual ~DofManExportModule();

    void initializeFrom(InputRecord &ir) override;
    void doOutput(TimeStep *tStep, bool forcedOutput = false) override;
    const char *giveClassName() const override { return "DofManExportModuleClass"; }
    const char *giveInputRecordName() const { return _IFT_DofManExportModule_Name; }

protected:
    FILE *giveOutputStream(TimeStep *tStep);
};
} // end namespace
#endif
