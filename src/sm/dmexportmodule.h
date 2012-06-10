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

#ifndef dmexportmodule_h
#define dmexportmodule_h

#include "exportmodule.h"
#include "domain.h"
#include "engngm.h"

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

public:
    /// Constructor
    DofManExportModule(int n, EngngModel *e);

    /// Destructor
    virtual ~DofManExportModule();

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void doOutput(TimeStep *tStep);
    virtual const char *giveClassName() const { return "DofManExportModuleClass"; }

protected:
    FILE *giveOutputStream(TimeStep *tStep);
};
} // end namespace
#endif



