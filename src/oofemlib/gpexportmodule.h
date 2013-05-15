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

#ifndef gpexportmodule_h_
#define gpexportmodule_h_

#include "exportmodule.h"

#include <cstdio>

///@name Input fields for Gausspoint export module
//@{
#define _IFT_GPExportModule_Name "gp"
#define _IFT_GPExportModule_vartypes "vars"
#define _IFT_GPExportModule_ncoords "ncoords"
//@}

namespace oofem {
/**
 * Represents GP (Gauss point) export module.
 * This module writes the coordinates of all Gauss points
 * along with the values of certain internal variables
 * for further processing.
 */
class GPExportModule : public ExportModule
{
protected:
    /// Identification numbers of variables to be exported
    IntArray vartypes;
    /// Number of coordinates to be exported (at each Gauss point)
    int ncoords;

public:
    /// Constructor. Creates empty Output Manager. By default all components are selected.
    GPExportModule(int n, EngngModel *e);
    /// Destructor
    virtual ~GPExportModule();

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void doOutput(TimeStep *tStep);
    virtual void initialize();
    virtual void terminate();
    virtual const char *giveClassName() const { return "GPExportModuleClass"; }

protected:
    /// Returns the output stream for given solution step
    FILE *giveOutputStream(TimeStep *tStep);
};

} // namespace oofem

#endif // gpexportmodule_h_
