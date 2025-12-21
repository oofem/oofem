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

#ifndef vtkxmlperiodicexportmodule_h
#define vtkxmlperiodicexportmodule_h

#include "vtkxmlexportmodule.h"
#include "floatmatrix.h"


///@name Input fields for VTK XML export module
//@{
#define _IFT_VTKXMLPeriodicExportModule_Name "vtkxmlperiodic"
//@}

namespace oofem {
/**
 *
 *
 * Represents VTK (Visualization Toolkit) export module for periodic cells. It uses VTK (.vtu) file format
 * It extends vtkxmlexportmodule by introducing methods to deal with boundary elements.
 * @author: Ignatios Athanasiadis, Adam Sciegaj, Peter Grassl
 *
 */
class OOFEM_EXPORT VTKXMLPeriodicExportModule : public VTKXMLExportModule
{
protected:
    FloatMatrix uniqueNodeTable;
    IntArray periodicMap, regionToUniqueMap;
    IntArray locationMap;
    IntArray elementNodePeriodicMap;
    int elemNodes;

    void giveSwitches(IntArray &answer, int location);

public:
    /// Constructor. Creates empty Output Manager. By default all components are selected.
    VTKXMLPeriodicExportModule(int n, EngngModel *e);
    /// Destructor
    virtual ~VTKXMLPeriodicExportModule();

    void initializeFrom(InputRecord &ir) override;

    const char *giveClassName() const override { return "VTKXMLPeriodicExportModule"; }

    void setupVTKPiece(ExportRegion &vtkPiece, TimeStep *tStep, Set& region) override;

    int initRegionNodeNumbering(ExportRegion& piece,
                                Domain *domain, TimeStep *tStep, Set& region) override;

    void exportPrimaryVars(ExportRegion &vtkPiece, Set& region, IntArray& primaryVarsToExport, NodalRecoveryModel& smoother, TimeStep *tStep) override;
    void exportIntVars(ExportRegion &vtkPiece, Set& region, IntArray& internalVarsToExport, NodalRecoveryModel& smoother, TimeStep *tStep) override;
};
} // end namespace oofem
#endif // vtkxmlperiodicexportmodule_h
