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
 *               Copyright (C) 1993 - 2019   Borek Patzak
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

    void setupVTKPiece(VTKPiece &vtkPiece, TimeStep *tStep, int region) override;

    int initRegionNodeNumbering(IntArray &mapG2L, IntArray &mapL2G,
                                int &regionDofMans,
                                int &totalcells,
                                Domain *domain, TimeStep *tStep, int reg) override;

    void exportPrimaryVars(VTKPiece &vtkPiece, IntArray &mapG2L, IntArray &mapL2G, int region, TimeStep *tStep) override;
    void exportIntVars(VTKPiece &vtkPiece, IntArray &mapG2L, IntArray &mapL2G, int region, TimeStep *tStep) override;
};
} // end namespace oofem
#endif // vtkxmlperiodicexportmodule_h
