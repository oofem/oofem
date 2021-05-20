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

#ifndef vtkxmllatticeexportmodule_h
#define vtkxmllatticeexportmodule_h

#include "vtkxmlexportmodule.h"
#include "floatmatrix.h"


///@name Input fields for VTK XML Lattice export module
//@{
#define _IFT_VTKXMLLatticeExportModule_Name "vtkxmllattice"
#define _IFT_VTKXMLLatticeExportModule_cross "cross"
//@}

namespace oofem {
/**
 *
 * Represents VTK (Visualization Toolkit) export module for lattice modes. It uses VTK (.vtu) file format
 * It extends vtkxmlexportmodule and uses ideas of vtkxmlperiodicmodule by introducing methods to deal with boundary elements of lattice elements.
 * It also provides the option to plot mid cross-sections for visualisation of crack patterns
 * @author: Peter Grassl
 *
 */
class OOFEM_EXPORT VTKXMLLatticeExportModule : public VTKXMLExportModule
{
protected:

    FloatMatrix uniqueNodeTable;
    FloatMatrix uniqueNodeTableCross;
    IntArray periodicMap, regionToUniqueMap;
    IntArray locationMap;
    IntArray elementNodePeriodicMap;
    int elemNodes;

    bool crossSectionExportFlag;

    void giveSwitches(IntArray &answer, int location);

public:

    /// Constructor. Creates empty Output Manager. By default all components are selected.
    VTKXMLLatticeExportModule(int n, EngngModel *e);
    /// Destructor
    virtual ~VTKXMLLatticeExportModule();

    std::ofstream fileStreamCross;

    VTKPiece defaultVTKPieceCross;

    void initializeFrom(InputRecord &ir) override;

    std::string giveOutputFileNameCross(TimeStep *tStep);

    std::ofstream giveOutputStreamCross(TimeStep *tStep);

    void doOutput(TimeStep *tStep, bool forcedOutput = false) override;

    void doOutputNormal(TimeStep *tStep, bool forcedOutput = false);

    void doOutputCross(TimeStep *tStep, bool forcedOutput = false);

    const char *giveClassName() const override { return "VTKXMLLatticeExportModule"; }

    bool writeVTKPieceCross(VTKPiece &vtkPiece, TimeStep *tStep);

    void setupVTKPiece(VTKPiece &vtkPiece, TimeStep *tStep, Set& region) override;

    void writeCellVarsCross(VTKPiece &vtkPiece);

    void setupVTKPieceCross(VTKPiece &vtkPiece, TimeStep *tStep, Set& region);

    int initRegionNodeNumbering(VTKPiece& piece, Domain *domain, TimeStep *tStep, Set& region) override;

    void exportPrimaryVars(VTKPiece &piece, Set& region, IntArray& primaryVarsToExport, NodalRecoveryModel& smoother, TimeStep *tStep) override;
    void exportIntVars(VTKPiece &piece, Set& region, IntArray& internalVarsToExport, NodalRecoveryModel& smoother, TimeStep *tStep) override;
};
} // end namespace oofem
#endif // vtkxmllatticeexportmodule_h
